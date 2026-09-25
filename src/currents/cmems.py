"""Copernicus Marine Service (CMEMS) via the official toolbox.

CMEMS needs a free account (https://marine.copernicus.eu/register) and the
``copernicusmarine`` toolbox. The engine works fully without it for the
NOAA path; requesting CMEMS without credentials raises a clear,
actionable error instead of hanging on an interactive login prompt.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .models import CurrentField
from .provenance import write_provenance

REGISTER_URL = "https://marine.copernicus.eu/register"


@dataclass(frozen=True)
class CmemsPreset:
    """A named CMEMS dataset + variable selection."""

    name: str
    dataset_id: str
    variables: Tuple[str, ...]      # (u, v, temperature) variable names
    description: str = ""


#: Curated presets. Any other dataset id the toolbox recognises can be
#: passed directly to :func:`subset_cmems` via a custom :class:`CmemsPreset`.
CMEMS_PRESETS: Dict[str, CmemsPreset] = {
    "global-physics-daily": CmemsPreset(
        "global-physics-daily",
        "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
        ("uo", "vo", "thetao"),
        "Global ocean physics analysis+forecast, 1/12 deg (~9 km), daily."),
    "global-physics-hourly": CmemsPreset(
        "global-physics-hourly",
        "cmems_mod_glo_phy_anfc_0.083deg_PT1H-m",
        ("uo", "vo", "thetao"),
        "Global ocean physics analysis+forecast, 1/12 deg (~9 km), hourly."),
}


def get_preset(name: str) -> CmemsPreset:
    try:
        return CMEMS_PRESETS[name]
    except KeyError as exc:
        known = ", ".join(sorted(CMEMS_PRESETS))
        raise KeyError(f"unknown CMEMS preset {name!r}; known: {known}") from exc


def require_toolbox():
    """Import copernicusmarine or raise an actionable error.

    Never triggers an interactive login: if credentials are missing the
    toolbox would prompt, so we check for them first and fail fast with
    setup instructions.
    """
    try:
        import copernicusmarine
    except ImportError as exc:
        raise ImportError(
            "CMEMS access needs the copernicusmarine toolbox "
            "(pip install 'survey-currents[cmems]')"
        ) from exc
    cred_file = os.path.expanduser("~/.copernicusmarine/.copernicusmarine-credentials")
    env_ok = bool(os.environ.get("COPERNICUSMARINE_SERVICE_USERNAME"))
    if not env_ok and not os.path.exists(cred_file):
        raise RuntimeError(
            "CMEMS credentials not found. CMEMS is free: register at "
            f"{REGISTER_URL}, then run "
            "`copernicusmarine login` (or set COPERNICUSMARINE_SERVICE_USERNAME "
            "and COPERNICUSMARINE_SERVICE_PASSWORD). The NOAA S3 path needs "
            "no account at all."
        )
    return copernicusmarine


def subset_cmems(preset: str,
                 bbox: Sequence[float],
                 start: str,
                 end: str,
                 work_dir: str,
                 variables: Optional[Sequence[str]] = None,
                 minimum_depth: float = 0.0,
                 maximum_depth: float = 0.5) -> str:
    """Subset a CMEMS dataset to bbox/time and download one NetCDF.

    Returns the local NetCDF path; a provenance sidecar is written
    alongside. ``bbox`` is ``(min_lon, min_lat, max_lon, max_lat)``;
    ``start``/``end`` are ISO datetimes.
    """
    cm = require_toolbox()
    pre = get_preset(preset)
    varlist = list(variables) if variables else list(pre.variables)
    minx, miny, maxx, maxy = (float(x) for x in bbox)
    os.makedirs(work_dir, exist_ok=True)
    out_name = (f"cmems_{pre.name}_{start[:10]}_{end[:10]}"
                f"_{minx}_{miny}_{maxx}_{maxy}.nc").replace(" ", "_")
    out_path = os.path.join(work_dir, out_name)
    cm.subset(
        dataset_id=pre.dataset_id,
        variables=varlist,
        minimum_longitude=minx, maximum_longitude=maxx,
        minimum_latitude=miny, maximum_latitude=maxy,
        minimum_depth=minimum_depth, maximum_depth=maximum_depth,
        start_datetime=start, end_datetime=end,
        output_filename=os.path.basename(out_path),
        output_directory=work_dir,
    )
    write_provenance(out_path, {
        "source": f"cmems:{pre.dataset_id}",
        "preset": pre.name,
        "variables": varlist,
        "bbox": [minx, miny, maxx, maxy],
        "time_window": [start, end],
        "depth_window": [minimum_depth, maximum_depth],
    })
    return out_path


def parse_cmems_netcdf(path: str,
                       u_name: str = "uo", v_name: str = "vo",
                       temp_name: str = "thetao") -> CurrentField:
    """Parse a CMEMS subset NetCDF (surface level) into a CurrentField."""
    try:
        import xarray as xr
    except ImportError as exc:
        raise ImportError(
            "parsing CMEMS NetCDF needs xarray/netCDF4 "
            "(pip install 'survey-currents[noaa]')"
        ) from exc
    ds = xr.open_dataset(path)
    try:
        def surface(da):
            for dim in ("depth", "z", "level"):
                if dim in da.dims:
                    return da.isel({dim: 0})
            return da

        u = surface(ds[u_name]).astype("float64")
        v = surface(ds[v_name]).astype("float64")
        temp = surface(ds[temp_name]).astype("float64") if temp_name in ds else None

        lat = np.asarray(ds["latitude"].values, dtype=float).ravel()
        lon = np.asarray(ds["longitude"].values, dtype=float).ravel()
        lat_o = np.argsort(lat)
        lon_o = np.argsort(lon)

        def grid(da):
            arr = np.squeeze(np.asarray(da.values, dtype=float))
            if arr.ndim == 2:
                arr = arr[None, :, :]
            if arr.ndim != 3:
                raise ValueError(
                    f"cannot parse {da.name!r} with shape {arr.shape}; "
                    "expected (time, lat, lon)")
            return arr[:, lat_o[:, None], lon_o]

        times = [str(t) for t in
                 np.atleast_1d(ds["time"].values)]
        return CurrentField(
            u=grid(u), v=grid(v),
            temperature=None if temp is None else grid(temp),
            times=times, lats=np.sort(lat), lons=np.sort(lon),
            source=f"cmems:{ds.attrs.get('dataset_id', 'subset')}",
            provenance={"variables": {"u": u_name, "v": v_name,
                                      "temperature": temp_name},
                        "netcdf_path": os.path.basename(path)})
    finally:
        ds.close()

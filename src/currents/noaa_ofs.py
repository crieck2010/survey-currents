"""NOAA Operational Forecast Systems via anonymous AWS S3 (no signup).

The NOAA OFS model NetCDFs live in two public buckets (verified
2026-09-25; see docs/DATA_SOURCES.md):

* ``noaa-ofs-pds`` — NOMADS production runs, last 30 days, prefix
  ``OFS.YYYYMMDD``.
* ``noaa-nos-ofs-pds`` — CO-OPS operational archive, prefix
  ``OFS/netcdf/YYYYMM/``.

Both allow unsigned requests. This module lists and downloads with the
stdlib (``urllib`` + S3's ListObjectsV2 XML API) and only uses ``boto3``
when it is installed. Filename conventions differ per OFS, so the
engine lists the date prefix and *filters keys by OFS-code, date,
cycle and forecast-hour tokens* instead of assuming an exact filename.
"""

from __future__ import annotations

import datetime as _dt
import os
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .models import CurrentField, get_ofs_model
from .provenance import write_provenance

BUCKET_RECENT = "noaa-ofs-pds"
BUCKET_ARCHIVE = "noaa-nos-ofs-pds"

ARCHIVE_PREFIX_TEMPLATE = "OFS/netcdf/{yyyymm}/"   # in noaa-nos-ofs-pds
RECENT_PREFIX_TEMPLATE = "OFS.{yyyymmdd}"          # in noaa-ofs-pds

# Candidate NetCDF variable names per physical quantity. NOS OFS files
# use different spellings across models, so resolution is tolerant and
# the chosen names are recorded in provenance.
U_CANDIDATES = ("u", "u_velocity", "eastward_velocity", "water_u", "uo")
V_CANDIDATES = ("v", "v_velocity", "northward_velocity", "water_v", "vo")
TEMP_CANDIDATES = ("temp", "temperature", "water_temp", "water_temperature",
                   "thetao", "sst", "surface_temperature")
LAT_CANDIDATES = ("lat", "latitude", "y")
LON_CANDIDATES = ("lon", "longitude", "x")
TIME_CANDIDATES = ("time", "date", "forecast_time")


# ---------------------------------------------------------------------------
# Path / key construction
# ---------------------------------------------------------------------------


def archive_prefix(date: str) -> str:
    """S3 prefix for ``noaa-nos-ofs-pds`` given an ISO date."""
    yyyymm = _dt.date.fromisoformat(date).strftime("%Y%m")
    return ARCHIVE_PREFIX_TEMPLATE.format(yyyymm=yyyymm)


def recent_prefix(date: str) -> str:
    """S3 prefix for ``noaa-ofs-pds`` given an ISO date."""
    yyyymmdd = _dt.date.fromisoformat(date).strftime("%Y%m%d")
    return RECENT_PREFIX_TEMPLATE.format(yyyymmdd=yyyymmdd)


def _tokens(ofs_code: str, date: str, cycle: str,
            forecast_hour: Optional[int]) -> List[str]:
    """Lowercase tokens a key must contain to match the request."""
    compact = date.replace("-", "")
    toks = [ofs_code.lower(), compact, f"t{cycle}z"]
    if forecast_hour is not None:
        toks.append(f"f{forecast_hour:03d}")
    return toks


def filter_keys(keys: Sequence[str], ofs_code: str, date: str, cycle: str,
                forecast_hour: Optional[int] = None) -> List[str]:
    """Keep keys matching the OFS/date/cycle/(hour) tokens, sorted."""
    toks = _tokens(ofs_code, date, cycle, forecast_hour)
    out = [k for k in keys
           if k.lower().endswith(".nc") and all(t in k.lower() for t in toks)]
    return sorted(out)


# ---------------------------------------------------------------------------
# Anonymous S3 access (stdlib; boto3 only if installed)
# ---------------------------------------------------------------------------


def _boto3_client():
    try:
        import boto3
        from botocore import UNSIGNED
        from botocore.config import Config
    except ImportError:
        return None
    return boto3.client("s3", config=Config(signature_version=UNSIGNED),
                        region_name="us-east-1")


def s3_list(bucket: str, prefix: str) -> List[str]:
    """List keys under ``prefix`` with unsigned requests.

    Uses boto3 when installed, otherwise the S3 ListObjectsV2 REST API
    over plain HTTPS (stdlib only).
    """
    client = _boto3_client()
    if client is not None:
        keys: List[str] = []
        token: Optional[str] = None
        while True:
            kwargs: Dict[str, Any] = {"Bucket": bucket, "Prefix": prefix}
            if token:
                kwargs["ContinuationToken"] = token
            resp = client.list_objects_v2(**kwargs)
            keys.extend(o["Key"] for o in resp.get("Contents", ()))
            token = resp.get("NextContinuationToken")
            if not token:
                return keys
    # stdlib fallback: ListObjectsV2 XML
    keys = []
    token = ""
    while True:
        query = urllib.parse.urlencode(
            {"list-type": "2", "prefix": prefix,
             **({"continuation-token": token} if token else {})})
        url = f"https://{bucket}.s3.amazonaws.com/?{query}"
        with urllib.request.urlopen(url, timeout=60) as resp:
            root = ET.fromstring(resp.read())
        ns = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
        keys.extend(e.text or "" for e in root.findall("s3:Contents/s3:Key", ns))
        trunc = root.findtext("s3:IsTruncated", namespaces=ns)
        if trunc != "true":
            return keys
        token = root.findtext("s3:NextContinuationToken", namespaces=ns) or ""


def s3_download(bucket: str, key: str, dest_path: str) -> str:
    """Download one key with unsigned requests; return ``dest_path``."""
    client = _boto3_client()
    os.makedirs(os.path.dirname(os.path.abspath(dest_path)), exist_ok=True)
    if client is not None:
        client.download_file(bucket, key, dest_path)
        return dest_path
    url = f"https://{bucket}.s3.amazonaws.com/{urllib.parse.quote(key)}"
    with urllib.request.urlopen(url, timeout=300) as resp, \
            open(dest_path, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    return dest_path


# ---------------------------------------------------------------------------
# High-level fetch
# ---------------------------------------------------------------------------


def list_field_files(ofs_code: str, date: str, cycle: str = "00",
                     bucket: str = BUCKET_ARCHIVE) -> List[str]:
    """List candidate NetCDF keys for one OFS model run (network)."""
    model = get_ofs_model(ofs_code)  # validates the code
    if cycle not in model.cycles:
        raise ValueError(
            f"cycle {cycle!r} not in {model.code} cycles {list(model.cycles)}")
    prefix = archive_prefix(date) if bucket == BUCKET_ARCHIVE else recent_prefix(date)
    return s3_list(bucket, prefix)


def select_files(keys: Sequence[str], ofs_code: str, date: str, cycle: str,
                 hours: Sequence[int]) -> Dict[int, str]:
    """Map each requested forecast hour to its best-matching key."""
    model = get_ofs_model(ofs_code)
    out: Dict[int, str] = {}
    for hour in hours:
        if not (0 <= hour <= model.horizon_hours):
            raise ValueError(
                f"forecast hour {hour} outside {model.code} horizon "
                f"({model.horizon_hours} h)")
        matches = filter_keys(keys, ofs_code, date, cycle, hour)
        if not matches:
            raise FileNotFoundError(
                f"no NetCDF for {ofs_code} {date} cycle {cycle} "
                f"forecast hour f{hour:03d}")
        out[hour] = matches[0]
    return out


def download_fields(ofs_code: str, date: str, cycle: str,
                    hours: Sequence[int], work_dir: str,
                    bucket: str = BUCKET_ARCHIVE) -> List[str]:
    """Download one NetCDF per forecast hour + provenance sidecars.

    Returns the local file paths, ordered by forecast hour.
    """
    keys = list_field_files(ofs_code, date, cycle, bucket)
    chosen = select_files(keys, ofs_code, date, cycle, hours)
    os.makedirs(work_dir, exist_ok=True)
    paths = []
    for hour in sorted(chosen):
        key = chosen[hour]
        dest = os.path.join(work_dir, os.path.basename(key))
        s3_download(bucket, key, dest)
        write_provenance(dest, {
            "source": f"s3://{bucket}/{key}",
            "ofs_code": ofs_code.upper(),
            "date": date, "cycle": cycle, "forecast_hour": hour,
            "bbox": None, "time_window": None,
        })
        paths.append(dest)
    return paths


# ---------------------------------------------------------------------------
# NetCDF -> CurrentField (xarray, lazy import)
# ---------------------------------------------------------------------------


def _require_xarray():
    try:
        import xarray as xr
    except ImportError as exc:
        raise ImportError(
            "parsing OFS NetCDF needs xarray/netCDF4 "
            "(pip install 'survey-currents[noaa]')"
        ) from exc
    return xr


def _resolve(ds, candidates: Sequence[str]) -> str:
    for name in candidates:
        if name in ds.variables or name in ds.coords:
            return name
    raise KeyError(
        f"none of {list(candidates)} found in dataset variables "
        f"{sorted(ds.variables)}")


def parse_ofs_netcdf(path: str,
                     u_name: Optional[str] = None,
                     v_name: Optional[str] = None,
                     temp_name: Optional[str] = None) -> CurrentField:
    """Parse one OFS NetCDF into a single-timestep :class:`CurrentField`.

    Takes the surface level when a depth dimension is present. Variable
    names are resolved tolerantly (see ``*_CANDIDATES``); explicit names
    override the search.

    Only regular lat/lon grids are supported: unstructured (e.g. FVCOM
    triangular) OFS grids raise a clear ``ValueError`` — regridding those
    is deferred (see README limitations).
    """
    xr = _require_xarray()
    ds = xr.open_dataset(path)
    try:
        u_var = u_name or _resolve(ds, U_CANDIDATES)
        v_var = v_name or _resolve(ds, V_CANDIDATES)
        t_var = temp_name
        if t_var is None:
            try:
                t_var = _resolve(ds, TEMP_CANDIDATES)
            except KeyError:
                t_var = None
        lat_var = _resolve(ds, LAT_CANDIDATES)
        lon_var = _resolve(ds, LON_CANDIDATES)

        # Regular-grid guard: u/v must be indexed by the lat/lon dims.
        u_dims = set(ds[u_var].dims)
        if lat_var not in u_dims or lon_var not in u_dims:
            raise ValueError(
                "unstructured OFS grid detected (u/v are not on a regular "
                "lat/lon grid, e.g. FVCOM triangular meshes). survey-currents "
                "v0.1.0 only parses regular-grid NetCDFs; regridding is "
                "planned for a later release.")

        def surface(da):
            # depth-like dims: take the first (surface) level
            for dim in ("depth", "z", "level", "s_rho", "siglay"):
                if dim in da.dims:
                    return da.isel({dim: 0})
            return da

        u = surface(ds[u_var]).astype("float64")
        v = surface(ds[v_var]).astype("float64")
        temp = surface(ds[t_var]).astype("float64") if t_var else None

        lat = np.asarray(ds[lat_var].values, dtype=float).ravel()
        lon = np.asarray(ds[lon_var].values, dtype=float).ravel()
        # OFS grids are sometimes descending; CurrentField wants increasing.
        lat_order = np.argsort(lat); lon_order = np.argsort(lon)
        lat, lon = lat[lat_order], lon[lon_order]

        def grid(da):
            arr = np.asarray(da.values, dtype=float)
            # squeeze singleton dims, then index the last two dims as (lat, lon)
            arr = np.squeeze(arr)
            if arr.ndim > 2:  # leftover time dim of length 1
                arr = arr[0]
            return arr[np.ix_(lat_order, lon_order)][None, :, :]

        times = ["unknown"]
        try:
            t_varname = _resolve(ds, TIME_CANDIDATES)
            vals = ds[t_varname].values
            vals = np.atleast_1d(vals)
            times = [str(np.datetime_as_string(np.datetime64(v), unit="s"))
                     for v in vals[:1]]
        except KeyError:
            pass

        prov = {
            "variables": {"u": u_var, "v": v_var, "temperature": t_var},
            "netcdf_path": os.path.basename(path),
        }
        return CurrentField(
            u=grid(u), v=grid(v),
            temperature=None if temp is None else grid(temp),
            times=times, lats=lat, lons=lon,
            source="noaa-ofs", provenance=prov)
    finally:
        ds.close()


def fetch_currents(ofs_code: str, date: str, cycle: str = "00",
                   hours: Sequence[int] = (0,),
                   work_dir: str = "currents_data",
                   bucket: str = BUCKET_ARCHIVE,
                   bbox: Optional[Sequence[float]] = None) -> CurrentField:
    """End-to-end: download OFS fields and stack them into a CurrentField.

    ``bbox`` optionally clips to ``(min_lon, min_lat, max_lon, max_lat)``.
    """
    paths = download_fields(ofs_code, date, cycle, hours, work_dir, bucket)
    steps = []
    for path, hour in zip(paths, sorted(hours)):
        step = parse_ofs_netcdf(path)
        step.forecast_hours = [hour]
        step.times = [f"{date}T{cycle}:00:00+00:00+f{hour:03d}"]
        step.model_run = f"{date}T{cycle}:00:00+00:00"
        step.source = f"noaa-ofs/{ofs_code.upper()}"
        step.provenance.update(read_provenance_sidecar(path))
        steps.append(step)
    field = stack_steps(steps)
    if bbox is not None:
        field = field.select_bbox(bbox)
    return field


def read_provenance_sidecar(path: str) -> Dict[str, Any]:
    from .provenance import read_provenance
    try:
        return read_provenance(path)
    except FileNotFoundError:
        return {}


def stack_steps(steps: Sequence[CurrentField]) -> CurrentField:
    """Stack single-timestep fields along the time axis (same grid)."""
    if not steps:
        raise ValueError("no steps to stack")
    first = steps[0]
    for s in steps[1:]:
        if s.u.shape[1:] != first.u.shape[1:]:
            raise ValueError("grids differ between timesteps; cannot stack")
        if not np.array_equal(s.lats, first.lats) or not np.array_equal(s.lons, first.lons):
            raise ValueError("coordinates differ between timesteps; cannot stack")
    has_temp = all(s.temperature is not None for s in steps)
    return CurrentField(
        u=np.concatenate([s.u for s in steps], axis=0),
        v=np.concatenate([s.v for s in steps], axis=0),
        temperature=np.concatenate([s.temperature for s in steps], axis=0)
        if has_temp else None,
        times=[t for s in steps for t in s.times],
        lats=first.lats, lons=first.lons, crs=first.crs,
        source=first.source, model_run=first.model_run,
        forecast_hours=[h for s in steps for h in s.forecast_hours],
        provenance={"steps": [s.provenance for s in steps]},
    )

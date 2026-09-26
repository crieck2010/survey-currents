"""Copernicus ERA5 reanalysis acquisition via the Climate Data Store (CDS).

:func:`fetch_era5` pulls hourly single-level reanalysis fields
(1940-present, 0.25-degree global grid) through the ``cdsapi`` package:

* ``"wind"`` — 10-m u/v wind components (rendered as wind speed)
* ``"msl"`` — mean sea-level pressure (Pa on the wire, hPa here)
* ``"t2m"`` — 2-m air temperature (Kelvin on the wire, °C here)
* ``"tp"`` — total precipitation (m on the wire, mm here)

Access needs a free Copernicus CDS account: without credentials a
:class:`CredentialsMissing` error explains the exact setup. Requests are
chunked by calendar month (CDS request-size limits) and one request is
issued per month per longitude window. ``cdsapi`` is a lazy import — the
module imports cleanly without it.

Both the field model (:class:`Era5Field`) and the provenance shape follow
the :mod:`currents.sst_global` conventions, and :class:`Era5Field`
duck-types into what ``survey-viz``'s ``render_viz`` consumes
(``.times``/``.lats``/``.lons`` plus a 3-D ``.values`` grid, plus
``.overlay_grids`` for contour overlays such as isobars).
"""

from __future__ import annotations

import calendar
import datetime as _dt
import hashlib
import os
import tempfile
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_datetime_utc
from .sst_global import CredentialsMissing as _BaseCredentialsMissing
from .sst_global import validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]

# ---------------------------------------------------------------------------
# ERA5 / CDS constants
# ---------------------------------------------------------------------------

#: CDS dataset id for the hourly single-level reanalysis.
ERA5_DATASET = "reanalysis-era5-single-levels"
#: ERA5 record start (hourly single levels, 1940-present).
ERA5_START = _dt.date(1940, 1, 1)
#: Native grid spacing, degrees.
ERA5_RES = 0.25
#: CDS request-size hygiene: one request per calendar month.
ERA5_CHUNK = "monthly"

#: Short variable key -> CDS request names, on-wire grid names, units, and
#: the linear conversion applied on ingest (value * scale + offset).
ERA5_VARIABLES: Dict[str, Dict[str, Any]] = {
    "wind": {
        "cds": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
        "grids": ["u10", "v10"],
        "units": "m/s",
        "label": "10-m wind",
        "scale": 1.0,
        "offset": 0.0,
    },
    "msl": {
        "cds": ["mean_sea_level_pressure"],
        "grids": ["msl"],
        "units": "hPa",
        "label": "Mean sea-level pressure",
        "scale": 0.01,   # Pa -> hPa
        "offset": 0.0,
    },
    "t2m": {
        "cds": ["2m_temperature"],
        "grids": ["t2m"],
        "units": "\u00b0C",
        "label": "2-m air temperature",
        "scale": 1.0,
        "offset": -273.15,  # K -> degC
    },
    "tp": {
        "cds": ["total_precipitation"],
        "grids": ["tp"],
        "units": "mm",
        "label": "Total precipitation",
        "scale": 1000.0,  # m -> mm (per hourly step)
        "offset": 0.0,
    },
}

#: Grid keys that make up the rendered base field for each variable key.
_BASE_GRID_KEYS: Dict[str, List[str]] = {
    "wind": ["u10", "v10"],
    "msl": ["msl"],
    "t2m": ["t2m"],
    "tp": ["tp"],
}


class CredentialsMissing(_BaseCredentialsMissing):
    """Copernicus CDS credentials are required but were not found.

    ERA5 downloads need a free Copernicus Climate Data Store account:

    1. Register (free) at https://cds.climate.copernicus.eu/ and accept
       the "ERA5 hourly data on single levels" Terms of Use.
    2. Create ``~/.cdsapirc`` with your personal API token
       (https://cds.climate.copernicus.eu/how-to-api)::

           url: https://cds.climate.copernicus.eu/api
           key: <your-personal-api-token>

       — or set the ``CDSAPI_URL`` / ``CDSAPI_KEY`` environment variables
       instead (best for CI / servers).

    Also install the client: ``pip install cdsapi``.
    """


# ---------------------------------------------------------------------------
# Era5Field — canonical multi-variable ERA5 field model
# ---------------------------------------------------------------------------


@dataclass
class Era5Field:
    """Time-indexed ERA5 reanalysis grids.

    ``grids`` maps on-wire grid names (``"u10"``, ``"v10"``, ``"msl"``,
    ``"t2m"``, ``"tp"``) to ``(nt, ny, nx)`` numpy masked arrays in the
    converted units from :data:`ERA5_VARIABLES`. ``variables`` records the
    requested short keys (``"wind"``, ``"msl"``, ``"t2m"``, ``"tp"``);
    ``base_variable`` is the one the renderer draws as the base map
    (``"wind"`` renders as wind speed). Follows the :class:`SstField`
    conventions in ``sst_global.py``.
    """

    grids: Dict[str, np.ma.MaskedArray]  # grid name -> (nt, ny, nx)
    times: List[str]                     # ISO-8601 timestamps, one per step
    lats: np.ndarray                     # (ny,) degrees north, increasing
    lons: np.ndarray                     # (nx,) degrees east, -180..180, increasing
    variables: List[str]                 # requested short keys
    base_variable: str                   # rendered base variable
    units: Dict[str, str] = field(default_factory=dict)
    crs: str = "EPSG:4326"
    source: str = "copernicus-cds/era5"
    provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        ny = self.lats.shape[0]
        nx = self.lons.shape[0]
        for name, grid in self.grids.items():
            arr = np.ma.asarray(grid, dtype=float)
            if arr.shape != (nt, ny, nx):
                raise ValueError(
                    f"grid {name!r} shape {arr.shape} does not match "
                    f"(nt, ny, nx)=({nt}, {ny}, {nx})")
            self.grids[name] = arr
        self.lats = np.asarray(self.lats, dtype=float)
        self.lons = np.asarray(self.lons, dtype=float)
        unknown = [v for v in self.variables if v not in ERA5_VARIABLES]
        if unknown:
            raise ValueError(f"unknown ERA5 variable keys: {unknown}")
        if self.base_variable not in self.variables:
            raise ValueError(
                f"base_variable {self.base_variable!r} not in variables "
                f"{self.variables}")
        missing = [g for g in _BASE_GRID_KEYS[self.base_variable]
                   if g not in self.grids]
        if missing:
            raise ValueError(
                f"base_variable {self.base_variable!r} needs grids {missing}")
        if not self.units:
            self.units = {g: ERA5_VARIABLES[v]["units"]
                          for v in self.variables
                          for g in ERA5_VARIABLES[v]["grids"]}
        self.bounds = (float(self.lons[0]), float(self.lats[0]),
                       float(self.lons[-1]), float(self.lats[-1]))

    # -- derived quantities --------------------------------------------------

    @property
    def wind_speed(self) -> np.ma.MaskedArray:
        """Wind speed magnitude ``sqrt(u10^2 + v10^2)`` in m/s, (nt, ny, nx)."""
        if "u10" not in self.grids or "v10" not in self.grids:
            raise KeyError("wind_speed needs the 'wind' variable (u10/v10 grids)")
        return np.ma.sqrt(self.grids["u10"] ** 2 + self.grids["v10"] ** 2)

    @property
    def values(self) -> np.ma.MaskedArray:
        """The rendered base grid, (nt, ny, nx).

        Wind speed for ``base_variable="wind"``, else the single grid.
        This is what ``survey-viz``'s ``render_viz`` picks up via its
        documented ``values`` attribute.
        """
        if self.base_variable == "wind":
            return self.wind_speed
        return self.grids[_BASE_GRID_KEYS[self.base_variable][0]]

    @property
    def overlay_grids(self) -> Dict[str, np.ma.MaskedArray]:
        """Non-base variables as render-ready 3-D grids, keyed by variable.

        ``"wind"`` maps to wind speed; the rest map to their single grid.
        ``survey-viz`` draws ``spec.overlays`` entries from this mapping
        (e.g. ``"msl"`` as contour isobars over a wind base map).
        """
        out: Dict[str, np.ma.MaskedArray] = {}
        for var in self.variables:
            if var == self.base_variable:
                continue
            out[var] = (self.wind_speed if var == "wind"
                        else self.grids[_BASE_GRID_KEYS[var][0]])
        return out

    def overlay_grid(self, name: str, index: int) -> np.ma.MaskedArray:
        """2-D overlay grid for variable ``name`` at timestep ``index``."""
        grids = self.overlay_grids
        if name not in grids:
            raise KeyError(
                f"no overlay {name!r} on this field; available: "
                f"{sorted(grids)}")
        return grids[name][index]

    def spatial_mean(self, variable: str, index: int) -> float:
        """Masked/NaN-aware spatial mean of ``variable`` at timestep ``index``."""
        grid = (self.wind_speed if variable == "wind"
                else self.grids[_BASE_GRID_KEYS[variable][0]])
        return float(np.ma.masked_invalid(grid[index]).mean())

    # -- selection -----------------------------------------------------------

    def select_time(self, index: int) -> "Era5Field":
        """Return the single-timestep field at ``index``."""
        return Era5Field(
            grids={k: v[index:index + 1] for k, v in self.grids.items()},
            times=[self.times[index]], lats=self.lats, lons=self.lons,
            variables=list(self.variables), base_variable=self.base_variable,
            units=dict(self.units), crs=self.crs, source=self.source,
            provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "Era5Field":
        """Spatial subset to ``(min_lon, min_lat, max_lon, max_lat)``."""
        minx, miny, maxx, maxy = (float(x) for x in bbox)
        ix = np.where((self.lons >= minx) & (self.lons <= maxx))[0]
        iy = np.where((self.lats >= miny) & (self.lats <= maxy))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(
                f"bbox {list(bbox)} does not overlap field bounds {self.bounds}")
        sl = (slice(None), slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1))
        return Era5Field(
            grids={k: v[sl] for k, v in self.grids.items()},
            times=list(self.times),
            lats=self.lats[iy[0]:iy[-1] + 1], lons=self.lons[ix[0]:ix[-1] + 1],
            variables=list(self.variables), base_variable=self.base_variable,
            units=dict(self.units), crs=self.crs, source=self.source,
            provenance=dict(self.provenance),
        )

    # -- serialization -------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dict (masked cells become NaN)."""
        return {
            "grids": {k: np.ma.filled(v, np.nan).tolist()
                      for k, v in self.grids.items()},
            "times": list(self.times),
            "lats": self.lats.tolist(),
            "lons": self.lons.tolist(),
            "variables": list(self.variables),
            "base_variable": self.base_variable,
            "units": dict(self.units),
            "crs": self.crs,
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Era5Field":
        """Rebuild from :meth:`to_dict` (NaN -> masked). Raises on missing keys."""
        required = ("grids", "times", "lats", "lons", "variables",
                    "base_variable")
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(f"Era5Field dict missing keys: {missing}")
        return cls(
            grids={k: np.ma.masked_invalid(np.asarray(v, dtype=float))
                   for k, v in data["grids"].items()},
            times=list(data["times"]),
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            variables=list(data["variables"]),
            base_variable=data["base_variable"],
            units=dict(data.get("units", {})),
            crs=data.get("crs", "EPSG:4326"),
            source=data.get("source", "copernicus-cds/era5"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "Era5Field":
        import json
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture (tests, demos, offline use) -----------------------

    @classmethod
    def synthetic(
        cls,
        variables: Sequence[str] = ("wind", "msl"),
        nt: int = 4, ny: int = 6, nx: int = 8,
        lats: Sequence[float] = (-30.0, 30.0),
        lons: Sequence[float] = (-60.0, 60.0),
        start: str = "2024-01-01T12:00:00+00:00",
        step_hours: int = 6,
        seed: int = 7,
        source: str = "synthetic",
    ) -> "Era5Field":
        """Deterministic synthetic ERA5 field: westerlies + a low-pressure center.

        Stdlib+numpy only. Used by the test suite and the offline demo.
        """
        variables = list(variables)
        unknown = [v for v in variables if v not in ERA5_VARIABLES]
        if unknown:
            raise ValueError(f"unknown ERA5 variable keys: {unknown}")
        rng = np.random.default_rng(seed)
        lat = np.linspace(lats[0], lats[1], ny)
        lon = np.linspace(lons[0], lons[1], nx)
        lon2d, lat2d = np.meshgrid(lon, lat)
        grids: Dict[str, np.ma.MaskedArray] = {}
        t0 = _dt.datetime.fromisoformat(start)
        times = [(t0 + _dt.timedelta(hours=k * step_hours)).isoformat()
                 for k in range(nt)]
        if "wind" in variables:
            # Zonal westerlies strengthening poleward + a cyclonic swirl.
            u = 8.0 + 0.25 * np.abs(lat2d) + rng.normal(0, 1.5, (ny, nx))
            v = 4.0 * np.exp(-((lon2d / 40.0) ** 2 + (lat2d / 25.0) ** 2))
            v = v * np.sign(lon2d + 1e-9) + rng.normal(0, 1.0, (ny, nx))
            grids["u10"] = np.ma.array(
                np.repeat(u[None, :, :], nt, axis=0)
                + rng.normal(0, 0.5, (nt, ny, nx)),
                mask=np.zeros((nt, ny, nx), dtype=bool))
            grids["v10"] = np.ma.array(
                np.repeat(v[None, :, :], nt, axis=0)
                + rng.normal(0, 0.5, (nt, ny, nx)),
                mask=np.zeros((nt, ny, nx), dtype=bool))
        if "msl" in variables:
            # A low-pressure center drifting east, hPa.
            center_lon = -20.0 + 5.0 * np.arange(nt)[:, None, None]
            r2 = ((lon2d[None, :, :] - center_lon) / 25.0) ** 2 \
                + (lat2d[None, :, :] / 20.0) ** 2
            p = 1013.0 - 24.0 * np.exp(-r2) + rng.normal(0, 0.8, (nt, ny, nx))
            grids["msl"] = np.ma.array(
                p, mask=np.zeros((nt, ny, nx), dtype=bool))
        if "t2m" in variables:
            base = 28.0 - 0.6 * np.abs(lat2d)
            grids["t2m"] = np.ma.array(
                np.repeat(base[None, :, :], nt, axis=0)
                + rng.normal(0, 0.5, (nt, ny, nx)),
                mask=np.zeros((nt, ny, nx), dtype=bool))
        if "tp" in variables:
            rain = np.maximum(
                0.0, rng.normal(0.4, 1.2, (nt, ny, nx)))
            grids["tp"] = np.ma.array(
                rain, mask=np.zeros((nt, ny, nx), dtype=bool))
        base_variable = variables[0]
        return cls(grids=grids, times=times, lats=lat, lons=lon,
                   variables=variables, base_variable=base_variable,
                   source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# bbox / variable / date validation
# ---------------------------------------------------------------------------


def validate_era5_bbox(bbox: Sequence[float]) -> Tuple[float, float, float, float]:
    """Validate ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180 degrees.

    ERA5 is a global grid, so this is the shared ascending-bbox check from
    ``sst_global.validate_sst_bbox`` (antimeridian-crossing boxes are valid
    and split into two CDS requests internally).
    """
    return validate_sst_bbox(bbox)


def normalize_era5_variables(
        variables: Union[str, Sequence[str]]) -> List[str]:
    """Normalize the ``variables`` argument to a deduped list of short keys.

    Accepts a single key (``"wind"``) or a sequence (``["wind", "msl"]``).
    Raises :class:`ValueError` on unknown or empty input.
    """
    if isinstance(variables, str):
        variables = [variables]
    keys = [str(v).strip().lower() for v in variables]
    if not keys:
        raise ValueError("fetch_era5 needs at least one variable; "
                         f"choose from {sorted(ERA5_VARIABLES)}")
    unknown = [k for k in keys if k not in ERA5_VARIABLES]
    if unknown:
        raise ValueError(
            f"unknown ERA5 variable keys {unknown}; "
            f"choose from {sorted(ERA5_VARIABLES)}")
    seen: List[str] = []
    for key in keys:
        if key not in seen:
            seen.append(key)
    return seen


def _validate_era5_dates(d0: _dt.datetime, d1: _dt.datetime) -> None:
    if d0 > d1:
        raise ValueError(f"start {d0.date()} is after end {d1.date()}")
    if d0.date() < ERA5_START:
        raise ValueError(
            f"ERA5 starts {ERA5_START.isoformat()}; "
            f"start {d0.date().isoformat()} is before the record")
    # ERA5 final release lags ~5 days behind real time; the CDS rejects
    # requests past the available data, so refuse future dates outright.
    if d1.date() > _dt.date.today():
        raise ValueError(
            f"end {d1.date().isoformat()} is in the future; "
            "ERA5 is a reanalysis with a ~5-day release lag")


# ---------------------------------------------------------------------------
# time sampling: stride_hours -> per-day (date, [hours]) plan
# ---------------------------------------------------------------------------


def era5_sample_plan(d0: _dt.datetime, d1: _dt.datetime,
                     stride_hours: int) -> List[Tuple[_dt.date, List[str]]]:
    """Sample ``[d0, d1]`` every ``stride_hours`` hours as ``(date, [hours])``.

    Sub-daily strides sample every day with hours ``00:00``, ``06:00``,
    ... for ``stride_hours=6``. Daily-or-longer strides sample one day in
    ``stride_hours // 24`` at ``12:00`` UTC.
    """
    if stride_hours < 1:
        raise ValueError(f"stride_hours must be >= 1, got {stride_hours}")
    plan: List[Tuple[_dt.date, List[str]]] = []
    if stride_hours >= 24:
        step_days = max(1, stride_hours // 24)
        cursor = d0.date()
        while cursor <= d1.date():
            plan.append((cursor, ["12:00"]))
            cursor += _dt.timedelta(days=step_days)
    else:
        hours = [f"{h:02d}:00" for h in range(0, 24, stride_hours)]
        cursor = d0.date()
        while cursor <= d1.date():
            plan.append((cursor, hours))
            cursor += _dt.timedelta(days=1)
    return plan


def era5_month_chunks(
        plan: Sequence[Tuple[_dt.date, List[str]]]
) -> List[Tuple[int, int, List[str], List[str]]]:
    """Group a sample plan into ``(year, month, [days], [hours])`` chunks.

    One chunk per calendar month (CDS request-size hygiene); ``days`` are
    zero-padded ``"DD"`` strings, ``hours`` the sorted unique ``"HH:00"``
    values across the month.
    """
    buckets: Dict[Tuple[int, int], List[Tuple[_dt.date, List[str]]]] = {}
    for day, hours in plan:
        buckets.setdefault((day.year, day.month), []).append((day, hours))
    chunks = []
    for (year, month) in sorted(buckets):
        entries = buckets[(year, month)]
        days = [f"{d.day:02d}" for d, _ in entries]
        hours = sorted({h for _, hs in entries for h in hs})
        chunks.append((year, month, days, hours))
    return chunks


# ---------------------------------------------------------------------------
# CDS area windows (antimeridian handling mirrors oisst_lon_windows)
# ---------------------------------------------------------------------------


def era5_area_windows(bbox: Sequence[float]) -> List[Tuple[float, float, float, float]]:
    """Convert a -180..180 bbox to CDS ``area=[N, W, S, E]`` request windows.

    * Antimeridian-crossing bboxes (``max_lon < min_lon``) become two
      windows: ``[N, min_lon, S, 180]`` and ``[N, -180, S, max_lon]``.
    * Full-globe bboxes (span >= 359.9°) request the whole grid.
    * Otherwise a single window.
    """
    minx, miny, maxx, maxy = validate_era5_bbox(bbox)
    span = maxx - minx
    if span < 0:
        span += 360.0
    if span >= 359.9:
        return [(90.0, -180.0, -90.0, 180.0)]
    if maxx < minx:  # antimeridian crossing
        return [(maxy, minx, miny, 180.0), (maxy, -180.0, miny, maxx)]
    return [(maxy, minx, miny, maxx)]


def era5_cds_names(variables: Sequence[str]) -> List[str]:
    """CDS long variable names for short keys (deduped, order-stable)."""
    names: List[str] = []
    for key in normalize_era5_variables(variables):
        for name in ERA5_VARIABLES[key]["cds"]:
            if name not in names:
                names.append(name)
    return names


def era5_request(variables: Sequence[str], year: int, month: int,
                 days: Sequence[str], hours: Sequence[str],
                 area: Sequence[float]) -> Dict[str, Any]:
    """Build one CDS API request dict for ``reanalysis-era5-single-levels``.

    ``area`` is ``[N, W, S, E]`` in degrees; the grid is pinned to the
    native 0.25° and the format to unarchived NetCDF.
    """
    n, w, s, e = (float(x) for x in area)
    return {
        "product_type": "reanalysis",
        "variable": era5_cds_names(variables),
        "year": f"{year:04d}",
        "month": f"{month:02d}",
        "day": [str(d) for d in days],
        "time": [str(h) for h in hours],
        "area": [n, w, s, e],
        "grid": "0.25/0.25",
        "data_format": "netcdf",
        "download_format": "unarchived",
    }


# ---------------------------------------------------------------------------
# cdsapi access (lazy import; the module imports cleanly without it)
# ---------------------------------------------------------------------------


def _require_cdsapi():
    """Lazy import of cdsapi with an actionable error."""
    try:
        import cdsapi  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "fetch_era5 needs the cdsapi package (pip install cdsapi); "
            "the survey-currents engine itself has no hard dependency on it."
        ) from exc
    return cdsapi


def _cds_client():
    """Build a ``cdsapi.Client``; raise :class:`CredentialsMissing`."""
    cdsapi = _require_cdsapi()
    try:
        return cdsapi.Client()
    except Exception as exc:
        raise CredentialsMissing(
            "Could not create the CDS API client "
            f"({type(exc).__name__}: {exc}). "
            "ERA5 downloads need a free Copernicus CDS account: register at "
            "https://cds.climate.copernicus.eu/, accept the ERA5 Terms of "
            "Use, and create ~/.cdsapirc (url + key) or set CDSAPI_URL / "
            "CDSAPI_KEY. See https://cds.climate.copernicus.eu/how-to-api."
        ) from exc


def _cds_retrieve_bytes(client: Any, request: Dict[str, Any]) -> bytes:
    """Run one CDS retrieve request and return the NetCDF payload bytes."""
    fd, path = tempfile.mkstemp(suffix=".nc")
    os.close(fd)
    try:
        try:
            client.retrieve(ERA5_DATASET, request).download(path)
        except Exception as exc:
            raise RuntimeError(
                f"CDS retrieve failed ({type(exc).__name__}: {exc}); "
                "check the CDS credentials, the accepted Terms of Use, and "
                "that the requested dates are within the ERA5 record."
            ) from exc
        with open(path, "rb") as fh:
            return fh.read()
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# NetCDF parsing (native ERA5 layout -> canonical field)
# ---------------------------------------------------------------------------


def _parse_era5_bytes(payload: bytes, grid_names: Sequence[str]
                      ) -> Tuple[np.ndarray, np.ndarray, List[str],
                                 Dict[str, np.ma.MaskedArray]]:
    """Parse one CDS ERA5 NetCDF payload.

    Returns ``(lats, lons, times, grids)``. The CDS latitude axis is
    descending (N->S); it is flipped to increasing here. Longitudes are
    normalized to -180..180 and sorted increasing (the CDS serves 0..360
    even for -180..180 area requests). Times come from ``valid_time`` (new
    CDS) or ``time`` via netCDF4's ``num2date``.
    """
    from .glsea import _require_netcdf4
    nc = _require_netcdf4()
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
        tmp.write(payload)
        tmp_path = tmp.name
    try:
        ds = nc.Dataset(tmp_path, mode="r")
        try:
            lat = np.asarray(ds.variables["latitude"][:], dtype=float)
            lon360 = np.asarray(ds.variables["longitude"][:], dtype=float)
            time_var = ds.variables.get("valid_time",
                                        ds.variables.get("time"))
            if time_var is None:
                raise ValueError(
                    "ERA5 payload has neither 'valid_time' nor 'time'")
            try:
                dts = nc.num2date(time_var[:], time_var.units)
            except Exception as exc:
                raise ValueError(
                    f"could not decode ERA5 time axis: {exc}") from exc
            grids: Dict[str, np.ma.MaskedArray] = {}
            for name in grid_names:
                if name not in ds.variables:
                    raise ValueError(
                        f"ERA5 payload missing variable {name!r}")
                grids[name] = np.ma.asarray(
                    ds.variables[name][:], dtype=float)
        finally:
            ds.close()
    finally:
        os.unlink(tmp_path)
    times = [_coerce_datetime_utc(d).isoformat() for d in dts]
    # Latitude: descending N->S on the wire -> flip to increasing.
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        lat_asc = True
    else:
        lat_asc = False
    lon = ((lon360 + 180.0) % 360.0) - 180.0
    lon_order = np.argsort(lon)
    lon = lon[lon_order]
    out: Dict[str, np.ma.MaskedArray] = {}
    for name, grid in grids.items():
        if grid.ndim != 3:
            raise ValueError(
                f"ERA5 variable {name!r} has shape {grid.shape}; "
                "expected (time, latitude, longitude)")
        g = grid[:, ::-1, :] if lat_asc else grid
        out[name] = np.ma.masked_invalid(g[:, :, lon_order])
    return lat, lon, times, out


def _concat_lon_parts(
        parts: List[Tuple[np.ndarray, np.ndarray, List[str],
                          Dict[str, np.ma.MaskedArray]]]
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, np.ma.MaskedArray]]:
    """Concatenate parsed ERA5 payloads along longitude, sorted increasing."""
    lats = parts[0][0]
    times = parts[0][2]
    names = list(parts[0][3])
    for part in parts[1:]:
        if part[2] != times:
            raise RuntimeError(
                "ERA5 lon-window parts returned different time axes; "
                "cannot concatenate")
    lons = np.concatenate([p[1] for p in parts])
    order = np.argsort(lons, kind="stable")
    lons_sorted = lons[order]
    # Antimeridian splits can share the seam meridian (180 == -180);
    # drop exact-duplicate longitudes so the grid stays strictly increasing.
    _, unique_idx = np.unique(lons_sorted, return_index=True)
    keep = order[np.sort(unique_idx)]
    merged = {name: np.ma.concatenate([p[3][name] for p in parts], axis=2)[:, :, keep]
              for name in names}
    return lats, lons[keep], times, merged


# ---------------------------------------------------------------------------
# public fetch
# ---------------------------------------------------------------------------


def fetch_era5(variables: Union[str, Sequence[str]],
               bbox: Sequence[float],
               start: DateLike, end: DateLike,
               stride_hours: int = 6) -> Era5Field:
    """Fetch Copernicus ERA5 reanalysis for ``bbox`` over [start, end].

    Args:
        variables: short key(s) — ``"wind"``, ``"msl"``, ``"t2m"``,
            ``"tp"`` (see :data:`ERA5_VARIABLES`). ``"wind"`` fetches the
            u/v component pair.
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180 degrees
            (ascending; antimeridian-crossing boxes split into two CDS
            requests internally).
        start, end: dates/datetimes/ISO strings within 1940-01-01..today.
        stride_hours: time-axis sampling stride, >= 1. ``6`` (default)
            gives 00/06/12/18 UTC; ``24`` gives daily 12:00 UTC;
            ``168`` gives weekly, etc.

    Returns:
        :class:`Era5Field` with converted units (m/s, hPa, °C, mm) and
        provenance (CDS dataset, the exact request dicts, SHA-256
        digests, byte counts, retrieval time).

    Raises:
        ImportError: ``cdsapi`` is not installed.
        CredentialsMissing: no CDS credentials / client creation failed.
        ValueError: invalid variables / bbox / dates / stride.
        RuntimeError: CDS request failures or NetCDF parse failures.
    """
    keys = normalize_era5_variables(variables)
    minx, miny, maxx, maxy = validate_era5_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_era5_dates(d0, d1)
    if stride_hours < 1:
        raise ValueError(f"stride_hours must be >= 1, got {stride_hours}")

    grid_names: List[str] = []
    for key in keys:
        for grid in ERA5_VARIABLES[key]["grids"]:
            if grid not in grid_names:
                grid_names.append(grid)

    client = _cds_client()  # raises ImportError / CredentialsMissing early
    plan = era5_sample_plan(d0, d1, stride_hours)
    chunks = era5_month_chunks(plan)
    windows = era5_area_windows((minx, miny, maxx, maxy))

    time_parts: List[Tuple[np.ndarray, np.ndarray, List[str],
                           Dict[str, np.ma.MaskedArray]]] = []
    requests: List[Dict[str, Any]] = []
    n_payload_bytes = 0
    combined = hashlib.sha256()
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()
    for year, month, days, hours in chunks:
        lon_parts = []
        for area in windows:
            request = era5_request(keys, year, month, days, hours, area)
            requests.append(request)
            payload = _cds_retrieve_bytes(client, request)
            n_payload_bytes += len(payload)
            combined.update(hashlib.sha256(payload).digest())
            try:
                lon_parts.append(_parse_era5_bytes(payload, grid_names))
            except ValueError:
                raise
            except Exception as exc:
                raise RuntimeError(
                    "ERA5 payload parse failed "
                    f"({type(exc).__name__}: {exc}). netCDF4 is required "
                    "(pip install netCDF4)."
                ) from exc
        time_parts.append(_concat_lon_parts(lon_parts))

    lats = time_parts[0][0]
    lons = time_parts[0][1]
    times = [t for part in time_parts for t in part[2]]
    merged = {name: np.ma.concatenate([p[3][name] for p in time_parts], axis=0)
              for name in grid_names}

    # Unit conversions (value * scale + offset per ERA5_VARIABLES).
    converted: Dict[str, np.ma.MaskedArray] = {}
    for key in keys:
        spec = ERA5_VARIABLES[key]
        for grid in spec["grids"]:
            converted[grid] = merged[grid] * spec["scale"] + spec["offset"]

    return Era5Field(
        grids=converted, times=times, lats=lats, lons=lons,
        variables=keys, base_variable=keys[0],
        source=f"copernicus-cds/{ERA5_DATASET}",
        provenance={
            "dataset": "Copernicus ERA5 hourly data on single levels "
                       "(reanalysis)",
            "cds_dataset": ERA5_DATASET,
            "requests": requests,
            "n_requests": len(requests),
            "sha256": combined.hexdigest(),
            "n_bytes": n_payload_bytes,
            "retrieved_at": retrieved_at,
            "bbox_requested": [minx, miny, maxx, maxy],
            "time_requested": [d0.date().isoformat(), d1.date().isoformat()],
            "variables_requested": keys,
            "stride_hours": stride_hours,
            "time_chunks": len(chunks),
            "lon_windows": [list(w) for w in windows],
            "units": {g: ERA5_VARIABLES[k]["units"]
                      for k in keys for g in ERA5_VARIABLES[k]["grids"]},
            "unit_notes": {
                "msl": "Pa on the wire, converted to hPa",
                "t2m": "Kelvin on the wire, converted to degC",
                "tp": "m on the wire, converted to mm per hourly step",
            },
            "grid": f"{ERA5_RES}-degree global (hourly single levels)",
        },
    )


# ---------------------------------------------------------------------------
# offline demo
# ---------------------------------------------------------------------------


def main_demo() -> None:
    """Print a small offline summary (no network). Mirrors sst_global.main_demo."""
    f = Era5Field.synthetic(variables=("wind", "msl"), nt=3)
    print(f"[era5] synthetic grids { {k: v.shape for k, v in f.grids.items()} }, "
          f"t0={f.times[0]}, mean wind={f.spatial_mean('wind', 0):.2f} m/s, "
          f"mean msl={f.spatial_mean('msl', 0):.1f} hPa")


__all__ = [
    "CredentialsMissing",
    "Era5Field",
    "DateLike",
    "ERA5_DATASET",
    "ERA5_START",
    "ERA5_RES",
    "ERA5_VARIABLES",
    "era5_area_windows",
    "era5_cds_names",
    "era5_month_chunks",
    "era5_request",
    "era5_sample_plan",
    "fetch_era5",
    "normalize_era5_variables",
    "validate_era5_bbox",
    "main_demo",
]

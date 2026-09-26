"""NOAA GLSEA (Great Lakes Surface Environmental Analysis) acquisition.

Two ERDDAP endpoints on NOAA's Great Lakes Environmental Research
Laboratory server (verified live 2026-09-26; see docs/DATA_SOURCES.md):

* ``GLSEA_ACSPO_GCS`` griddap — daily sea-surface-temperature analysis
  on a 0.014-deg (~1.5 km) grid, 2006-present, variable ``sst`` in
  degrees Celsius. The longitude axis is clipped to the lakes region:
  its minimum is exactly **-92.4199507342304** — any requested bbox west
  of that floor raises a clear ``ValueError`` before any download.
* ``glsea_avgtemps_3`` tabledap — lake-wide average surface
  temperature per day-of-year (columns ``Year, Day, Sup, Mich, Huron,
  Erie, Ont``), parsed with stdlib ``csv``.

Downloads use stdlib ``urllib`` only. ``netCDF4`` is a lazy import:
the engine imports fine without it and raises an informative error
naming ``pip install netCDF4`` only when a NetCDF actually has to be
parsed.
"""

from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import io
import urllib.request
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np


# ---------------------------------------------------------------------------
# Dataset constants (verified live against ERDDAP 2026-09-26)
# ---------------------------------------------------------------------------

ERDDAP_BASE = "https://apps.glerl.noaa.gov/erddap"
GRID_DATASET = "GLSEA_ACSPO_GCS"
AVG_DATASET = "glsea_avgtemps_3"

#: Longitude axis minimum of GLSEA_ACSPO_GCS — the grid is clipped to the
#: lakes region, so requests west of this floor are rejected up front.
GLSEA_LON_MIN = -92.4199507342304
GLSEA_LON_MAX = -75.8816402880531
GLSEA_LAT_MIN = 38.8749871947297
GLSEA_LAT_MAX = 50.6059751976539
#: Daily timesteps are stamped at 12:00 UTC.
GLSEA_TIME_OF_DAY = _dt.time(12, 0, 0)
#: Analysis coverage per the dataset title.
GLSEA_START_YEAR = 2006

#: Canonical lake name -> glsea_avgtemps_3 CSV column.
LAKE_COLUMNS: Dict[str, str] = {
    "superior": "Sup",
    "michigan": "Mich",
    "huron": "Huron",
    "erie": "Erie",
    "ontario": "Ont",
}

DateLike = Union[_dt.date, _dt.datetime, str]


# ---------------------------------------------------------------------------
# Date / time coercion
# ---------------------------------------------------------------------------


def _coerce_date(value: DateLike) -> _dt.date:
    """Accept date/datetime/ISO-string; return a :class:`date`."""
    if isinstance(value, _dt.datetime):
        return value.date()
    if isinstance(value, _dt.date):
        return value
    text = str(value).strip()
    if text.endswith(("Z", "z")):
        text = text[:-1] + "+00:00"
    try:
        return _dt.datetime.fromisoformat(text).date()
    except ValueError:
        pass
    try:
        return _dt.date.fromisoformat(text[:10])
    except ValueError as exc:
        raise ValueError(
            f"cannot parse {value!r} as a date (want YYYY-MM-DD or ISO datetime)"
        ) from exc


def _coerce_datetime_utc(value: DateLike) -> _dt.datetime:
    """Accept date/datetime/ISO-string; return a UTC-aware datetime.

    Date-only inputs default to 12:00 UTC, the GLSEA daily timestamp.
    """
    if isinstance(value, _dt.date) and not isinstance(value, _dt.datetime):
        return _dt.datetime.combine(value, GLSEA_TIME_OF_DAY,
                                    tzinfo=_dt.timezone.utc)
    if isinstance(value, _dt.datetime):
        dt = value if value.tzinfo else value.replace(tzinfo=_dt.timezone.utc)
        return dt.astimezone(_dt.timezone.utc)
    text = str(value).strip()
    date_only = "T" not in text.upper() and ":" not in text
    if text.endswith(("Z", "z")):
        text = text[:-1] + "+00:00"
    try:
        dt = _dt.datetime.fromisoformat(text)
    except ValueError as exc:
        raise ValueError(
            f"cannot parse {value!r} as a datetime (want YYYY-MM-DD or ISO datetime)"
        ) from exc
    if date_only:  # "2016-01-01" -> noon, the GLSEA daily timestamp
        dt = _dt.datetime.combine(dt.date(), GLSEA_TIME_OF_DAY)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=_dt.timezone.utc)
    return dt.astimezone(_dt.timezone.utc)


def _erddap_iso(value: DateLike) -> str:
    """``2016-01-01T12:00:00Z``-style ISO for ERDDAP constraints."""
    return _coerce_datetime_utc(value).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# Gridded SST: URL construction + bbox validation
# ---------------------------------------------------------------------------


def validate_glsea_bbox(bbox: Sequence[float]) -> Tuple[float, float, float, float]:
    """Check ``(lon_min, lat_min, lon_max, lat_max)`` against the GLSEA grid.

    Constraints must be ascending, and the box must lie inside the
    dataset's clipped lakes-region grid (notably the longitude floor
    ``-92.4199507342304``). Raises ``ValueError`` with an actionable
    message otherwise.
    """
    parts = [float(x) for x in bbox]
    if len(parts) != 4:
        raise ValueError(
            f"bbox needs (lon_min, lat_min, lon_max, lat_max); got {list(bbox)!r}")
    lon_min, lat_min, lon_max, lat_max = parts
    if not lon_min < lon_max:
        raise ValueError(
            f"longitude constraint must be ascending: lon_min={lon_min} "
            f"is not less than lon_max={lon_max}")
    if not lat_min < lat_max:
        raise ValueError(
            f"latitude constraint must be ascending: lat_min={lat_min} "
            f"is not less than lat_max={lat_max}")
    if lon_min < GLSEA_LON_MIN:
        raise ValueError(
            f"lon_min={lon_min} is west of the GLSEA grid's longitude floor "
            f"{GLSEA_LON_MIN} — the GLSEA_ACSPO_GCS grid is clipped to the "
            "Great Lakes region and has no data west of that line")
    if lon_max > GLSEA_LON_MAX:
        raise ValueError(
            f"lon_max={lon_max} is east of the GLSEA grid's eastern edge "
            f"{GLSEA_LON_MAX}")
    if lat_min < GLSEA_LAT_MIN or lat_max > GLSEA_LAT_MAX:
        raise ValueError(
            f"latitude [{lat_min}, {lat_max}] is outside the GLSEA grid's "
            f"latitude range [{GLSEA_LAT_MIN}, {GLSEA_LAT_MAX}]")
    return lon_min, lat_min, lon_max, lat_max


def glsea_sst_url(bbox: Sequence[float], start: DateLike, end: DateLike,
                  stride_days: int = 30) -> str:
    """ERDDAP griddap URL for the daily SST subset.

    Time constraint uses daily-index stride ``stride_days``; lat/lon use
    stride 1 (full resolution).
    """
    lon_min, lat_min, lon_max, lat_max = validate_glsea_bbox(bbox)
    t0, t1 = _coerce_datetime_utc(start), _coerce_datetime_utc(end)
    if t0 > t1:
        raise ValueError(f"start {t0.isoformat()} is after end {t1.isoformat()}")
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1; got {stride_days}")
    query = (
        f"sst[({_erddap_iso(t0)}):{stride_days}:({_erddap_iso(t1)})]"
        f"[({lat_min}):1:({lat_max})][({lon_min}):1:({lon_max})]"
    )
    return f"{ERDDAP_BASE}/griddap/{GRID_DATASET}.nc?{query}"


def _require_netcdf4():
    try:
        import netCDF4  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "parsing GLSEA NetCDF needs netCDF4 (pip install netCDF4)"
        ) from exc
    return netCDF4


def _download_bytes(url: str, timeout: int = 300) -> bytes:
    """Download ``url`` with stdlib urllib; retry transient connection drops."""
    import http.client
    import time
    import urllib.error
    last: Optional[Exception] = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, http.client.RemoteDisconnected,
                TimeoutError, ConnectionError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"download failed after 3 attempts: {url}") from last


# ---------------------------------------------------------------------------
# GlseaField — canonical SST field model
# ---------------------------------------------------------------------------


@dataclass
class GlseaField:
    """Time-indexed sea-surface-temperature grids from GLSEA.

    ``sst`` is ``(nt, ny, nx)`` in degrees Celsius as a numpy masked
    array (land / missing cells masked). ``lats``/``lons`` are 1-D
    coordinate vectors (increasing). Follows the :class:`CurrentField`
    conventions in ``models.py``.
    """

    sst: np.ma.MaskedArray               # (nt, ny, nx) sea-surface temperature, degC
    times: List[str]                     # ISO-8601 timestamps, one per step
    lats: np.ndarray                     # (ny,) degrees north, increasing
    lons: np.ndarray                     # (nx,) degrees east, increasing
    crs: str = "EPSG:4326"
    source: str = "noaa-glsea/GLSEA_ACSPO_GCS"
    provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        self.sst = np.ma.asarray(self.sst, dtype=float)
        if self.sst.shape != (nt, self.lats.shape[0], self.lons.shape[0]):
            raise ValueError(
                f"sst shape {self.sst.shape} does not match "
                f"(nt, ny, nx)=({nt}, {self.lats.shape[0]}, {self.lons.shape[0]})"
            )
        self.lats = np.asarray(self.lats, dtype=float)
        self.lons = np.asarray(self.lons, dtype=float)
        self.bounds = self._bounds()

    def _bounds(self) -> Tuple[float, float, float, float]:
        return (float(self.lons[0]), float(self.lats[0]),
                float(self.lons[-1]), float(self.lats[-1]))

    # -- derived quantities --------------------------------------------------

    def spatial_mean(self, index: int) -> float:
        """Masked/NaN-aware mean SST in degC at timestep ``index``."""
        arr = np.ma.masked_invalid(self.sst[index])
        return float(arr.mean())

    # -- selection -----------------------------------------------------------

    def select_time(self, index: int) -> "GlseaField":
        """Return the single-timestep field at ``index``."""
        return GlseaField(
            sst=self.sst[index:index + 1], times=[self.times[index]],
            lats=self.lats, lons=self.lons, crs=self.crs,
            source=self.source, provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "GlseaField":
        """Spatial subset to ``(min_lon, min_lat, max_lon, max_lat)``."""
        minx, miny, maxx, maxy = (float(x) for x in bbox)
        ix = np.where((self.lons >= minx) & (self.lons <= maxx))[0]
        iy = np.where((self.lats >= miny) & (self.lats <= maxy))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(
                f"bbox {list(bbox)} does not overlap field bounds {self.bounds}")
        return GlseaField(
            sst=self.sst[:, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1],
            times=list(self.times),
            lats=self.lats[iy[0]:iy[-1] + 1], lons=self.lons[ix[0]:ix[-1] + 1],
            crs=self.crs, source=self.source,
            provenance=dict(self.provenance),
        )

    # -- serialization -------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dict (masked cells become NaN)."""
        return {
            "sst": np.ma.filled(self.sst, np.nan).tolist(),
            "times": list(self.times),
            "lats": self.lats.tolist(),
            "lons": self.lons.tolist(),
            "crs": self.crs,
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GlseaField":
        """Rebuild from :meth:`to_dict` (NaN -> masked). Raises on missing keys."""
        required = ("sst", "times", "lats", "lons")
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(f"GlseaField dict missing keys: {missing}")
        return cls(
            sst=np.ma.masked_invalid(np.asarray(data["sst"], dtype=float)),
            times=list(data["times"]),
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            crs=data.get("crs", "EPSG:4326"),
            source=data.get("source", "noaa-glsea/GLSEA_ACSPO_GCS"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "GlseaField":
        import json
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture (tests, demos, offline use) -----------------------

    @classmethod
    def synthetic(
        cls,
        nt: int = 4, ny: int = 6, nx: int = 8,
        lats: Sequence[float] = (41.5, 46.5),
        lons: Sequence[float] = (-92.4, -84.5),
        start: str = "2026-09-16T12:00:00+00:00",
        step_days: int = 1,
        seed: int = 7,
        source: str = "synthetic",
    ) -> "GlseaField":
        """Deterministic synthetic SST field: warm core + cooling trend.

        Stdlib+numpy only. Cells outside the central ellipse are masked
        (irregular-shoreline analog). Used by the test suite and the
        offline demo.
        """
        rng = np.random.default_rng(seed)
        lat = np.linspace(lats[0], lats[1], ny)
        lon = np.linspace(lons[0], lons[1], nx)
        yy, xx = np.meshgrid(np.linspace(-1, 1, ny), np.linspace(-1, 1, nx),
                             indexing="ij")
        base = 18.0 - 4.0 * (xx ** 2 + yy ** 2)   # warm core, degC
        mask = (xx ** 2 + yy ** 2) > 1.0          # masked "land" corners
        sst = np.ma.empty((nt, ny, nx))
        t0 = _dt.datetime.fromisoformat(start)
        times = []
        for k in range(nt):
            sst[k] = np.ma.array(
                base + rng.normal(0.0, 0.15, size=(ny, nx)) - 0.3 * k,
                mask=mask)
            times.append((t0 + _dt.timedelta(days=k * step_days)).isoformat())
        return cls(sst=sst, times=times, lats=lat, lons=lon,
                   source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# NetCDF parsing (netCDF4, lazy import)
# ---------------------------------------------------------------------------


def _parse_glsea_dataset(ds) -> Tuple[np.ndarray, np.ndarray, List[str],
                                      np.ma.MaskedArray]:
    """Read (lats, lons, times, sst) from an open netCDF4 Dataset.

    Masks the ``_FillValue`` sentinel (-99999.0) as missing.
    """
    var = ds.variables["sst"]
    fill = var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else None
    raw = np.asarray(var[:], dtype=float)
    if raw.ndim != 3:
        raise ValueError(f"expected sst(time, lat, lon); got shape {raw.shape}")
    lats = np.asarray(ds.variables["latitude"][:], dtype=float).ravel()
    lons = np.asarray(ds.variables["longitude"][:], dtype=float).ravel()

    time_var = ds.variables["time"]
    units = time_var.getncattr("units") if "units" in time_var.ncattrs() else None
    calendar = (time_var.getncattr("calendar")
                if "calendar" in time_var.ncattrs() else "standard")
    times: List[str] = []
    if units:
        netCDF4 = _require_netcdf4()
        for val in np.atleast_1d(np.asarray(time_var[:])).ravel():
            d = netCDF4.num2date(val, units, calendar=calendar)
            times.append(d.isoformat())
    else:
        times = [str(v) for v in np.atleast_1d(np.asarray(time_var[:])).ravel()]

    if fill is not None:
        sst = np.ma.masked_values(raw, fill)
    else:
        sst = np.ma.asarray(raw)
    return lats, lons, times, sst


def parse_glsea_bytes(payload: bytes) -> Tuple[np.ndarray, np.ndarray, List[str],
                                              np.ma.MaskedArray]:
    """Parse a downloaded GLSEA_ACSPO_GCS NetCDF (in-memory)."""
    netCDF4 = _require_netcdf4()
    ds = netCDF4.Dataset("glsea-in-memory", memory=payload)
    try:
        return _parse_glsea_dataset(ds)
    finally:
        ds.close()


def fetch_glsea_sst(bbox: Sequence[float], start: DateLike, end: DateLike,
                    stride_days: int = 30) -> GlseaField:
    """End-to-end: download a GLSEA daily SST subset -> :class:`GlseaField`.

    ``bbox`` is ``(lon_min, lat_min, lon_max, lat_max)`` and is validated
    against the dataset's clipped lakes-region grid (see
    :func:`validate_glsea_bbox`). ``start``/``end`` accept
    date/datetime/ISO strings; the time axis is sampled every
    ``stride_days`` days.
    """
    url = glsea_sst_url(bbox, start, end, stride_days=stride_days)
    payload = _download_bytes(url)
    lats, lons, times, sst = parse_glsea_bytes(payload)
    retrieved = _dt.datetime.now(_dt.timezone.utc).isoformat()
    provenance = {
        "url": url,
        "dataset": GRID_DATASET,
        "bbox": [float(x) for x in bbox],
        "time_window": [_erddap_iso(start), _erddap_iso(end)],
        "stride_days": stride_days,
        "bytes_downloaded": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "retrieved_at": retrieved,
        "units": "degree_C",
    }
    return GlseaField(sst=sst, times=times, lats=lats, lons=lons,
                      provenance=provenance)


# ---------------------------------------------------------------------------
# Lake-average series (tabledap CSV, stdlib only)
# ---------------------------------------------------------------------------


@dataclass
class LakeSeries:
    """Lake-wide average surface temperature series from GLSEA.

    ``dates`` are ISO ``YYYY-MM-DD`` strings, ``temps`` are degrees
    Celsius; both lists are aligned and filtered to the requested window.
    """

    lake: str                          # canonical lowercase name, e.g. "superior"
    dates: List[str]                   # ISO dates, ascending
    temps: List[float]                 # degC, one per date
    provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.dates) != len(self.temps):
            raise ValueError(
                f"dates ({len(self.dates)}) and temps ({len(self.temps)}) "
                "lengths differ")
        if self.lake not in LAKE_COLUMNS:
            raise ValueError(
                f"unknown lake {self.lake!r}; known lakes: "
                f"{sorted(LAKE_COLUMNS)}")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "lake": self.lake,
            "lake_column": LAKE_COLUMNS[self.lake],
            "dates": list(self.dates),
            "temps": list(self.temps),
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "LakeSeries":
        return cls(
            lake=data["lake"],
            dates=list(data["dates"]),
            temps=[float(t) for t in data["temps"]],
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "LakeSeries":
        import json
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    @property
    def n(self) -> int:
        return len(self.dates)

    def mean(self) -> float:
        """Mean temperature over the series in degC."""
        if not self.temps:
            raise ValueError("empty series")
        return float(np.mean(self.temps))


def glsea_averages_url(lake: str, year0: int) -> str:
    """ERDDAP tabledap CSV URL for one lake's average temps since ``year0``."""
    key = lake.strip().lower()
    if key not in LAKE_COLUMNS:
        raise ValueError(
            f"unknown lake {lake!r}; known lakes: {sorted(LAKE_COLUMNS)}")
    column = LAKE_COLUMNS[key]
    # The ">=" is percent-encoded: some proxies/middleboxes drop the
    # connection when a raw ">" appears in the request line (verified
    # live 2026-09-26); ERDDAP decodes %3E identically.
    return (f"{ERDDAP_BASE}/tabledap/{AVG_DATASET}.csv"
            f"?Year,Day,{column}&Year%3E={year0}")


def parse_glsea_averages_csv(text: str, lake: str,
                             start: _dt.date, end: _dt.date) -> "LakeSeries":
    """Parse a ``glsea_avgtemps_3`` CSV payload into a filtered LakeSeries.

    Skips the ERDDAP units row and any rows with missing/non-numeric
    temperatures; keeps rows whose date falls in ``[start, end]``.
    """
    key = lake.strip().lower()
    if key not in LAKE_COLUMNS:
        raise ValueError(
            f"unknown lake {lake!r}; known lakes: {sorted(LAKE_COLUMNS)}")
    dates: List[str] = []
    temps: List[float] = []
    reader = csv.reader(io.StringIO(text))
    for row in reader:
        if len(row) < 3:
            continue
        year_txt, day_txt, temp_txt = (row[0].strip(), row[1].strip(),
                                       row[2].strip())
        if not (year_txt.isdigit() and day_txt.isdigit()):
            continue  # header / units row
        try:
            temp = float(temp_txt)
        except ValueError:
            continue  # missing temperature (NaN / empty)
        year, doy = int(year_txt), int(day_txt)
        try:
            day = _dt.date(year, 1, 1) + _dt.timedelta(days=doy - 1)
        except (ValueError, OverflowError):
            continue
        if day.year != year:  # day-of-year ran past Dec 31
            continue
        if start <= day <= end:
            dates.append(day.isoformat())
            temps.append(temp)
    return LakeSeries(lake=key, dates=dates, temps=temps)


def fetch_glsea_lake_averages(lake: str, start: DateLike,
                             end: DateLike) -> LakeSeries:
    """End-to-end: download one lake's average-temperature series.

    ``lake`` is case-insensitive (``superior``, ``michigan``, ``huron``,
    ``erie``, ``ontario``); anything else raises ``ValueError``.
    ``start``/``end`` accept date/datetime/ISO strings.
    """
    key = lake.strip().lower()
    d0, d1 = _coerce_date(start), _coerce_date(end)
    if d0 > d1:
        raise ValueError(f"start {d0.isoformat()} is after end {d1.isoformat()}")
    url = glsea_averages_url(key, d0.year)
    payload = _download_bytes(url)
    text = payload.decode("utf-8")
    series = parse_glsea_averages_csv(text, key, d0, d1)
    series.provenance = {
        "url": url,
        "dataset": AVG_DATASET,
        "lake_column": LAKE_COLUMNS[key],
        "window": [d0.isoformat(), d1.isoformat()],
        "rows_returned": series.n,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "retrieved_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
    }
    return series

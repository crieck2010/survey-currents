"""CSR GRACE / GRACE-FO RL06.3 terrestrial water storage (TWS) anomalies.

:func:`fetch_grace` downloads the CSR mascon "all-corrections" NetCDF
(keyless anonymous HTTPS, ~107 MB) once into a SHA-256-verified local
cache and returns a :class:`WaterField`: monthly liquid-water-equivalent
thickness anomalies (cm) on a 0.25° lat/lon grid, land-only (the CSR
land mask is applied — ocean cells are NaN), with gap months carried
as all-NaN frames so the monthly timeline stays complete and honest.

Access truth (verified live 2026-09-27 — see docs/DATA_SOURCES.md):

* ``https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/``
  — keyless anonymous HTTPS (HTTP 200, ``Accept-Ranges: bytes``).
* ``CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc``
  (112,611,569 bytes): monthly TWS anomaly grids, 2002-04 – 2026-06,
  258 months, already gridded at 0.25° (no mascon mapping needed).
* ``CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc`` (4,174,189
  bytes): 0/1 land mask on the same grid (``LO_val`` variable).
* Units are **cm** of liquid-water-equivalent thickness; the
  2004.000–2009.999 time-mean has been removed (``time_mean_removed``
  global attribute) — every value is an anomaly against that baseline.
* The archive's ``months_missing`` global attribute lists every
  missing month explicitly, including the 2017-07 … 2018-05
  GRACE/GRACE-FO gap. Missing months are absent from the time axis —
  they are never interpolated; :func:`fetch_grace` inserts them as
  all-NaN frames and records them in ``WaterField.gap_months`` and in
  provenance.

Cache discipline (mirrors :mod:`currents.basemaps` and
:mod:`currents.storms`): both files are downloaded once into
``$SURVEY_CURRENTS_CACHE/grace`` (else ``~/.cache/survey-currents/grace``)
with atomic writes and ``.sha256`` sidecars; a corrupt cache entry is
re-downloaded. CSR re-releases the file monthly, so entries older than
``max_cache_age_days`` (default 30) are revalidated against the remote
``Last-Modified`` header and refreshed when the server copy is newer.

Subsetting: months within ``[start, end]`` are selected (missing
months become all-NaN gap frames); the spatial window is cut to
``bbox``. Longitudes in the archive run 0..360 ascending; the field
normalizes them to the conventional -180..180 axis.
"""

from __future__ import annotations

import calendar as _cal
import datetime as _dt
import hashlib
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from dataclasses import field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date
from .sst_global import validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]

__all__ = [
    "GRACE_PRODUCT",
    "GRACE_PRODUCT_VERSION",
    "GRACE_SOLUTION_FILE",
    "GRACE_SOLUTION_URL",
    "GRACE_MASK_FILE",
    "GRACE_MASK_URL",
    "GRACE_EPOCH",
    "GRACE_RECORD_START",
    "GRACE_BASELINE",
    "GRACE_UNITS",
    "GRACE_MAX_CACHE_AGE_DAYS",
    "WaterField",
    "grace_cache_dir",
    "grace_file_url",
    "fetch_grace",
]

#: Product identity (from the file's global attributes, verified live
#: 2026-09-27).
GRACE_PRODUCT = "CSR GRACE and GRACE-FO MASCON RL0603M"
GRACE_PRODUCT_VERSION = "RL06.3"
#: Keyless CSR download directory (verified live 2026-09-27).
_GRACE_DIR = "https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/"
#: Monthly TWS anomaly grids: 258 months, 2002-04 – 2026-06, 107 MB.
GRACE_SOLUTION_FILE = "CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc"
GRACE_SOLUTION_URL = _GRACE_DIR + GRACE_SOLUTION_FILE
#: 0/1 land mask (``LO_val``) on the same 0.25° grid, ~4 MB.
GRACE_MASK_FILE = "CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc"
GRACE_MASK_URL = _GRACE_DIR + GRACE_MASK_FILE
#: NetCDF ``time`` axis unit: days since this epoch (mid-month values).
GRACE_EPOCH = _dt.date(2002, 1, 1)
#: First month with data (``time_coverage_start`` 2002-04-05).
GRACE_RECORD_START = _dt.date(2002, 4, 1)
#: Anomaly baseline (``time_mean_removed`` global attribute): the
#: 2004–2009 time-mean has been removed from every grid cell.
GRACE_BASELINE = ("anomaly vs 2004–2009 time-mean removed "
                  "(CSR RL06.3 time_mean_removed 2004.000–2009.999)")
GRACE_UNITS = "cm"
#: Default cache freshness: revalidate files older than this many days.
GRACE_MAX_CACHE_AGE_DAYS = 30


def grace_cache_dir() -> str:
    """User cache root for GRACE downloads.

    ``$SURVEY_CURRENTS_CACHE/grace`` when set, else
    ``~/.cache/survey-currents/grace``. Created on demand.
    """
    root = os.environ.get("SURVEY_CURRENTS_CACHE") or os.path.join(
        os.path.expanduser("~"), ".cache", "survey-currents")
    path = os.path.join(root, "grace")
    os.makedirs(path, exist_ok=True)
    return path


def grace_file_url(which: str = "solution") -> str:
    """Keyless CSR download URL for ``"solution"`` or ``"mask"``."""
    if which == "solution":
        return GRACE_SOLUTION_URL
    if which == "mask":
        return GRACE_MASK_URL
    raise ValueError(f"grace_file_url: expected 'solution' or 'mask', got {which!r}")


def _tool_version() -> str:
    from . import __version__
    return __version__


# ---------------------------------------------------------------------------
# Download + cache (same discipline as currents.storms)
# ---------------------------------------------------------------------------

def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _sidecar_digest(path: str) -> Optional[str]:
    sidecar = path + ".sha256"
    if not (os.path.isfile(path) and os.path.isfile(sidecar)):
        return None
    try:
        with open(sidecar, "r", encoding="utf-8") as fh:
            expected = fh.read().strip().split()[0]
    except OSError:
        return None
    if _sha256_file(path) == expected:
        return expected
    return None


def _write_atomic(path: str, data: bytes) -> None:
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, path)


def _http_head_last_modified(url: str, timeout: int = 60) -> Optional[str]:
    req = urllib.request.Request(
        url, method="HEAD",
        headers={"User-Agent": f"survey-currents/{_tool_version()}"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.headers.get("Last-Modified")
    except (urllib.error.URLError, OSError):
        return None


def _download_bytes(url: str, timeout: int = 900) -> bytes:
    import http.client
    import time
    last: Optional[Exception] = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent":
                              f"survey-currents/{_tool_version()}"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, http.client.RemoteDisconnected,
                http.client.IncompleteRead, TimeoutError,
                ConnectionError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GRACE download failed after 3 attempts: {url}") from last


def _ensure_one_file(filename: str, url: str, cache: str,
                     max_cache_age_days: int,
                     refresh: bool) -> Tuple[str, Dict[str, Any]]:
    """Download-once helper shared by the solution and mask files."""
    path = os.path.join(cache, filename)
    provenance: Dict[str, Any] = {
        "file": filename, "file_url": url, "cache_path": path,
        "downloaded": False, "cache_hit": False,
    }

    def _fresh_download() -> str:
        data = _download_bytes(url)
        _write_atomic(path, data)
        digest = _sha256_file(path)
        with open(path + ".sha256", "w", encoding="utf-8") as fh:
            fh.write(digest + "  " + filename + "\n")
        return digest

    digest = None if refresh else _sidecar_digest(path)
    if digest is not None and not refresh:
        try:
            age_days = ((_dt.datetime.now(_dt.timezone.utc)
                         - _dt.datetime.fromtimestamp(
                             os.path.getmtime(path),
                             tz=_dt.timezone.utc)).days)
        except OSError:
            age_days = max_cache_age_days + 1
        if age_days <= max_cache_age_days:
            provenance.update(cache_hit=True, sha256=digest,
                              cache_age_days=age_days)
            return path, provenance
        remote_lm = _http_head_last_modified(url)
        lm_path = path + ".last_modified"
        local_lm = None
        if os.path.isfile(lm_path):
            with open(lm_path, "r", encoding="utf-8") as fh:
                local_lm = fh.read().strip() or None
        if remote_lm and remote_lm == local_lm:
            provenance.update(cache_hit=True, sha256=digest,
                              cache_age_days=age_days,
                              revalidated=True)
            return path, provenance
        digest = _fresh_download()
        if remote_lm:
            with open(lm_path, "w", encoding="utf-8") as fh:
                fh.write(remote_lm)
        provenance.update(downloaded=True, sha256=digest,
                          refreshed=True)
        return path, provenance

    digest = _fresh_download()
    remote_lm = _http_head_last_modified(url)
    if remote_lm:
        with open(path + ".last_modified", "w", encoding="utf-8") as fh:
            fh.write(remote_lm)
    provenance.update(downloaded=True, sha256=digest)
    return path, provenance


def ensure_grace_files(cache_dir: Optional[str] = None,
                       max_cache_age_days: int = GRACE_MAX_CACHE_AGE_DAYS,
                       refresh: bool = False) -> Tuple[str, str, Dict[str, Any]]:
    """Return ``(solution_path, mask_path, download_provenance)``.

    Downloads the CSR solution + land-mask NetCDFs once into the GRACE
    cache; reuses verified copies, refreshes corrupt or stale ones —
    see :func:`_ensure_one_file`. Because CSR re-releases the solution
    monthly, copies older than ``max_cache_age_days`` are revalidated
    against the remote ``Last-Modified`` header.
    """
    cache = cache_dir or grace_cache_dir()
    os.makedirs(cache, exist_ok=True)
    sol_path, sol_prov = _ensure_one_file(
        GRACE_SOLUTION_FILE, GRACE_SOLUTION_URL, cache,
        max_cache_age_days, refresh)
    mask_path, mask_prov = _ensure_one_file(
        GRACE_MASK_FILE, GRACE_MASK_URL, cache,
        max_cache_age_days, refresh)
    provenance = {
        "solution": sol_prov,
        "mask": mask_prov,
    }
    return sol_path, mask_path, provenance


# ---------------------------------------------------------------------------
# NetCDF parsing (pure over duck-typed datasets — offline-testable)
# ---------------------------------------------------------------------------

def _require_netcdf4():
    """Lazily import netCDF4 (optional dependency)."""
    try:
        import netCDF4
    except ImportError as exc:
        raise ImportError(
            "reading GRACE needs the optional 'netCDF4' package "
            "(the CSR mascon archive is NetCDF-4/HDF5).\nInstall it with:\n\n"
            "    pip install netCDF4\n\n"
            "WaterField.synthetic() and the parsing helpers keep working "
            "without it."
        ) from exc
    return netCDF4


def _month_label(time_days: float,
                 bound_lo: float, bound_hi: float) -> Tuple[int, int]:
    """(year, month) containing the midpoint of a month's time bounds.

    ``time_days``/bounds are days since :data:`GRACE_EPOCH`. The bounds
    midpoint is the documented month extent (verified: bounds
    ``[94, 120]`` -> 2002-04).
    """
    mid = GRACE_EPOCH + _dt.timedelta(days=(bound_lo + bound_hi) / 2.0)
    return mid.year, mid.month


def _file_months(ds: Any) -> List[Tuple[int, int]]:
    """(year, month) labels for every frame on the file's time axis."""
    times = np.asarray(ds.variables["time"][:]).astype(float).ravel()
    bounds = np.asarray(ds.variables["time_bounds"][:]).astype(float)
    if bounds.shape != (times.shape[0], 2):
        # Fall back to the mid-month time value itself.
        return [((GRACE_EPOCH + _dt.timedelta(days=float(t))).year,
                 (GRACE_EPOCH + _dt.timedelta(days=float(t))).month)
                for t in times]
    return [_month_label(float(t), float(b[0]), float(b[1]))
            for t, b in zip(times, bounds)]


def _window_months(d0: _dt.date, d1: _dt.date) -> List[Tuple[int, int]]:
    """Every (year, month) in ``[d0, d1]`` (month granularity)."""
    months: List[Tuple[int, int]] = []
    y, m = d0.year, d0.month
    while (y, m) <= (d1.year, d1.month):
        months.append((y, m))
        m += 1
        if m > 12:
            m, y = 1, y + 1
    return months


def _norm_lon_360(lon: float) -> float:
    return float(lon) % 360.0


def _lon_windows_360(bbox: Tuple[float, float, float, float]
                     ) -> List[Tuple[float, float]]:
    """bbox (-180..180) -> one or two [west, east] windows on 0..360."""
    lon_min, lat_min, lon_max, lat_max = bbox
    west, east = _norm_lon_360(lon_min), _norm_lon_360(lon_max)
    if west == east:
        # Full-globe bbox (-180..180): both ends normalize to 180.
        return [(0.0, 360.0)]
    if west <= east:
        return [(west, east)]
    # A bbox straddling the prime meridian in -180..180 terms (e.g.
    # -10..10) wraps on the 0..360 axis: two windows.
    return [(west, 360.0), (0.0, east)]


def _parse_grace_dataset(ds: Any, land_mask: np.ndarray,
                         bbox: Tuple[float, float, float, float],
                         d0: _dt.date, d1: _dt.date) -> Dict[str, Any]:
    """Subset the CSR solution into arrays (pure, offline-testable).

    ``ds`` is duck-typed (``ds.variables[name][:]``) so tests can pass
    fake datasets without netCDF4. ``land_mask`` is a (nlat, nlon)
    0/1 array on the solution grid. Returns a dict with ``times``
    (month-start dates, gap months included), ``lats``, ``lons``
    (-180..180 ascending), ``values`` (ntime, nlat, nlon) float with
    ocean cells and gap months as NaN, ``gap_months`` (month-start
    dates), and ``file_months``.
    """
    import warnings
    with warnings.catch_warnings():
        # The CSR file carries valid_min/valid_max attributes that
        # netCDF4 cannot safely cast for these scaled coordinate
        # variables — cosmetic warnings only; the values are fine.
        warnings.simplefilter("ignore")
        file_lats = np.asarray(ds.variables["lat"][:]).astype(float).ravel()
        file_lons = np.asarray(ds.variables["lon"][:]).astype(float).ravel()
        months = _file_months(ds)
    var = ds.variables["lwe_thickness"]

    lon_min, lat_min, lon_max, lat_max = bbox
    lat_idx = np.nonzero((file_lats >= lat_min) & (file_lats <= lat_max))[0]
    if lat_idx.size == 0:
        raise ValueError(
            f"GRACE: bbox latitude range [{lat_min}, {lat_max}] selects no "
            "grid cells")

    # Longitude selection on the 0..360 axis (possibly two windows).
    lon_sel = np.zeros(file_lons.shape[0], dtype=bool)
    for west, east in _lon_windows_360(bbox):
        lon_sel |= (file_lons >= west) & (file_lons <= east)
    lon_idx = np.nonzero(lon_sel)[0]
    if lon_idx.size == 0:
        raise ValueError(
            f"GRACE: bbox longitude range [{lon_min}, {lon_max}] selects no "
            "grid cells")

    lats = file_lats[lat_idx]
    # Normalize selected longitudes to -180..180 and keep ascending order.
    lons_360 = file_lons[lon_idx]
    lons = np.array([x - 360.0 if x > 180.0 else x for x in lons_360])
    order = np.argsort(lons)
    lons = lons[order]
    lon_idx = lon_idx[order]

    mask = np.asarray(land_mask, dtype=float)
    if mask.shape != (file_lats.shape[0], file_lons.shape[0]):
        raise ValueError(
            f"GRACE: land mask shape {mask.shape} does not match solution "
            f"grid ({file_lats.shape[0]}, {file_lons.shape[0]})")
    mask_sub = mask[np.ix_(lat_idx, lon_idx)]

    month_index = {ym: i for i, ym in enumerate(months)}
    want = _window_months(d0, d1)
    times: List[_dt.date] = []
    gap_months: List[_dt.date] = []
    frames: List[np.ndarray] = []
    for ym in want:
        t0 = _dt.date(ym[0], ym[1], 1)
        times.append(t0)
        fi = month_index.get(ym)
        if fi is None:
            # Gap month: all-NaN frame — never interpolated, never
            # silently skipped (the renderer draws a NO-DATA panel).
            gap_months.append(t0)
            frames.append(np.full((lat_idx.shape[0], lon_idx.shape[0]),
                                  np.nan))
            continue
        raw = np.asanyarray(var[fi])  # asanyarray: keep netCDF4 masks
        if np.ma.isMaskedArray(raw):
            raw = np.ma.filled(raw, np.nan)
        grid = np.asarray(raw, dtype=float)[np.ix_(lat_idx, lon_idx)]
        grid = np.where(mask_sub == 1.0, grid, np.nan)  # land-only
        frames.append(grid)

    values = np.stack(frames) if frames else np.zeros(
        (0, lat_idx.shape[0], lon_idx.shape[0]))
    return {
        "times": times,
        "lats": lats,
        "lons": lons,
        "values": values,
        "gap_months": gap_months,
        "file_months": months,
    }


# ---------------------------------------------------------------------------
# WaterField
# ---------------------------------------------------------------------------

@dataclass
class WaterField:
    """Monthly CSR GRACE/GRACE-FO terrestrial water storage anomalies.

    ``values`` is ``(ntime, nlat, nlon)`` float, cm of
    liquid-water-equivalent thickness anomaly (:data:`GRACE_BASELINE`),
    on a regular lat/lon grid (ascending axes, longitudes -180..180).
    Ocean cells are NaN (CSR land mask applied). Gap months — months
    with no GRACE/GRACE-FO solution, e.g. the 2017-07 … 2018-05
    inter-mission gap — are present in ``times`` as all-NaN frames and
    listed in ``gap_months``: never interpolated, never silently
    dropped.
    """

    times: List[_dt.date]
    lats: np.ndarray
    lons: np.ndarray
    values: np.ndarray
    gap_months: List[_dt.date] = _dc_field(default_factory=list)
    units: str = GRACE_UNITS
    anomaly_baseline: str = GRACE_BASELINE
    bbox: Tuple[float, float, float, float] = (-180.0, -90.0, 180.0, 90.0)
    source: str = "grace"
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        self.times = [_coerce_date(t) for t in self.times]
        self.lats = np.asarray(self.lats, dtype=float).reshape(-1)
        self.lons = np.asarray(self.lons, dtype=float).reshape(-1)
        self.values = np.asarray(self.values, dtype=float).reshape(
            nt, self.lats.shape[0], self.lons.shape[0])
        self.gap_months = [_coerce_date(t) for t in self.gap_months]
        self.bbox = validate_sst_bbox(self.bbox)

    def __len__(self) -> int:
        return len(self.times)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    @property
    def shape(self) -> Tuple[int, int, int]:
        """(ntime, nlat, nlon)."""
        return self.values.shape  # type: ignore[return-value]

    @property
    def time_range(self) -> Tuple[Optional[_dt.date], Optional[_dt.date]]:
        """(earliest, latest) month, or (None, None) when empty."""
        if not self.times:
            return None, None
        return min(self.times), max(self.times)

    @property
    def n_gap_months(self) -> int:
        """Number of gap (all-NaN) months in the field."""
        return len(self.gap_months)

    def is_gap(self, index: int) -> bool:
        """True when frame ``index`` is a gap month (all-NaN)."""
        return bool(np.all(np.isnan(self.values[index])))

    def spatial_mean(self, index: int) -> float:
        """NaN-aware spatial mean of frame ``index`` (cm); NaN for gaps."""
        frame = self.values[index]
        if not np.isfinite(frame).any():
            return float("nan")  # avoids the all-NaN np.nanmean warning
        return float(np.nanmean(frame))

    # -- filters ---------------------------------------------------------

    def select_time(self, start: DateLike, end: DateLike) -> "WaterField":
        """Months with ``start <= month <= end`` (month granularity)."""
        d0, d1 = _coerce_date(start), _coerce_date(end)
        keep = [i for i, t in enumerate(self.times) if d0 <= t <= d1]
        kept_gaps = [t for t in self.gap_months if d0 <= t <= d1]
        return WaterField(
            times=[self.times[i] for i in keep],
            lats=self.lats.copy(), lons=self.lons.copy(),
            values=self.values[keep] if keep else np.zeros(
                (0, self.lats.shape[0], self.lons.shape[0])),
            gap_months=kept_gaps, units=self.units,
            anomaly_baseline=self.anomaly_baseline,
            bbox=self.bbox, source=self.source,
            provenance=dict(self.provenance))

    def select_bbox(self, bbox: Sequence[float]) -> "WaterField":
        """Spatial subset to ``bbox`` (gap months preserved)."""
        box = validate_sst_bbox(bbox)
        lon_min, lat_min, lon_max, lat_max = box
        lat_sel = (self.lats >= lat_min) & (self.lats <= lat_max)
        if lon_min <= lon_max:
            lon_sel = (self.lons >= lon_min) & (self.lons <= lon_max)
        else:  # antimeridian-crossing bbox on the -180..180 axis
            lon_sel = (self.lons >= lon_min) | (self.lons <= lon_max)
        return WaterField(
            times=list(self.times),
            lats=self.lats[lat_sel], lons=self.lons[lon_sel],
            values=self.values[:, lat_sel, :][:, :, lon_sel],
            gap_months=list(self.gap_months), units=self.units,
            anomaly_baseline=self.anomaly_baseline,
            bbox=box, source=self.source,
            provenance=dict(self.provenance))

    # -- serialization ---------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "times": [t.isoformat() for t in self.times],
            "lats": [float(x) for x in self.lats],
            "lons": [float(x) for x in self.lons],
            "values": self.values.tolist(),
            "gap_months": [t.isoformat() for t in self.gap_months],
            "units": self.units,
            "anomaly_baseline": self.anomaly_baseline,
            "bbox": list(self.bbox),
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WaterField":
        return cls(
            times=[_coerce_date(t) for t in data["times"]],
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            values=np.asarray(data["values"], dtype=float),
            gap_months=[_coerce_date(t)
                        for t in data.get("gap_months", [])],
            units=str(data.get("units", GRACE_UNITS)),
            anomaly_baseline=str(data.get("anomaly_baseline",
                                          GRACE_BASELINE)),
            bbox=tuple(data.get("bbox", (-180.0, -90.0, 180.0, 90.0))),
            source=str(data.get("source", "grace")),
            provenance=dict(data.get("provenance", {})))

    def to_json(self, path: str) -> str:
        """Write the field as JSON; return ``path``."""
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "WaterField":
        """Read a field written by :meth:`to_json`."""
        import json
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture -------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-125.0, 30.0, -110.0, 45.0),
                  start: DateLike = "2020-01-01",
                  end: DateLike = "2020-12-01",
                  resolution: float = 2.5,
                  gap: Sequence[str] = ("2020-06-01",),
                  seed: int = 11,
                  source: str = "synthetic") -> "WaterField":
        """Deterministic synthetic TWS anomalies (offline tests / demos).

        Builds monthly anomaly grids with a declining (drying) trend
        plus seasonal cycle and noise; ``gap`` lists month-start ISO
        dates rendered as all-NaN gap frames. ~35% of cells are masked
        as "ocean".
        """
        rng = np.random.default_rng(seed)
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        d0 = _coerce_date(start)
        d1 = _coerce_date(end)
        d0 = _dt.date(d0.year, d0.month, 1)
        d1 = _dt.date(d1.year, d1.month, 1)
        if d1 < d0:
            raise ValueError(f"end {d1} is before start {d0}")
        months = _window_months(d0, d1)
        lons = np.arange(lon_min, lon_max + resolution / 2, resolution)
        lats = np.arange(lat_min, lat_max + resolution / 2, resolution)
        ny, nx = lats.shape[0], lons.shape[0]
        ocean = rng.random((ny, nx)) < 0.35
        gap_set = {_coerce_date(g) for g in gap}
        times: List[_dt.date] = []
        frames: List[np.ndarray] = []
        gap_months: List[_dt.date] = []
        lon2d, lat2d = np.meshgrid(lons, lats)
        for k, (y, m) in enumerate(months):
            t = _dt.date(y, m, 1)
            times.append(t)
            if t in gap_set:
                gap_months.append(t)
                frames.append(np.full((ny, nx), np.nan))
                continue
            seasonal = 6.0 * np.sin(2 * np.pi * (m - 1) / 12.0)
            trend = -0.9 * k  # drying trend (cm)
            field = (trend + seasonal
                     + 3.0 * np.sin(np.radians(lon2d * 2.0))
                     * np.cos(np.radians(lat2d * 3.0))
                     + rng.normal(0.0, 1.5, (ny, nx)))
            field = np.where(ocean, np.nan, field)
            frames.append(field)
        return cls(times=times, lats=lats, lons=lons,
                   values=np.stack(frames),
                   gap_months=gap_months, units=GRACE_UNITS,
                   anomaly_baseline=GRACE_BASELINE,
                   bbox=(lon_min, lat_min, lon_max, lat_max),
                   source=source,
                   provenance={"synthetic": True, "seed": seed,
                               "resolution": resolution})


# ---------------------------------------------------------------------------
# Public fetch
# ---------------------------------------------------------------------------

def fetch_grace(bbox: Sequence[float], start: DateLike, end: DateLike,
                cache_dir: Optional[str] = None,
                max_cache_age_days: int = GRACE_MAX_CACHE_AGE_DAYS,
                refresh: bool = False) -> WaterField:
    """Fetch CSR GRACE/GRACE-FO RL06.3 TWS anomalies for ``bbox`` x ``[start, end]``.

    Downloads the solution + land-mask NetCDFs once (cached,
    SHA-256-verified, freshness-revalidated — see
    :func:`ensure_grace_files`), then returns a :class:`WaterField`
    with monthly anomaly frames (cm LWE, land-only). Months with no
    GRACE/GRACE-FO solution in the window — including the 2017-07 …
    2018-05 inter-mission gap — are present as all-NaN frames and
    listed in ``field.gap_months``: never interpolated, never silently
    dropped.
    """
    box = validate_sst_bbox(bbox)
    d0 = _coerce_date(start)
    d1 = _coerce_date(end)
    if d1 < d0:
        raise ValueError(f"end {d1.isoformat()} is before start {d0.isoformat()}")
    if _dt.date(d0.year, d0.month, 1) < GRACE_RECORD_START:
        raise ValueError(
            f"GRACE record starts {GRACE_RECORD_START.isoformat()}; "
            f"start {d0.isoformat()} is before the record")
    if d1 > _dt.date.today():
        raise ValueError(f"end {d1.isoformat()} is in the future")

    sol_path, mask_path, dl_prov = ensure_grace_files(
        cache_dir=cache_dir, max_cache_age_days=max_cache_age_days,
        refresh=refresh)
    netCDF4 = _require_netcdf4()
    with netCDF4.Dataset(sol_path, "r") as ds, \
            netCDF4.Dataset(mask_path, "r") as mds:
        mask_var = mds.variables.get("LO_val")
        if mask_var is None:
            # Fall back to the first 2D variable (the mask file's
            # layout is a documented constant, but be liberal).
            names = [n for n in mds.variables
                     if len(mds.variables[n].dimensions) == 2]
            if not names:
                raise RuntimeError(
                    "GRACE land-mask file has no 2D variable")
            mask_var = mds.variables[names[0]]
        land_mask = np.asarray(mask_var[:])
        if np.ma.isMaskedArray(land_mask):
            land_mask = np.ma.filled(land_mask, 0.0)
        parsed = _parse_grace_dataset(ds, np.asarray(land_mask, dtype=float),
                                      box, d0, d1)

    provenance = {
        "source": "grace",
        "product": GRACE_PRODUCT,
        "product_version": GRACE_PRODUCT_VERSION,
        "solution_url": dl_prov["solution"]["file_url"],
        "solution_sha256": dl_prov["solution"].get("sha256"),
        "mask_url": dl_prov["mask"]["file_url"],
        "mask_sha256": dl_prov["mask"].get("sha256"),
        "retrieved_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "bbox": list(box),
        "start": d0.isoformat(),
        "end": d1.isoformat(),
        "units": GRACE_UNITS,
        "anomaly_baseline": GRACE_BASELINE,
        "n_months": len(parsed["times"]),
        "gap_months": [t.isoformat() for t in parsed["gap_months"]],
        "n_gap_months": len(parsed["gap_months"]),
        "ocean_masked": True,
        "cache_hit": bool(dl_prov["solution"].get("cache_hit", False)
                          and dl_prov["mask"].get("cache_hit", False)),
        "tool": f"survey-currents {_tool_version()}",
    }
    return WaterField(times=parsed["times"], lats=parsed["lats"],
                      lons=parsed["lons"], values=parsed["values"],
                      gap_months=parsed["gap_months"],
                      units=GRACE_UNITS, anomaly_baseline=GRACE_BASELINE,
                      bbox=box, source="grace", provenance=provenance)


def main_demo() -> None:
    """Offline demo: synthetic field summary (no network)."""
    field = WaterField.synthetic()
    print(f"synthetic WaterField: {len(field)} months, "
          f"{field.n_gap_months} gap months, shape={field.shape}")
    print(f"units={field.units} baseline={field.anomaly_baseline}")

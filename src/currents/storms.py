"""NOAA NCEI IBTrACS v4 tropical-cyclone best-track acquisition.

:func:`fetch_ibtracs` downloads the IBTrACS best-track archive (keyless
NOAA NCEI HTTPS, NetCDF-4) once into a SHA-256-verified local cache and
returns a :class:`StormField`: per-storm tracks of 3–6-hourly positions
with max sustained wind (kt), min central pressure (hPa), storm name,
SID, season, and basin.

Wind priority is USA-agency first (``usa_wind``, 1-minute sustained,
the NHC/JTWC operational value), falling back to ``wmo_wind`` per
observation; pressure is ``usa_pres`` then ``wmo_pres``. Longitudes are
normalized to the conventional -180..180 axis (the archive mixes
-180..180 and 0..360 longitudes); tracks that cross the dateline keep
continuous longitudes in the field and renderers break polyline
segments on |Δlon| > 180° — documented, never silently wrapped.

Access truth (verified live 2026-09-27 — see docs/DATA_SOURCES.md):

* ``https://www.ncei.noaa.gov/data/international-best-track-archive-
  for-climate-stewardship-ibtracs/v04r01/access/netcdf/`` — keyless
  anonymous HTTPS (HTTP 200, ``Accept-Ranges: bytes``).
* ``IBTrACS.since1980.v04r01.nc`` (10,811,444 bytes, 1980–present) is
  the default: small enough to cache whole.
* ``IBTrACS.ALL.v04r01.nc`` (23,386,955 bytes, 1842–present) is used
  when ``full_archive=True``.
* Both files are NetCDF-4/HDF5, so reading needs the ``netCDF4``
  package — a LAZY optional import (like ``cdsapi`` for ERA5): the
  module imports cleanly without it and only the read path raises an
  actionable ImportError. The pure parsing logic
  (:func:`_parse_ibtracs_dataset`) works over a duck-typed dataset, so
  the whole test suite runs offline with fake datasets.

Cache discipline (mirrors :mod:`currents.basemaps`): the archive file
is downloaded once into ``$SURVEY_CURRENTS_CACHE/storms`` (else
``~/.cache/survey-currents/storms``) with atomic writes and a
``.sha256`` sidecar; a corrupt cache entry is re-downloaded. IBTrACS
v04r01 is re-released as new storms are added, so entries older than
``max_cache_age_days`` (default 30) are revalidated against the
remote ``Last-Modified`` header and refreshed when the server copy is
newer — the cache never silently goes stale, and ``refresh=True``
forces a re-download.

Subsetting selects STORMS (any observation inside the bbox and the
time window keeps the storm) and clips observations to the time
window; spatial clipping is deliberately NOT applied — cutting a
polyline at the bbox edge would misrepresent the track. Optional
``min_wind`` keeps storms whose lifetime maximum sustained wind
(kt) reaches the threshold; optional ``storm_name`` selects a single
named storm (case-insensitive exact match on the stripped IBTrACS
name, e.g. ``"katrina"``).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import os
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from dataclasses import field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date, _coerce_datetime_utc
from .sst_global import validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]

__all__ = [
    "IBTRACS_VERSION",
    "IBTRACS_NC_DIR",
    "IBTRACS_SINCE1980_FILE",
    "IBTRACS_ALL_FILE",
    "IBTRACS_START_SINCE1980",
    "IBTRACS_START_ALL",
    "SSHS_CATEGORIES",
    "SSHS_COLORS",
    "SSHS_LABELS",
    "sshs_category",
    "StormTrack",
    "StormField",
    "storm_cache_dir",
    "ibtracs_filename",
    "ibtracs_url",
    "fetch_ibtracs",
]

#: Pinned IBTrACS release.
IBTRACS_VERSION = "v04r01"
#: Keyless NCEI HTTPS directory (verified live 2026-09-27).
IBTRACS_NC_DIR = (
    "https://www.ncei.noaa.gov/data/international-best-track-archive-"
    "for-climate-stewardship-ibtracs/v04r01/access/netcdf/"
)
#: Default file: 1980–present, ~10.8 MB.
IBTRACS_SINCE1980_FILE = "IBTrACS.since1980.v04r01.nc"
#: Full archive: 1842–present, ~23.4 MB.
IBTRACS_ALL_FILE = "IBTrACS.ALL.v04r01.nc"
#: Record starts.
IBTRACS_START_SINCE1980 = _dt.date(1980, 1, 1)
IBTRACS_START_ALL = _dt.date(1842, 1, 1)
#: NetCDF fill sentinels (verified against the live file 2026-09-27).
IBTRACS_FILL_INT = -9999
IBTRACS_SSHS_FILL = -15

#: Saffir-Simpson hurricane wind scale (1-minute sustained wind, kt):
#: (lower_bound_inclusive, upper_bound_exclusive, code). Below 34 kt is
#: a tropical depression; 34–63 kt a tropical storm.
SSHS_CATEGORIES: Tuple[Tuple[float, float, str], ...] = (
    (0.0, 34.0, "TD"),
    (34.0, 64.0, "TS"),
    (64.0, 83.0, "C1"),
    (83.0, 96.0, "C2"),
    (96.0, 113.0, "C3"),
    (113.0, 137.0, "C4"),
    (137.0, float("inf"), "C5"),
)
#: Track colors per category (NHC-style palette, documented in
#: survey-viz docs/STORMS.md).
SSHS_COLORS: Dict[str, str] = {
    "TD": "#5ebaff",
    "TS": "#00faf4",
    "C1": "#ffffcc",
    "C2": "#ffe775",
    "C3": "#ffc140",
    "C4": "#ff8f20",
    "C5": "#ff6060",
    "unknown": "#9aa0a6",
}
SSHS_LABELS: Dict[str, str] = {
    "TD": "Tropical Depression (<34 kt)",
    "TS": "Tropical Storm (34–63 kt)",
    "C1": "Category 1 (64–82 kt)",
    "C2": "Category 2 (83–95 kt)",
    "C3": "Category 3 (96–112 kt)",
    "C4": "Category 4 (113–136 kt)",
    "C5": "Category 5 (≥137 kt)",
    "unknown": "Unknown intensity",
}

#: Default cache freshness: revalidate files older than this many days.
IBTRACS_MAX_CACHE_AGE_DAYS = 30


def sshs_category(wind_kt: float) -> str:
    """Saffir-Simpson code for ``wind_kt`` (``"unknown"`` when NaN)."""
    try:
        w = float(wind_kt)
    except (TypeError, ValueError):
        return "unknown"
    if not (w == w) or w < 0:  # NaN or negative
        return "unknown"
    for lo, hi, code in SSHS_CATEGORIES:
        if lo <= w < hi:
            return code
    return "unknown"


def storm_cache_dir() -> str:
    """User cache root for IBTrACS downloads.

    ``$SURVEY_CURRENTS_CACHE/storms`` when set, else
    ``~/.cache/survey-currents/storms``. Created on demand.
    """
    root = os.environ.get("SURVEY_CURRENTS_CACHE") or os.path.join(
        os.path.expanduser("~"), ".cache", "survey-currents")
    path = os.path.join(root, "storms")
    os.makedirs(path, exist_ok=True)
    return path


def ibtracs_filename(full_archive: bool = False) -> str:
    """Archive file name for the requested coverage."""
    return IBTRACS_ALL_FILE if full_archive else IBTRACS_SINCE1980_FILE


def ibtracs_url(full_archive: bool = False) -> str:
    """Keyless NCEI HTTPS URL of the archive file."""
    return IBTRACS_NC_DIR + ibtracs_filename(full_archive)


def _require_netcdf4():
    """Lazily import netCDF4 (optional dependency)."""
    try:
        import netCDF4
    except ImportError as exc:
        raise ImportError(
            "reading IBTrACS needs the optional 'netCDF4' package "
            "(the archive is NetCDF-4/HDF5).\nInstall it with:\n\n"
            "    pip install netCDF4\n\n"
            "StormField.synthetic() and the parsing helpers keep working "
            "without it."
        ) from exc
    return netCDF4


# ---------------------------------------------------------------------------
# Download + cache
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


def _tool_version() -> str:
    from . import __version__
    return __version__


def _download_bytes(url: str, timeout: int = 600) -> bytes:
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
    raise RuntimeError(f"IBTrACS download failed after 3 attempts: {url}") from last


def ensure_ibtracs_file(full_archive: bool = False,
                        cache_dir: Optional[str] = None,
                        max_cache_age_days: int = IBTRACS_MAX_CACHE_AGE_DAYS,
                        refresh: bool = False) -> Tuple[str, Dict[str, Any]]:
    """Return ``(local_path, download_provenance)`` for the archive file.

    Downloads once; reuses the cached copy while its ``.sha256``
    sidecar verifies. Because IBTrACS v04r01 is re-released as new
    storms are added, a cached copy older than ``max_cache_age_days``
    is revalidated against the remote ``Last-Modified`` header and
    refreshed when the server copy is newer; ``refresh=True`` forces a
    fresh download. A corrupt cache entry (sidecar mismatch) is always
    re-downloaded.
    """
    cache = cache_dir or storm_cache_dir()
    os.makedirs(cache, exist_ok=True)
    filename = ibtracs_filename(full_archive)
    path = os.path.join(cache, filename)
    url = ibtracs_url(full_archive)
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
        # Revalidate freshness: IBTrACS is re-released with new storms.
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
        local_lm = None
        lm_path = path + ".last_modified"
        if os.path.isfile(lm_path):
            with open(lm_path, "r", encoding="utf-8") as fh:
                local_lm = fh.read().strip() or None
        if remote_lm and remote_lm == local_lm:
            provenance.update(cache_hit=True, sha256=digest,
                              cache_age_days=age_days,
                              revalidated=True)
            return path, provenance
        # Server copy is newer (or revalidation failed) — refresh.
        digest = _fresh_download()
        if remote_lm:
            with open(lm_path, "w", encoding="utf-8") as fh:
                fh.write(remote_lm)
        provenance.update(downloaded=True, sha256=digest,
                          refreshed=True)
        return path, provenance

    # No usable cache (or forced refresh): download.
    digest = _fresh_download()
    remote_lm = _http_head_last_modified(url)
    if remote_lm:
        with open(path + ".last_modified", "w", encoding="utf-8") as fh:
            fh.write(remote_lm)
    provenance.update(downloaded=True, sha256=digest)
    return path, provenance


# ---------------------------------------------------------------------------
# NetCDF parsing (pure over a duck-typed dataset — offline-testable)
# ---------------------------------------------------------------------------

def _chars_to_str(arr: Any) -> str:
    """Decode a trailing-dimension char array (``|S1``) to a clean str."""
    raw = np.asarray(arr)
    if raw.dtype.kind == "S":
        blob = b"".join(bytes(x) for x in raw.reshape(-1))
    else:
        blob = b"".join(bytes(str(x), "ascii", "replace")
                         for x in raw.reshape(-1))
    return blob.decode("ascii", "replace").replace("\x00", " ").strip()


def _read_float(var: Any, fill: Optional[float] = None) -> np.ndarray:
    """Read a numeric variable as float64 with missing -> NaN."""
    data = np.asarray(var[:])
    if np.ma.isMaskedArray(data):
        data = np.ma.filled(data, np.nan)
    out = np.asarray(data, dtype=float)
    sentinel = float(fill if fill is not None
                     else getattr(var, "_FillValue", "nan"))
    if sentinel == sentinel:  # not NaN
        out[out == sentinel] = np.nan
    return out


def _parse_iso_time(raw: bytes) -> Optional[_dt.datetime]:
    try:
        text = raw.decode("ascii", "replace").replace("\x00", " ").strip()
    except Exception:
        return None
    if len(text) < 10 or text.startswith(" "):
        return None
    try:
        return _dt.datetime.strptime(text[:19], "%Y-%m-%d %H:%M:%S").replace(
            tzinfo=_dt.timezone.utc)
    except ValueError:
        return None


def _norm_lon(lon: float) -> float:
    """Normalize any longitude to the conventional [-180, 180) axis."""
    return ((float(lon) + 180.0) % 360.0) - 180.0


def _in_bbox(lon: float, lat: float,
             bbox: Tuple[float, float, float, float]) -> bool:
    lon_min, lat_min, lon_max, lat_max = bbox
    in_lat = lat_min <= lat <= lat_max
    if lon_min <= lon_max:
        return bool(in_lat and lon_min <= lon <= lon_max)
    # Antimeridian-crossing bbox.
    return bool(in_lat and (lon >= lon_min or lon <= lon_max))


def _parse_ibtracs_dataset(ds: Any, bbox: Tuple[float, float, float, float],
                           d0: _dt.date, d1: _dt.date,
                           min_wind: Optional[float] = None,
                           storm_name: Optional[str] = None) -> List["StormTrack"]:
    """Parse matching storms out of an open IBTrACS NetCDF dataset.

    ``ds`` is duck-typed (``ds.variables[name][:]``) so tests can pass
    fake datasets without netCDF4. Wind priority per observation is
    ``usa_wind`` then ``wmo_wind``; pressure ``usa_pres`` then
    ``wmo_pres``. Observations are clipped to ``[d0, d1]``; a storm is
    kept when at least one in-window observation falls inside
    ``bbox``. Longitudes are normalized to [-180, 180).
    """
    v = ds.variables
    numobs = np.asarray(v["numobs"][:]).astype(int).ravel()
    n_storms = numobs.shape[0]
    lats = _read_float(v["lat"])
    lons = _read_float(v["lon"])
    usa_wind = _read_float(v["usa_wind"], IBTRACS_FILL_INT)
    wmo_wind = _read_float(v["wmo_wind"], IBTRACS_FILL_INT)
    usa_pres = _read_float(v["usa_pres"], IBTRACS_FILL_INT)
    wmo_pres = _read_float(v["wmo_pres"], IBTRACS_FILL_INT)
    sid_raw = v["sid"][:]
    name_raw = v["name"][:]
    season = np.asarray(v["season"][:]).astype(int).ravel()
    basin_raw = v["basin"][:]
    iso_raw = v["iso_time"][:]

    want_name = (storm_name or "").strip().upper()
    tracks: List[StormTrack] = []
    for i in range(n_storms):
        n = int(numobs[i])
        if n <= 0:
            continue
        name = _chars_to_str(name_raw[i])
        if want_name and name.upper() != want_name:
            continue
        sid = _chars_to_str(sid_raw[i])
        obs_times: List[_dt.datetime] = []
        obs_lats: List[float] = []
        obs_lons: List[float] = []
        obs_wind: List[float] = []
        obs_pres: List[float] = []
        for j in range(n):
            t = _parse_iso_time(bytes(bytearray(
                np.asarray(iso_raw[i, j]).astype("S1").tobytes())))
            if t is None:
                continue
            if not (d0 <= t.date() <= d1):
                continue
            la, lo = lats[i, j], lons[i, j]
            if not (la == la and lo == lo):
                continue
            w = usa_wind[i, j]
            if not (w == w):
                w = wmo_wind[i, j]
            p = usa_pres[i, j]
            if not (p == p):
                p = wmo_pres[i, j]
            obs_times.append(t)
            obs_lats.append(float(la))
            obs_lons.append(_norm_lon(float(lo)))
            obs_wind.append(float(w) if w == w else float("nan"))
            obs_pres.append(float(p) if p == p else float("nan"))
        if not obs_times:
            continue
        if not any(_in_bbox(lo, la, bbox)
                   for lo, la in zip(obs_lons, obs_lats)):
            continue
        max_wind = max((w for w in obs_wind if w == w), default=float("nan"))
        if min_wind is not None and not (max_wind >= float(min_wind)):
            continue
        basin = _chars_to_str(basin_raw[i, 0])
        tracks.append(StormTrack(
            sid=sid, name=name, season=int(season[i]), basin=basin,
            times=obs_times, lats=np.array(obs_lats),
            lons=np.array(obs_lons),
            winds=np.array(obs_wind), press=np.array(obs_pres)))
    # Deterministic order: season, then first-observation time.
    tracks.sort(key=lambda t: (t.season, t.times[0].timestamp()))
    return tracks


# ---------------------------------------------------------------------------
# StormTrack / StormField
# ---------------------------------------------------------------------------

@dataclass
class StormTrack:
    """One tropical-cyclone track: time-ordered best-track fixes.

    ``winds`` are max sustained winds in kt (USA-agency priority, WMO
    fallback; NaN when neither agency reported), ``press`` min central
    pressure in hPa (same priority). Longitudes are on [-180, 180);
    dateline-crossing tracks keep continuous values here — renderers
    break polyline segments on |Δlon| > 180°.
    """

    sid: str
    name: str
    season: int
    basin: str
    times: List[_dt.datetime]
    lats: np.ndarray
    lons: np.ndarray
    winds: np.ndarray
    press: np.ndarray

    def __post_init__(self) -> None:
        n = len(self.times)
        self.lats = np.asarray(self.lats, dtype=float).reshape(n)
        self.lons = np.asarray(self.lons, dtype=float).reshape(n)
        self.winds = np.asarray(self.winds, dtype=float).reshape(n)
        self.press = np.asarray(self.press, dtype=float).reshape(n)

    def __len__(self) -> int:
        return len(self.times)

    @property
    def max_wind(self) -> float:
        """Lifetime maximum sustained wind (kt); NaN when unreported."""
        valid = self.winds[self.winds == self.winds]
        return float(np.max(valid)) if valid.size else float("nan")

    @property
    def min_pres(self) -> float:
        """Lifetime minimum central pressure (hPa); NaN when unreported."""
        valid = self.press[self.press == self.press]
        return float(np.min(valid)) if valid.size else float("nan")

    @property
    def max_category(self) -> str:
        """Saffir-Simpson code of the lifetime maximum wind."""
        return sshs_category(self.max_wind)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sid": self.sid,
            "name": self.name,
            "season": self.season,
            "basin": self.basin,
            "times": [t.isoformat() for t in self.times],
            "lats": [float(x) for x in self.lats],
            "lons": [float(x) for x in self.lons],
            "winds": [float(x) for x in self.winds],
            "press": [float(x) for x in self.press],
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StormTrack":
        return cls(
            sid=str(data["sid"]), name=str(data["name"]),
            season=int(data["season"]), basin=str(data.get("basin", "")),
            times=[_dt.datetime.fromisoformat(t) for t in data["times"]],
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            winds=np.asarray(data["winds"], dtype=float),
            press=np.asarray(data["press"], dtype=float))


@dataclass
class StormField:
    """A set of IBTrACS storm tracks for one query.

    ``provenance`` records the IBTrACS version, file URL, SHA-256, and
    the subset parameters — following the :mod:`currents.sst_global`
    conventions.
    """

    tracks: List[StormTrack]
    bbox: Tuple[float, float, float, float]
    start: _dt.date
    end: _dt.date
    min_wind: Optional[float] = None
    storm_name: str = ""
    source: str = "ibtracs"
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        self.bbox = validate_sst_bbox(self.bbox)
        self.start = _coerce_date(self.start)
        self.end = _coerce_date(self.end)
        if self.start > self.end:
            raise ValueError(
                f"StormField: start {self.start} is after end {self.end}")
        self.tracks = list(self.tracks)
        self.storm_name = str(self.storm_name or "")

    def __len__(self) -> int:
        return len(self.tracks)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    @property
    def n_obs(self) -> int:
        """Total best-track fixes across all tracks."""
        return sum(len(t) for t in self.tracks)

    @property
    def time_range(self) -> Tuple[Optional[_dt.datetime],
                                  Optional[_dt.datetime]]:
        """(earliest, latest) fix time, or (None, None) when empty."""
        if not self.tracks or not any(len(t) for t in self.tracks):
            return None, None
        lo = min(t.times[0] for t in self.tracks if len(t))
        hi = max(t.times[-1] for t in self.tracks if len(t))
        return lo, hi

    def rank_by_intensity(self) -> List[StormTrack]:
        """Tracks sorted by lifetime max wind, descending.

        The documented aggregation rule behind "strongest hurricanes"
        requests: intensity ranking = lifetime maximum sustained wind
        (kt, USA-agency priority per fix), NaN-max tracks last.
        """
        return sorted(
            self.tracks,
            key=lambda t: (t.max_wind != t.max_wind,  # NaN last
                           -(t.max_wind if t.max_wind == t.max_wind else 0.0)))

    # -- filters ---------------------------------------------------------

    def select_time(self, start: DateLike, end: DateLike) -> "StormField":
        """Fixes with ``start <= time <= end`` (inclusive); empty tracks dropped."""
        d0, d1 = _coerce_date(start), _coerce_date(end)
        kept: List[StormTrack] = []
        for t in self.tracks:
            idx = [k for k, tm in enumerate(t.times) if d0 <= tm.date() <= d1]
            if not idx:
                continue
            kept.append(StormTrack(
                sid=t.sid, name=t.name, season=t.season, basin=t.basin,
                times=[t.times[k] for k in idx],
                lats=t.lats[idx], lons=t.lons[idx],
                winds=t.winds[idx], press=t.press[idx]))
        return StormField(tracks=kept, bbox=self.bbox, start=d0, end=d1,
                          min_wind=self.min_wind, storm_name=self.storm_name,
                          source=self.source,
                          provenance=dict(self.provenance))

    def select_bbox(self, bbox: Sequence[float]) -> "StormField":
        """Storms with at least one fix inside ``bbox`` (full tracks kept)."""
        box = validate_sst_bbox(bbox)
        kept = [t for t in self.tracks
                if any(_in_bbox(lo, la, box)
                       for lo, la in zip(t.lons, t.lats))]
        return StormField(tracks=kept, bbox=box, start=self.start,
                          end=self.end, min_wind=self.min_wind,
                          storm_name=self.storm_name, source=self.source,
                          provenance=dict(self.provenance))

    # -- serialization ---------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "storm_tracks": [t.to_dict() for t in self.tracks],
            "bbox": list(self.bbox),
            "start": self.start.isoformat(),
            "end": self.end.isoformat(),
            "min_wind": self.min_wind,
            "storm_name": self.storm_name,
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StormField":
        return cls(
            tracks=[StormTrack.from_dict(t)
                    for t in data.get("storm_tracks", [])],
            bbox=tuple(data["bbox"]),
            start=_coerce_date(data["start"]),
            end=_coerce_date(data["end"]),
            min_wind=data.get("min_wind"),
            storm_name=data.get("storm_name", ""),
            source=data.get("source", "ibtracs"),
            provenance=dict(data.get("provenance", {})))

    def to_json(self, path: str) -> str:
        """Write the field as JSON; return ``path``."""
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)
        return path

    @classmethod
    def from_json(cls, path: str) -> "StormField":
        """Read a field written by :meth:`to_json`."""
        import json
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture -------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-100.0, 10.0, -60.0, 40.0),
                  start: DateLike = "2024-08-01",
                  end: DateLike = "2024-09-15",
                  n_storms: int = 3, seed: int = 7,
                  source: str = "synthetic") -> "StormField":
        """Deterministic synthetic tracks (offline tests / demos).

        Builds ``n_storms`` curved tracks inside ``bbox`` with rising
        then falling winds (genesis -> peak -> decay). The last track
        crosses the dateline when the bbox touches it — otherwise it
        stays inside the bbox. One track is always named ``"TESTALPHA"``
        so name-selection tests are deterministic.
        """
        rng = np.random.default_rng(seed)
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        d0, d1 = _coerce_date(start), _coerce_date(end)
        if d1 < d0:
            raise ValueError(f"end {d1} is before start {d0}")
        ndays = (d1 - d0).days + 1
        names = ["TESTALPHA", "TESTBRAVO", "TESTCHARLIE", "TESTDELTA",
                 "TESTECHO", "TESTFOXTROT"]
        tracks: List[StormTrack] = []
        for s in range(n_storms):
            n = int(rng.integers(12, 30))
            t0 = _dt.datetime(d0.year, d0.month, d0.day,
                              tzinfo=_dt.timezone.utc) + _dt.timedelta(
                days=int(rng.integers(0, max(1, ndays - n // 4))))
            times = [t0 + _dt.timedelta(hours=6 * k) for k in range(n)]
            # Genesis in the lower-left, curving poleward and eastward.
            lon0 = lon_min + 0.15 * (lon_max - lon_min) * rng.random()
            lat0 = lat_min + 0.15 * (lat_max - lat_min) * rng.random()
            k = np.arange(n, dtype=float)
            lons = lon0 + 0.55 * k + 0.4 * np.sin(k / 3.0)
            lats = lat0 + 0.35 * k + 0.2 * np.cos(k / 4.0)
            if lon_min <= 170.0 or lon_max >= -170.0:
                # Keep inside the bbox (wrap only for dateline bboxes).
                span = (lon_max - lon_min) % 360.0 or 360.0
                lons = lon_min + (lons - lon_min) % span
                lons = np.where(lons > 180.0, lons - 360.0, lons)
            peak = float(rng.uniform(65.0, 145.0))
            winds = peak * np.sin(np.pi * (k + 1) / (n + 1)) ** 0.7
            winds = np.clip(winds, 20.0, None)
            press = 1010.0 - 0.55 * (winds - 20.0)
            tracks.append(StormTrack(
                sid=f"2099{s:03d}N00000",
                name=names[s % len(names)],
                season=2099, basin="NA",
                times=times, lats=np.clip(lats, -90.0, 90.0), lons=lons,
                winds=winds, press=press))
        return cls(tracks=tracks, bbox=(lon_min, lat_min, lon_max, lat_max),
                   start=d0, end=d1, source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# Public fetch
# ---------------------------------------------------------------------------

def fetch_ibtracs(bbox: Sequence[float], start: DateLike, end: DateLike,
                  min_wind: Optional[float] = None,
                  storm_name: Optional[str] = None,
                  full_archive: bool = False,
                  cache_dir: Optional[str] = None,
                  max_cache_age_days: int = IBTRACS_MAX_CACHE_AGE_DAYS,
                  refresh: bool = False) -> StormField:
    """Fetch IBTrACS v4 best tracks for ``bbox`` x ``[start, end]``.

    Downloads the archive NetCDF once (cached, SHA-256-verified,
    freshness-revalidated — see :func:`ensure_ibtracs_file`), then
    selects storms with at least one fix inside ``bbox`` and clips
    fixes to the date window. ``min_wind`` (kt) keeps storms whose
    lifetime maximum sustained wind reaches the threshold;
    ``storm_name`` selects one named storm (case-insensitive);
    ``full_archive=True`` uses the 1842–present file instead of the
    default 1980–present one.
    """
    box = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    if d1 < d0:
        raise ValueError(f"end {d1.date()} is before start {d0.date()}")
    earliest = IBTRACS_START_ALL if full_archive else IBTRACS_START_SINCE1980
    if d0.date() < earliest:
        raise ValueError(
            f"IBTrACS {'full archive' if full_archive else 'since-1980 file'} "
            f"starts {earliest.isoformat()}; start {d0.date().isoformat()} "
            "is before the record (pass full_archive=True for 1842–1979)")
    if d1.date() > _dt.date.today() + _dt.timedelta(days=2):
        raise ValueError(f"end {d1.date().isoformat()} is in the future")
    if min_wind is not None and float(min_wind) < 0:
        raise ValueError(f"min_wind must be >= 0, got {min_wind}")

    path, dl_prov = ensure_ibtracs_file(
        full_archive=full_archive, cache_dir=cache_dir,
        max_cache_age_days=max_cache_age_days, refresh=refresh)
    netCDF4 = _require_netcdf4()
    with netCDF4.Dataset(path, "r") as ds:
        tracks = _parse_ibtracs_dataset(
            ds, box, d0.date(), d1.date(),
            min_wind=float(min_wind) if min_wind is not None else None,
            storm_name=storm_name)

    provenance = {
        "source": "ibtracs",
        "ibtracs_version": IBTRACS_VERSION,
        "coverage": "1842-present" if full_archive else "1980-present",
        "file_url": dl_prov["file_url"],
        "sha256": dl_prov.get("sha256"),
        "retrieved_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "bbox": list(box),
        "start": d0.date().isoformat(),
        "end": d1.date().isoformat(),
        "min_wind": float(min_wind) if min_wind is not None else None,
        "storm_name": (storm_name or "").strip(),
        "full_archive": bool(full_archive),
        "n_storms": len(tracks),
        "n_obs": sum(len(t) for t in tracks),
        "cache_hit": dl_prov.get("cache_hit", False),
        "wind_priority": "usa_wind then wmo_wind (kt, 1-min sustained)",
        "pressure_priority": "usa_pres then wmo_pres (hPa)",
        "tool": f"survey-currents {_tool_version()}",
    }
    return StormField(tracks=tracks, bbox=box, start=d0.date(),
                      end=d1.date(),
                      min_wind=float(min_wind) if min_wind is not None else None,
                      storm_name=(storm_name or "").strip(),
                      source="ibtracs", provenance=provenance)

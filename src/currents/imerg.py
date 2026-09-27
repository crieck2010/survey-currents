"""NASA GPM IMERG precipitation acquisition via GES DISC HTTPS.

:func:`fetch_imerg` pulls half-hourly observed precipitation
(2000-present, 0.1-degree global grid) from NASA's Goddard Earth
Sciences Data and Information Services Center (GES DISC):

* ``run="early"`` — GPM_3IMERGHHE (latency ~4 hours)
* ``run="late"`` — GPM_3IMERGHHL (latency ~14 hours; the sane default)
* ``run="final"`` — GPM_3IMERGHH (latency ~3.5 months; gauge-adjusted)

The on-wire ``precipitation`` variable is a mm/hr *rate* (gauge-adjusted
merged microwave-infrared estimate, formerly ``precipitationCal``);
``accumulate="daily"`` (default) sums the 48 half-hourly rates x 0.5 h
into daily totals (mm/day), while ``accumulate="native"`` passes the
half-hourly rates through unchanged.

Access truth (verified live 2026-09-26 — see docs/DATA_SOURCES.md):

* Catalog browsing and OPeNDAP metadata (``.dds``/``.das``) at
  ``https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/...`` are
  keyless (HTTP 200).
* Actual data downloads (``/data/GPM_L3/...``) require a free Earthdata
  Login account: unauthenticated requests 302 to
  ``urs.earthdata.nasa.gov`` and then 401. Without credentials a
  :class:`CredentialsMissing` error explains the exact setup (mirrors
  the OSCAR / MUR adapters).
* File naming (verified against the live catalog and NASA CMR):
  ``3B-HHR[-E|-L].MS.MRG.3IMERG.YYYYMMDD-SHHMMSS-EHHMMSS.HHMM.V07B.HDF5``
  under ``{collection}.07/{YYYY}/{DOY}/`` — ``-E``/``-L`` infix for the
  Early/Late runs, and the 4th field is the start half-hour (``HHMM``).
* Grid (from the live ``.das``): 0.1-degree, lon -179.95..179.95 (3600,
  ascending), lat -89.95..89.95 (1800, ascending); ``precipitation``
  units mm/hr, ``_FillValue`` -9999.900391 (masked to NaN); time is
  seconds since 1980-01-06 (GPS epoch, no leap seconds) — per-granule
  timestamps are derived from the filename half-hour slot instead
  (deterministic, no leap-second table needed).
* Record starts (per the NASA V07 product documentation): 2000-06-01
  for all three runs — TRMM era June 2000 - May 2014, GPM era June
  2014 - present. (The live CMR archive shows V07 reprocessing reaching
  1998-01-01 for Early/Final; the adapter keeps the documented floor.)

Downloads use stdlib ``urllib`` with an Earthdata-authenticated opener
(the URS OAuth redirect dance, mirroring the OSCAR adapter); parsing
the HDF5 payloads needs ``h5py`` (lazy import — the engine core stays
stdlib+numpy, like the netCDF4/rasterio precedents).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import os
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass, field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date
from .sst_global import CredentialsMissing as _BaseCredentialsMissing
from .sst_global import validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]


# ---------------------------------------------------------------------------
# IMERG / GES DISC constants (verified live 2026-09-26)
# ---------------------------------------------------------------------------

#: GES DISC host serving GPM IMERG.
IMERG_BASE = "https://gpm1.gesdisc.eosdis.nasa.gov"
#: GES DISC OPeNDAP catalog root (keyless browsing + .dds/.das metadata).
IMERG_OPENDAP_BASE = f"{IMERG_BASE}/opendap/GPM_L3"
#: HTTPS data root (Earthdata Login required).
IMERG_DATA_BASE = f"{IMERG_BASE}/data/GPM_L3"
#: Product version in directory names / file names.
IMERG_VERSION = "07"

#: run key -> GES DISC collection id (verified via NASA CMR 2026-09-26).
IMERG_COLLECTIONS: Dict[str, str] = {
    "early": "GPM_3IMERGHHE",
    "late": "GPM_3IMERGHHL",
    "final": "GPM_3IMERGHH",
}
#: run key -> file-name infix ("" for the Final run).
IMERG_RUN_INFIX: Dict[str, str] = {
    "early": "-E",
    "late": "-L",
    "final": "",
}
# TRMM era (June 2000 - May 2014) and the GPM era (June 2014 - present).
#: First day with granules per run. The documented V07 record starts
#: 2000-06-01 for all three runs (TRMM June 2000 - May 2014, GPM June
#: 2014 - present; Early ~4 h, Late ~14 h, Final ~3.5 months latency).
#: Note: the live CMR archive shows V07 reprocessing reaching back to
#: 1998-01-01 for the Early and Final collections, but the adapter
#: enforces the documented 2000-06-01 floor — pre-2000 TRMM-era
#: estimates come from a much smaller satellite constellation and are
#: of lower quality, especially at high latitudes.
IMERG_START: Dict[str, _dt.date] = {
    "early": _dt.date(2000, 6, 1),
    "late": _dt.date(2000, 6, 1),
    "final": _dt.date(2000, 6, 1),
}
#: Documented product latencies (NASA IMERG documentation).
IMERG_LATENCY: Dict[str, str] = {
    "early": "~4 hours",
    "late": "~14 hours",
    "final": "~3.5 months",
}
#: Native grid spacing, degrees.
IMERG_RES = 0.1
#: On-wire precipitation units (mm/hr rates).
IMERG_UNITS_NATIVE = "mm/hr"
#: Daily-accumulated precipitation units (mm/day totals).
IMERG_UNITS_DAILY = "mm/day"
#: Fallback fill value for the precipitation variable (the live .das
#: reports -9999.900391; the parser reads the dataset attribute first).
IMERG_FILL = -9999.9
#: Half-hour slots per day.
IMERG_SLOTS_PER_DAY = 48
#: Valid accumulation modes.
IMERG_ACCUMULATIONS = ("daily", "native")

#: HDF5 paths inside each granule (from the live .das fullnamepath).
IMERG_H5_PRECIP = "/Grid/precipitation"
IMERG_H5_LAT = "/Grid/lat"
IMERG_H5_LON = "/Grid/lon"


class CredentialsMissing(_BaseCredentialsMissing):
    """Earthdata Login credentials are required but were not found.

    GPM IMERG downloads from NASA GES DISC require a free Earthdata
    Login account (catalog browsing and OPeNDAP metadata are keyless;
    the data files are not). Provide credentials one of two ways:

    1. Environment variables (best for CI / servers)::

           export EARTHDATA_USERNAME="your_username"
           export EARTHDATA_PASSWORD="your_password"

    2. A ``~/.netrc`` entry (best for interactive use)::

           machine urs.earthdata.nasa.gov
           login your_username
           password your_password

    Register for free at https://urs.earthdata.nasa.gov/users/new.
    """


# ---------------------------------------------------------------------------
# pure, offline-testable helpers
# ---------------------------------------------------------------------------


def normalize_imerg_run(run: str) -> str:
    """Normalize the ``run`` argument to ``"early"``/``"late"``/``"final"``.

    Raises :class:`ValueError` on anything else.
    """
    key = str(run).strip().lower()
    if key not in IMERG_COLLECTIONS:
        raise ValueError(
            f"run must be one of {sorted(IMERG_COLLECTIONS)}, got {run!r} "
            "(early: ~4 h latency, late: ~14 h latency [default], "
            "final: ~3.5 month latency)")
    return key


def normalize_imerg_accumulate(accumulate: str) -> str:
    """Normalize the ``accumulate`` argument to ``"daily"``/``"native"``.

    ``"daily"`` (default) sums the 48 half-hourly mm/hr rates x 0.5 h
    into daily totals (mm/day); ``"native"`` passes the half-hourly
    rates through unchanged (mm/hr).
    """
    key = str(accumulate).strip().lower()
    if key not in IMERG_ACCUMULATIONS:
        raise ValueError(
            f"accumulate must be one of {sorted(IMERG_ACCUMULATIONS)}, "
            f"got {accumulate!r}")
    return key


def imerg_granule_name(day: _dt.date, half_hour: int, run: str = "late") -> str:
    """Build the GPM IMERG half-hourly granule file name.

    ``half_hour`` is the 0-47 slot index (00:00-00:30 UTC -> 0, ...,
    23:30-24:00 UTC -> 47). The 4th dot-field is the slot start as
    ``HHMM`` — verified against the live catalog (``...-S000000-E002959
    .0000.V07B.HDF5``, ``...-S003000-E005959.0030.V07B.HDF5``).
    """
    key = normalize_imerg_run(run)
    if not 0 <= half_hour <= 47:
        raise ValueError(f"half_hour must be 0..47, got {half_hour}")
    start_min = half_hour * 30
    end_min = start_min + 30
    infix = IMERG_RUN_INFIX[key]
    return (
        f"3B-HHR{infix}.MS.MRG.3IMERG."
        f"{day.year:04d}{day.month:02d}{day.day:02d}"
        f"-S{start_min // 60:02d}{start_min % 60:02d}00"
        f"-E{(end_min - 1) // 60:02d}{(end_min - 1) % 60:02d}59"
        f".{start_min // 60:02d}{start_min % 60:02d}"
        f".V{IMERG_VERSION}B.HDF5"
    )


def imerg_file_url(day: _dt.date, half_hour: int, run: str = "late") -> str:
    """Build the GES DISC HTTPS data URL for one IMERG granule.

    Pure and offline-testable. Layout ``{collection}.07/{YYYY}/{DOY}/``
    verified against the live OPeNDAP catalog 2026-09-26.
    """
    key = normalize_imerg_run(run)
    collection = IMERG_COLLECTIONS[key]
    doy = day.timetuple().tm_yday
    return (f"{IMERG_DATA_BASE}/{collection}.{IMERG_VERSION}/"
            f"{day.year:04d}/{doy:03d}/{imerg_granule_name(day, half_hour, key)}")


def imerg_sample_days(d0: _dt.date, d1: _dt.date,
                      stride_days: int = 1) -> List[_dt.date]:
    """Sample ``[d0, d1]`` every ``stride_days`` days (pure, offline-testable)."""
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")
    days: List[_dt.date] = []
    cur = d0
    while cur <= d1:
        days.append(cur)
        cur += _dt.timedelta(days=stride_days)
    return days


def _validate_imerg_dates(d0: _dt.date, d1: _dt.date, run: str) -> None:
    if d1 < d0:
        raise ValueError(f"end {d1} is before start {d0}")
    start = IMERG_START[run]
    if d0 < start:
        raise ValueError(
            f"IMERG {run}-run granules start {start.isoformat()} "
            f"(requested {d0.isoformat()})")
    if d1 > _dt.date.today():
        raise ValueError(
            f"end {d1.isoformat()} is in the future; IMERG is an observed "
            f"product with a {IMERG_LATENCY[run]} latency for the "
            f"{run} run")


def _slot_timestamp(day: _dt.date, half_hour: int) -> _dt.datetime:
    """UTC timestamp for the start of a half-hour slot (filename-derived).

    The granule's internal ``Grid/time`` is seconds since the GPS epoch
    1980-01-06 *without* leap seconds, so converting it to UTC needs a
    leap-second table — the filename slot is exact and deterministic.
    """
    return (_dt.datetime(day.year, day.month, day.day,
                         tzinfo=_dt.timezone.utc)
            + _dt.timedelta(minutes=30 * half_hour))


# ---------------------------------------------------------------------------
# HDF5 parsing (h5py is a lazy import — actionable error without it)
# ---------------------------------------------------------------------------


def _require_h5py():
    """Lazy import of h5py with an actionable error."""
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "fetch_imerg needs the h5py package (pip install h5py); "
            "the survey-currents engine itself has no hard dependency on it."
        ) from exc
    return h5py


def _imerg_axis_permutation(shape: Sequence[int], nlat: int,
                            nlon: int) -> Tuple[int, int, int]:
    """Map an IMERG precipitation dataspace to (time, lat, lon) axis order.

    The live ``.das`` reports the dimensions as ``time=1, lon=3600,
    lat=1800`` — i.e. (time, lon, lat) on the wire. The axis order is
    resolved from the actual sizes of the ``Grid/lat`` (``nlat``) and
    ``Grid/lon`` (``nlon``) vectors — an honest :class:`ValueError` if a
    future reprocessing changes the layout instead of silently
    transposing the grid.
    """
    sizes = tuple(int(s) for s in shape)
    if len(sizes) != 3:
        raise ValueError(
            f"IMERG precipitation: expected a 3-D dataspace, got shape "
            f"{sizes}")
    axes: Dict[str, int] = {}
    for i, size in enumerate(sizes):
        if size == 1:
            if "time" in axes:
                raise ValueError(
                    f"IMERG precipitation: ambiguous time axis in shape "
                    f"{sizes}")
            axes["time"] = i
        elif size == nlat:
            if "lat" in axes:
                raise ValueError(
                    f"IMERG precipitation: ambiguous latitude axis in shape "
                    f"{sizes}")
            axes["lat"] = i
        elif size == nlon:
            if "lon" in axes:
                raise ValueError(
                    f"IMERG precipitation: ambiguous longitude axis in shape "
                    f"{sizes}")
            axes["lon"] = i
    missing = {"time", "lat", "lon"} - set(axes)
    if missing:
        raise ValueError(
            f"IMERG precipitation: shape {sizes} does not match the "
            f"documented (time=1, lat={nlat}, lon={nlon}) layout "
            f"(missing: {sorted(missing)})")
    return (axes["time"], axes["lat"], axes["lon"])


def _parse_imerg_bytes(payload: bytes) -> Tuple[np.ndarray, np.ndarray,
                                               np.ndarray]:
    """Parse one IMERG half-hourly HDF5 granule payload.

    Returns ``(lats, lons, values)`` where ``values`` is the 2-D
    precipitation *rate* (mm/hr) shaped (nlat, nlon), ascending axes,
    with fill/missing cells as NaN. The lat/lon axes are read from the
    file itself (``/Grid/lat``, ``/Grid/lon``) and asserted ascending;
    the precipitation axis order is resolved from the actual dataspace
    sizes by :func:`_imerg_axis_permutation`.

    Precipitation rates are physically non-negative, so any negative
    value (the ``_FillValue`` -9999.9 and friends) becomes NaN rather
    than depending on the exact fill constant.
    """
    h5py = _require_h5py()
    if payload[:4] != b"\x89HDF":
        raise ValueError(
            "not an HDF5 payload (bad magic); the GES DISC URL may have "
            "returned an error page")
    fd, path = tempfile.mkstemp(suffix=".HDF5")
    os.close(fd)
    try:
        with open(path, "wb") as fh:
            fh.write(payload)
        with h5py.File(path, "r") as ds:
            for key in (IMERG_H5_LAT, IMERG_H5_LON, IMERG_H5_PRECIP):
                if key not in ds:
                    raise ValueError(
                        f"IMERG payload missing dataset {key!r}")
            lats = np.asarray(ds[IMERG_H5_LAT][:], dtype=float).reshape(-1)
            lons = np.asarray(ds[IMERG_H5_LON][:], dtype=float).reshape(-1)
            precip_ds = ds[IMERG_H5_PRECIP]
            perm = _imerg_axis_permutation(precip_ds.shape,
                                           len(lats), len(lons))
            arr = np.asarray(precip_ds[:], dtype=float)
            # perm = (time_axis, lat_axis, lon_axis) -> reorder to
            # (time, lat, lon), then take the single time step.
            arr = np.moveaxis(arr, perm, (0, 1, 2))[0]
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass
    if lats[0] > lats[-1]:
        raise ValueError("IMERG latitude axis is not ascending")
    if lons[0] > lons[-1]:
        raise ValueError("IMERG longitude axis is not ascending")
    # Precipitation rates are non-negative; the _FillValue (-9999.9)
    # and any other negative placeholder become NaN.
    values = np.where(arr < 0.0, np.nan, arr)
    return lats, lons, values


def imerg_index_windows(bbox: Sequence[float], lats: np.ndarray,
                        lons: np.ndarray) -> List[Tuple[int, int, int, int]]:
    """Map a -180..180 bbox onto IMERG grid index windows.

    Returns ``(lon0, lon1, lat0, lat1)`` (inclusive) — one window, or
    two when the bbox crosses the antimeridian (the lon axis is
    -180..180 ascending). Mirrors :func:`currents.era5.era5_area_windows`
    in index space.
    """
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    lats = np.asarray(lats, dtype=float)
    lons = np.asarray(lons, dtype=float)

    def lat_idx(v: float) -> int:
        return int(np.clip(np.searchsorted(lats, v, side="left"), 0, len(lats) - 1))

    def lon_idx(v: float) -> int:
        return int(np.clip(np.searchsorted(lons, v, side="left"), 0, len(lons) - 1))

    j0 = lat_idx(max(lat_min, float(lats[0])))
    j1 = lat_idx(min(lat_max, float(lats[-1])))
    if j0 > j1:
        raise ValueError(
            f"bbox latitudes {lat_min}..{lat_max} do not overlap the IMERG "
            f"grid {float(lats[0])}..{float(lats[-1])}")
    lon_lo, lon_hi = float(lons[0]), float(lons[-1])
    windows = []
    if lon_max < lon_min:  # antimeridian crossing -> two index windows
        i0 = lon_idx(max(lon_min, lon_lo))
        i1 = lon_idx(lon_hi)
        windows.append((i0, i1, j0, j1))
        i0 = lon_idx(lon_lo)
        i1 = lon_idx(min(lon_max, lon_hi))
        windows.append((i0, i1, j0, j1))
    else:
        i0 = lon_idx(max(lon_min, lon_lo))
        i1 = lon_idx(min(lon_max, lon_hi))
        if i0 > i1:
            raise ValueError(
                f"bbox longitudes {lon_min}..{lon_max} do not overlap the "
                f"IMERG grid {lon_lo}..{lon_hi}")
        windows.append((i0, i1, j0, j1))
    return windows


# ---------------------------------------------------------------------------
# authenticated download (Earthdata Login; mirrors the OSCAR adapter)
# ---------------------------------------------------------------------------


def _imerg_opener() -> urllib.request.OpenerDirector:
    """Build an Earthdata-authenticated urllib opener.

    GES DISC redirects unauthenticated data requests to URS
    (``urs.earthdata.nasa.gov``) for an OAuth handshake, so the password
    manager carries both hosts and a cookie processor keeps the session.

    Raises :class:`CredentialsMissing` when no credentials are found.
    """
    from .sst_global import earthdata_credentials
    creds = earthdata_credentials()
    if creds is None:
        raise CredentialsMissing(
            "GPM IMERG needs Earthdata Login credentials, but none were "
            "found. Set EARTHDATA_USERNAME/EARTHDATA_PASSWORD or add a "
            "~/.netrc entry for urs.earthdata.nasa.gov (register free at "
            "https://urs.earthdata.nasa.gov/users/new).")
    user, pw = creds
    mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    mgr.add_password(None, IMERG_BASE, user, pw)
    mgr.add_password(None, "https://urs.earthdata.nasa.gov", user, pw)
    auth = urllib.request.HTTPBasicAuthHandler(mgr)
    return urllib.request.build_opener(
        auth, urllib.request.HTTPCookieProcessor())


def _imerg_download(url: str, opener: urllib.request.OpenerDirector,
                    timeout: int = 300) -> bytes:
    """Download one IMERG granule; 401 -> CredentialsMissing, else RuntimeError."""
    import http.client
    import time
    req = urllib.request.Request(
        url, headers={"User-Agent": f"survey-currents/{_tool_version()}"})
    last: Optional[Exception] = None
    for attempt in range(3):
        try:
            with opener.open(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 401:
                raise CredentialsMissing(
                    "Earthdata Login rejected the credentials for GES DISC "
                    "(HTTP 401). Check EARTHDATA_USERNAME/EARTHDATA_PASSWORD "
                    "or the ~/.netrc entry for urs.earthdata.nasa.gov.") from exc
            if exc.code == 404:
                raise FileNotFoundError(
                    f"IMERG granule not found (HTTP 404): {url}. The "
                    f"{url.split('/')[5] if len(url.split('/')) > 5 else ''} "
                    "run may not cover this date yet (early ~4 h, late "
                    "~14 h, final ~3.5 month latency).") from exc
            raise RuntimeError(
                f"IMERG download failed (HTTP {exc.code}: {exc.reason}): "
                f"{url}") from exc
        except (urllib.error.URLError, http.client.RemoteDisconnected,
                http.client.IncompleteRead, TimeoutError,
                ConnectionError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"IMERG download failed after 3 attempts: {url}"
                       ) from last


# ---------------------------------------------------------------------------
# RainField — canonical precipitation field model
# ---------------------------------------------------------------------------


@dataclass
class RainField:
    """Time-indexed GPM IMERG precipitation grids.

    ``values`` is ``(nt, ny, nx)`` float on a regular lat/lon grid
    (ascending axes, longitudes -180..180); fill / missing cells are
    NaN. ``units`` is ``"mm/day"`` for ``accumulate="daily"`` (daily
    totals) or ``"mm/hr"`` for ``accumulate="native"`` (half-hourly
    rates). ``times`` are UTC datetimes (one per day, or one per
    half-hour slot for ``"native"``).

    ``provenance`` records the exact file URLs, per-file SHA-256
    digests, the run/accumulate choices, per-day granule coverage, the
    retrieval timestamp, and the tool version — following the
    :mod:`currents.sea_ice` conventions.
    """

    times: List[_dt.datetime]
    lats: np.ndarray
    lons: np.ndarray
    values: np.ndarray
    units: str = IMERG_UNITS_DAILY
    run: str = "late"
    accumulate: str = "daily"
    bbox: Tuple[float, float, float, float] = (-180.0, -90.0, 180.0, 90.0)
    source: str = "nasa-gpm/imerg"
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        self.lats = np.asarray(self.lats, dtype=float).reshape(-1)
        self.lons = np.asarray(self.lons, dtype=float).reshape(-1)
        self.values = np.asarray(self.values, dtype=float).reshape(
            nt, self.lats.shape[0], self.lons.shape[0])
        self.run = normalize_imerg_run(self.run)
        self.accumulate = normalize_imerg_accumulate(self.accumulate)
        self.bbox = validate_sst_bbox(self.bbox)

    def __len__(self) -> int:
        return len(self.times)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    @property
    def time_range(self) -> Tuple[Optional[_dt.datetime], Optional[_dt.datetime]]:
        """(earliest, latest) timestep, or (None, None) when empty."""
        if not self.times:
            return None, None
        return min(self.times), max(self.times)

    @property
    def shape(self) -> Tuple[int, int, int]:
        """(ntime, nlat, nlon)."""
        return self.values.shape  # type: ignore[return-value]

    def total(self, index: int) -> float:
        """NaN-aware spatial total of timestep ``index`` (mm or mm/day)."""
        return float(np.nansum(self.values[index]))

    # -- filters ---------------------------------------------------------

    def select_time(self, start: DateLike, end: DateLike) -> "RainField":
        """Timesteps with ``start <= date <= end`` (inclusive)."""
        d0 = _coerce_date(start)
        d1 = _coerce_date(end)
        keep = [i for i, t in enumerate(self.times) if d0 <= t.date() <= d1]
        return RainField(
            times=[self.times[i] for i in keep],
            lats=self.lats, lons=self.lons,
            values=self.values[keep],
            units=self.units, run=self.run, accumulate=self.accumulate,
            bbox=self.bbox, source=self.source,
            provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "RainField":
        """Spatial subset to ``bbox`` (must lie inside the field grid)."""
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        if lon_max < lon_min:
            raise ValueError(
                f"RainField.select_bbox does not wrap antimeridian-crossing "
                f"bboxes: {tuple(bbox)!r}")
        iy = np.flatnonzero((self.lats >= lat_min) & (self.lats <= lat_max))
        ix = np.flatnonzero((self.lons >= lon_min) & (self.lons <= lon_max))
        if len(iy) == 0 or len(ix) == 0:
            raise ValueError(
                f"bbox {tuple(bbox)!r} has no overlap with the field grid")
        return RainField(
            times=list(self.times),
            lats=self.lats[iy], lons=self.lons[ix],
            values=self.values[:, iy[:, None], ix],
            units=self.units, run=self.run, accumulate=self.accumulate,
            bbox=(lon_min, lat_min, lon_max, lat_max),
            source=self.source, provenance=dict(self.provenance),
        )

    # -- serialization ---------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "times": [t.isoformat() for t in self.times],
            "lats": [float(x) for x in self.lats],
            "lons": [float(x) for x in self.lons],
            "values": self.values.tolist(),
            "units": self.units,
            "run": self.run,
            "accumulate": self.accumulate,
            "bbox": list(self.bbox),
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RainField":
        required = ("times", "lats", "lons", "values", "bbox")
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(f"RainField dict missing keys: {missing}")
        return cls(
            times=[_dt.datetime.fromisoformat(t) for t in data["times"]],
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            values=np.asarray(data["values"], dtype=float),
            units=data.get("units", IMERG_UNITS_DAILY),
            run=data.get("run", "late"),
            accumulate=data.get("accumulate", "daily"),
            bbox=tuple(data["bbox"]),
            source=data.get("source", "nasa-gpm/imerg"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        """Write the field as JSON; return ``path``."""
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "RainField":
        """Read a field written by :meth:`to_json`."""
        import json
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture -------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-125.0, 25.0, -66.0, 49.0),
                  start: DateLike = "2024-01-01", end: DateLike = "2024-01-05",
                  resolution: float = 1.0, seed: int = 7,
                  accumulate: str = "daily",
                  source: str = "synthetic") -> "RainField":
        """Deterministic synthetic precipitation (offline tests / demos).

        A few Gaussian rain cells drift eastward over a light drizzle
        background; values are mm/day totals (``accumulate="daily"``)
        or mm/hr rates (``accumulate="native"``, 48 slots/day). A NaN
        corner exercises missing-data handling.
        """
        acc = normalize_imerg_accumulate(accumulate)
        rng = np.random.default_rng(seed)
        d0, d1 = _coerce_date(start), _coerce_date(end)
        if d1 < d0:
            raise ValueError(f"end {d1} is before start {d0}")
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        if lon_max < lon_min:
            raise ValueError("synthetic() needs a non-wrapping bbox")
        nx = max(1, int(round((lon_max - lon_min) / resolution)) + 1)
        ny = max(1, int(round((lat_max - lat_min) / resolution)) + 1)
        lons = np.linspace(lon_min, lon_max, nx)
        lats = np.linspace(lat_min, lat_max, ny)
        ndays = (d1 - d0).days + 1
        if acc == "daily":
            nsteps = ndays
            times = [_dt.datetime(d0.year, d0.month, d0.day,
                                  tzinfo=_dt.timezone.utc)
                     + _dt.timedelta(days=i) for i in range(ndays)]
            units = IMERG_UNITS_DAILY
            scale = 1.0
        else:
            nsteps = ndays * IMERG_SLOTS_PER_DAY
            t0 = _dt.datetime(d0.year, d0.month, d0.day,
                              tzinfo=_dt.timezone.utc)
            times = [t0 + _dt.timedelta(minutes=30 * i)
                     for i in range(nsteps)]
            units = IMERG_UNITS_NATIVE
            scale = 1.0 / IMERG_SLOTS_PER_DAY
        LAT, LON = np.meshgrid(lats, lons, indexing="ij")
        values = np.empty((nsteps, ny, nx), dtype=float)
        for i in range(nsteps):
            field = 0.3 + 0.2 * rng.standard_normal((ny, nx))
            for c in range(3):
                cx = lon_min + (lon_max - lon_min) * (
                    (0.2 + 0.25 * c + 0.02 * i) % 1.0)
                cy = lat_min + (lat_max - lat_min) * (0.3 + 0.2 * c)
                field += (18.0 * scale
                          * np.exp(-(((LON - cx) / 6.0) ** 2
                                     + ((LAT - cy) / 4.0) ** 2)))
            values[i] = np.clip(field, 0.0, None)
        values[:, :1, :1] = np.nan  # missing-data corner
        return cls(times=times, lats=lats, lons=lons, values=values,
                   units=units, run="late", accumulate=acc,
                   bbox=(lon_min, lat_min, lon_max, lat_max),
                   source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# daily accumulation
# ---------------------------------------------------------------------------


def accumulate_daily(times: List[_dt.datetime],
                     frames: List[np.ndarray]
                     ) -> Tuple[List[_dt.datetime], np.ndarray]:
    """Accumulate half-hourly mm/hr rate frames into daily mm/day totals.

    ``times``/``frames`` are parallel (one 2-D rate grid per half-hour
    slot, NaN for fill/missing). Each day's total is
    ``sum(rate_i * 0.5 h)`` over that day's available slots; days with
    no slots at all are dropped. Returns ``(day_times, daily_values)``
    with ``daily_values`` shaped (ndays, ny, nx).

    Pure and offline-testable.
    """
    if len(times) != len(frames):
        raise ValueError(
            f"times ({len(times)}) and frames ({len(frames)}) differ in length")
    buckets: Dict[_dt.date, List[np.ndarray]] = {}
    for t, frame in zip(times, frames):
        buckets.setdefault(t.date(), []).append(np.asarray(frame, dtype=float))
    day_times: List[_dt.datetime] = []
    daily: List[np.ndarray] = []
    for day in sorted(buckets):
        day_frames = buckets[day]
        if not day_frames:
            continue
        total = np.zeros_like(day_frames[0])
        for frame in day_frames:
            total = total + np.where(np.isnan(frame), 0.0, frame) * 0.5
        # Cells NaN in every slot stay NaN (all-fill grid points).
        all_nan = np.ones_like(day_frames[0], dtype=bool)
        for frame in day_frames:
            all_nan = all_nan & np.isnan(frame)
        total[all_nan] = np.nan
        day_times.append(_dt.datetime(day.year, day.month, day.day,
                                      tzinfo=_dt.timezone.utc))
        daily.append(total)
    if not day_times:
        raise ValueError("accumulate_daily: no frames to accumulate")
    return day_times, np.stack(daily)


# ---------------------------------------------------------------------------
# public fetch
# ---------------------------------------------------------------------------


def fetch_imerg(bbox: Sequence[float], start: DateLike, end: DateLike,
                accumulate: str = "daily", run: str = "late",
                stride_days: int = 1, timeout: int = 300) -> RainField:
    """Fetch NASA GPM IMERG half-hourly precipitation for ``bbox``.

    Args:
        bbox: (lon_min, lat_min, lon_max, lat_max) in -180..180 degrees;
            antimeridian-crossing boxes are subset as two index windows.
        start/end: inclusive date range (dates, datetimes, or ISO strings).
        accumulate: ``"daily"`` (default) — sum the 48 half-hourly
            mm/hr rates x 0.5 h into daily totals (mm/day);
            ``"native"`` — pass the half-hourly rates through (mm/hr).
        run: ``"early"`` (GPM_3IMERGHHE, ~4 h latency), ``"late"``
            (GPM_3IMERGHHL, ~14 h latency — the sane default), or
            ``"final"`` (GPM_3IMERGHH, ~3.5 month latency,
            gauge-adjusted).
        stride_days: keep every Nth day (default 1).
        timeout: per-request HTTP timeout in seconds.

    Returns:
        A :class:`RainField` with per-file provenance (exact URLs,
        SHA-256 digests, run/accumulate choices, per-day granule
        coverage, retrieval time). Granules that 404 (latency window,
        processing gaps) are skipped with a provenance note; days with
        no granules at all are dropped from the field.

    Raises:
        CredentialsMissing: no Earthdata Login credentials found.
        ImportError: ``h5py`` is not installed.
        ValueError: bad bbox / dates / run / accumulate, or no granules
            retrieved at all.
    """
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    d0, d1 = _coerce_date(start), _coerce_date(end)
    key = normalize_imerg_run(run)
    acc = normalize_imerg_accumulate(accumulate)
    _validate_imerg_dates(d0, d1, key)
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")

    _require_h5py()  # fail fast with the actionable ImportError
    opener = _imerg_opener()  # raises CredentialsMissing early

    days = imerg_sample_days(d0, d1, stride_days)
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()

    slot_times: List[_dt.datetime] = []
    slot_frames: List[np.ndarray] = []
    file_records: List[Dict[str, Any]] = []
    skipped: List[str] = []
    day_coverage: Dict[str, Dict[str, int]] = {}
    grid_lats: Optional[np.ndarray] = None
    grid_lons: Optional[np.ndarray] = None

    for day in days:
        day_key = day.isoformat()
        n_expected = IMERG_SLOTS_PER_DAY
        n_got = 0
        for half_hour in range(IMERG_SLOTS_PER_DAY):
            url = imerg_file_url(day, half_hour, key)
            try:
                payload = _imerg_download(url, opener, timeout=timeout)
            except FileNotFoundError as exc:
                skipped.append(f"{day_key} slot {half_hour:02d}: 404")
                continue
            except Exception as exc:
                raise RuntimeError(
                    f"IMERG download failed for {url} "
                    f"({type(exc).__name__}: {exc})") from exc
            try:
                lats, lons, frame = _parse_imerg_bytes(payload)
            except ValueError as exc:
                skipped.append(f"{day_key} slot {half_hour:02d}: {exc}")
                continue
            if grid_lats is None:
                grid_lats, grid_lons = lats, lons
            windows = imerg_index_windows(
                (lon_min, lat_min, lon_max, lat_max), grid_lats, grid_lons)
            parts = [frame[j0:j1 + 1, i0:i1 + 1]
                     for i0, i1, j0, j1 in windows]
            sub = np.concatenate(parts, axis=1)
            # Normalize longitudes to -180..180, sorted increasing
            # (antimeridian windows concatenate out of order).
            sub_lons = np.concatenate(
                [grid_lons[i0:i1 + 1] for i0, i1, _j0, _j1 in windows])
            order = np.argsort(sub_lons, kind="stable")
            sub = sub[:, order]
            slot_times.append(_slot_timestamp(day, half_hour))
            slot_frames.append(sub)
            n_got += 1
            file_records.append({
                "url": url,
                "day": day_key,
                "half_hour": half_hour,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            })
        day_coverage[day_key] = {"expected": n_expected, "retrieved": n_got}

    if not slot_frames:
        raise ValueError(
            f"no IMERG granules retrieved for {d0.isoformat()}.."
            f"{d1.isoformat()} (run={key}; {len(skipped)} slots skipped: "
            f"{skipped[:3]}{'...' if len(skipped) > 3 else ''})")

    # Rebuild the canonical lon axis from the last frame's ordering.
    first_windows = imerg_index_windows(
        (lon_min, lat_min, lon_max, lat_max), grid_lats, grid_lons)
    # Derive the latitude vector from the same inclusive window slice used
    # for the frames above (never an independent mask).
    j0, j1 = first_windows[0][2], first_windows[0][3]
    sub_lats = grid_lats[j0:j1 + 1]
    canon_lons = np.concatenate(
        [grid_lons[i0:i1 + 1] for i0, i1, _j0, _j1 in first_windows])
    canon_lons = canon_lons[np.argsort(canon_lons, kind="stable")]

    if acc == "daily":
        times, values = accumulate_daily(slot_times, slot_frames)
        units = IMERG_UNITS_DAILY
    else:
        times = slot_times
        values = np.stack(slot_frames)
        units = IMERG_UNITS_NATIVE

    provenance = {
        "source": ("NASA GPM IMERG precipitation "
                   f"({IMERG_COLLECTIONS[key]} V{IMERG_VERSION}, {key} run)"),
        "collection": IMERG_COLLECTIONS[key],
        "version": IMERG_VERSION,
        "run": key,
        "run_latency": IMERG_LATENCY[key],
        "accumulate": acc,
        "units": units,
        "on_wire": ("precipitation is a mm/hr rate (gauge-adjusted merged "
                    "microwave-infrared estimate, formerly precipitationCal); "
                    "daily totals are sum(rate_i x 0.5 h) over each day's "
                    "retrieved half-hour slots"),
        "bbox": [lon_min, lat_min, lon_max, lat_max],
        "time_window": [d0.isoformat(), d1.isoformat()],
        "stride_days": stride_days,
        "grid": (f"{IMERG_RES}-degree global (lon -179.95..179.95, lat "
                 f"-89.95..89.95, ascending)"),
        "files": file_records,
        "n_files": len(file_records),
        "day_coverage": day_coverage,
        "skipped_slots": skipped,
        "retrieved_utc": retrieved_at,
        "authentication": ("Earthdata Login (free at "
                           "https://urs.earthdata.nasa.gov/users/new)"),
        "tool": f"survey-currents/{_tool_version()}",
    }
    return RainField(
        times=times, lats=sub_lats, lons=canon_lons, values=values,
        units=units, run=key, accumulate=acc,
        bbox=(lon_min, lat_min, lon_max, lat_max),
        provenance=provenance,
    )


# ---------------------------------------------------------------------------
# offline demo
# ---------------------------------------------------------------------------


def main_demo() -> None:
    """Print a small offline summary (no network). Mirrors sea_ice.main_demo."""
    f = RainField.synthetic()
    print(f"[imerg] synthetic field {f.shape}, units={f.units}, "
          f"t0={f.times[0].date()}, day-0 total={f.total(0):.1f} mm")


def _tool_version() -> str:
    from . import __version__
    return __version__


__all__ = [
    "CredentialsMissing",
    "RainField",
    "DateLike",
    "IMERG_BASE",
    "IMERG_VERSION",
    "IMERG_COLLECTIONS",
    "IMERG_RUN_INFIX",
    "IMERG_START",
    "IMERG_LATENCY",
    "IMERG_RES",
    "IMERG_UNITS_NATIVE",
    "IMERG_UNITS_DAILY",
    "IMERG_FILL",
    "IMERG_ACCUMULATIONS",
    "accumulate_daily",
    "fetch_imerg",
    "imerg_file_url",
    "imerg_granule_name",
    "imerg_index_windows",
    "imerg_sample_days",
    "normalize_imerg_accumulate",
    "normalize_imerg_run",
    "main_demo",
]

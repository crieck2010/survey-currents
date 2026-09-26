"""Global sea-surface temperature acquisition: NOAA OISST v2.1 + NASA JPL MUR v4.1.

Two complementary global SST sources (both verified live 2026-09-26;
see docs/DATA_SOURCES.md):

* :func:`fetch_oisst` — NOAA OISST v2.1 (``ncdcOisst21Agg`` griddap on
  NOAA CoastWatch ERDDAP): daily SST, 1981-present, 0.25-degree global
  grid, variable ``sst`` in degrees Celsius. The longitude axis runs
  0.125..359.875 (0-360 convention); requested bboxes in the usual
  -180..180 convention are converted, including antimeridian-crossing
  boxes. Keyless.
* :func:`fetch_mur` — NASA JPL MUR v4.1 (GHRSST L4, collection
  ``C1996881146-POCLOUD``) via the Earthdata OPeNDAP endpoint: daily
  SST, 2002-present, 0.01-degree (~1 km) global grid, variable
  ``analysed_sst`` (Kelvin on the wire, converted to Celsius here).
  Needs free Earthdata Login credentials (env vars or netrc); without
  them a :class:`CredentialsMissing` error explains the setup.

Both return :class:`SstField`, which follows the :class:`GlseaField`
conventions in ``glsea.py`` (and duck-types into what
``survey-viz``'s ``render_viz`` consumes: ``.times``/``.lats``/``.lons``
plus a 3-D ``.sst`` grid).

Downloads use stdlib ``urllib`` only. ``netCDF4`` is a lazy import, via
:func:`glsea._require_netcdf4`.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import (_coerce_date, _coerce_datetime_utc, _download_bytes,
                    _erddap_iso, _require_netcdf4)


# ---------------------------------------------------------------------------
# NOAA OISST v2.1 constants (verified live against CoastWatch ERDDAP 2026-09-26)
# ---------------------------------------------------------------------------

OISST_ERDDAP_BASE = "https://coastwatch.pfeg.noaa.gov/erddap"
OISST_DATASET = "ncdcOisst21Agg"
#: Grid axes per the dataset .das (verified live).
OISST_LON_MIN, OISST_LON_MAX = 0.125, 359.875   # 0-360 convention
OISST_LAT_MIN, OISST_LAT_MAX = -89.875, 89.875
OISST_RES = 0.25
OISST_FILL = -9.99
OISST_START = _dt.date(1981, 9, 1)
#: Daily timesteps are stamped at 12:00 UTC (per time_coverage_end).
OISST_TIME_OF_DAY = _dt.time(12, 0, 0)
#: Long time ranges make ERDDAP drop the connection; chunk requests into
#: windows no longer than this (the heatwaveR docs report trouble past
#: ~9 years — we stay well under that).
OISST_MAX_YEARS_PER_REQUEST = 5

# ---------------------------------------------------------------------------
# NASA JPL MUR v4.1 constants (verified via NASA CMR 2026-09-26)
# ---------------------------------------------------------------------------

#: CMR collection concept id for "GHRSST Level 4 MUR Global Foundation
#: Sea Surface Temperature Analysis (v4.1)".
MUR_COLLECTION_ID = "C1996881146-POCLOUD"
MUR_OPENDAP_BASE = (
    "https://opendap.earthdata.nasa.gov/collections/"
    f"{MUR_COLLECTION_ID}/granules"
)
#: GHRSST L4 variable name for the foundation SST analysis field.
MUR_VAR = "analysed_sst"
#: MUR granule titles embed the analysis time, always 09:00 UTC, e.g.
#: "20260925090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1".
MUR_GRANULE_TIME = _dt.time(9, 0, 0)
MUR_START = _dt.date(2002, 6, 1)

DateLike = Union[_dt.date, _dt.datetime, str]


class CredentialsMissing(RuntimeError):
    """Earthdata Login credentials are required but were not found.

    NASA's OPeNDAP endpoint requires a free Earthdata Login account.
    Provide credentials one of two ways:

    1. Environment variables (best for CI / servers)::

           export EARTHDATA_USERNAME="your_username"
           export EARTHDATA_PASSWORD="your_password"

    2. A ``~/.netrc`` entry (best for interactive use)::

           machine opendap.earthdata.nasa.gov
           login your_username
           password your_password

    Register for free at https://urs.earthdata.nasa.gov/users/new.
    The keyless NOAA OISST source (:func:`fetch_oisst`) keeps working
    without any credentials.
    """


# ---------------------------------------------------------------------------
# SstField — canonical global-SST field model
# ---------------------------------------------------------------------------


@dataclass
class SstField:
    """Time-indexed sea-surface-temperature grids (global sources).

    ``sst`` is ``(nt, ny, nx)`` in degrees Celsius as a numpy masked
    array (land / missing cells masked). ``lats``/``lons`` are 1-D
    coordinate vectors, increasing; longitudes are normalized to the
    -180..180 convention. Follows the :class:`GlseaField` conventions
    in ``glsea.py``.
    """

    sst: np.ma.MaskedArray               # (nt, ny, nx) sea-surface temperature, degC
    times: List[str]                     # ISO-8601 timestamps, one per step
    lats: np.ndarray                     # (ny,) degrees north, increasing
    lons: np.ndarray                     # (nx,) degrees east, -180..180, increasing
    crs: str = "EPSG:4326"
    source: str = "noaa-oisst/ncdcOisst21Agg"
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

    def select_time(self, index: int) -> "SstField":
        """Return the single-timestep field at ``index``."""
        return SstField(
            sst=self.sst[index:index + 1], times=[self.times[index]],
            lats=self.lats, lons=self.lons, crs=self.crs,
            source=self.source, provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "SstField":
        """Spatial subset to ``(min_lon, min_lat, max_lon, max_lat)``."""
        minx, miny, maxx, maxy = (float(x) for x in bbox)
        ix = np.where((self.lons >= minx) & (self.lons <= maxx))[0]
        iy = np.where((self.lats >= miny) & (self.lats <= maxy))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(
                f"bbox {list(bbox)} does not overlap field bounds {self.bounds}")
        return SstField(
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
    def from_dict(cls, data: Dict[str, Any]) -> "SstField":
        """Rebuild from :meth:`to_dict` (NaN -> masked). Raises on missing keys."""
        required = ("sst", "times", "lats", "lons")
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(f"SstField dict missing keys: {missing}")
        return cls(
            sst=np.ma.masked_invalid(np.asarray(data["sst"], dtype=float)),
            times=list(data["times"]),
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            crs=data.get("crs", "EPSG:4326"),
            source=data.get("source", "noaa-oisst/ncdcOisst21Agg"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "SstField":
        import json
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture (tests, demos, offline use) -----------------------

    @classmethod
    def synthetic(
        cls,
        nt: int = 4, ny: int = 6, nx: int = 8,
        lats: Sequence[float] = (-30.0, 30.0),
        lons: Sequence[float] = (-60.0, 60.0),
        start: str = "2020-01-01T12:00:00+00:00",
        step_days: int = 30,
        seed: int = 11,
        source: str = "synthetic",
    ) -> "SstField":
        """Deterministic synthetic global-SST field: warm equator, cool poles.

        Stdlib+numpy only. Used by the test suite and the offline demo.
        """
        rng = np.random.default_rng(seed)
        lat = np.linspace(lats[0], lats[1], ny)
        lon = np.linspace(lons[0], lons[1], nx)
        yy, _ = np.meshgrid(np.linspace(-1, 1, ny), np.linspace(-1, 1, nx),
                            indexing="ij")
        base = 27.0 - 18.0 * (yy ** 2)   # warm equator, cool poles, degC
        sst = np.ma.empty((nt, ny, nx))
        t0 = _dt.datetime.fromisoformat(start)
        times = []
        for k in range(nt):
            noise = rng.normal(0.0, 0.2, size=(ny, nx))
            sst[k] = np.ma.array(base + noise + 0.05 * k,
                                 mask=np.zeros((ny, nx), dtype=bool))
            times.append((t0 + _dt.timedelta(days=k * step_days)).isoformat())
        return cls(sst=sst, times=times, lats=lat, lons=lon,
                   source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# shared bbox validation (conventional -180..180 longitudes)
# ---------------------------------------------------------------------------


def validate_sst_bbox(bbox: Sequence[float]) -> Tuple[float, float, float, float]:
    """Validate ``(min_lon, min_lat, max_lon, max_lat)`` in the -180..180
    convention (ascending; latitudes within -90..90).

    Antimeridian-crossing bboxes such as (170, -10, -170, 10) are VALID —
    they are converted to a wrapped request window internally.
    """
    if len(bbox) != 4:
        raise ValueError(f"bbox must have 4 numbers, got {list(bbox)}")
    minx, miny, maxx, maxy = (float(x) for x in bbox)
    if not (-90.0 <= miny <= 90.0 and -90.0 <= maxy <= 90.0):
        raise ValueError(f"bbox latitudes must be within [-90, 90], got {list(bbox)}")
    if not (-180.0 <= minx <= 180.0 and -180.0 <= maxx <= 180.0):
        raise ValueError(f"bbox longitudes must be within [-180, 180], got {list(bbox)}")
    if miny >= maxy:
        raise ValueError(f"bbox latitudes must ascend, got {list(bbox)}")
    return (minx, miny, maxx, maxy)


def _validate_sst_dates(d0: _dt.datetime, d1: _dt.datetime,
                        earliest: _dt.date, source_name: str) -> None:
    if d0 > d1:
        raise ValueError(f"start {d0.date()} is after end {d1.date()}")
    if d0.date() < earliest:
        raise ValueError(
            f"{source_name} starts {earliest.isoformat()}; "
            f"start {d0.date().isoformat()} is before the record")
    if d1.date() > _dt.date.today() + _dt.timedelta(days=2):
        raise ValueError(
            f"end {d1.date().isoformat()} is in the future")


# ---------------------------------------------------------------------------
# NOAA OISST v2.1
# ---------------------------------------------------------------------------


def _lon_to_360(lon: float) -> float:
    return lon % 360.0


def oisst_lon_windows(bbox: Sequence[float]) -> List[Tuple[float, float]]:
    """Convert a -180..180 bbox to OISST 0-360 longitude request windows.

    Returns one window for the common case. Two cases need special care:

    * **Antimeridian-crossing bboxes** (e.g. ``(170, ., -170, .)``) wrap to
      a single 0-360 window (``(170.0, 190.0)``).
    * **Windows that cross 360° in 0-360 space** (e.g. ``(-170, ., 170, .)``
      becomes ``(190.0, 530.0)``) cannot be expressed in one ERDDAP query,
      so they are split into two windows: ``(190.0, 359.875)`` and
      ``(0.125, 170.0)``. :func:`fetch_oisst` issues both requests and
      concatenates along longitude, so no data is silently dropped.
    * **Full-globe bboxes** (span >= 359.9°) request the entire grid.

    All windows stay within the OISST grid floor/ceiling (0.125..359.875).
    """
    minx, _, maxx, _ = validate_sst_bbox(bbox)
    span = maxx - minx
    if span < 0:
        span += 360.0
    if span >= 359.9:
        return [(OISST_LON_MIN, OISST_LON_MAX)]
    lo, hi = _lon_to_360(minx), _lon_to_360(maxx)
    if hi <= lo:
        hi += 360.0
    if hi <= OISST_LON_MAX:
        return [(max(lo, OISST_LON_MIN), hi)]
    return [(max(lo, OISST_LON_MIN), OISST_LON_MAX),
            (OISST_LON_MIN, hi - 360.0)]


def oisst_sst_urls(bbox: Sequence[float], start: DateLike, end: DateLike,
                   stride_days: int = 30) -> List[str]:
    """Build the OISST griddap NetCDF URL(s) (exact bytes land in provenance).

    Usually a single URL; two URLs when the 0-360 longitude window crosses
    360° (see :func:`oisst_lon_windows`). Each query indexes the singleton
    ``zlev`` axis explicitly with ``[(0.0)]`` — ERDDAP 404s without it.
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_sst_dates(d0, d1, OISST_START, "OISST v2.1")
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")
    t0, t1 = _erddap_iso(d0), _erddap_iso(d1)
    step = f"{stride_days}"
    urls = []
    for lo360, hi360 in oisst_lon_windows(bbox):
        query = (
            f"sst[({t0}):{step}:({t1})][(0.0)][({miny}):({maxy})]"
            f"[({lo360}):({hi360})]"
        )
        urls.append(f"{OISST_ERDDAP_BASE}/griddap/{OISST_DATASET}.nc?{query}")
    return urls


def _chunk_date_range(d0: _dt.datetime, d1: _dt.datetime,
                      max_years: int) -> List[Tuple[_dt.datetime, _dt.datetime]]:
    """Split [d0, d1] into windows no longer than ``max_years`` years."""
    chunks: List[Tuple[_dt.datetime, _dt.datetime]] = []
    cursor = d0
    while cursor <= d1:
        try:
            next_end = cursor.replace(year=cursor.year + max_years)
        except ValueError:  # Feb 29 edge
            next_end = cursor.replace(year=cursor.year + max_years, day=28)
        chunk_end = min(d1, next_end - _dt.timedelta(seconds=1))
        chunks.append((cursor, chunk_end))
        cursor = chunk_end + _dt.timedelta(seconds=1)
    return chunks


import contextlib


@contextlib.contextmanager
def _open_nc_bytes(payload: bytes):
    """Open NetCDF ``payload`` bytes as a read-only dataset.

    netCDF4's in-memory open only handles NetCDF-4/HDF5 payloads —
    classic NetCDF-3 (what CoastWatch ERDDAP serves for OISST) fails
    there, so fall back to a temp file in that case. Yields the open
    dataset; callers use ``with _open_nc_bytes(payload) as ds:``.
    """
    nc = _require_netcdf4()
    try:
        ds = nc.Dataset("sst-global-in-memory", memory=payload, mode="r")
    except Exception:
        ds = None
    if ds is not None:
        try:
            yield ds
        finally:
            ds.close()
        return
    import os
    import tempfile
    fd, path = tempfile.mkstemp(suffix=".nc")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(payload)
        ds = nc.Dataset(path, mode="r")
        try:
            yield ds
        finally:
            ds.close()
    finally:
        os.unlink(path)


def _parse_oisst_bytes(payload: bytes) -> Tuple[np.ndarray, np.ndarray,
                                               List[str], np.ma.MaskedArray]:
    """Parse one OISST griddap NetCDF payload.

    Returns (lats, lons, times, sst). Lons are normalized to -180..180 and
    sorted increasing; the sst array has shape (nt, ny, nx) in degC with
    fill values masked.
    """
    nc = _require_netcdf4()
    with _open_nc_bytes(payload) as ds:
        if "sst" not in ds.variables:
            raise ValueError("OISST payload missing variable 'sst'")
        svar = ds.variables["sst"]
        if svar.dimensions != ("time", "zlev", "latitude", "longitude"):
            raise ValueError(
                f"unexpected OISST dimensions {svar.dimensions}; "
                "expected ('time', 'zlev', 'latitude', 'longitude')")
        lats = np.asarray(ds.variables["latitude"][:], dtype=float)
        lons360 = np.asarray(ds.variables["longitude"][:], dtype=float)
        # np.ma.asarray keeps netCDF4's automatic _FillValue mask; the extra
        # masked_values with tolerance catches fill stored as float32
        # (-9.99f != -9.99 in float64).
        raw = np.ma.asarray(svar[:, 0, :, :], dtype=float)  # drop singleton zlev
        time_var = ds.variables["time"]
        try:
            dts = nc.num2date(time_var[:], time_var.units)
        except Exception as exc:
            raise ValueError(f"could not decode OISST time axis: {exc}") from exc
    times = [_coerce_datetime_utc(d).isoformat() for d in dts]
    lons = ((lons360 + 180.0) % 360.0) - 180.0
    order = np.argsort(lons)
    lons = lons[order]
    sst = np.ma.masked_values(raw[:, :, order], OISST_FILL,
                              rtol=1e-5, atol=1e-3)
    sst = np.ma.masked_invalid(sst)
    return lats, lons, times, sst


def _concat_lon_parts(
        parts: List[Tuple[np.ndarray, np.ndarray, List[str],
                          np.ma.MaskedArray]]
) -> Tuple[np.ndarray, np.ndarray, List[str], np.ma.MaskedArray]:
    """Concatenate parsed OISST payloads along longitude, sorted increasing.

    Each part is ``(lats, lons, times, sst)`` from :func:`_parse_oisst_bytes`
    (lons already in -180..180, increasing within each part).
    """
    lats = parts[0][0]
    times = parts[0][2]
    for p in parts[1:]:
        if p[2] != times:
            raise RuntimeError(
                "OISST lon-window parts returned different time axes; "
                "cannot concatenate")
    lons = np.concatenate([p[1] for p in parts])
    sst = np.ma.concatenate([p[3] for p in parts], axis=2)
    order = np.argsort(lons)
    return lats, lons[order], times, sst[:, :, order]


def fetch_oisst(bbox: Sequence[float], start: DateLike, end: DateLike,
                stride_days: int = 30) -> SstField:
    """Fetch NOAA OISST v2.1 daily SST for ``bbox`` over [start, end].

    Args:
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180 degrees
            (ascending; antimeridian-crossing boxes wrap internally).
        start, end: dates/datetimes/ISO strings; sampled every
            ``stride_days`` along the time axis (same griddap stride
            semantics as the GLSEA adapter).
        stride_days: time-axis stride, >= 1.

    Returns:
        :class:`SstField` in degC with provenance (exact URLs, SHA-256,
        byte counts, retrieval time, dataset, bbox, time window, units).

    Raises:
        ValueError: invalid bbox / dates / stride.
        RuntimeError: network or NetCDF failures (with context).
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_sst_dates(d0, d1, OISST_START, "OISST v2.1")
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")

    chunks = _chunk_date_range(d0, d1, OISST_MAX_YEARS_PER_REQUEST)
    urls: List[str] = []
    time_parts: List[Tuple[np.ndarray, np.ndarray, List[str],
                           np.ma.MaskedArray]] = []
    n_payload_bytes = 0
    payload_for_hash = hashlib.sha256()
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()
    for c0, c1 in chunks:
        lon_parts = []
        for url in oisst_sst_urls((minx, miny, maxx, maxy), c0, c1,
                                  stride_days):
            urls.append(url)
            try:
                payload = _download_bytes(url)
            except Exception as exc:
                raise RuntimeError(
                    f"OISST request failed ({type(exc).__name__}: {exc}) — "
                    f"URL: {url}"
                ) from exc
            try:
                lon_parts.append(_parse_oisst_bytes(payload))
            except ValueError:
                raise
            except Exception as exc:
                raise RuntimeError(
                    f"OISST payload parse failed ({type(exc).__name__}: {exc}). "
                    "netCDF4 is required (pip install netCDF4)."
                ) from exc
            n_payload_bytes += len(payload)
            payload_for_hash.update(payload)
        time_parts.append(_concat_lon_parts(lon_parts))
    lats = time_parts[0][0]
    lons = time_parts[0][1]
    times = [t for part in time_parts for t in part[2]]
    sst = np.ma.concatenate([p[3] for p in time_parts], axis=0)
    return SstField(
        sst=sst, times=times, lats=lats, lons=lons,
        source=f"noaa-oisst/{OISST_DATASET}",
        provenance={
            "dataset": f"NOAA OISST v2.1 ({OISST_DATASET})",
            "urls": urls,
            "sha256": payload_for_hash.hexdigest(),
            "n_bytes": n_payload_bytes,
            "retrieved_at": retrieved_at,
            "bbox_requested": [minx, miny, maxx, maxy],
            "time_requested": [_erddap_iso(d0), _erddap_iso(d1)],
            "stride_days": stride_days,
            "time_chunks": len(chunks),
            "lon_windows": [list(w) for w in oisst_lon_windows(
                (minx, miny, maxx, maxy))],
            "units": "degree_C",
            "grid": f"{OISST_RES}-degree global (0-360 longitude, "
                    f"{OISST_LON_MIN}..{OISST_LON_MAX})",
        },
    )


# ---------------------------------------------------------------------------
# NASA JPL MUR v4.1 (Earthdata OPeNDAP)
# ---------------------------------------------------------------------------


def earthdata_credentials() -> Optional[Tuple[str, str]]:
    """Return ``(username, password)`` for Earthdata Login, or None.

    Checks ``EARTHDATA_USERNAME``/``EARTHDATA_PASSWORD`` first, then a
    ``~/.netrc`` entry for ``opendap.earthdata.nasa.gov``.
    """
    user = os.environ.get("EARTHDATA_USERNAME", "").strip()
    pw = os.environ.get("EARTHDATA_PASSWORD", "").strip()
    if user and pw:
        return (user, pw)
    try:
        import netrc
        hosts = netrc.netrc().hosts
    except Exception:
        return None
    for machine in ("opendap.earthdata.nasa.gov", "urs.earthdata.nasa.gov"):
        if machine in hosts:
            login, _acct, password = hosts[machine]
            if login and password:
                return (login, password)
    return None


def _earthdata_opener() -> urllib.request.OpenerDirector:
    """Build an authenticated urllib opener; raise :class:`CredentialsMissing`."""
    creds = earthdata_credentials()
    if creds is None:
        raise CredentialsMissing(
            "MUR SST needs Earthdata Login credentials, but none were found. "
            "Set EARTHDATA_USERNAME/EARTHDATA_PASSWORD or add a ~/.netrc entry "
            "for opendap.earthdata.nasa.gov (register free at "
            "https://urs.earthdata.nasa.gov/users/new). "
            "The keyless NOAA OISST source (fetch_oisst) needs no credentials.")
    user, pw = creds
    mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    mgr.add_password(None, MUR_OPENDAP_BASE, user, pw)
    auth = urllib.request.HTTPBasicAuthHandler(mgr)
    return urllib.request.build_opener(auth)


def mur_granule_title(day: DateLike) -> str:
    """Granule title for the MUR analysis stamped on ``day``.

    e.g. ``20260925090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1``
    (verified live via NASA CMR 2026-09-26).
    """
    d = _coerce_date(day)
    return (f"{d.strftime('%Y%m%d')}{MUR_GRANULE_TIME.strftime('%H%M%S')}"
            "-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1")


def mur_opendap_url(day: DateLike) -> str:
    """Base OPeNDAP service URL for one MUR granule (no constraint)."""
    return f"{MUR_OPENDAP_BASE}/{mur_granule_title(day)}"


#: Public NASA CMR granule search endpoint (keyless; discovery needs no
#: Earthdata credentials — only the data download does).
CMR_GRANULE_SEARCH = "https://cmr.earthdata.nasa.gov/search/granules.json"


def cmr_search_mur_granules(start: DateLike, end: DateLike,
                            page_size: int = 200) -> List[Dict[str, Any]]:
    """Discover real MUR v4.1 granules via the public NASA CMR API.

    Queries ``cmr.earthdata.nasa.gov`` (keyless) for granules of
    ``MUR_COLLECTION_ID`` overlapping ``[start, end]``. Each result is a
    dict with ``title`` (the authoritative granule name), ``time_start`` /
    ``time_end`` (ISO), ``opendap_url`` (the CMR-advertised OPeNDAP data
    link, if any), and ``links`` (the raw CMR link list).

    This is the authoritative discovery step: :func:`fetch_mur` uses the
    titles returned here rather than constructing them, so a rename or
    reprocessing of the MUR collection fails loudly instead of fetching
    the wrong granule.
    """
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    if page_size < 1:
        raise ValueError(f"page_size must be >= 1, got {page_size}")
    params = urllib.parse.urlencode({
        "collection_concept_id": MUR_COLLECTION_ID,
        "temporal": f"{d0.date().isoformat()}T00:00:00Z,"
                    f"{d1.date().isoformat()}T23:59:59Z",
        "page_size": page_size,
        "sort_key": "start_date",
    })
    url = f"{CMR_GRANULE_SEARCH}?{params}"
    payload = _download_bytes(url)
    try:
        data = json.loads(payload.decode("utf-8"))
    except ValueError as exc:
        raise RuntimeError(
            f"CMR granule search returned non-JSON ({exc}).") from exc
    granules: List[Dict[str, Any]] = []
    for entry in data.get("feed", {}).get("entry", []):
        links = entry.get("links") or []
        opendap_url = None
        for link in links:
            href = str(link.get("href", ""))
            if ("opendap" in href.lower()
                    or "OPENDAP" in str(link.get("rel", ""))):
                opendap_url = href
                break
        granules.append({
            "title": str(entry.get("title", "")),
            "time_start": entry.get("time_start"),
            "time_end": entry.get("time_end"),
            "opendap_url": opendap_url,
            "links": links,
        })
    return granules


def mur_match_granules(days: Sequence[_dt.date],
                       discovered: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Match each sampled ``day`` to one CMR-discovered granule.

    MUR granules are one daily analysis (09:00 UTC); the match is made on
    the granule title's leading ``YYYYMMDD`` stamp, falling back to the
    ``time_start`` date. Raises :class:`RuntimeError` if any day has no
    matching granule — an honest gap rather than a silently substituted
    granule.
    """
    matched: List[Dict[str, Any]] = []
    for day in days:
        stamp = day.strftime("%Y%m%d")
        hit = None
        for granule in discovered:
            if str(granule.get("title", "")).startswith(stamp):
                hit = granule
                break
        if hit is None:
            for granule in discovered:
                time_start = str(granule.get("time_start", ""))
                if time_start[:10] == day.isoformat():
                    hit = granule
                    break
        if hit is None:
            raise RuntimeError(
                f"CMR granule search returned no MUR granule for {day.isoformat()}; "
                "cannot fetch that date without inventing a granule name. "
                "Try a narrower date range or check the MUR collection status.")
        matched.append(hit)
    return matched


def _mur_sample_dates(d0: _dt.datetime, d1: _dt.datetime,
                      stride_days: int) -> List[_dt.date]:
    days: List[_dt.date] = []
    cursor = d0.date()
    while cursor <= d1.date():
        days.append(cursor)
        cursor += _dt.timedelta(days=stride_days)
    return days


def _opendap_get(opener: urllib.request.OpenerDirector, url: str,
                 context: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "survey-currents/0.3.0"})
    try:
        with opener.open(req, timeout=120) as resp:
            return resp.read()
    except urllib.error.HTTPError as exc:
        if exc.code == 401:
            raise CredentialsMissing(
                f"Earthdata Login rejected the credentials for {context} "
                "(HTTP 401). Check EARTHDATA_USERNAME/EARTHDATA_PASSWORD or "
                "the ~/.netrc entry for opendap.earthdata.nasa.gov.") from exc
        raise RuntimeError(
            f"MUR {context} failed (HTTP {exc.code}: {exc.reason}).") from exc
    except Exception as exc:
        raise RuntimeError(
            f"MUR {context} failed ({type(exc).__name__}: {exc}).") from exc


def _parse_mur_grid(das: str, dds: str) -> Dict[str, Tuple[int, float, float]]:
    """Recover the MUR lat/lon grid geometry from granule metadata.

    Axis sizes come from the DAP2 ``.dds`` (``Float64 lat[lat = N]``);
    coordinate ranges come from the ``.das`` ``actual_range`` attributes.
    Returns ``{"lat": (size, min, max), "lon": (size, min, max)}``.
    """
    out: Dict[str, Tuple[int, float, float]] = {}
    for name in ("lat", "lon"):
        m = re.search(rf"\b{name}\s*\[\s*{name}\s*=\s*(\d+)\s*\]", dds)
        if not m:
            raise ValueError(f"MUR .dds missing axis declaration for '{name}'")
        size = int(m.group(1))
        # actual_range lives inside the axis's container block in the .das:
        #   lat { String long_name "latitude"; Float64 actual_range -89.99, 89.99; }
        block = re.search(rf"(?s)\b{name}\s*\{{(.*?)\}}", das)
        if not block:
            raise ValueError(f"MUR .das missing attribute block for '{name}'")
        r = re.search(r"actual_range\s+\"?([^\n;\"]+)\"?", block.group(1))
        if not r:
            raise ValueError(f"MUR .das missing actual_range for '{name}'")
        lo_s, hi_s = [x.strip() for x in r.group(1).split(",")[:2]]
        out[name] = (size, float(lo_s), float(hi_s))
    return out


def _mur_grid_windows(grid: Dict[str, Tuple[int, float, float]],
                      bbox: Sequence[float]
                      ) -> List[Tuple[int, int, int, int]]:
    """Map a -180..180 bbox onto MUR grid index windows ``(i0, i1, j0, j1)``.

    Antimeridian-crossing bboxes (e.g. ``(170, ., -170, .)``) become two
    windows; everything else is a single window. MUR's axis is -180..180,
    so only clamping to the grid floor/ceiling is needed otherwise.
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    nlat, lat0, lat1 = grid["lat"]
    nlon, lon0, lon1 = grid["lon"]
    step_lat = (lat1 - lat0) / (nlat - 1)
    step_lon = (lon1 - lon0) / (nlon - 1)

    def idx_of(val: float, v0: float, step: float, n: int) -> int:
        return max(0, min(n - 1, int(round((val - v0) / step))))

    j0 = idx_of(miny, lat0, step_lat, nlat)
    j1 = idx_of(maxy, lat0, step_lat, nlat)
    if j0 > j1:
        j0, j1 = j1, j0

    spans: List[Tuple[float, float]] = []
    if maxx < minx:  # antimeridian crossing: split into two lon ranges
        spans = [(minx, lon1), (lon0, maxx)]
    else:
        spans = [(minx, maxx)]
    windows = []
    for lo, hi in spans:
        i0 = idx_of(lo, lon0, step_lon, nlon)
        i1 = idx_of(hi, lon0, step_lon, nlon)
        if i0 <= i1:
            windows.append((i0, i1, j0, j1))
    if not windows:
        raise ValueError(
            f"bbox {list(bbox)} does not overlap the MUR grid")
    return windows


def _mur_service_url(granule: Dict[str, Any]) -> Tuple[str, str]:
    """Base OPeNDAP URL for a CMR-discovered granule -> ``(url, source)``.

    Prefers the OPeNDAP link CMR advertises on the granule; otherwise
    falls back to the documented Earthdata OPeNDAP URL pattern for the
    collection (verified live via CMR 2026-09-26). ``source`` is
    ``"cmr-link"`` or ``"constructed"`` and is recorded in provenance.
    """
    link = granule.get("opendap_url")
    if link:
        return str(link).rstrip("/"), "cmr-link"
    return f"{MUR_OPENDAP_BASE}/{granule['title']}", "constructed"


def mur_subset_urls(day: DateLike, bbox: Sequence[float],
                    grid: Dict[str, Tuple[int, float, float]],
                    base_url: Optional[str] = None) -> List[str]:
    """Build the constrained OPeNDAP URL(s) for one MUR granule subset.

    ``grid`` comes from :func:`_parse_mur_grid` on the granule's metadata.
    The time axis is indexed ``[0:1:0]`` (one analysis per granule).
    Antimeridian-crossing bboxes yield two URLs (a single DAP2 constraint
    cannot wrap); :func:`fetch_mur` concatenates them along longitude.

    ``base_url`` (from :func:`_mur_service_url` on a CMR-discovered
    granule) overrides the constructed :func:`mur_opendap_url`; when
    omitted the legacy constructed URL is used (kept for tests).
    """
    base = base_url or mur_opendap_url(day)
    urls = []
    for i0, i1, j0, j1 in _mur_grid_windows(grid, bbox):
        constraint = f"{MUR_VAR}[0:1:0][{j0}:1:{j1}][{i0}:1:{i1}]"
        urls.append(f"{base}.nc?{constraint}")
    return urls


def _parse_mur_bytes(payload: bytes) -> Tuple[np.ndarray, np.ndarray,
                                              np.ma.MaskedArray]:
    """Parse one MUR subset NetCDF payload -> ``(lats, lons, sst[1, ny, nx])``.

    ``analysed_sst`` is Kelvin on the wire (per the GHRSST spec and the
    variable's ``units`` attribute, usually packed as scaled shorts);
    netCDF4 applies the packing and fill mask automatically, and the
    Kelvin -> Celsius conversion happens here.
    """
    nc = _require_netcdf4()
    with _open_nc_bytes(payload) as ds:
        if MUR_VAR not in ds.variables:
            raise ValueError(f"MUR payload missing variable '{MUR_VAR}'")
        var = ds.variables[MUR_VAR]
        units = str(getattr(var, "units", "")).strip().lower()
        data = np.ma.asarray(var[:], dtype=float)  # keeps auto fill-mask
        if data.ndim == 4:      # (time, lat, lon)
            data = data[0]
        elif data.ndim == 3 and data.shape[0] == 1:
            data = data[0]
        lats = np.asarray(ds.variables["lat"][:], dtype=float)
        lons = np.asarray(ds.variables["lon"][:], dtype=float)
    if units in ("kelvin", "k", "degrees_kelvin", "degree_kelvin",
                 "degreesk", "degree_k"):
        data = data - 273.15
    return lats, lons, np.ma.masked_invalid(data)[None, :, :]


def fetch_mur(bbox: Sequence[float], start: DateLike, end: DateLike,
              stride_days: int = 30) -> SstField:
    """Fetch NASA JPL MUR v4.1 SST for ``bbox`` over [start, end].

    Granules are discovered first through the public NASA CMR granule
    search API (keyless) — the OPeNDAP service URLs come from real
    granule titles, never from constructed names. Each sampled day is
    then subset through the Earthdata OPeNDAP endpoint; the spatial grid
    comes from the granule's own ``.das``.

    Args:
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180 degrees
            (ascending). MUR's axis is -180..180; no wrap conversion.
        start, end: dates/datetimes/ISO strings; one granule per
            ``stride_days``.
        stride_days: granule sampling stride, >= 1.

    Returns:
        :class:`SstField` in degC with provenance (CMR search record,
        discovered granule titles, service URLs, SHA-256 digests,
        retrieval time, dataset, units).

    Raises:
        CredentialsMissing: no Earthdata Login credentials found (or the
            server rejected them with HTTP 401).
        ValueError: invalid bbox / dates / stride.
        RuntimeError: CMR search failures, a sampled day with no
            discovered granule, or network/NetCDF failures (with context).
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_sst_dates(d0, d1, MUR_START, "MUR v4.1")
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")

    opener = _earthdata_opener()  # raises CredentialsMissing early
    days = _mur_sample_dates(d0, d1, stride_days)

    # Authoritative granule discovery via the public CMR API (keyless):
    # the OPeNDAP service URLs below come from real granule titles, never
    # from constructed names.
    discovered = cmr_search_mur_granules(d0, d1)
    matched = mur_match_granules(days, discovered)
    service = [_mur_service_url(g) for g in matched]
    base_url = service[0][0]

    # Grid geometry from the first granule's .das/.dds (stable across granules).
    das_text = _opendap_get(opener, f"{base_url}.das",
                            "metadata (.das)").decode("utf-8", errors="replace")
    dds_text = _opendap_get(opener, f"{base_url}.dds",
                            "metadata (.dds)").decode("utf-8", errors="replace")
    try:
        grid = _parse_mur_grid(das_text, dds_text)
    except ValueError as exc:
        raise RuntimeError(
            f"MUR metadata parse failed ({exc}); "
            f"granule {matched[0]['title']}"
        ) from exc

    lats_all = lons_all = None
    frames: List[np.ma.MaskedArray] = []
    times: List[str] = []
    urls: List[str] = []
    digests: List[str] = []
    titles: List[str] = []
    url_sources: List[str] = []
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()
    for day, granule, (granule_url, url_source) in zip(days, matched, service):
        title = granule["title"]
        titles.append(title)
        url_sources.append(url_source)
        # Frame timestamp: the analysis time embedded in the granule title
        # (YYYYMMDDHHMMSS, always 09:00 UTC for MUR) when it parses, else
        # the CMR-advertised time_start, else the documented 09:00 UTC.
        m = re.match(r"^(\d{14})", title)
        if m:
            stamp = _dt.datetime.strptime(
                m.group(1), "%Y%m%d%H%M%S").replace(
                    tzinfo=_dt.timezone.utc).isoformat()
        elif granule.get("time_start"):
            stamp = str(granule["time_start"])
        else:
            stamp = _dt.datetime.combine(
                _coerce_date(day), MUR_GRANULE_TIME,
                tzinfo=_dt.timezone.utc).isoformat()
        day_frames: List[np.ma.MaskedArray] = []
        day_lons: List[np.ndarray] = []
        for url in mur_subset_urls(day, (minx, miny, maxx, maxy), grid,
                                   base_url=granule_url):
            urls.append(url)
            payload = _opendap_get(opener, url, f"granule {title}")
            digests.append(hashlib.sha256(payload).hexdigest())
            try:
                lats, lons, sst = _parse_mur_bytes(payload)
            except ValueError:
                raise
            except Exception as exc:
                raise RuntimeError(
                    f"MUR granule {title} parse failed "
                    f"({type(exc).__name__}: {exc}). netCDF4 is required "
                    "(pip install netCDF4)."
                ) from exc
            if lats_all is None:
                lats_all = lats
            day_frames.append(sst)
            day_lons.append(lons)
        order = np.argsort(np.concatenate(day_lons))
        day_cube = np.ma.concatenate(day_frames, axis=2)[:, :, order]
        if lons_all is None:
            lons_all = np.concatenate(day_lons)[order]
        frames.append(day_cube)
        times.append(stamp)

    sst_cube = np.ma.concatenate(frames, axis=0)
    combined = hashlib.sha256()
    for d in digests:
        combined.update(d.encode())
    return SstField(
        sst=sst_cube, times=times, lats=lats_all, lons=lons_all,
        source="nasa-jpl-mur/MUR-JPL-L4-GLOB-v4.1",
        provenance={
            "dataset": "NASA JPL MUR v4.1 (GHRSST Level 4 MUR Global "
                       "Foundation SST Analysis)",
            "collection": MUR_COLLECTION_ID,
            "cmr_search": {
                "endpoint": CMR_GRANULE_SEARCH,
                "collection_concept_id": MUR_COLLECTION_ID,
                "temporal": [d0.date().isoformat(), d1.date().isoformat()],
                "n_results": len(discovered),
            },
            "granule_titles": titles,
            "granule_time_start": [g.get("time_start") for g in matched],
            "opendap_url_source": url_sources,
            "service_urls": urls,
            "granule_sha256": digests,
            "combined_sha256": combined.hexdigest(),
            "retrieved_at": retrieved_at,
            "bbox_requested": [minx, miny, maxx, maxy],
            "time_requested": [d0.date().isoformat(), d1.date().isoformat()],
            "stride_days": stride_days,
            "n_granules": len(days),
            "units": "degree_C (converted from Kelvin analysed_sst)",
            "grid": "~0.01-degree global MUR analysis grid",
        },
    )


# ---------------------------------------------------------------------------
# offline demo
# ---------------------------------------------------------------------------


def main_demo() -> None:
    """Print a small offline summary (no network). Mirrors glsea.main_demo."""
    for label, src in (("oisst", "noaa-oisst/ncdcOisst21Agg"),
                       ("mur", "nasa-jpl-mur/synthetic")):
        f = SstField.synthetic(source=src, nt=3)
        print(f"[{label}] synthetic {f.sst.shape} degC, "
              f"t0={f.times[0]}, mean={f.spatial_mean(0):.2f} degC")


__all__ = [
    "CredentialsMissing",
    "SstField",
    "DateLike",
    "OISST_ERDDAP_BASE",
    "OISST_DATASET",
    "OISST_LON_MIN",
    "OISST_LON_MAX",
    "OISST_LAT_MIN",
    "OISST_LAT_MAX",
    "OISST_START",
    "MUR_COLLECTION_ID",
    "MUR_OPENDAP_BASE",
    "MUR_VAR",
    "MUR_START",
    "CMR_GRANULE_SEARCH",
    "cmr_search_mur_granules",
    "mur_match_granules",
    "validate_sst_bbox",
    "oisst_lon_windows",
    "oisst_sst_urls",
    "fetch_oisst",
    "earthdata_credentials",
    "mur_granule_title",
    "mur_opendap_url",
    "mur_subset_urls",
    "fetch_mur",
    "main_demo",
]

"""NASA Black Marble night-lights acquisition via LAADS DAAC.

:func:`fetch_blackmarble` pulls daily moonlight-adjusted nighttime lights
(2012-01-19-present) from NASA's Level-1 and Atmosphere Archive &
Distribution System (LAADS DAAC):

* ``product="daily"`` — VNP46A2 V002, the gap-filled lunar BRDF-adjusted
  nighttime-lights product (the sane default, and the only product wired
  in this version).

The science dataset is ``Gap_Filled_DNB_BRDF_Corrected_NTL`` in
nW·cm⁻²·sr⁻¹ (DNB radiance, moonlight- and atmosphere-corrected,
cloud gaps filled from neighboring nights).

Access truth (verified live 2026-09-26 — see docs/DATA_SOURCES.md):

* Granule *discovery* is keyless: NASA CMR's granule search
  (``https://cmr.earthdata.nasa.gov/search/granules.json``) returns the
  exact download URLs for any (day, tile) — no login needed.
* Actual *downloads* require a free Earthdata Login account: both the
  LAADS on-prem archive (``ladsweb.modaps.eosdis.nasa.gov``) and the
  Earthdata Cloud copy
  (``data.laadsdaac.earthdatacloud.nasa.gov``) redirect unauthenticated
  requests to ``urs.earthdata.nasa.gov``. Without credentials a
  :class:`CredentialsMissing` error explains the exact setup (mirrors
  the OSCAR / MUR / IMERG adapters).
* Tile layout: VNP46A2 V002 keeps the ``hHHvVV`` tile naming but the
  tiles are 10°×10° *lat/lon* tiles (the V2 "15 arc-second linear lat/lon
  grid"), not sinusoidal: ``h = floor((lon + 180) / 10)``,
  ``v = floor((90 - lat_top) / 10)`` — verified against live CMR
  granule footprints (``h07v10`` = lon −110..−100, lat −20..−10).
  Each tile is natively 2400×2400 (15 arc-second); the adapter mosaics
  the tiles covering the bbox and resamples to the requested
  ``resolution`` (default 0.05°) with NaN-aware block averaging.
* CMR collection: ``C3365931269-LAADS`` (VNP46A2, version 2, record
  starts 2012-01-19).

Parsing the HDF5 payloads needs ``h5py`` (lazy import — the engine core
stays stdlib+numpy, like the netCDF4/rasterio precedents).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from dataclasses import field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date
from .sst_global import CredentialsMissing as _BaseCredentialsMissing
from .sst_global import earthdata_credentials, validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]


# ---------------------------------------------------------------------------
# Black Marble / LAADS constants (verified live 2026-09-26)
# ---------------------------------------------------------------------------

#: CMR collection concept id for VNP46A2 V002 (daily gap-filled NTL).
BM_COLLECTION_CONCEPT_ID = "C3365931269-LAADS"
#: Product short name / version.
BM_SHORT_NAME = "VNP46A2"
BM_VERSION = "2"
#: First day of the VNP46A2 record (per the CMR collection metadata).
BM_START = _dt.date(2012, 1, 19)
#: Primary science dataset (DNB radiance, nW/cm^2/sr).
BM_SDS = "Gap_Filled_DNB_BRDF_Corrected_NTL"
#: Fallback SDS when the gap-filled layer is absent.
BM_SDS_FALLBACK = "DNB_BRDF_Corrected_NTL"
#: On-wire radiance units.
BM_UNITS = "nW/cm^2/sr"
#: V2 tile span, degrees (10°×10° lat/lon tiles).
BM_TILE_DEG = 10.0
#: Native V2 grid spacing: 15 arc-seconds.
BM_NATIVE_RES = 15.0 / 3600.0
#: CMR granule search endpoint (keyless).
CMR_GRANULES_URL = "https://cmr.earthdata.nasa.gov/search/granules.json"
#: LAADS Earthdata Cloud data host (Earthdata Login required).
BM_CLOUD_HOST = "data.laadsdaac.earthdatacloud.nasa.gov"
#: LAADS on-prem archive host (Earthdata Login required).
BM_ARCHIVE_HOST = "ladsweb.modaps.eosdis.nasa.gov"

#: Filename tile token, e.g. ``h18v11`` in
#: ``VNP46A2.A2026262.h18v11.002.2026263100808.h5``.
_TILE_TOKEN_RE = re.compile(r"\.h(\d{2})v(\d{2})\.")


class CredentialsMissing(_BaseCredentialsMissing):
    """Earthdata Login credentials are required but were not found.

    Black Marble downloads need a free Earthdata Login account:

    1. Register (free) at https://urs.earthdata.nasa.gov/users/new.
    2. Provide credentials one of two ways:

       Environment variables (best for CI / servers)::

           export EARTHDATA_USERNAME="your_username"
           export EARTHDATA_PASSWORD=<redacted>

       A ``~/.netrc`` entry (best for interactive use)::

           machine urs.earthdata.nasa.gov
           login your_username
           password your_password

    Granule *discovery* (which tiles cover a bbox on a given day) stays
    keyless through NASA CMR — only the file downloads need the login.
    """


# ---------------------------------------------------------------------------
# Tile math — VNP46A2 V002 10°×10° lat/lon tiles
# ---------------------------------------------------------------------------


def normalize_bm_product(product: str) -> str:
    """Validate the product key; only ``"daily"`` (VNP46A2) is wired."""
    key = str(product).strip().lower()
    if key == "daily":
        return key
    raise ValueError(
        f"unknown Black Marble product {product!r}: only 'daily' "
        f"(VNP46A2 V{BM_VERSION}) is implemented — 'monthly' (VNP46A3) "
        "and 'annual' (VNP46A4) are documented future work.")


def tile_for_lon_lat(lon: float, lat: float) -> Tuple[int, int]:
    """Return the ``(h, v)`` VNP46A2 V002 tile containing ``(lon, lat)``.

    V2 tiles are 10°×10° lat/lon tiles: ``h = floor((lon + 180) / 10)``
    (0..35), ``v = floor((90 - lat_top) / 10)`` (0..17) where ``lat_top``
    is the tile's northern edge. Verified against live CMR footprints
    (``h07v10`` = lon −110..−100, lat −20..−10).
    """
    h = int(np.floor((float(lon) + 180.0) / BM_TILE_DEG))
    v = int(np.floor((90.0 - float(lat)) / BM_TILE_DEG))
    return max(0, min(35, h)), max(0, min(17, v))


def tile_bounds(h: int, v: int) -> Tuple[float, float, float, float]:
    """Return ``(lon_min, lat_min, lon_max, lat_max)`` for tile ``(h, v)``."""
    if not (0 <= h <= 35 and 0 <= v <= 17):
        raise ValueError(f"tile (h={h}, v={v}) out of range "
                         "(h: 0..35, v: 0..17)")
    lon_min = -180.0 + BM_TILE_DEG * h
    lat_max = 90.0 - BM_TILE_DEG * v
    return (lon_min, lat_max - BM_TILE_DEG, lon_min + BM_TILE_DEG, lat_max)


def tile_token(h: int, v: int) -> str:
    """Filename tile token, e.g. ``h18v11``."""
    return f"h{h:02d}v{v:02d}"


def tiles_for_bbox(bbox: Sequence[float]) -> List[Tuple[int, int]]:
    """All ``(h, v)`` tiles intersecting ``bbox`` (sorted, deterministic)."""
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    # Antimeridian-crossing boxes wrap: cover both longitude windows.
    lon_windows = ([(lon_min, 180.0), (-180.0, lon_max)]
                   if lon_min > lon_max else [(lon_min, lon_max)])
    tiles = set()
    for wmin, wmax in lon_windows:
        h0, _ = tile_for_lon_lat(wmin, lat_min + 1e-9)
        h1, _ = tile_for_lon_lat(wmax - 1e-9, lat_min + 1e-9)
        _, v0 = tile_for_lon_lat(wmin, lat_max - 1e-9)
        _, v1 = tile_for_lon_lat(wmin, lat_min + 1e-9)
        for h in range(min(h0, h1), max(h0, h1) + 1):
            for v in range(min(v0, v1), max(v0, v1) + 1):
                tiles.add((h, v))
    return sorted(tiles)


# ---------------------------------------------------------------------------
# CMR granule discovery (keyless) + authenticated download
# ---------------------------------------------------------------------------


def _validate_bm_dates(d0: _dt.date, d1: _dt.date) -> None:
    if d0 > d1:
        raise ValueError(f"start ({d0}) must not be after end ({d1})")
    if d0 < BM_START:
        raise ValueError(
            f"Black Marble VNP46A2 starts {BM_START.isoformat()}; "
            f"requested start {d0.isoformat()} predates the record.")


def bm_sample_days(d0: _dt.date, d1: _dt.date,
                   stride_days: int = 1) -> List[_dt.date]:
    """Sample dates from ``d0``..``d1`` inclusive, every ``stride_days``."""
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")
    n = (d1 - d0).days
    return [d0 + _dt.timedelta(days=k)
            for k in range(0, n + 1, stride_days)]


def cmr_granule_query(day: _dt.date, h: int, v: int,
                      page_size: int = 20) -> str:
    """CMR granule-search URL discovering the VNP46A2 granule for a tile.

    The search is keyless HTTPS; the response's ``GET DATA`` link is the
    authenticated download URL. The tile's 10°×10° footprint is passed as
    the search bounding box.
    """
    lon_min, lat_min, lon_max, lat_max = tile_bounds(h, v)
    params = {
        "short_name": BM_SHORT_NAME,
        "version": BM_VERSION,
        "temporal": (f"{day.isoformat()}T00:00:00Z,"
                     f"{day.isoformat()}T23:59:59Z"),
        "bounding_box": f"{lon_min},{lat_min},{lon_max},{lat_max}",
        "page_size": str(page_size),
    }
    return CMR_GRANULES_URL + "?" + urllib.parse.urlencode(params)


def _cmr_get_json(url: str, timeout: int = 60) -> Dict[str, Any]:
    """GET a CMR JSON endpoint (keyless); RuntimeError on failure."""
    import json
    req = urllib.request.Request(
        url, headers={"User-Agent": f"survey-currents/{_tool_version()}",
                      "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"CMR granule search failed (HTTP {exc.code}: {exc.reason}): "
            f"{url}") from exc
    except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
        raise RuntimeError(
            f"CMR granule search failed ({type(exc).__name__}: {exc}): "
            f"{url}") from exc


def discover_granule_url(day: _dt.date, h: int, v: int,
                         timeout: int = 60) -> Optional[str]:
    """Return the authenticated download URL for one (day, tile).

    Queries NASA CMR (keyless) and picks the ``GET DATA`` link whose
    filename carries this tile's ``hHHvVV`` token. Returns ``None``
    when CMR lists no granule for the tile/day (processing gap) —
    callers record the skip in provenance.
    """
    token = tile_token(h, v)
    doc = _cmr_get_json(cmr_granule_query(day, h, v), timeout=timeout)
    for entry in doc.get("feed", {}).get("entry", []):
        for link in entry.get("links", []):
            rel = str(link.get("rel", ""))
            href = str(link.get("href", ""))
            if "data#" not in rel or not href:
                continue
            if re.search(r"\." + re.escape(token) + r"\.", href):
                return href
    return None


def _blackmarble_opener() -> urllib.request.OpenerDirector:
    """Build an Earthdata-authenticated urllib opener.

    LAADS redirects unauthenticated data requests to URS
    (``urs.earthdata.nasa.gov``) for an OAuth handshake, so the password
    manager carries the data hosts plus URS, and a cookie processor
    keeps the session.

    Raises :class:`CredentialsMissing` when no credentials are found.
    """
    creds = earthdata_credentials()
    if creds is None:
        raise CredentialsMissing(
            "Black Marble downloads need Earthdata Login credentials, but "
            "none were found. Set EARTHDATA_USERNAME/EARTHDATA_PASSWORD or "
            "add a ~/.netrc entry for urs.earthdata.nasa.gov (register free "
            "at https://urs.earthdata.nasa.gov/users/new).")
    user, pw = creds
    mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    for host in (f"https://{BM_CLOUD_HOST}",
                 f"https://{BM_ARCHIVE_HOST}",
                 "https://urs.earthdata.nasa.gov"):
        mgr.add_password(None, host, user, pw)
    auth = urllib.request.HTTPBasicAuthHandler(mgr)
    return urllib.request.build_opener(
        auth, urllib.request.HTTPCookieProcessor())


def _bm_download(url: str, opener: urllib.request.OpenerDirector,
                 timeout: int = 300) -> bytes:
    """Download one Black Marble granule; 401 -> CredentialsMissing."""
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
                    "Earthdata Login rejected the credentials for LAADS "
                    "(HTTP 401). Check EARTHDATA_USERNAME/EARTHDATA_PASSWORD "
                    "or the ~/.netrc entry for urs.earthdata.nasa.gov."
                ) from exc
            if exc.code == 404:
                raise FileNotFoundError(
                    f"Black Marble granule not found (HTTP 404): {url}. "
                    "The day may sit inside the processing latency window "
                    "or a reprocessing gap.") from exc
            raise RuntimeError(
                f"Black Marble download failed (HTTP {exc.code}: "
                f"{exc.reason}): {url}") from exc
        except (urllib.error.URLError, http.client.RemoteDisconnected,
                http.client.IncompleteRead, TimeoutError,
                ConnectionError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Black Marble download failed after 3 attempts: "
                       f"{url}") from last


# ---------------------------------------------------------------------------
# HDF5 parsing (h5py is a lazy import)
# ---------------------------------------------------------------------------


def _require_h5py():
    try:
        import h5py  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Reading Black Marble HDF5 granules needs h5py: "
            "pip install h5py (or pip install \"survey-currents[blackmarble]\")."
        ) from exc
    import h5py
    return h5py


def _find_sds(h5py_file) -> Tuple[str, Any]:
    """Locate the NTL science dataset, searching the file tree.

    Prefers the gap-filled layer, falls back to the plain BRDF-corrected
    layer; raises ValueError listing the datasets found when neither
    exists (defensive against PGE layout changes).
    """
    found: Dict[str, Any] = {}

    def _visit(name, obj):
        if getattr(obj, "shape", None) and len(obj.shape) == 2:
            found[name.split("/")[-1]] = obj

    h5py_file.visititems(_visit)
    for key in (BM_SDS, BM_SDS_FALLBACK):
        if key in found:
            return key, found[key]
    raise ValueError(
        f"no NTL science dataset found (looked for {BM_SDS!r} and "
        f"{BM_SDS_FALLBACK!r}); datasets present: {sorted(found)}")


def _sds_to_radiance(ds) -> np.ndarray:
    """Read an SDS into float nW/cm²/sr, masking fill values to NaN.

    Applies ``scale_factor``/``add_offset`` attributes when present and
    masks ``_FillValue`` (checked under both common spellings).
    """
    raw = np.asarray(ds[...], dtype=float)
    attrs = dict(ds.attrs)
    fill = attrs.get("_FillValue", attrs.get("_Fillvalue", None))
    scale = float(attrs.get("scale_factor", 1.0))
    offset = float(attrs.get("add_offset", 0.0))
    out = raw * scale + offset
    if fill is not None:
        out = np.where(raw == float(fill), np.nan, out)
    return out


def _parse_blackmarble_bytes(payload: bytes, h: int, v: int,
                             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                        str]:
    """Parse one VNP46A2 tile granule.

    Returns ``(lats, lons, radiance, sds_name)`` where ``lats`` descends
    (row 0 = northern edge) and ``lons`` ascends; radiance is float
    nW/cm²/sr with missing cells as NaN. The tile spans exactly
    10°×10° (V2 linear lat/lon grid); pixel centers are derived from
    the array shape, so the parser tolerates PGE resolution changes.
    """
    import io
    h5py = _require_h5py()
    lon_min, lat_min, lon_max, lat_max = tile_bounds(h, v)
    with h5py.File(io.BytesIO(payload), "r") as fh:
        sds_name, ds = _find_sds(fh)
        radiance = _sds_to_radiance(ds)
    if radiance.ndim != 2:
        raise ValueError(
            f"expected a 2-D NTL grid, got shape {radiance.shape}")
    ny, nx = radiance.shape
    if ny < 2 or nx < 2:
        raise ValueError(
            f"NTL grid too small to geolocate: shape {radiance.shape}")
    lons = lon_min + (np.arange(nx) + 0.5) * (BM_TILE_DEG / nx)
    lats = lat_max - (np.arange(ny) + 0.5) * (BM_TILE_DEG / ny)
    return lats, lons, radiance, sds_name


def _block_average_to_grid(t_lats: np.ndarray, t_lons: np.ndarray,
                           data: np.ndarray,
                           lon_min: float, lat_max: float,
                           res: float, ny: int, nx: int,
                           ) -> Tuple[np.ndarray, np.ndarray]:
    """NaN-aware block average of one tile into the output accumulator.

    Each native pixel center is binned to its output cell
    (``row = floor((lat_max - lat) / res)``); ``(sum, count)`` arrays are
    returned for the caller to accumulate across tiles before dividing.
    Pixels falling outside the output grid are clipped away.
    """
    rows = np.floor((lat_max - t_lats[:, None]) / res).astype(np.int64)
    cols = np.floor((t_lons[None, :] - lon_min) / res).astype(np.int64)
    rows, cols = np.broadcast_arrays(rows, cols)
    valid = (~np.isnan(data)) & (rows >= 0) & (rows < ny) \
        & (cols >= 0) & (cols < nx)
    r = rows[valid]
    c = cols[valid]
    acc_sum = np.zeros((ny, nx), dtype=float)
    acc_count = np.zeros((ny, nx), dtype=float)
    np.add.at(acc_sum, (r, c), data[valid])
    np.add.at(acc_count, (r, c), 1.0)
    return acc_sum, acc_count


# ---------------------------------------------------------------------------
# LightsField — canonical night-lights field model
# ---------------------------------------------------------------------------


@dataclass
class LightsField:
    """Daily Black Marble night-lights grids.

    ``values`` is float DNB radiance (nW/cm²/sr) shaped (ntime, nlat,
    nlon) on a regular lat/lon grid; unlit / missing cells are NaN.
    ``times`` are UTC datetimes (one per day). Follows the
    :class:`IceField` / :class:`RainField` conventions.

    ``provenance`` records the exact file URLs, per-file SHA-256
    digests, skipped tiles/days, the retrieval timestamp, and the tool
    version — following the :mod:`currents.sst_global` conventions.
    """

    times: List[_dt.datetime]
    lats: np.ndarray
    lons: np.ndarray
    values: np.ndarray
    bbox: Tuple[float, float, float, float] = (-125.0, 25.0, -66.0, 49.0)
    resolution: float = 0.05
    source: str = "blackmarble"
    units: str = BM_UNITS
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        self.lats = np.asarray(self.lats, dtype=float).reshape(-1)
        self.lons = np.asarray(self.lons, dtype=float).reshape(-1)
        self.values = np.asarray(self.values, dtype=float).reshape(
            nt, self.lats.shape[0], self.lons.shape[0])
        self.bbox = validate_sst_bbox(self.bbox)
        if not self.resolution > 0:
            raise ValueError(
                f"LightsField.resolution must be > 0, got {self.resolution}")

    def __len__(self) -> int:
        return len(self.times)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    @property
    def time_range(self) -> Tuple[Optional[_dt.datetime],
                                  Optional[_dt.datetime]]:
        """(earliest, latest) timestep, or (None, None) when empty."""
        if not self.times:
            return None, None
        return min(self.times), max(self.times)

    @property
    def shape(self) -> Tuple[int, int, int]:
        """(ntime, nlat, nlon)."""
        return self.values.shape  # type: ignore[return-value]

    @property
    def total_radiance(self) -> np.ndarray:
        """Per-timestep NaN-aware spatial sum of DNB radiance."""
        return np.nansum(self.values, axis=(1, 2))

    # -- filters ---------------------------------------------------------

    def select_time(self, start: DateLike, end: DateLike) -> "LightsField":
        """Timesteps with ``start <= date <= end`` (inclusive)."""
        d0 = _coerce_date(start)
        d1 = _coerce_date(end)
        keep = [i for i, t in enumerate(self.times) if d0 <= t.date() <= d1]
        return LightsField(
            times=[self.times[i] for i in keep],
            lats=self.lats, lons=self.lons,
            values=self.values[keep],
            bbox=self.bbox, resolution=self.resolution,
            source=self.source, units=self.units,
            provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "LightsField":
        """Spatial subset to ``bbox`` (must lie inside the field grid)."""
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        iy = np.flatnonzero((self.lats >= lat_min) & (self.lats <= lat_max))
        ix = np.flatnonzero((self.lons >= lon_min) & (self.lons <= lon_max))
        if len(iy) == 0 or len(ix) == 0:
            raise ValueError(
                f"bbox {list(bbox)} does not overlap field lats/lons")
        return LightsField(
            times=list(self.times),
            lats=self.lats[iy], lons=self.lons[ix],
            values=self.values[:, iy][:, :, ix],
            bbox=(lon_min, lat_min, lon_max, lat_max),
            resolution=self.resolution,
            source=self.source, units=self.units,
            provenance=dict(self.provenance),
        )

    # -- serialization ----------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dict (values as nested lists, NaN -> None)."""
        return {
            "times": [t.isoformat() for t in self.times],
            "lats": [float(x) for x in self.lats],
            "lons": [float(x) for x in self.lons],
            "values": np.where(np.isnan(self.values), None,
                               self.values).tolist(),
            "bbox": list(self.bbox),
            "resolution": self.resolution,
            "source": self.source,
            "units": self.units,
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "LightsField":
        """Rebuild from :meth:`to_dict` output."""
        import math
        times = [_dt.datetime.fromisoformat(t) for t in data["times"]]
        values = np.array(
            [[[math.nan if v is None else float(v) for v in row]
              for row in frame]
             for frame in data["values"]],
            dtype=float)
        return cls(
            times=times,
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            values=values,
            bbox=tuple(data["bbox"]),  # type: ignore[arg-type]
            resolution=float(data.get("resolution", 0.05)),
            source=str(data.get("source", "blackmarble")),
            units=str(data.get("units", BM_UNITS)),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        """Write JSON to ``path``; returns ``path``."""
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "LightsField":
        """Read a field written by :meth:`to_json`."""
        import json
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture --------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-125.0, 25.0, -66.0, 49.0),
                  start: DateLike = "2024-01-01",
                  end: DateLike = "2024-01-05",
                  resolution: float = 1.0,
                  seed: int = 7,
                  source: str = "synthetic") -> "LightsField":
        """Deterministic offline fixture: city-like light blobs on darkness.

        A few Gaussian "urban" blobs (bright cores, dim sprawl) over a
        near-zero background with slight day-to-day flicker — enough
        structure to exercise renderers and mosaicking without a
        network. Fully offline and deterministic for ``seed``.
        """
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        d0 = _coerce_date(start)
        d1 = _coerce_date(end)
        days = [(d0 + _dt.timedelta(days=k))
                for k in range((d1 - d0).days + 1)]
        rng = np.random.default_rng(seed)
        lons = np.arange(lon_min + resolution / 2, lon_max, resolution)
        lats = np.arange(lat_max - resolution / 2, lat_min, -resolution)
        # Deterministic city centers inside the bbox.
        n_cities = 4
        cx = rng.uniform(lon_min, lon_max, n_cities)
        cy = rng.uniform(lat_min, lat_max, n_cities)
        amp = rng.uniform(20.0, 120.0, n_cities)
        spread = rng.uniform(0.4, 1.2, n_cities)
        lon2, lat2 = np.meshgrid(lons, lats)
        base = np.zeros_like(lon2)
        for x, y, a, s in zip(cx, cy, amp, spread):
            base += a * np.exp(-(((lon2 - x) ** 2 + (lat2 - y) ** 2)
                                 / (2 * s * s)))
        frames = []
        for k in range(len(days)):
            flicker = rng.normal(1.0, 0.03, size=base.shape)
            frame = np.clip(base * flicker, 0.0, None)
            frame[frame < 0.5] = np.nan  # unlit background
            frames.append(frame)
        return cls(
            times=[_dt.datetime(d.year, d.month, d.day,
                                tzinfo=_dt.timezone.utc) for d in days],
            lats=lats, lons=lons, values=np.asarray(frames),
            bbox=(lon_min, lat_min, lon_max, lat_max),
            resolution=resolution, source=source, units=BM_UNITS,
            provenance={"synthetic": True, "seed": seed,
                        "generator": "LightsField.synthetic"},
        )


# ---------------------------------------------------------------------------
# fetch_blackmarble
# ---------------------------------------------------------------------------


def fetch_blackmarble(bbox: Sequence[float], start: DateLike, end: DateLike,
                      product: str = "daily",
                      stride_days: int = 1,
                      resolution: float = 0.05,
                      timeout: int = 300) -> LightsField:
    """Fetch NASA Black Marble daily night lights for ``bbox``.

    Args:
        bbox: (lon_min, lat_min, lon_max, lat_max) in -180..180 degrees;
            antimeridian-crossing boxes wrap.
        start/end: inclusive date range (dates, datetimes, or ISO
            strings); must be on/after 2012-01-19.
        product: ``"daily"`` (VNP46A2 V002 — the only wired product).
        stride_days: keep every Nth day (default 1).
        resolution: output grid spacing in degrees (default 0.05);
            native 15 arc-second tiles are NaN-aware block-averaged.
        timeout: per-request HTTP timeout in seconds.

    Returns:
        A :class:`LightsField` of daily DNB radiance (nW/cm²/sr) with
        per-file provenance (exact URLs, SHA-256 digests, tile/day,
        CMR discovery queries, skipped tiles/days, retrieval time).
        Tiles/days CMR lists no granule for are skipped with a
        provenance note; days with no tiles at all are dropped.

    Raises:
        CredentialsMissing: no Earthdata Login credentials found.
        ImportError: ``h5py`` is not installed.
        ValueError: bad bbox / dates / product / resolution, dates
            before 2012-01-19, or no granules retrieved at all.
    """
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    prod = normalize_bm_product(product)
    d0, d1 = _coerce_date(start), _coerce_date(end)
    _validate_bm_dates(d0, d1)
    if not resolution > 0:
        raise ValueError(f"resolution must be > 0, got {resolution}")
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")

    _require_h5py()  # fail fast with the actionable ImportError
    opener = _blackmarble_opener()  # raises CredentialsMissing early

    tiles = tiles_for_bbox((lon_min, lat_min, lon_max, lat_max))
    days = bm_sample_days(d0, d1, stride_days)
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()

    # Output grid (pixel centers), lats descending.
    out_lons = np.arange(lon_min + resolution / 2, lon_max, resolution)
    out_lats = np.arange(lat_max - resolution / 2, lat_min, -resolution)
    if out_lons.size == 0 or out_lats.size == 0:
        raise ValueError(
            f"bbox {list(bbox)} is smaller than resolution {resolution}")
    ny, nx = out_lats.size, out_lons.size

    frames: List[np.ndarray] = []
    frame_times: List[_dt.datetime] = []
    file_records: List[Dict[str, Any]] = []
    skipped: List[str] = []
    day_tiles: Dict[str, Dict[str, int]] = {}

    for day in days:
        day_key = day.isoformat()
        acc_sum = np.zeros((ny, nx), dtype=float)
        acc_count = np.zeros((ny, nx), dtype=float)
        n_expected = len(tiles)
        n_got = 0
        for h, v in tiles:
            url = discover_granule_url(day, h, v)
            if url is None:
                skipped.append(f"{day_key} tile {tile_token(h, v)}: "
                               "no CMR granule")
                continue
            try:
                payload = _bm_download(url, opener, timeout=timeout)
            except FileNotFoundError as exc:
                skipped.append(f"{day_key} tile {tile_token(h, v)}: 404")
                continue
            except Exception as exc:
                raise RuntimeError(
                    f"Black Marble download failed for {url} "
                    f"({type(exc).__name__}: {exc})") from exc
            try:
                t_lats, t_lons, radiance, sds_name = _parse_blackmarble_bytes(
                    payload, h, v)
            except ValueError as exc:
                skipped.append(f"{day_key} tile {tile_token(h, v)}: {exc}")
                continue
            part_sum, part_count = _block_average_to_grid(
                t_lats, t_lons, radiance,
                lon_min, lat_max, resolution, ny, nx)
            acc_sum += part_sum
            acc_count += part_count
            n_got += 1
            file_records.append({
                "url": url,
                "day": day_key,
                "tile": tile_token(h, v),
                "sds": sds_name,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            })
        day_tiles[day_key] = {"expected": n_expected, "retrieved": n_got}
        if n_got == 0:
            skipped.append(f"{day_key}: no tiles retrieved; day dropped")
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            frame = np.where(acc_count > 0, acc_sum / acc_count, np.nan)
        frames.append(frame)
        frame_times.append(_dt.datetime(day.year, day.month, day.day,
                                        tzinfo=_dt.timezone.utc))

    if not frames:
        raise ValueError(
            f"no Black Marble granules retrieved for {d0.isoformat()}.."
            f"{d1.isoformat()} ({len(skipped)} tiles/days skipped: "
            f"{skipped[:3]}{'...' if len(skipped) > 3 else ''})")

    provenance: Dict[str, Any] = {
        "product": BM_SHORT_NAME,
        "version": BM_VERSION,
        "collection_concept_id": BM_COLLECTION_CONCEPT_ID,
        "sds": BM_SDS,
        "units": BM_UNITS,
        "resolution_deg": resolution,
        "tiles": [tile_token(h, v) for h, v in tiles],
        "n_files": len(file_records),
        "files": file_records,
        "day_tiles": day_tiles,
        "skipped": skipped,
        "cmr_endpoint": CMR_GRANULES_URL,
        "retrieved_at": retrieved_at,
        "tool": f"survey-currents/{_tool_version()}",
    }
    return LightsField(
        times=frame_times, lats=out_lats, lons=out_lons,
        values=np.asarray(frames),
        bbox=(lon_min, lat_min, lon_max, lat_max),
        resolution=resolution, source="blackmarble", units=BM_UNITS,
        provenance=provenance,
    )


def main_demo() -> None:
    """Print a small deterministic synthetic field summary (offline)."""
    field = LightsField.synthetic()
    t0, t1 = field.time_range
    print(f"synthetic LightsField: {len(field)} days "
          f"{t0.date()}..{t1.date()} shape={field.shape} "
          f"units={field.units}")


def _tool_version() -> str:
    try:
        from . import __version__
        return str(__version__)
    except Exception:
        return "0.0.0"

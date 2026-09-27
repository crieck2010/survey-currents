"""NSIDC Sea Ice Index (G02135 v4.0) acquisition: concentration GeoTIFFs -> IceField.

Daily passive-microwave sea-ice concentration, November 1978-present, from
the NOAA/NSIDC Sea Ice Index (G02135 v4.0). This module:

* :func:`fetch_nsidc_sic` — downloads the daily *concentration* GeoTIFFs
  from NSIDC's public HTTPS archive (verified live 2026-09-26; see
  docs/DATA_SOURCES.md), reprojects the native 25 km polar-stereographic
  grids to a regular lat/lon grid (nearest-neighbor, documented below),
  and returns an :class:`IceField`.
* :class:`IceField` — the canonical sea-ice concentration model
  (``times``/``lats``/``lons``/``values`` on a regular lat/lon grid),
  with provenance, JSON round-trip, ``select_time``/``select_bbox``,
  and a deterministic ``synthetic()`` fixture. The dict form is the
  exact shape ``survey-viz``'s renderer consumes, so ice maps render
  with zero renderer changes.

Access truth (verified live 2026-09-26 — no Earthdata login needed):

* The G02135 archive at ``https://noaadata.apps.nsidc.org/NOAA/G02135/``
  serves daily concentration GeoTIFFs over plain keyless HTTPS:
  ``{north,south}/daily/geotiff/YYYY/MM_Mon/{N,S}_YYYYMMDD_concentration_v4.0.tif``.
* GeoTIFF value encoding (G02135 v4.0 user guide, Table 5): unsigned
  16-bit, concentration scaled x10 (divide by 10 -> percent 0-100);
  2510 = Arctic pole hole, 2530 = coast line, 2540 = land, 2550 = missing
  data, 0 = open ocean. Flag cells become NaN in :class:`IceField`.
* Native grids: NSIDC polar stereographic, Hughes 1980 ellipsoid,
  25 km — North (EPSG:3411): 304x448, South (EPSG:3412): 316x332.
  Reprojection is nearest-neighbor: each target lat/lon cell takes the
  value of the nearest native cell (documented choice — it preserves the
  15%-threshold ice edge exactly, which bilinear interpolation would
  smear).

Temporal notes: the record starts 1978-11-01; the SMMR portion
(1978-10-26 to 1987-08-20) is every other day; there are no data
1987-12-03 to 1988-01-13 (satellite problems). Missing days are skipped
with a provenance note. Values 1-150 (1-15%) are statistically
irrelevant per the user guide (passive-microwave uncertainty below
15%); they are kept as-is and documented rather than silently dropped,
but the guide's 15% cutoff defines the ice *extent* edge.

Downloads use stdlib ``urllib`` only (via :func:`glsea._download_bytes`)
and the GeoTIFFs are uncompressed single-band 16-bit, so decoding needs
nothing but ``numpy`` — the keyless NSIDC path has zero extra
dependencies. (The optional ``raster`` extra is not required here.)
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import struct
from dataclasses import dataclass, field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date, _download_bytes

DateLike = Union[_dt.date, _dt.datetime, str]


# ---------------------------------------------------------------------------
# NSIDC G02135 v4.0 constants (verified live 2026-09-26)
# ---------------------------------------------------------------------------

#: Public HTTPS archive root (no account, no key).
NSIDC_BASE = "https://noaadata.apps.nsidc.org"
#: G02135 product path under the archive root.
NSIDC_G02135_PATH = "/NOAA/G02135"
#: Product version in file names.
NSIDC_VERSION = "v4.0"
#: First day of the Sea Ice Index record.
NSIDC_START = _dt.date(1978, 11, 1)
#: Concentration GeoTIFFs are scaled x10 -> divide by 10 for percent.
NSIDC_CONC_SCALE = 10.0
#: Flag values (G02135 v4.0 user guide, Table 5) -> NaN in IceField.
NSIDC_FLAGS = frozenset({2510, 2530, 2540, 2550})
#: Human labels for the flags (docs/provenance).
NSIDC_FLAG_LABELS = {
    2510: "arctic pole hole",
    2530: "coast line",
    2540: "land",
    2550: "missing data",
}
#: Hughes 1980 ellipsoid used by the NSIDC polar stereographic grids.
HUGHES_A = 6378273.0
HUGHES_E2 = 0.006693883
#: Nominal native cell size (m).
NSIDC_CELL_M = 25000.0

#: Native grid definitions. x0/y0_top are the map coordinates (m) of the
#: upper-left corner of pixel (0, 0); dx is the cell size; lon0/lat_ts
#: are the projection center/standard parallel. Verified against the
#: GeoTIFF georeferencing tags of a live file (2026-02-01, north):
#: pixel scale (25000, 25000), tiepoint (-3850000, 5850000), and the
#: pole hole (flag 2510) lands on pixel (154, 234) as predicted.
NSIDC_GRIDS: Dict[str, Dict[str, Any]] = {
    "north": {
        "epsg": 3411, "nx": 304, "ny": 448,
        "x0": -3850000.0, "y0_top": 5850000.0, "dx": NSIDC_CELL_M,
        "lon0": -45.0, "lat_ts": 70.0, "prefix": "N",
    },
    "south": {
        "epsg": 3412, "nx": 316, "ny": 332,
        "x0": -3950000.0, "y0_top": 4350000.0, "dx": NSIDC_CELL_M,
        "lon0": 0.0, "lat_ts": -70.0, "prefix": "S",
    },
}

#: Ice-domain latitude limits (user guide §4.2.1): the grids cover
#: north of 30.98°N / south of 39.23°S.
NSIDC_NORTH_LAT_MIN = 30.98
NSIDC_SOUTH_LAT_MAX = -39.23

_MONTH_ABBR = ("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


# ---------------------------------------------------------------------------
# Pure, offline-testable helpers
# ---------------------------------------------------------------------------

def validate_ice_bbox(bbox: Sequence[float]) -> Tuple[float, float, float, float]:
    """Validate a (lon_min, lat_min, lon_max, lat_max) bbox.

    Antimeridian-crossing boxes (lon_min > lon_max) are allowed — the
    polar grids wrap the pole, so Arctic requests routinely span -180..180.
    """
    try:
        lon_min, lat_min, lon_max, lat_max = (float(x) for x in bbox)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"bbox must be 4 numbers, got {bbox!r}") from exc
    if not (lat_min < lat_max and -90.0 <= lat_min <= 90.0
            and -90.0 <= lat_max <= 90.0):
        raise ValueError(
            f"bbox latitudes must satisfy -90 <= lat_min < lat_max <= 90, "
            f"got {bbox!r}")
    if not (-180.0 <= lon_min <= 180.0 and -180.0 <= lon_max <= 180.0):
        raise ValueError(
            f"bbox longitudes must be within [-180, 180], got {bbox!r}")
    if lon_min == lon_max:
        raise ValueError(f"bbox has zero longitude span: {bbox!r}")
    return lon_min, lat_min, lon_max, lat_max


def pick_hemisphere(bbox: Sequence[float],
                    hemisphere: str = "auto") -> str:
    """Pick ``"north"``/``"south"`` for ``bbox``.

    ``"auto"`` (default) chooses from the bbox latitudes: entirely
    non-negative -> north, entirely non-positive -> south. Bboxes
    straddling the equator need an explicit hemisphere (sea ice never
    does, so this is almost certainly a caller bug).
    """
    hemi = str(hemisphere).strip().lower()
    if hemi in ("north", "south"):
        return hemi
    if hemi != "auto":
        raise ValueError(
            f"hemisphere must be 'north', 'south', or 'auto', got "
            f"{hemisphere!r}")
    _, lat_min, _, lat_max = validate_ice_bbox(bbox)
    if lat_min >= 0.0:
        return "north"
    if lat_max <= 0.0:
        return "south"
    raise ValueError(
        f"bbox {tuple(bbox)!r} straddles the equator: pass "
        "hemisphere='north' or 'south' explicitly")


def nsidc_conc_url(day: _dt.date, hemisphere: str = "north") -> str:
    """Build the G02135 daily concentration GeoTIFF URL for ``day``.

    Pure and offline-testable. Verified live 2026-09-26, e.g.
    ``.../north/daily/geotiff/2026/02_Feb/N_20260201_concentration_v4.0.tif``.
    """
    hemi = str(hemisphere).strip().lower()
    if hemi not in NSIDC_GRIDS:
        raise ValueError(
            f"hemisphere must be one of {sorted(NSIDC_GRIDS)}, got "
            f"{hemisphere!r}")
    prefix = NSIDC_GRIDS[hemi]["prefix"]
    month_dir = f"{day.month:02d}_{_MONTH_ABBR[day.month - 1]}"
    name = (f"{prefix}_{day.year:04d}{day.month:02d}{day.day:02d}"
            f"_concentration_{NSIDC_VERSION}.tif")
    return (f"{NSIDC_BASE}{NSIDC_G02135_PATH}/{hemi}/daily/geotiff/"
            f"{day.year:04d}/{month_dir}/{name}")


def decode_concentration(raw: np.ndarray) -> np.ndarray:
    """Decode a raw uint16 concentration array to percent (0-100).

    Raw values are scaled x10 (803 -> 80.3%); flag values 2510/2530/
    2540/2550 (pole hole, coast, land, missing) become NaN, as do any
    out-of-range values above 1000 that are not recognized flags.
    """
    arr = np.asarray(raw)
    out = np.full(arr.shape, np.nan, dtype=float)
    valid = (arr <= 1000) & ~np.isin(arr, tuple(NSIDC_FLAGS))
    out[valid] = arr[valid].astype(float) / NSIDC_CONC_SCALE
    return out


def ps_forward(lat_deg: float, lon_deg: float,
               hemisphere: str = "north") -> Tuple[float, float]:
    """Forward NSIDC polar-stereographic projection: (lat, lon) -> (x, y) m.

    Ellipsoidal (Hughes 1980), true scale at |70°| (EPSG:3411/3412
    variant B), north-up map coordinates. Used to locate native-grid
    pixels for target lat/lon cells during nearest-neighbor
    reprojection. Verified against a live file: the pole maps to (0, 0),
    i.e. pixel (154, 234) on the north grid, where the pole-hole flag
    sits.
    """
    hemi = str(hemisphere).strip().lower()
    if hemi not in NSIDC_GRIDS:
        raise ValueError(f"unknown hemisphere {hemisphere!r}")
    g = NSIDC_GRIDS[hemi]
    e = math.sqrt(HUGHES_E2)
    phi = math.radians(abs(float(lat_deg)))
    lam = math.radians(float(lon_deg))
    lam0 = math.radians(g["lon0"])
    phi1 = math.radians(abs(g["lat_ts"]))

    def _t(p: float) -> float:
        s = math.sin(p)
        return (math.tan(math.pi / 4 - p / 2)
                / ((1 - e * s) / (1 + e * s)) ** (e / 2))

    m1 = math.cos(phi1) / math.sqrt(1 - HUGHES_E2 * math.sin(phi1) ** 2)
    rho = HUGHES_A * m1 * _t(phi) / _t(phi1)
    dlam = lam - lam0
    x = rho * math.sin(dlam)
    y = (-rho * math.cos(dlam)) if hemi == "north" else (rho * math.cos(dlam))
    return x, y


# ---------------------------------------------------------------------------
# Minimal stdlib GeoTIFF reader (uncompressed 16-bit single-band)
# ---------------------------------------------------------------------------

_TIFF_MAGIC = b"II*\x00"


def read_concentration_geotiff(payload: bytes) -> Tuple[np.ndarray, Dict[str, float]]:
    """Read one G02135 concentration GeoTIFF payload.

    The files are uncompressed single-band 16-bit little-endian TIFFs,
    so decoding needs nothing but the stdlib + numpy. Returns
    ``(raw, grid)`` where ``raw`` is the uint16 array shaped (ny, nx)
    and ``grid`` carries ``x0``/``y0_top``/``dx``/``nx``/``ny`` read
    from the GeoTIFF's own georeferencing tags (ModelPixelScaleTag /
    ModelTiepointTag), falling back to :data:`NSIDC_GRIDS` when the
    tags are absent.

    Raises:
        ValueError: not a TIFF, or not the expected uncompressed
            16-bit single-band layout.
    """
    if payload[:4] != _TIFF_MAGIC:
        raise ValueError(
            "not a little-endian TIFF payload "
            f"(magic={payload[:4]!r}); the NSIDC URL may have returned "
            "an error page")
    (ifd_off,) = struct.unpack("<I", payload[4:8])
    (n_tags,) = struct.unpack("<H", payload[ifd_off:ifd_off + 2])
    tags: Dict[int, Tuple[int, int, int]] = {}
    for i in range(n_tags):
        entry = payload[ifd_off + 2 + i * 12: ifd_off + 2 + (i + 1) * 12]
        tag, typ, cnt, val = struct.unpack("<HHI4s", entry)
        tags[tag] = (typ, cnt, struct.unpack("<I", val)[0])

    def _tag(tag: int, name: str) -> Tuple[int, int, int]:
        if tag not in tags:
            raise ValueError(f"TIFF is missing {name} (tag {tag})")
        return tags[tag]

    width = _tag(256, "ImageWidth")[2]
    height = _tag(257, "ImageLength")[2]
    if _tag(258, "BitsPerSample")[2] != 16:
        raise ValueError("expected 16-bit samples")
    if _tag(259, "Compression")[2] != 1:
        raise ValueError("expected uncompressed TIFF (Compression=1)")
    if _tag(277, "SamplesPerPixel")[2] != 1:
        raise ValueError("expected single-band TIFF")
    rows_per_strip = _tag(278, "RowsPerStrip")[2]
    typ, cnt, off = _tag(273, "StripOffsets")
    strip_offsets = struct.unpack(f"<{cnt}{'I' if typ == 4 else 'H'}",
                                  payload[off:off + cnt * (4 if typ == 4 else 2)])
    typ, cnt, off = _tag(279, "StripByteCounts")
    strip_counts = struct.unpack(f"<{cnt}{'I' if typ == 4 else 'H'}",
                                  payload[off:off + cnt * (4 if typ == 4 else 2)])

    buf = bytearray()
    rows_left = height
    for s_off, s_cnt in zip(strip_offsets, strip_counts):
        buf.extend(payload[s_off:s_off + s_cnt])
        rows_left -= rows_per_strip
        if rows_left <= 0:
            break
    raw = np.frombuffer(bytes(buf[:width * height * 2]),
                        dtype="<u2").reshape(height, width)

    # Georeferencing from the file itself; fall back to the documented
    # grid constants if the tags are missing.
    grid: Dict[str, float] = {"nx": float(width), "ny": float(height)}
    try:
        typ, cnt, off = _tag(33550, "ModelPixelScaleTag")
        dx, dy, _ = struct.unpack("<3d", payload[off:off + 24])
        typ, cnt, off = _tag(33922, "ModelTiepointTag")
        _, _, _, x0, y0, _ = struct.unpack("<6d", payload[off:off + 48])
        grid.update({"x0": x0, "y0_top": y0, "dx": dx, "dy": dy})
    except ValueError:
        grid.update({"x0": None, "y0_top": None, "dx": None, "dy": None})  # type: ignore[dict-item]
    return raw, grid


def resolve_native_grid(grid: Dict[str, float],
                        hemisphere: str) -> Dict[str, float]:
    """Fill missing georeferencing from :data:`NSIDC_GRIDS`.

    ``grid`` comes from :func:`read_concentration_geotiff`; when the
    file lacked georeferencing tags, the documented constants for
    ``hemisphere`` are used (and the pixel dimensions must match).
    """
    g = NSIDC_GRIDS[str(hemisphere).strip().lower()]
    out = dict(grid)
    if out.get("x0") is None:
        if int(out["nx"]) != g["nx"] or int(out["ny"]) != g["ny"]:
            raise ValueError(
                f"TIFF is {int(out['nx'])}x{int(out['ny'])} but the "
                f"{hemisphere} grid is {g['nx']}x{g['ny']} and carries no "
                "georeferencing tags")
        out.update({"x0": g["x0"], "y0_top": g["y0_top"],
                    "dx": g["dx"], "dy": g["dx"]})
    return out


def reproject_to_latlon(raw: np.ndarray, native: Dict[str, float],
                        lats: np.ndarray, lons: np.ndarray,
                        hemisphere: str) -> np.ndarray:
    """Nearest-neighbor reproject of a decoded grid onto a lat/lon mesh.

    ``raw`` is the uint16 native array (ny, nx); ``native`` its grid
    dict from :func:`resolve_native_grid`. ``lats``/``lons`` are 1-D
    target axes (cell centers). Each target cell takes the decoded
    value of the nearest native pixel (documented choice — it preserves
    the 15%-threshold ice edge instead of smearing it); cells with no
    native pixel (outside the grid, or flag/NaN source cells) are NaN.
    Returns float percent (nt-less) shaped (len(lats), len(lons)).
    """
    hemi = str(hemisphere).strip().lower()
    decoded = decode_concentration(np.asarray(raw))
    ny_t, nx_t = len(lats), len(lons)
    out = np.full((ny_t, nx_t), np.nan, dtype=float)

    lat_a = np.asarray(lats, dtype=float)
    lon_a = np.asarray(lons, dtype=float)
    LAT, LON = np.meshgrid(lat_a, lon_a, indexing="ij")

    # Vectorized forward projection.
    e = math.sqrt(HUGHES_E2)
    g = NSIDC_GRIDS[hemi]
    phi = np.radians(np.abs(LAT))
    dlam = np.radians(LON - g["lon0"])
    phi1 = math.radians(abs(g["lat_ts"]))
    sp = np.sin(phi)
    t = (np.tan(math.pi / 4 - phi / 2)
         / ((1 - e * sp) / (1 + e * sp)) ** (e / 2))
    s1 = math.sin(phi1)
    t1 = (math.tan(math.pi / 4 - phi1 / 2)
          / ((1 - e * s1) / (1 + e * s1)) ** (e / 2))
    m1 = math.cos(phi1) / math.sqrt(1 - HUGHES_E2 * s1 ** 2)
    rho = HUGHES_A * m1 * t / t1
    xs = rho * np.sin(dlam)
    ys = np.where(hemi == "north",
                  -rho * np.cos(dlam), rho * np.cos(dlam))

    x0 = float(native["x0"])
    y0 = float(native["y0_top"])
    dx = float(native["dx"])
    dy = float(native.get("dy") or dx)
    cols = np.rint((xs - x0) / dx).astype(int)
    rows = np.rint((y0 - ys) / dy).astype(int)
    nx_n, ny_n = int(native["nx"]), int(native["ny"])
    ok = (cols >= 0) & (cols < nx_n) & (rows >= 0) & (rows < ny_n)
    # The polar stereographic formulas are symmetric in |lat|: without
    # this mask, a southern-hemisphere target cell would sample the
    # northern grid's ring of the same |latitude|. Each grid only covers
    # its own hemisphere's ice domain.
    ok = ok & (LAT >= 0.0) if hemi == "north" else ok & (LAT <= 0.0)
    src = np.full((ny_t, nx_t), np.nan, dtype=float)
    src[ok] = decoded[rows[ok], cols[ok]]
    out[:, :] = src
    return out


def target_grid(bbox: Sequence[float],
                resolution: float = 0.25) -> Tuple[np.ndarray, np.ndarray]:
    """Regular lat/lon cell-center axes spanning ``bbox``.

    Antimeridian-crossing bboxes (lon_min > lon_max) produce a grid on
    the caller's wrapped axis (lons may exceed 180°), matching the
    convention in :mod:`currents.fires`.
    """
    if resolution <= 0:
        raise ValueError(f"resolution must be > 0, got {resolution}")
    lon_min, lat_min, lon_max, lat_max = validate_ice_bbox(bbox)
    span = (lon_max - lon_min) if lon_min < lon_max else (lon_max - lon_min) % 360.0
    nx = max(1, int(math.ceil(span / resolution)))
    ny = max(1, int(math.ceil((lat_max - lat_min) / resolution)))
    lons = lon_min + (np.arange(nx) + 0.5) * resolution
    lats = lat_min + (np.arange(ny) + 0.5) * resolution
    return lats, lons


# ---------------------------------------------------------------------------
# IceField — canonical sea-ice concentration model
# ---------------------------------------------------------------------------

@dataclass
class IceField:
    """Daily sea-ice concentration from NSIDC G02135 v4.0.

    ``values`` is float percent (0-100) shaped (ntime, nlat, nlon) on a
    regular lat/lon grid; land, coast, pole-hole, and missing-data cells
    are NaN. ``times`` are UTC datetimes (one per day), ``hemisphere``
    is ``"north"`` or ``"south"``.

    ``provenance`` records the exact file URLs, per-file SHA-256
    digests, skipped days, the retrieval timestamp, and the tool
    version — following the :mod:`currents.sst_global` conventions.
    """

    times: List[_dt.datetime]
    lats: np.ndarray
    lons: np.ndarray
    values: np.ndarray
    hemisphere: str = "north"
    bbox: Tuple[float, float, float, float] = (-180.0, 66.0, 180.0, 90.0)
    source: str = "nsidc"
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        self.lats = np.asarray(self.lats, dtype=float).reshape(-1)
        self.lons = np.asarray(self.lons, dtype=float).reshape(-1)
        self.values = np.asarray(self.values, dtype=float).reshape(
            nt, self.lats.shape[0], self.lons.shape[0])
        self.hemisphere = str(self.hemisphere).strip().lower()
        if self.hemisphere not in NSIDC_GRIDS:
            raise ValueError(
                f"IceField.hemisphere: {self.hemisphere!r} not in "
                f"{sorted(NSIDC_GRIDS)}")
        self.bbox = validate_ice_bbox(self.bbox)

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

    # -- filters ---------------------------------------------------------

    def select_time(self, start: DateLike, end: DateLike) -> "IceField":
        """Timesteps with ``start <= date <= end`` (inclusive)."""
        d0 = _coerce_date(start)
        d1 = _coerce_date(end)
        keep = [i for i, t in enumerate(self.times) if d0 <= t.date() <= d1]
        return IceField(
            times=[self.times[i] for i in keep],
            lats=self.lats, lons=self.lons,
            values=self.values[keep],
            hemisphere=self.hemisphere, bbox=self.bbox,
            source=self.source, provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "IceField":
        """Spatial subset to ``bbox`` (must lie inside the field grid)."""
        lon_min, lat_min, lon_max, lat_max = validate_ice_bbox(bbox)
        iy = np.flatnonzero((self.lats >= lat_min) & (self.lats <= lat_max))
        ix = np.flatnonzero((self.lons >= lon_min) & (self.lons <= lon_max))
        if len(iy) == 0 or len(ix) == 0:
            raise ValueError(f"bbox {tuple(bbox)!r} has no overlap with the field grid")
        return IceField(
            times=list(self.times),
            lats=self.lats[iy], lons=self.lons[ix],
            values=self.values[:, iy[:, None], ix],
            hemisphere=self.hemisphere,
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
            "hemisphere": self.hemisphere,
            "bbox": list(self.bbox),
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "IceField":
        return cls(
            times=[_dt.datetime.fromisoformat(t) for t in data["times"]],
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            values=np.asarray(data["values"], dtype=float),
            hemisphere=data.get("hemisphere", "north"),
            bbox=tuple(data["bbox"]),
            source=data.get("source", "nsidc"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        """Write the field as JSON; return ``path``."""
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "IceField":
        """Read a field written by :meth:`to_json`."""
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture -------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-180.0, 66.0, 180.0, 90.0),
                  start: DateLike = "2024-02-01", end: DateLike = "2024-02-05",
                  resolution: float = 1.0, seed: int = 7,
                  hemisphere: str = "auto",
                  source: str = "synthetic") -> "IceField":
        """Deterministic synthetic concentration (offline tests / demos).

        Concentration falls off smoothly from 100% at the pole to 0%
        below an ice edge whose latitude wobbles with longitude and
        time, so frames show a plausible retreating/advancing edge.
        """
        rng = np.random.default_rng(seed)
        d0, d1 = _coerce_date(start), _coerce_date(end)
        if d1 < d0:
            raise ValueError(f"end {d1} is before start {d0}")
        hemi = pick_hemisphere(bbox, hemisphere)
        lats, lons = target_grid(bbox, resolution)
        ndays = (d1 - d0).days + 1
        times = [_dt.datetime(d0.year, d0.month, d0.day,
                              tzinfo=_dt.timezone.utc)
                 + _dt.timedelta(days=i) for i in range(ndays)]
        LAT, LON = np.meshgrid(lats, lons, indexing="ij")
        pole = 90.0 if hemi == "north" else -90.0
        sign = 1.0 if hemi == "north" else -1.0
        values = np.empty((ndays, len(lats), len(lons)), dtype=float)
        for i in range(ndays):
            edge = 72.0 + 3.0 * np.sin(np.radians(LON * 3.0)) + 0.8 * i
            dist = sign * (pole - LAT)  # 0 at the pole, grows equatorward
            edge_dist = sign * (pole - edge)
            conc = 100.0 * np.clip(1.0 - (dist - edge_dist * 0.55)
                                   / (edge_dist * 0.45 + 1e-9), 0.0, 1.0)
            conc += rng.normal(0.0, 2.0, conc.shape)
            values[i] = np.clip(conc, 0.0, 100.0)
        # A land-like NaN cap so NaN handling is exercised.
        values[:, :1, :] = np.nan
        return cls(times=times, lats=lats, lons=lons, values=values,
                   hemisphere=hemi, bbox=validate_ice_bbox(bbox),
                   source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# fetch_nsidc_sic
# ---------------------------------------------------------------------------

def fetch_nsidc_sic(bbox: Sequence[float], start: DateLike, end: DateLike,
                    hemisphere: str = "auto", stride_days: int = 1,
                    resolution: float = 0.25,
                    timeout: int = 300) -> IceField:
    """Fetch NSIDC G02135 v4.0 daily sea-ice concentration.

    Args:
        bbox: (lon_min, lat_min, lon_max, lat_max); antimeridian-crossing
            boxes are allowed.
        start/end: inclusive date range (dates, datetimes, or ISO strings).
        hemisphere: ``"north"``/``"south"``/``"auto"`` (default; picked
            from the bbox latitudes).
        stride_days: keep every Nth day (default 1).
        resolution: target lat/lon grid resolution in degrees (default
            0.25); the 25 km native grid is nearest-neighbor reprojected
            onto it.
        timeout: per-request HTTP timeout in seconds.

    Returns:
        An :class:`IceField` with per-file provenance. Days whose file
        is missing (SMMR every-other-day era, the 1987-12-03 to
        1988-01-13 outage) are skipped and listed in
        ``provenance["skipped_days"]``.

    Raises:
        ValueError: bad bbox / dates / hemisphere, or every requested
            day missing.
    """
    lon_min, lat_min, lon_max, lat_max = validate_ice_bbox(bbox)
    d0, d1 = _coerce_date(start), _coerce_date(end)
    if d1 < d0:
        raise ValueError(f"end {d1} is before start {d0}")
    if d0 < NSIDC_START:
        raise ValueError(
            f"NSIDC G02135 coverage starts {NSIDC_START} (requested {d0})")
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")
    hemi = pick_hemisphere((lon_min, lat_min, lon_max, lat_max), hemisphere)

    lats, lons = target_grid((lon_min, lat_min, lon_max, lat_max), resolution)
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()

    days: List[_dt.date] = []
    cur = d0
    while cur <= d1:
        days.append(cur)
        cur += _dt.timedelta(days=stride_days)

    frames: List[np.ndarray] = []
    times: List[_dt.datetime] = []
    file_records: List[Dict[str, Any]] = []
    skipped: List[str] = []
    native_grid: Optional[Dict[str, float]] = None

    for day in days:
        url = nsidc_conc_url(day, hemi)
        try:
            payload = _download_bytes(url, timeout=timeout)
        except Exception as exc:  # 404s for SMMR-era off days, outages
            skipped.append(f"{day.isoformat()}: {type(exc).__name__}")
            continue
        try:
            raw, grid = read_concentration_geotiff(payload)
        except ValueError as exc:
            skipped.append(f"{day.isoformat()}: {exc}")
            continue
        if native_grid is None:
            native_grid = resolve_native_grid(grid, hemi)
        frame = reproject_to_latlon(raw, native_grid, lats, lons, hemi)
        frames.append(frame)
        times.append(_dt.datetime(day.year, day.month, day.day,
                                  tzinfo=_dt.timezone.utc))
        file_records.append({
            "url": url,
            "date": day.isoformat(),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "bytes": len(payload),
        })

    if not frames:
        raise ValueError(
            f"no NSIDC sea-ice files retrieved for "
            f"{d0.isoformat()}..{d1.isoformat()} ({len(skipped)} days "
            f"skipped: {skipped[:3]}{'...' if len(skipped) > 3 else ''})")

    provenance = {
        "source": "NSIDC G02135 v4.0 Sea Ice Index (daily concentration)",
        "archive": f"{NSIDC_BASE}{NSIDC_G02135_PATH}",
        "hemisphere": hemi,
        "bbox": [lon_min, lat_min, lon_max, lat_max],
        "time_window": [d0.isoformat(), d1.isoformat()],
        "stride_days": stride_days,
        "resolution_deg": resolution,
        "reprojection": ("nearest-neighbor from the native 25 km NSIDC "
                         "polar stereographic grid (EPSG "
                         f"{NSIDC_GRIDS[hemi]['epsg']})"),
        "value_encoding": ("uint16 scaled x10 -> percent 0-100; flags "
                           "2510/2530/2540/2550 (pole hole/coast/land/"
                           "missing) -> NaN"),
        "files": file_records,
        "n_files": len(file_records),
        "skipped_days": skipped,
        "retrieved_utc": retrieved_at,
        "authentication": "none (public keyless HTTPS)",
        "tool": f"survey-currents/{_tool_version()}",
    }
    return IceField(
        times=times, lats=lats, lons=lons,
        values=np.stack(frames), hemisphere=hemi,
        bbox=(lon_min, lat_min, lon_max, lat_max),
        provenance=provenance,
    )


def _tool_version() -> str:
    from . import __version__
    return __version__

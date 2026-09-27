"""GEBCO bathymetry/topography + Natural Earth vector basemaps (stdlib-only).

Two basemap families, both keyless HTTPS, both cached locally on first
use so repeated calls never re-download:

* **GEBCO 2024** (``fetch_gebco``) — the 15-arc-second global terrain
  model for ocean and land (elevation in meters, positive up, WGS 84 /
  mean sea level). Distributed as eight 90°×90° GeoTIFF tiles inside one
  ~4.3 GB zip on CEDA (``dap.ceda.ac.uk``). Downloading the whole zip
  for a regional subset would be wasteful, so the adapter extracts only
  the intersecting tile(s) with HTTP Range requests: it reads the zip's
  central directory from the tail, locates the tile entry, and pulls just
  that entry's byte range (~0.5 GB per tile, cached forever after). The
  cached tile is then window-read with a minimal stdlib GeoTIFF reader
  (tiled deflate int16, or stripped uncompressed int16 as shipped by
  the official GEBCO 2024 tiles) and block-averaged to the requested
  ``resolution``.
* **Natural Earth vectors** (``fetch_naturalearth``) — coastline,
  country, land, and lake vectors as GeoJSON FeatureCollections, parsed
  from the canonical Natural Earth S3 zips with a minimal stdlib
  shapefile reader (no fiona/geopandas needed).

Access truth (verified live 2026-09-27 — see docs/DATA_SOURCES.md):

* ``https://www.bodc.ac.uk/data/open_download/gebco/gebco_2024/geotiff/``
  301-redirects (keyless, no login) to
  ``https://dap.ceda.ac.uk/bodc/gebco/global/gebco_2024/ice_surface_elevation/geotiff/gebco_2024_geotiff.zip``
  (4,257,768,030 bytes, ``Accept-Ranges: bytes``). The zip holds the
  eight tiles named
  ``gebco_2024_n{north}_s{south}_w{west}_e{east}.tif`` plus the grid
  documentation PDFs. GEBCO 2025 exists (released Aug 2025) but is not
  on the keyless open_download path yet, so this adapter pins 2024.
* Natural Earth zips (``https://naturalearth.s3.amazonaws.com/...``)
  return HTTP 200 anonymously, e.g.
  ``110m_physical/ne_110m_coastline.zip`` and
  ``110m_cultural/ne_110m_admin_0_countries.zip``.

The GEBCO tiles are 21600×21600 int16 (15 arc-second, pixel-center
registered per the GEBCO grid documentation). ``fetch_gebco`` refuses
native-resolution reads whose pixel window would exceed 8192×8192
(request a coarser ``resolution`` instead) so a careless call cannot
try to materialize a continent at 15 arc-seconds.

The engine core is stdlib+numpy (the GeoTIFF and shapefile readers are
hand-rolled with ``struct``/``zlib``/``zipfile``); heavy work is lazy.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import io
import math
import os
import re
import struct
import urllib.error
import urllib.request
import zipfile
import zlib
from dataclasses import dataclass
from dataclasses import field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .sst_global import validate_sst_bbox

__all__ = [
    "GEBCO_VERSION",
    "GEBCO_ZIP_URL",
    "GEBCO_NATIVE_ARCSEC",
    "TopoField",
    "fetch_gebco",
    "synthetic_topography",
    "fetch_naturalearth",
    "synthetic_naturalearth",
    "gebco_tiles_for_bbox",
    "basemap_cache_dir",
    "NATURAL_EARTH_SCALES",
    "NATURAL_EARTH_LAYERS",
]

#: Pinned GEBCO grid release (2025 is not on the keyless path yet).
GEBCO_VERSION = "2024"

#: BODC open_download directory URL (301 -> CEDA zip, keyless).
GEBCO_ZIP_URL = (
    "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2024/geotiff/"
)

#: Native GEBCO grid spacing, arc-seconds.
GEBCO_NATIVE_ARCSEC = 15.0
GEBCO_NATIVE_DEG = GEBCO_NATIVE_ARCSEC / 3600.0

#: Safety cap: native-resolution window reads may not exceed this many
#: pixels on a side (~67 M cells).
GEBCO_MAX_NATIVE_WINDOW = 8192

#: Natural Earth S3 root (keyless, verified live 2026-09-27).
NE_BASE_URL = "https://naturalearth.s3.amazonaws.com"

#: Supported Natural Earth scales and their vector layers. Each value is
#: (subdirectory, zip stem) with ``{scale}`` interpolated.
NE_LAYER_FILES: Dict[str, Tuple[str, str]] = {
    "coastline": ("{scale}_physical",
                  "ne_{scale}_coastline.zip"),
    "countries": ("{scale}_cultural",
                  "ne_{scale}_admin_0_countries.zip"),
    "land": ("{scale}_physical",
             "ne_{scale}_land.zip"),
    "lakes": ("{scale}_physical",
              "ne_{scale}_lakes.zip"),
}
NATURAL_EARTH_SCALES = ("110m", "50m", "10m")
NATURAL_EARTH_LAYERS = tuple(NE_LAYER_FILES)

#: The single structural timestamp carried by TopoField. GEBCO is a
#: static compilation, not a time series; the field model (and the viz
#: renderer) need a ``times`` axis, so it carries one element. Renderers
#: should label it "static compilation", never a data date.
GEBCO_STATIC_TIME = _dt.datetime(2024, 1, 1, tzinfo=_dt.timezone.utc)

_TILE_NAME_RE = re.compile(
    r"gebco_2024_n(?P<north>-?\d+(?:\.\d+)?)"
    r"_s(?P<south>-?\d+(?:\.\d+)?)"
    r"_w(?P<west>-?\d+(?:\.\d+)?)"
    r"_e(?P<east>-?\d+(?:\.\d+)?)\.tif$"
)


# ---------------------------------------------------------------------------
# Cache discipline
# ---------------------------------------------------------------------------

def basemap_cache_dir() -> str:
    """User cache root for basemap downloads.

    ``$SURVEY_CURRENTS_CACHE/basemap`` when set, else
    ``~/.cache/survey-currents/basemap``. Created on demand.
    """
    root = os.environ.get("SURVEY_CURRENTS_CACHE") or os.path.join(
        os.path.expanduser("~"), ".cache", "survey-currents")
    path = os.path.join(root, "basemap")
    os.makedirs(path, exist_ok=True)
    return path


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_atomic(path: str, data: bytes) -> None:
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, path)


def _sidecar_ok(path: str) -> bool:
    """True when ``path`` exists and matches its ``.sha256`` sidecar."""
    sidecar = path + ".sha256"
    if not (os.path.isfile(path) and os.path.isfile(sidecar)):
        return False
    try:
        with open(sidecar, "r", encoding="utf-8") as fh:
            expected = fh.read().strip().split()[0]
    except OSError:
        return False
    return _sha256_file(path) == expected


def _write_sidecar(path: str) -> str:
    digest = _sha256_file(path)
    with open(path + ".sha256", "w", encoding="utf-8") as fh:
        fh.write(digest + "  " + os.path.basename(path) + "\n")
    return digest


def _http_get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent":
                                               "survey-currents/0.10.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"basemap download failed: HTTP {exc.code} for {url}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"basemap download failed for {url}: {exc.reason}") from exc


def _http_head_size(url: str, timeout: int = 60) -> Tuple[str, int]:
    """Follow redirects with HEAD; return (final URL, Content-Length)."""
    req = urllib.request.Request(url, method="HEAD",
                                 headers={"User-Agent":
                                          "survey-currents/0.10.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            size = resp.headers.get("Content-Length")
            if size is None:
                raise RuntimeError(
                    f"no Content-Length for {resp.url}; cannot range-read")
            return resp.url, int(size)
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"basemap HEAD failed: HTTP {exc.code} for {url}") from exc


def _http_range(url: str, start: int, end: int,
                timeout: int = 120) -> bytes:
    """GET bytes [start, end) — ``end`` exclusive, like slicing."""
    req = urllib.request.Request(
        url, headers={"User-Agent": "survey-currents/0.10.0",
                      "Range": f"bytes={start}-{end - 1}"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if resp.status != 206:
                raise RuntimeError(
                    f"range request not honored (HTTP {resp.status}) for "
                    f"{url}; server must support Accept-Ranges")
            return resp.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"range download failed: HTTP {exc.code} for {url}") from exc


# ---------------------------------------------------------------------------
# GEBCO tile layout + ranged zip extraction
# ---------------------------------------------------------------------------

def _tile_bounds(name: str) -> Tuple[float, float, float, float]:
    """(lon_min, lat_min, lon_max, lat_max) parsed from a tile file name."""
    m = _TILE_NAME_RE.match(name)
    if not m:
        raise ValueError(f"not a GEBCO 2024 tile name: {name!r}")
    south = float(m.group("south"))
    north = float(m.group("north"))
    west = float(m.group("west"))
    east = float(m.group("east"))
    return (west, south, east, north)


def _all_gebco_tiles() -> List[Tuple[str, Tuple[float, float, float, float]]]:
    """The eight 90°×90° GEBCO 2024 tiles: (file name, bounds)."""
    tiles = []
    for north, south in ((90.0, 0.0), (0.0, -90.0)):
        for west, east in ((-180.0, -90.0), (-90.0, 0.0),
                           (0.0, 90.0), (90.0, 180.0)):
            name = (f"gebco_2024_n{north:.1f}_s{south:.1f}_"
                    f"w{west:.1f}_e{east:.1f}.tif")
            tiles.append((name, (west, south, east, north)))
    return tiles


def _bbox_intersects(a: Tuple[float, float, float, float],
                     b: Tuple[float, float, float, float]) -> bool:
    return not (a[2] <= b[0] or a[0] >= b[2] or
                a[3] <= b[1] or a[1] >= b[3])


def gebco_tiles_for_bbox(
        bbox: Sequence[float]) -> List[Tuple[str, Tuple[float, float,
                                                        float, float]]]:
    """GEBCO 2024 tiles intersecting ``bbox`` (antimeridian-safe).

    Returns ``[(tile file name, tile bounds), ...]`` in a deterministic
    north-to-south, west-to-east order.
    """
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    boxes = [(lon_min, lat_min, lon_max, lat_max)]
    if lon_max < lon_min:  # antimeridian wrap: split into two boxes
        boxes = [(lon_min, lat_min, 180.0, lat_max),
                 (-180.0, lat_min, lon_max, lat_max)]
    out = []
    for name, tb in _all_gebco_tiles():
        if any(_bbox_intersects(tb, b) for b in boxes):
            out.append((name, tb))
    return out


def _zip_entries(url: str, timeout: int = 120) -> Tuple[str, List[Dict[str, Any]]]:
    """Read a remote zip's central directory via HEAD + ranged GETs.

    Returns (final URL, entries) where each entry has ``name``,
    ``method``, ``compressed_size``, ``uncompressed_size``, and
    ``local_header_offset``. Only the zip tail (~central directory) is
    transferred — the multi-GB payload is never downloaded whole.
    """
    final_url, size = _http_head_size(url, timeout=timeout)
    tail_len = min(size, 1 << 20)
    tail = _http_range(final_url, size - tail_len, size, timeout=timeout)
    eocd_at = tail.rfind(b"PK\x05\x06")
    if eocd_at < 0:  # EOCD not in the last MiB: widen the tail search
        tail = _http_range(final_url, max(0, size - (1 << 22)), size,
                           timeout=timeout)
        eocd_at = tail.rfind(b"PK\x05\x06")
        if eocd_at < 0:
            raise RuntimeError(
                f"could not find zip end-of-central-directory in {final_url}")
    cd_size, cd_offset = struct.unpack("<II", tail[eocd_at + 12:eocd_at + 20])
    cd = _http_range(final_url, cd_offset, cd_offset + cd_size,
                     timeout=timeout)
    entries: List[Dict[str, Any]] = []
    pos = 0
    while pos + 46 <= len(cd):
        if cd[pos:pos + 4] != b"PK\x01\x02":
            raise RuntimeError(
                f"bad central-directory signature at offset {pos}")
        method = struct.unpack("<H", cd[pos + 10:pos + 12])[0]
        comp_size, _uncomp_size = struct.unpack("<II", cd[pos + 20:pos + 28])
        nl, el, cl = struct.unpack("<HHH", cd[pos + 28:pos + 34])
        lho = struct.unpack("<I", cd[pos + 42:pos + 46])[0]
        name = cd[pos + 46:pos + 46 + nl].decode("utf-8", "replace")
        entries.append({"name": name, "method": method,
                        "compressed_size": comp_size,
                        "local_header_offset": lho})
        pos += 46 + nl + el + cl
    return final_url, entries


def _extract_zip_entry(url: str, entry: Dict[str, Any], dest_path: str,
                       timeout: int = 600) -> str:
    """Range-download one zip entry and write it atomically to ``dest_path``.

    Returns the SHA-256 of the extracted file. Supports stored (0) and
    deflated (8) entries.
    """
    method = entry["method"]
    if method not in (0, 8):
        raise RuntimeError(
            f"unsupported zip method {method} for {entry['name']!r}")
    final_url = url  # callers pass the already-resolved final URL
    lho = entry["local_header_offset"]
    lh = _http_range(final_url, lho, lho + 64, timeout=timeout)
    if lh[:4] != b"PK\x03\x04":
        raise RuntimeError(
            f"bad local file header for {entry['name']!r}")
    nl, el = struct.unpack("<HH", lh[26:30])
    data_off = lho + 30 + nl + el
    comp_size = entry["compressed_size"]
    tmp = dest_path + ".tmp"
    sha = hashlib.sha256()
    chunk = 8 * 1024 * 1024
    if method == 0:
        with open(tmp, "wb") as fh:
            got = 0
            while got < comp_size:
                buf = _http_range(final_url, data_off + got,
                                  min(data_off + got + chunk,
                                      data_off + comp_size),
                                  timeout=timeout)
                fh.write(buf)
                sha.update(buf)
                got += len(buf)
    else:
        dec = zlib.decompressobj(-15)
        with open(tmp, "wb") as fh:
            got = 0
            while got < comp_size:
                buf = _http_range(final_url, data_off + got,
                                  min(data_off + got + chunk,
                                      data_off + comp_size),
                                  timeout=timeout)
                raw = dec.decompress(buf)
                fh.write(raw)
                sha.update(raw)
                got += len(buf)
            tail = dec.flush()
            fh.write(tail)
            sha.update(tail)
    os.replace(tmp, dest_path)
    digest = _write_sidecar(dest_path)
    assert digest == sha.hexdigest()
    return digest


def ensure_gebco_tile(tile_name: str,
                      timeout: int = 600) -> Tuple[str, str]:
    """Ensure a GEBCO 2024 tile is in the local cache; return (path, sha256).

    Cache hits (valid ``.sha256`` sidecar) never touch the network. A
    corrupt cache entry is re-downloaded, never trusted.
    """
    cache = os.path.join(basemap_cache_dir(), "gebco", "gebco_2024")
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, tile_name)
    if _sidecar_ok(path):
        with open(path + ".sha256", "r", encoding="utf-8") as fh:
            return path, fh.read().strip().split()[0]
    final_url, entries = _zip_entries(GEBCO_ZIP_URL, timeout=timeout)
    match = [e for e in entries if e["name"] == tile_name]
    if not match:
        raise RuntimeError(
            f"tile {tile_name!r} not found in the GEBCO 2024 zip index")
    digest = _extract_zip_entry(final_url, match[0], path, timeout=timeout)
    return path, digest

# ---------------------------------------------------------------------------
# Minimal tiled-GeoTIFF reader (stdlib only)
# ---------------------------------------------------------------------------

class TiledGeoTIFF:
    """Window reader for tiled or stripped, deflate or uncompressed, int16
    GeoTIFFs.

    Parses the IFD with ``struct`` and decodes only the internal
    tiles/strips intersecting the requested pixel window (``zlib`` for
    deflate). Enough for the GEBCO 15-arc-second tiles (official 2024
    tiles: stripped, uncompressed); not a general TIFF library.
    """

    def __init__(self, path: str) -> None:
        # File-backed: only the TIFF header/IFD is read up front; internal
        # tiles/strips are fetched with seek+read on demand so a ~933 MB
        # GEBCO tile never has to sit in RAM.
        self._fh = open(path, "rb")
        try:
            self._parse_header(path)
        except Exception:
            self.close()
            raise

    def _parse_header(self, path: str) -> None:
        self._head = self._fh.read(65536)
        d = self._head
        if d[:2] == b"II":
            self._bo = "<"
        elif d[:2] == b"MM":
            self._bo = ">"
        else:
            raise ValueError(f"{path}: not a TIFF (bad byte order mark)")
        if struct.unpack(self._bo + "H", d[2:4])[0] != 42:
            raise ValueError(f"{path}: not a TIFF (bad magic)")
        ifd_off = struct.unpack(self._bo + "I", d[4:8])[0]
        n = struct.unpack(self._bo + "H", d[ifd_off:ifd_off + 2])[0]
        # The IFD entry table may extend past the initial header buffer
        # (big tile-offset tables); extend it if needed.
        need = ifd_off + 2 + 12 * n + 4
        if need > len(self._head):
            self._head += self._fh.read(need - len(self._head))
            d = self._head
        tags: Dict[int, Any] = {}
        for i in range(n):
            e = d[ifd_off + 2 + 12 * i:ifd_off + 2 + 12 * (i + 1)]
            tag, typ, cnt = struct.unpack(self._bo + "HHI", e[:8])
            raw = e[8:12]
            tags[tag] = self._tag_value(typ, cnt, raw, ifd_off)
        bo = self._bo

        def _req(tag: int, name: str) -> Any:
            if tag not in tags:
                raise ValueError(f"{path}: GeoTIFF missing {name} (tag {tag})")
            return tags[tag]

        self.width = int(_req(256, "ImageWidth"))
        self.height = int(_req(257, "ImageLength"))
        bps = _req(258, "BitsPerSample")
        bps = bps[0] if isinstance(bps, (list, tuple)) else bps
        if int(bps) != 16:
            raise ValueError(f"{path}: only 16-bit samples supported, "
                             f"got BitsPerSample={bps}")
        sfmt = tags.get(339, 1)
        sfmt = sfmt[0] if isinstance(sfmt, (list, tuple)) else sfmt
        if int(sfmt) != 2:
            raise ValueError(f"{path}: only signed samples supported, "
                             f"got SampleFormat={sfmt}")
        comp = int(_req(259, "Compression"))
        if comp not in (1, 8):
            raise ValueError(f"{path}: only uncompressed (1) or deflate (8) "
                             f"supported, got Compression={comp}")
        self._compression = comp
        self._tiled = 322 in tags
        if self._tiled:
            self.tile_w = int(_req(322, "TileWidth"))
            self.tile_h = int(_req(323, "TileLength"))
            offs = _req(324, "TileOffsets")
            bcs = _req(325, "TileByteCounts")
            if not isinstance(offs, (list, tuple)):
                offs = [offs]
            if not isinstance(bcs, (list, tuple)):
                bcs = [bcs]
            self._block_offsets = [int(v) for v in offs]
            self._block_bytecounts = [int(v) for v in bcs]
            ntx = (self.width + self.tile_w - 1) // self.tile_w
            nty = (self.height + self.tile_h - 1) // self.tile_h
            if len(self._block_offsets) != ntx * nty:
                raise ValueError(
                    f"{path}: tile table has {len(self._block_offsets)} entries, "
                    f"expected {ntx * nty}")
            self._nbx, self._nby = ntx, nty
        else:
            # Stripped layout (official GEBCO 2024 tiles: 21600 strips of
            # one row each). A "block" is one strip: full image width.
            self.tile_w = self.width
            self.tile_h = int(_req(278, "RowsPerStrip"))
            offs = _req(273, "StripOffsets")
            bcs = _req(279, "StripByteCounts")
            if not isinstance(offs, (list, tuple)):
                offs = [offs]
            if not isinstance(bcs, (list, tuple)):
                bcs = [bcs]
            self._block_offsets = [int(v) for v in offs]
            self._block_bytecounts = [int(v) for v in bcs]
            nstrips = ((self.height + self.tile_h - 1) // self.tile_h)
            if len(self._block_offsets) != nstrips:
                raise ValueError(
                    f"{path}: strip table has {len(self._block_offsets)} entries, "
                    f"expected {nstrips}")
            self._nbx, self._nby = 1, nstrips
        pixscale = _req(33550, "ModelPixelScaleTag")
        tiepoint = _req(33922, "ModelTiepointTag")
        self._sx = float(pixscale[0])
        # GeoTIFF convention: ScaleY is positive and latitude *decreases*
        # with row (verified on the official GEBCO 2024 tile, 2026-09-27).
        # Normalize to a signed "latitude change per row" so both
        # conventions (positive sy per spec, negative sy in older test
        # fixtures) behave identically.
        self._sy = -abs(float(pixscale[1]))
        # Tiepoint (I,J,K, X,Y,Z): pixel (I,J) sits at geo (X,Y).
        self._tie_i = float(tiepoint[0])
        self._tie_j = float(tiepoint[1])
        self._tie_x = float(tiepoint[3])
        self._tie_y = float(tiepoint[4])
        self.path = path

    def _read_range(self, off: int, n: int) -> bytes:
        """Read ``n`` bytes at ``off``, from the header buffer if possible."""
        end = off + n
        if end <= len(self._head):
            return self._head[off:end]
        self._fh.seek(off)
        data = self._fh.read(n)
        if len(data) != n:
            raise ValueError(f"{self.path}: short read at offset {off}")
        return data

    def close(self) -> None:
        fh, self._fh = self._fh, None  # type: ignore[assignment]
        if fh is not None:
            fh.close()

    def __enter__(self) -> "TiledGeoTIFF":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    def _tag_value(self, typ: int, cnt: int, raw: bytes,
                   ifd_off: int) -> Any:
        """Decode one IFD entry value (inline or offset)."""
        bo = self._bo
        sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 12: 8, 16: 8}
        if typ not in sizes:
            raise ValueError(f"unsupported TIFF field type {typ}")
        total = sizes[typ] * cnt
        buf = raw if total <= 4 else None
        if buf is None:
            off = struct.unpack(bo + "I", raw)[0]
            buf = self._read_range(off, total)
        if typ == 1:
            vals = list(buf[:cnt])
        elif typ == 2:
            return buf[:cnt].split(b"\x00")[0].decode("ascii", "replace")
        elif typ == 3:
            vals = list(struct.unpack(bo + f"{cnt}H", buf[:2 * cnt]))
        elif typ == 4:
            vals = list(struct.unpack(bo + f"{cnt}I", buf[:4 * cnt]))
        elif typ == 5:
            vals = [struct.unpack(bo + "II", buf[8 * i:8 * i + 8])
                    for i in range(cnt)]
            vals = [a / b if b else float("nan") for a, b in vals]
        elif typ == 12:
            vals = list(struct.unpack(bo + f"{cnt}d", buf[:8 * cnt]))
        elif typ == 16:
            vals = list(struct.unpack(bo + f"{cnt}q", buf[:8 * cnt]))
        else:  # pragma: no cover - guarded above
            raise ValueError(f"unsupported TIFF field type {typ}")
        return vals[0] if cnt == 1 else vals

    # -- georeferencing ----------------------------------------------------

    def pixel_to_geo(self, col: float, row: float) -> Tuple[float, float]:
        """Geo (lon, lat) of pixel ``(col, row)`` (fractional allowed)."""
        lon = self._tie_x + (col - self._tie_i) * self._sx
        lat = self._tie_y + (row - self._tie_j) * self._sy
        return lon, lat

    def geo_to_pixel(self, lon: float, lat: float) -> Tuple[float, float]:
        """Fractional pixel ``(col, row)`` of geo ``(lon, lat)``."""
        col = self._tie_i + (lon - self._tie_x) / self._sx
        row = self._tie_j + (lat - self._tie_y) / self._sy
        return col, row

    def cell_center_lons(self) -> np.ndarray:
        return np.array(
            [self.pixel_to_geo(c + 0.5, 0)[0] for c in range(self.width)])

    def cell_center_lats(self) -> np.ndarray:
        return np.array(
            [self.pixel_to_geo(0, r + 0.5)[1] for r in range(self.height)])

    # -- window reads -------------------------------------------------------

    def _decode_block(self, bx: int, by: int) -> np.ndarray:
        """Decode one internal block: tile (bx, by), or strip ``by``.

        Tiled fixtures pad edge tiles to the full tile size; strips
        hold exactly their rows. The caller slices to the image.
        """
        idx = by * self._nbx + bx
        off = self._block_offsets[idx]
        nbytes = self._block_bytecounts[idx]
        raw = self._read_range(off, nbytes)
        if self._compression == 8:
            # TIFF Compression=8 (Deflate) is the zlib (RFC 1950) stream.
            raw = zlib.decompress(raw)
        arr = np.frombuffer(raw, dtype=self._bo + "i2")
        if self._tiled:
            return arr.reshape(self.tile_h, self.tile_w)
        rows = min(self.tile_h, self.height - by * self.tile_h)
        return arr.reshape(rows, self.width)

    def _decode_tile(self, tx: int, ty: int) -> np.ndarray:
        """Decode tile (tx, ty) — kept for backward compatibility."""
        return self._decode_block(tx, ty)

    def read_window(self, col0: int, row0: int,
                    ncol: int, nrow: int) -> np.ndarray:
        """Read an int16 pixel window ``[col0, col0+ncol)``×``[row0, …)``.

        Only the internal tiles/strips intersecting the window are
        decoded. Out-of-range edges are clipped to the image.
        """
        col0 = max(0, col0)
        row0 = max(0, row0)
        col1 = min(self.width, col0 + ncol)
        row1 = min(self.height, row0 + nrow)
        if col1 <= col0 or row1 <= row0:
            raise ValueError("window does not intersect the image")
        out = np.empty((row1 - row0, col1 - col0), dtype=np.int16)
        for by in range(row0 // self.tile_h, (row1 - 1) // self.tile_h + 1):
            brow0 = by * self.tile_h
            brow1 = min(brow0 + self.tile_h, self.height)
            for bx in range(col0 // self.tile_w,
                             (col1 - 1) // self.tile_w + 1):
                block = self._decode_block(bx, by)
                tc0 = max(col0, bx * self.tile_w)
                tc1 = min(col1, (bx + 1) * self.tile_w)
                tr0 = max(row0, brow0)
                tr1 = min(row1, brow1)
                out[tr0 - row0:tr1 - row0, tc0 - col0:tc1 - col0] = \
                    block[tr0 - brow0:tr1 - brow0,
                          tc0 - bx * self.tile_w:tc1 - bx * self.tile_w]
        return out


def _block_average(a: np.ndarray, factor: int) -> np.ndarray:
    """NaN-aware block mean over ``factor``×``factor`` cells (2D)."""
    if factor == 1:
        return a.astype(float)
    ny, nx = a.shape
    ny_t, nx_t = (ny // factor) * factor, (nx // factor) * factor
    a = a[:ny_t, :nx_t].astype(float)
    with np.errstate(invalid="ignore"):
        sums = np.nansum(a.reshape(ny_t // factor, factor,
                                   nx_t // factor, factor), axis=(1, 3))
        counts = np.sum(~np.isnan(a).reshape(ny_t // factor, factor,
                                             nx_t // factor, factor),
                        axis=(1, 3))
    out = sums / np.where(counts == 0, np.nan, counts)
    return out

# ---------------------------------------------------------------------------
# TopoField + fetch_gebco
# ---------------------------------------------------------------------------

@dataclass
class TopoField:
    """Static GEBCO topography/bathymetry grid.

    ``values`` is float elevation in meters (positive up; ocean depths
    negative) shaped ``(1, nlat, nlon)`` on a regular lat/lon grid —
    the single time axis is structural (``times ==
    [2024-01-01T00:00:00+00:00]``): GEBCO is a static compilation, not a
    time series. Renderers must label it "static", never a data date.

    ``provenance`` records the grid version, the requested subset bbox,
    the tile file name(s) read, their SHA-256 digests, the cache paths,
    the output resolution, and the retrieval timestamp.
    """

    times: List[_dt.datetime]
    lats: np.ndarray
    lons: np.ndarray
    values: np.ndarray
    bbox: Tuple[float, float, float, float] = (-125.0, 25.0, -66.0, 49.0)
    resolution: float = 1.0
    source: str = "gebco-2024"
    units: str = "m"
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
                f"TopoField.resolution must be > 0, got {self.resolution}")

    def __len__(self) -> int:
        return len(self.times)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    @property
    def shape(self) -> Tuple[int, int, int]:
        """(ntime, nlat, nlon) — ntime is always 1 (static grid)."""
        return self.values.shape  # type: ignore[return-value]

    @property
    def elevation(self) -> np.ndarray:
        """The single static elevation frame, shape (nlat, nlon)."""
        return np.asarray(self.values[0])

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
    def from_dict(cls, data: Dict[str, Any]) -> "TopoField":
        """Rebuild from :meth:`to_dict` output."""
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
            resolution=float(data.get("resolution", 1.0)),
            source=str(data.get("source", "gebco-2024")),
            units=str(data.get("units", "m")),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        """Write JSON to ``path``; returns ``path``."""
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "TopoField":
        """Read a field written by :meth:`to_json`."""
        import json
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture --------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-125.0, 25.0, -66.0, 49.0),
                  resolution: float = 1.0,
                  seed: int = 7,
                  source: str = "synthetic") -> "TopoField":
        """Deterministic offline fixture: fake continents and ocean basins.

        A few Gaussian "continents" (positive elevation with ridged
        noise) over a negative "ocean" background — enough structure to
        exercise tinting, hillshading, and mosaicking without a network.
        Fully offline and deterministic for ``seed``.
        """
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        if not resolution > 0:
            raise ValueError(f"resolution must be > 0, got {resolution}")
        rng = np.random.default_rng(seed)
        lons = np.arange(lon_min + resolution / 2, lon_max, resolution)
        lats = np.arange(lat_max - resolution / 2, lat_min, -resolution)
        lon2, lat2 = np.meshgrid(lons, lats)
        elev = np.full(lon2.shape, -3800.0)
        n_land = 5
        cx = rng.uniform(lon_min, lon_max, n_land)
        cy = rng.uniform(lat_min, lat_max, n_land)
        for x, y in zip(cx, cy):
            dist2 = (lon2 - x) ** 2 + (lat2 - y) ** 2
            elev += 5200.0 * np.exp(-dist2 / (2 * 6.0 ** 2))
        elev += rng.normal(0.0, 180.0, size=elev.shape)
        # A trench: a deep linear gash somewhere in the ocean part.
        trench_x = lon_min + 0.3 * (lon_max - lon_min)
        elev -= 4500.0 * np.exp(-((lon2 - trench_x) ** 2) / (2 * 0.8 ** 2))
        return cls(
            times=[GEBCO_STATIC_TIME],
            lats=lats, lons=lons, values=elev[np.newaxis, :, :],
            bbox=(lon_min, lat_min, lon_max, lat_max),
            resolution=resolution, source=source, units="m",
            provenance={"synthetic": True, "seed": seed,
                        "generator": "TopoField.synthetic"},
        )


def _parse_resolution(resolution: Any) -> Tuple[bool, float]:
    """Normalize ``resolution`` to (is_native, degrees)."""
    if isinstance(resolution, str):
        s = resolution.strip().lower()
        if s in ("15s", "15arcsec", "15-arc-second", "native"):
            return True, GEBCO_NATIVE_DEG
        m = re.match(r"^(\d+(?:\.\d+)?)\s*(deg|°|arcmin|')$", s)
        if m:
            val = float(m.group(1))
            if m.group(2) in ("arcmin", "'"):
                val /= 60.0
            return False, val
        raise ValueError(
            f"resolution must be '15s' or a positive number of degrees "
            f"(optionally suffixed 'deg' or 'arcmin'), got {resolution!r}")
    val = float(resolution)
    if not val > 0:
        raise ValueError(f"resolution must be > 0, got {resolution!r}")
    return False, val


def _native_window_for_bbox(tif: TiledGeoTIFF,
                            bbox: Tuple[float, float, float, float]
                            ) -> Tuple[int, int, int, int]:
    """Pixel window ``(col0, row0, ncol, nrow)`` covering ``bbox``."""
    lon_min, lat_min, lon_max, lat_max = bbox
    c0f, r1f = tif.geo_to_pixel(lon_min, lat_min)
    c1f, r0f = tif.geo_to_pixel(lon_max, lat_max)
    # Row 0 is the northern edge (lat decreases with row in GEBCO tiles).
    col0 = int(math.floor(min(c0f, c1f)))
    col1 = int(math.ceil(max(c0f, c1f)))
    row0 = int(math.floor(min(r0f, r1f)))
    row1 = int(math.ceil(max(r0f, r1f)))
    col0 = max(0, col0)
    row0 = max(0, row0)
    ncol = min(tif.width, col1) - col0
    nrow = min(tif.height, row1) - row0
    if ncol <= 0 or nrow <= 0:
        raise ValueError(f"bbox {list(bbox)} does not intersect tile "
                         f"{tif.path}")
    return col0, row0, ncol, nrow


def fetch_gebco(bbox: Sequence[float],
                resolution: Any = "15s",
                timeout: int = 600) -> TopoField:
    """Fetch the GEBCO 2024 elevation grid subset for ``bbox``.

    Args:
        bbox: (lon_min, lat_min, lon_max, lat_max) in -180..180 degrees;
            antimeridian-crossing boxes wrap.
        resolution: ``"15s"`` (native 15 arc-second) or a positive
            number of degrees (e.g. ``0.25``); degree values must be an
            integer multiple of the native 15 arc-seconds and are
            produced by NaN-aware block averaging. The suffix forms
            ``"30arcmin"`` / ``"1deg"`` are also accepted.
        timeout: seconds per HTTP request.

    Tiles are downloaded once (ranged zip-entry extraction, ~0.5 GB per
    90°×90° tile) into the user cache with SHA-256 sidecars; repeated
    calls never re-download. Returns a :class:`TopoField` (static,
    single-timestep). ``provenance`` records the grid version, the
    requested subset bbox, the tile(s) read, and the output resolution.
    """
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    is_native, res_deg = _parse_resolution(resolution)
    factor = 1
    if not is_native:
        ratio = res_deg / GEBCO_NATIVE_DEG
        factor = int(round(ratio))
        if factor < 1 or abs(ratio - factor) > 1e-6:
            raise ValueError(
                f"resolution {res_deg}° is not an integer multiple of the "
                f"native 15 arc-seconds ({GEBCO_NATIVE_DEG}°); use e.g. "
                f"0.25, 0.5, or 1.0")
        res_deg = factor * GEBCO_NATIVE_DEG
    # The 90° tiles must land exactly on the output grid: res_deg has to
    # divide 90° (all of 15s, 0.25°, 0.5°, 1°, ... do).
    if abs(90.0 / res_deg - round(90.0 / res_deg)) > 1e-9:
        raise ValueError(
            f"resolution {res_deg}° does not divide the 90° GEBCO tiles "
            f"evenly; use a resolution whose 90° quotient is an integer "
            f"(e.g. 0.25, 0.5, 1.0)")

    # Output grid, snapped to the global res_deg lattice. Antimeridian
    # wraps are handled in unwrapped longitude space.
    wrap = lon_max < lon_min
    lon_max_u = lon_max + 360.0 if wrap else lon_max
    gx0 = math.floor(lon_min / res_deg + 1e-9) * res_deg
    gx1 = math.ceil(lon_max_u / res_deg - 1e-9) * res_deg
    gy1 = math.ceil(lat_max / res_deg - 1e-9) * res_deg
    gy0 = math.floor(lat_min / res_deg + 1e-9) * res_deg
    nlon = int(round((gx1 - gx0) / res_deg))
    nlat = int(round((gy1 - gy0) / res_deg))
    if nlon <= 0 or nlat <= 0:
        raise ValueError(f"bbox {list(bbox)} is smaller than one "
                         f"{res_deg}° output cell")
    if is_native and (nlon > GEBCO_MAX_NATIVE_WINDOW or
                      nlat > GEBCO_MAX_NATIVE_WINDOW):
        raise ValueError(
            f"native-resolution output for {list(bbox)} would be "
            f"{nlon}×{nlat} cells, over the {GEBCO_MAX_NATIVE_WINDOW} px "
            f"safety cap — request a coarser resolution (e.g. 0.25)")

    lons_u = gx0 + (np.arange(nlon) + 0.5) * res_deg
    lats = gy1 - (np.arange(nlat) + 0.5) * res_deg
    grid = np.full((nlat, nlon), np.nan)

    tiles = gebco_tiles_for_bbox(bbox)
    if not tiles:
        raise ValueError(f"bbox {list(bbox)} intersects no GEBCO tile")

    tile_records: List[Dict[str, Any]] = []
    for tile_name, (tw, ts, te, tn) in tiles:
        path, digest = ensure_gebco_tile(tile_name, timeout=timeout)
        with TiledGeoTIFF(path) as tif:
            # In the wrap case a tile west of lon_min is used shifted +360°.
            shifts = (360.0,) if (wrap and te <= lon_min) else (0.0,)
            for shift in shifts:
                tws, tes = tw + shift, te + shift
                if tes <= gx0 or tws >= gx1:
                    continue
                if abs((tws - gx0) / res_deg -
                       round((tws - gx0) / res_deg)) > 1e-6:
                    raise RuntimeError(
                        f"tile {tile_name} does not align with the "
                        f"{res_deg}° output grid")
                i_t = int(round((tws - gx0) / res_deg))
                ni = int(round((tes - tws) / res_deg))
                j_t = int(round((gy1 - tn) / res_deg))
                nj = int(round((tn - ts) / res_deg))
                a0, a1 = max(0, i_t), min(nlon, i_t + ni)
                b0, b1 = max(0, j_t), min(nlat, j_t + nj)
                if a1 <= a0 or b1 <= b0:
                    continue
                # Exact index math: output cell (j, i) covers native cells
                # [(i-i_t)*factor, (i-i_t+1)*factor) × same for rows, because
                # the tile origin sits on the global grid lattice.
                c0 = (a0 - i_t) * factor
                r0 = (b0 - j_t) * factor
                nc = (a1 - a0) * factor
                nr = (b1 - b0) * factor
                win = tif.read_window(c0, r0, nc, nr).astype(float)
                avg = _block_average(win, factor)
                if avg.shape != (b1 - b0, a1 - a0):
                    raise RuntimeError(
                        f"tile {tile_name}: averaged window shape "
                        f"{avg.shape} != placement {(b1 - b0, a1 - a0)}")
                grid[b0:b1, a0:a1] = avg
        tile_records.append({"tile": tile_name, "sha256": digest,
                             "cache_path": path,
                             "bounds": [tw, ts, te, tn]})

    # Wrap longitudes back to -180..180 for the output axis.
    lons = np.where(lons_u > 180.0, lons_u - 360.0, lons_u)

    provenance = {
        "grid": f"GEBCO_{GEBCO_VERSION}",
        "grid_resolution_arcsec": GEBCO_NATIVE_ARCSEC,
        "subset_bbox": [lon_min, lat_min, lon_max, lat_max],
        "output_grid": {"lon_min": float(gx0), "lon_max": float(gx1),
                        "lat_min": float(gy0), "lat_max": float(gy1),
                        "resolution_deg": res_deg,
                        "nlon": nlon, "nlat": nlat},
        "output_resolution_deg": res_deg,
        "tiles": tile_records,
        "retrieved_utc": _dt.datetime.now(
            _dt.timezone.utc).isoformat(),
        "static_compilation": True,
    }
    return TopoField(
        times=[GEBCO_STATIC_TIME],
        lats=lats, lons=lons, values=grid[np.newaxis, :, :],
        bbox=(lon_min, lat_min, lon_max, lat_max),
        resolution=res_deg, source=f"gebco-{GEBCO_VERSION}", units="m",
        provenance=provenance,
    )

# ---------------------------------------------------------------------------
# Natural Earth vectors (stdlib shapefile reader -> GeoJSON)
# ---------------------------------------------------------------------------

def _ne_zip_url(scale: str, layer: str) -> str:
    subdir, stem = NE_LAYER_FILES[layer]
    return (f"{NE_BASE_URL}/{subdir.format(scale=scale)}/"
            f"{stem.format(scale=scale)}")


def _ne_cache_path(scale: str, layer: str) -> str:
    subdir, stem = NE_LAYER_FILES[layer]
    cache = os.path.join(basemap_cache_dir(), "naturalearth", scale)
    os.makedirs(cache, exist_ok=True)
    return os.path.join(cache, stem.format(scale=scale))


def ensure_naturalearth_zip(scale: str, layer: str,
                            timeout: int = 120) -> Tuple[str, str]:
    """Ensure a Natural Earth zip is cached; return (path, sha256)."""
    if scale not in NATURAL_EARTH_SCALES:
        raise ValueError(
            f"scale must be one of {NATURAL_EARTH_SCALES}, got {scale!r}")
    if layer not in NE_LAYER_FILES:
        raise ValueError(
            f"layer must be one of {NATURAL_EARTH_LAYERS}, got {layer!r}")
    path = _ne_cache_path(scale, layer)
    if _sidecar_ok(path):
        with open(path + ".sha256", "r", encoding="utf-8") as fh:
            return path, fh.read().strip().split()[0]
    data = _http_get(_ne_zip_url(scale, layer), timeout=timeout)
    if not data.startswith(b"PK"):
        raise RuntimeError(
            f"{_ne_zip_url(scale, layer)} did not return a zip archive")
    _write_atomic(path, data)
    return path, _write_sidecar(path)


def _read_shp_polylines(shp: bytes) -> List[List[Tuple[float, float]]]:
    """Parse an ESRI shapefile of PolyLine (3) into part polylines."""
    return _read_shp_parts(shp, allowed=(3,))


def _read_shp_parts(shp: bytes, allowed: Tuple[int, ...]) -> List[Any]:
    """Parse shape records into ``[(shape_type, [parts])]``.

    Each part is a list of (x, y) tuples. Only 2D point arrays are read
    (Z/M ordinates, when present, are skipped).
    """
    if len(shp) < 100 or struct.unpack(">i", shp[0:4])[0] != 9994:
        raise ValueError("not an ESRI shapefile (bad header)")
    file_type = struct.unpack("<i", shp[32:36])[0]
    if file_type not in allowed:
        raise ValueError(
            f"shapefile shape type {file_type} not in {allowed}")
    out: List[Any] = []
    off = 100
    has_z = file_type in (13, 15)
    has_m = file_type in (23, 25)
    while off + 8 <= len(shp):
        _rec_no, content_len = struct.unpack(">ii", shp[off:off + 8])
        rec_end = off + 8 + content_len * 2
        if rec_end > len(shp):
            break
        st = struct.unpack("<i", shp[off + 8:off + 12])[0]
        if st in allowed:
            n_parts, n_points = struct.unpack(
                "<ii", shp[off + 44:off + 52])
            parts = struct.unpack(
                f"<{n_parts}i", shp[off + 52:off + 52 + 4 * n_parts])
            pts_off = off + 52 + 4 * n_parts
            pts: List[Tuple[float, float]] = []
            for i in range(n_points):
                x, y = struct.unpack(
                    "<dd", shp[pts_off + 16 * i:pts_off + 16 * i + 16])
                pts.append((x, y))
            # Skip Z/M blocks when present (keeps offsets honest).
            _ = has_z, has_m
            rings = []
            for pi, p0 in enumerate(parts):
                p1 = parts[pi + 1] if pi + 1 < n_parts else n_points
                rings.append(pts[p0:p1])
            out.append((st, rings))
        off = rec_end
    return out


def _signed_area(ring: List[Tuple[float, float]]) -> float:
    s = 0.0
    n = len(ring)
    for i in range(n):
        x0, y0 = ring[i]
        x1, y1 = ring[(i + 1) % n]
        s += x0 * y1 - x1 * y0
    return s / 2.0


def _point_in_ring(pt: Tuple[float, float],
                   ring: List[Tuple[float, float]]) -> bool:
    x, y = pt
    inside = False
    n = len(ring)
    for i in range(n):
        x0, y0 = ring[i]
        x1, y1 = ring[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            xin = x0 + (x1 - x0) * (y - y0) / (y1 - y0)
            if x < xin:
                inside = not inside
    return inside


def _rings_to_polygons(rings: List[List[Tuple[float, float]]]
                       ) -> List[List[List[Tuple[float, float]]]]:
    """Group shapefile rings into polygons: [exterior, hole, ...].

    Clockwise rings (negative signed area, x-right/y-up) are exteriors;
    counter-clockwise rings are holes, assigned to the smallest-area
    exterior containing them. Degenerate rings (< 4 points) are dropped.
    """
    rings = [r for r in rings if len(r) >= 4]
    exteriors = [r for r in rings if _signed_area(r) < 0]
    holes = [r for r in rings if _signed_area(r) >= 0]
    polys: List[List[List[Tuple[float, float]]]] = [[e] for e in exteriors]
    for h in holes:
        best, best_area = -1, math.inf
        for i, (e, *_) in enumerate(polys):
            area = abs(_signed_area(e))
            if area < best_area and _point_in_ring(h[0], e):
                best, best_area = i, area
        if best >= 0:
            polys[best].append(h)
        # Orphan holes (no containing exterior) are dropped: emitting
        # them as exteriors would be a lie about the source geometry.
    return polys


def _clip_polyline_bbox(line: List[Tuple[float, float]],
                        bbox: Tuple[float, float, float, float]
                        ) -> List[List[Tuple[float, float]]]:
    """Clip a polyline to ``bbox`` (Liang–Barsky per segment)."""
    xmin, ymin, xmax, ymax = bbox
    segs: List[List[Tuple[float, float]]] = []
    cur: List[Tuple[float, float]] = []

    def _clip_pt(px: float, py: float, dx: float, dy: float):
        t0, t1 = 0.0, 1.0
        for p, q in ((-dx, px - xmin), (dx, xmax - px),
                     (-dy, py - ymin), (dy, ymax - py)):
            if p == 0:
                if q < 0:
                    return None
            else:
                r = q / p
                if p < 0:
                    if r > t1:
                        return None
                    t0 = max(t0, r)
                else:
                    if r < t0:
                        return None
                    t1 = min(t1, r)
        return (px + t0 * dx, py + t0 * dy, px + t1 * dx, py + t1 * dy)

    for (x0, y0), (x1, y1) in zip(line, line[1:]):
        hit = _clip_pt(x0, y0, x1 - x0, y1 - y0)
        if hit is None:
            if cur:
                segs.append(cur)
                cur = []
            continue
        ax, ay, bx, by = hit
        if cur and cur[-1] == (ax, ay):
            cur.append((bx, by))
        else:
            if cur:
                segs.append(cur)
            cur = [(ax, ay), (bx, by)]
    if cur:
        segs.append(cur)
    return [s for s in segs if len(s) >= 2]


def _shp_to_geojson(shp_name: str, shp: bytes, layer: str,
                    bbox: Tuple[float, float, float, float]) -> Dict[str, Any]:
    """Convert shapefile bytes to a GeoJSON FeatureCollection.

    Polylines (coastline) are hard-clipped to ``bbox``; polygons
    (countries/land/lakes) are prefiltered to features whose bbox
    intersects (documented, not true-clipped — the renderer clips
    visually via axes limits).
    """
    is_line = layer == "coastline"
    recs = _read_shp_parts(shp, allowed=(3,) if is_line else (5,))
    features: List[Dict[str, Any]] = []
    for st, rings in recs:
        if is_line:
            for ring in rings:
                for seg in _clip_polyline_bbox(ring, bbox):
                    features.append({
                        "type": "Feature",
                        "properties": {"layer": layer},
                        "geometry": {"type": "LineString",
                                     "coordinates": [[x, y] for x, y in seg]},
                    })
        else:
            xs = [x for r in rings for x, _ in r]
            ys = [y for r in rings for _, y in r]
            if not xs or max(xs) < bbox[0] or min(xs) > bbox[2] or \
               max(ys) < bbox[1] or min(ys) > bbox[3]:
                continue
            polys = _rings_to_polygons(rings)
            if not polys:
                continue
            coords = [[[x, y] for x, y in ring] for poly in polys
                      for ring in [poly[0]]]
            # Build proper Polygon / MultiPolygon nesting.
            poly_coords = [[[ [x, y] for x, y in ring] for ring in poly]
                           for poly in polys]
            if len(poly_coords) == 1:
                geom: Dict[str, Any] = {"type": "Polygon",
                                        "coordinates": poly_coords[0]}
            else:
                geom = {"type": "MultiPolygon", "coordinates": poly_coords}
            features.append({
                "type": "Feature",
                "properties": {"layer": layer},
                "geometry": geom,
            })
    return {"type": "FeatureCollection", "features": features}


def fetch_naturalearth(bbox: Sequence[float],
                       scale: str = "110m",
                       layers: Sequence[str] = ("coastline", "countries"),
                       timeout: int = 120) -> Dict[str, Dict[str, Any]]:
    """Fetch Natural Earth vectors as GeoJSON FeatureCollections.

    Args:
        bbox: (lon_min, lat_min, lon_max, lat_max); antimeridian-crossing
            boxes wrap (clipped per 180° half).
        scale: ``"110m"`` (default), ``"50m"``, or ``"10m"``.
        layers: any of ``"coastline"``, ``"countries"``, ``"land"``,
            ``"lakes"``.
        timeout: seconds per HTTP request.

    Zips are downloaded once per (scale, layer) into the user cache with
    SHA-256 sidecars; repeated calls never re-download. Returns
    ``{layer: FeatureCollection}``.
    """
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    layers = tuple(layers)
    for layer in layers:
        if layer not in NE_LAYER_FILES:
            raise ValueError(
                f"layer must be one of {NATURAL_EARTH_LAYERS}, "
                f"got {layer!r}")
    if scale not in NATURAL_EARTH_SCALES:
        raise ValueError(
            f"scale must be one of {NATURAL_EARTH_SCALES}, got {scale!r}")
    boxes = [(lon_min, lat_min, lon_max, lat_max)]
    if lon_max < lon_min:
        boxes = [(lon_min, lat_min, 180.0, lat_max),
                 (-180.0, lat_min, lon_max, lat_max)]
    out: Dict[str, Dict[str, Any]] = {}
    for layer in layers:
        path, digest = ensure_naturalearth_zip(scale, layer, timeout=timeout)
        with zipfile.ZipFile(path) as zf:
            shp_names = [n for n in zf.namelist()
                         if n.lower().endswith(".shp")]
            if not shp_names:
                raise RuntimeError(f"no .shp found in {path}")
            shp = zf.read(shp_names[0])
        feats: List[Dict[str, Any]] = []
        for box in boxes:
            fc = _shp_to_geojson(shp_names[0], shp, layer, box)
            feats.extend(fc["features"])
        out[layer] = {"type": "FeatureCollection", "features": feats,
                      "properties": {
                          "scale": scale, "layer": layer,
                          "source": "Natural Earth",
                          "source_url": _ne_zip_url(scale, layer),
                          "zip_sha256": digest,
                          "subset_bbox": [lon_min, lat_min,
                                         lon_max, lat_max],
                          "retrieved_utc": _dt.datetime.now(
                              _dt.timezone.utc).isoformat(),
                      }}
    return out


def synthetic_naturalearth(
        bbox: Sequence[float] = (-125.0, 25.0, -66.0, 49.0),
        layers: Sequence[str] = ("coastline", "countries"),
        seed: int = 7) -> Dict[str, Dict[str, Any]]:
    """Deterministic offline Natural Earth stand-in (for tests/demos).

    One wavy "coastline" polyline and two rectangular "countries" —
    enough geometry to exercise the GeoJSON plumbing without a network.
    """
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    rng = np.random.default_rng(seed)
    layers = tuple(layers)
    out: Dict[str, Dict[str, Any]] = {}
    if "coastline" in layers:
        n = 40
        xs = np.linspace(lon_min, lon_max, n)
        mid = (lat_min + lat_max) / 2
        amp = (lat_max - lat_min) / 8
        ys = mid + amp * np.sin(np.linspace(0, 3 * np.pi, n) +
                                rng.uniform(0, 2 * np.pi))
        out["coastline"] = {
            "type": "FeatureCollection",
            "features": [{
                "type": "Feature",
                "properties": {"layer": "coastline", "synthetic": True},
                "geometry": {"type": "LineString",
                             "coordinates": [[float(x), float(y)]
                                             for x, y in zip(xs, ys)]},
            }],
        }
    if "countries" in layers:
        feats = []
        for i, name in enumerate(("Westland", "Eastland")):
            x0 = lon_min + i * (lon_max - lon_min) / 2
            x1 = lon_min + (i + 1) * (lon_max - lon_min) / 2
            feats.append({
                "type": "Feature",
                "properties": {"layer": "countries", "name": name,
                               "synthetic": True},
                "geometry": {"type": "Polygon", "coordinates": [[
                    [x0, lat_min], [x1, lat_min], [x1, lat_max],
                    [x0, lat_max], [x0, lat_min]]]},
            })
        out["countries"] = {"type": "FeatureCollection", "features": feats}
    for layer in layers:
        if layer not in out:
            raise ValueError(
                f"synthetic_naturalearth has no fixture for {layer!r}; "
                f"choose from ('coastline', 'countries')")
    return out

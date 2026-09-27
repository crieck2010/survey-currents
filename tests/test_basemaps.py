"""Tests for the GEBCO + Natural Earth basemap adapters (currents.basemaps).

Fully offline: the ranged zip machinery is exercised against synthetic
zip bytes served through monkeypatched HTTP helpers, the tiled-GeoTIFF
reader against hand-built deflate tiles, and the shapefile reader
against hand-built .shp bytes. ``TopoField.synthetic()`` and
``synthetic_naturalearth()`` cover the field-model paths.
"""

from __future__ import annotations

import hashlib
import io
import os
import struct
import zipfile
import zlib

import numpy as np
import pytest

import currents.basemaps as bm
from currents.basemaps import (
    GEBCO_VERSION,
    TopoField,
    _block_average,
    _clip_polyline_bbox,
    _rings_to_polygons,
    _read_shp_parts,
    basemap_cache_dir,
    ensure_gebco_tile,
    ensure_naturalearth_zip,
    fetch_gebco,
    fetch_naturalearth,
    gebco_tiles_for_bbox,
    synthetic_naturalearth,
)


# ---------------------------------------------------------------------------
# Synthetic fixtures: tiled deflate GeoTIFF + shapefile bytes
# ---------------------------------------------------------------------------

def _build_tiled_geotiff(width, height, tw, th, sx, sy, tie_x, tie_y,
                         values):
    """Minimal tiled deflate int16 GeoTIFF -> bytes.

    ``values``: (height, width) int16 array. Georeferencing:
    pixel (0,0) at (tie_x, tie_y), pixel size (sx, sy).
    """
    assert values.shape == (height, width)
    assert values.dtype == np.int16
    ntx = (width + tw - 1) // tw
    nty = (height + th - 1) // th
    blobs = []
    for ty in range(nty):
        for tx in range(ntx):
            blk = np.zeros((th, tw), dtype=np.int16)
            r0, c0 = ty * th, tx * tw
            blk[:min(th, height - r0), :min(tw, width - c0)] = \
                values[r0:min(r0 + th, height), c0:min(c0 + tw, width)]
            blobs.append(zlib.compress(blk.tobytes(), 9))

    def entry(tag, typ, cnt, val):
        return struct.pack("<HHI", tag, typ, cnt) + val

    tags = []
    tags.append(entry(256, 4, 1, struct.pack("<I", width)))
    tags.append(entry(257, 4, 1, struct.pack("<I", height)))
    tags.append(entry(258, 3, 1, struct.pack("<H", 16) + b"\x00\x00"))
    tags.append(entry(259, 3, 1, struct.pack("<H", 8) + b"\x00\x00"))
    tags.append(entry(262, 3, 1, struct.pack("<H", 1) + b"\x00\x00"))
    tags.append(entry(277, 3, 1, struct.pack("<H", 1) + b"\x00\x00"))
    tags.append(entry(339, 3, 1, struct.pack("<H", 2) + b"\x00\x00"))
    tags.append(entry(322, 4, 1, struct.pack("<I", tw)))
    tags.append(entry(323, 4, 1, struct.pack("<I", th)))
    ntags = len(tags) + 4
    ifd_len = 2 + 12 * ntags + 4
    data_off = 8 + ifd_len
    # Layout: [tile blobs][offsets][bytecounts][pixscale][tiepoint]
    offs, cursor = [], data_off
    for b in blobs:
        offs.append(cursor)
        cursor += len(b)
    bcs = [len(b) for b in blobs]
    offs_off = cursor
    cursor += 4 * len(offs)
    bcs_off = cursor
    cursor += 4 * len(bcs)
    ps_off = cursor
    cursor += 24
    tp_off = cursor

    tags.append(entry(324, 4, len(offs), struct.pack("<I", offs_off)))
    tags.append(entry(325, 4, len(bcs), struct.pack("<I", bcs_off)))
    tags.append(entry(33550, 12, 3, struct.pack("<I", ps_off)))
    tags.append(entry(33922, 12, 6, struct.pack("<I", tp_off)))

    out = bytearray()
    out += b"II" + struct.pack("<H", 42) + struct.pack("<I", 8)
    out += struct.pack("<H", len(tags))
    for t in tags:
        out += t
    out += struct.pack("<I", 0)
    for b in blobs:
        out += b
    for v in offs:
        out += struct.pack("<I", v)
    for v in bcs:
        out += struct.pack("<I", v)
    out += struct.pack("<ddd", sx, sy, 0.0)
    out += struct.pack("<dddddd", 0.0, 0.0, 0.0, tie_x, tie_y, 0.0)
    return bytes(out)


def _make_fake_tile(tmp_path, name, lon_min, lat_min, lon_max, lat_max,
                    px_deg=1 / 240, tile_px=256, seed=3):
    """Write a small fake GEBCO-like tile; return its path."""
    w = int(round((lon_max - lon_min) / px_deg))
    h = int(round((lat_max - lat_min) / px_deg))
    rng = np.random.default_rng(seed)
    vals = (rng.integers(-5000, 5000, size=(h, w))).astype(np.int16)
    blob = _build_tiled_geotiff(w, h, tile_px, tile_px, px_deg, -px_deg,
                                lon_min, lat_max, vals)
    p = tmp_path / name
    p.write_bytes(blob)
    return str(p), vals


def _build_shp(records, shape_type):
    """Minimal .shp bytes. records: list of list-of-rings (x, y)."""
    body = bytearray()
    for i, rings in enumerate(records, start=1):
        pts = [p for r in rings for p in r]
        parts = []
        k = 0
        for r in rings:
            parts.append(k)
            k += len(r)
        content = bytearray()
        content += struct.pack("<i", shape_type)
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        content += struct.pack("<4d", min(xs), min(ys), max(xs), max(ys))
        content += struct.pack("<ii", len(rings), len(pts))
        content += struct.pack(f"<{len(parts)}i", *parts)
        for x, y in pts:
            content += struct.pack("<dd", x, y)
        body += struct.pack(">ii", i, len(content) // 2) + content
    header = struct.pack(">i", 9994) + b"\x00" * 20
    header += struct.pack(">i", (100 + len(body)) // 2)
    header += struct.pack("<i", 1000) + struct.pack("<i", shape_type)
    header += struct.pack("<8d", -180, -90, 180, 90, 0, 0, 0, 0)
    return bytes(header) + bytes(body)


def _build_ne_zip(shp_bytes, stem="ne_110m_coastline"):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(stem + ".shp", shp_bytes)
        zf.writestr(stem + ".dbf", b"")
    return buf.getvalue()


@pytest.fixture()
def cache_env(tmp_path, monkeypatch):
    monkeypatch.setenv("SURVEY_CURRENTS_CACHE", str(tmp_path / "cache"))
    return tmp_path


# ---------------------------------------------------------------------------
# Tile layout
# ---------------------------------------------------------------------------

def test_tile_name_bounds_parse():
    assert bm._tile_bounds(
        "gebco_2024_n90.0_s0.0_w-180.0_e-90.0.tif") == (-180.0, 0.0, -90.0, 90.0)


def test_tile_name_regex_rejects():
    with pytest.raises(ValueError):
        bm._tile_bounds("gebco_2024_n90_s0.tif")


def test_all_tiles_cover_globe_no_overlap():
    tiles = bm._all_gebco_tiles()
    assert len(tiles) == 8
    for i, (n1, b1) in enumerate(tiles):
        for n2, b2 in tiles[i + 1:]:
            assert not bm._bbox_intersects(
                (b1[0] + 1e-9, b1[1] + 1e-9, b1[2] - 1e-9, b1[3] - 1e-9),
                b2), f"{n1} overlaps {n2}"
    # Union spans the globe.
    lons = sorted({b[0] for _, b in tiles} | {b[2] for _, b in tiles})
    lats = sorted({b[1] for _, b in tiles} | {b[3] for _, b in tiles})
    assert lons == [-180.0, -90.0, 0.0, 90.0, 180.0]
    assert lats == [-90.0, 0.0, 90.0]


def test_gebco_tiles_for_bbox_single():
    tiles = gebco_tiles_for_bbox((-80.0, 25.0, -70.0, 35.0))
    assert len(tiles) == 1
    assert tiles[0][0] == "gebco_2024_n90.0_s0.0_w-90.0_e0.0.tif"


def test_gebco_tiles_for_bbox_multi():
    tiles = gebco_tiles_for_bbox((-95.0, -5.0, 5.0, 5.0))
    assert len(tiles) == 6  # crosses lon -90, lon 0, and lat 0


def test_gebco_tiles_for_bbox_antimeridian():
    tiles = gebco_tiles_for_bbox((170.0, -10.0, -170.0, 10.0))
    names = {t[0] for t in tiles}
    assert names == {"gebco_2024_n90.0_s0.0_w90.0_e180.0.tif",
                     "gebco_2024_n90.0_s0.0_w-180.0_e-90.0.tif",
                     "gebco_2024_n0.0_s-90.0_w90.0_e180.0.tif",
                     "gebco_2024_n0.0_s-90.0_w-180.0_e-90.0.tif"}


# ---------------------------------------------------------------------------
# Ranged zip machinery (offline: monkeypatched HTTP)
# ---------------------------------------------------------------------------

def _serve(data):
    def fake_head(url, timeout=60):
        return "https://example.invalid/x.zip", len(data)

    def fake_range(url, start, end, timeout=120):
        assert 0 <= start <= end <= len(data)
        return data[start:end]

    return fake_head, fake_range


def test_zip_entries_roundtrip(monkeypatch):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("a.tif", b"0123456789" * 100)
        zf.writestr("b.tif", b"abcdefghij" * 50)
    data = buf.getvalue()
    fake_head, fake_range = _serve(data)
    monkeypatch.setattr(bm, "_http_head_size", fake_head)
    monkeypatch.setattr(bm, "_http_range", fake_range)
    final_url, entries = bm._zip_entries("https://example.invalid/x.zip")
    assert final_url == "https://example.invalid/x.zip"
    assert [e["name"] for e in entries] == ["a.tif", "b.tif"]
    assert all(e["method"] == 8 for e in entries)


def test_zip_entries_bad_signature(monkeypatch):
    fake_head, fake_range = _serve(b"not a zip at all" * 100)
    monkeypatch.setattr(bm, "_http_head_size", fake_head)
    monkeypatch.setattr(bm, "_http_range", fake_range)
    with pytest.raises(RuntimeError):
        bm._zip_entries("https://example.invalid/x.zip")


def test_extract_zip_entry_deflated(monkeypatch, tmp_path):
    payload = b"hello gebco" * 1000
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("tile.tif", payload)
    data = buf.getvalue()
    _, fake_range = _serve(data)
    monkeypatch.setattr(bm, "_http_range", fake_range)
    # Parse entries with the real central-directory reader.
    fake_head, _ = _serve(data)
    monkeypatch.setattr(bm, "_http_head_size", fake_head)
    _, entries = bm._zip_entries("https://example.invalid/x.zip")
    dest = str(tmp_path / "tile.tif")
    digest = bm._extract_zip_entry("https://example.invalid/x.zip",
                                   entries[0], dest)
    assert open(dest, "rb").read() == payload
    assert digest == hashlib.sha256(payload).hexdigest()
    assert bm._sidecar_ok(dest)


def test_extract_zip_entry_unsupported_method():
    with pytest.raises(RuntimeError):
        bm._extract_zip_entry("https://example.invalid/x.zip",
                              {"name": "x", "method": 99,
                               "compressed_size": 10,
                               "local_header_offset": 0},
                              "/tmp/never_written_xyz")


def test_ensure_gebco_tile_cache_hit_no_network(cache_env, tmp_path,
                                                monkeypatch):
    cache = tmp_path / "cache" / "basemap" / "gebco" / "gebco_2024"
    cache.mkdir(parents=True)
    name = "gebco_2024_n90.0_s0.0_w-90.0_e0.0.tif"
    p = cache / name
    p.write_bytes(b"fake-tile-bytes")
    digest = hashlib.sha256(b"fake-tile-bytes").hexdigest()
    (cache / (name + ".sha256")).write_text(digest + "\n")

    def _boom(*a, **k):
        raise AssertionError("network must not be touched on cache hit")

    monkeypatch.setattr(bm, "_zip_entries", _boom)
    path, got = ensure_gebco_tile(name)
    assert path == str(p)
    assert got == digest


def test_ensure_gebco_tile_corrupt_sidecar_redownloads(cache_env, tmp_path,
                                                       monkeypatch):
    cache = tmp_path / "cache" / "basemap" / "gebco" / "gebco_2024"
    cache.mkdir(parents=True)
    name = "gebco_2024_n90.0_s0.0_w-90.0_e0.0.tif"
    p = cache / name
    p.write_bytes(b"corrupted")
    (cache / (name + ".sha256")).write_text("0" * 64 + "\n")
    fresh = b"fresh-tile-bytes"
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as zf:
        zf.writestr(name, fresh)
    data = buf.getvalue()
    fake_head, fake_range = _serve(data)
    # Bypass _zip_entries with a canned index (the central-directory
    # reader itself is covered by test_zip_entries_roundtrip).
    entries = [{"name": name, "method": 0,
                "compressed_size": len(fresh),
                "local_header_offset": 0}]
    # Recompute the real local-header offset from the zip bytes.
    lho = data.find(fresh) - (30 + len(name))
    entries[0]["local_header_offset"] = lho
    monkeypatch.setattr(bm, "_zip_entries",
                        lambda url, timeout=600: ("https://example.invalid/x.zip",
                                                  entries))
    monkeypatch.setattr(bm, "_http_range", fake_range)
    path, got = ensure_gebco_tile(name)
    assert open(path, "rb").read() == fresh
    assert got == hashlib.sha256(fresh).hexdigest()


# ---------------------------------------------------------------------------
# TiledGeoTIFF reader
# ---------------------------------------------------------------------------

def test_tiff_reader_full_window(tmp_path):
    vals = (np.arange(48, dtype=np.int16).reshape(6, 8) * 7 - 100)
    blob = _build_tiled_geotiff(8, 6, 4, 3, 0.5, -0.5, -180.0, 90.0, vals)
    p = tmp_path / "t.tif"
    p.write_bytes(blob)
    tif = bm.TiledGeoTIFF(str(p))
    assert (tif.width, tif.height) == (8, 6)
    assert (tif.tile_w, tif.tile_h) == (4, 3)
    np.testing.assert_array_equal(tif.read_window(0, 0, 8, 6), vals)


def test_tiff_reader_partial_window_crossing_tiles(tmp_path):
    vals = (np.arange(64, dtype=np.int16).reshape(8, 8) * 3)
    blob = _build_tiled_geotiff(8, 8, 4, 4, 1.0, -1.0, 0.0, 8.0, vals)
    p = tmp_path / "t.tif"
    p.write_bytes(blob)
    tif = bm.TiledGeoTIFF(str(p))
    np.testing.assert_array_equal(tif.read_window(3, 2, 4, 5),
                                  vals[2:7, 3:7])


def test_tiff_reader_georeferencing(tmp_path):
    vals = np.zeros((4, 4), dtype=np.int16)
    blob = _build_tiled_geotiff(4, 4, 4, 4, 0.25, -0.25, 10.0, 20.0, vals)
    p = tmp_path / "t.tif"
    p.write_bytes(blob)
    tif = bm.TiledGeoTIFF(str(p))
    lon, lat = tif.pixel_to_geo(0.5, 0.5)
    assert (lon, lat) == pytest.approx((10.125, 19.875))
    c, r = tif.geo_to_pixel(10.125, 19.875)
    assert (c, r) == pytest.approx((0.5, 0.5))
    assert tif.cell_center_lons()[0] == pytest.approx(10.125)
    assert tif.cell_center_lats()[0] == pytest.approx(19.875)


def test_tiff_reader_bad_magic(tmp_path):
    p = tmp_path / "t.tif"
    p.write_bytes(b"NOPE" + b"\x00" * 100)
    with pytest.raises(ValueError):
        bm.TiledGeoTIFF(str(p))


def test_tiff_reader_uncompressed_rejected(tmp_path):
    # Compression=1 (uncompressed) is now SUPPORTED — the official GEBCO
    # 2024 tiles ship stripped and uncompressed (verified 2026-09-27).
    # An unsupported compression (e.g. 5 = LZW) is still rejected.
    vals = np.zeros((4, 4), dtype=np.int16)
    blob = bytearray(_build_tiled_geotiff(4, 4, 4, 4, 1.0, -1.0,
                                          0.0, 4.0, vals))
    # Patch Compression tag (259) from 8 to 5 (LZW) in the IFD.
    n = struct.unpack("<H", blob[8:10])[0]
    for i in range(n):
        off = 10 + 12 * i
        if struct.unpack("<H", blob[off:off + 2])[0] == 259:
            blob[off + 8:off + 10] = struct.pack("<H", 5)
    p = tmp_path / "t.tif"
    p.write_bytes(bytes(blob))
    with pytest.raises(ValueError, match="only uncompressed"):
        bm.TiledGeoTIFF(str(p))


def test_block_average_nan_aware():
    a = np.array([[1.0, np.nan], [3.0, 4.0]])
    out = _block_average(a, 2)
    assert out.shape == (1, 1)
    assert out[0, 0] == pytest.approx((1 + 3 + 4) / 3)
    out1 = _block_average(a, 1)
    assert out1.shape == (2, 2)


# ---------------------------------------------------------------------------
# TopoField
# ---------------------------------------------------------------------------

def test_topofield_synthetic_deterministic():
    f1 = TopoField.synthetic()
    f2 = TopoField.synthetic()
    np.testing.assert_array_equal(f1.values, f2.values)
    assert f1.values.shape[0] == 1  # static: single time axis
    assert len(f1.times) == 1
    assert np.nanmax(f1.values) > 0  # land exists
    assert np.nanmin(f1.values) < 0  # ocean exists
    assert f1.units == "m"


def test_topofield_json_roundtrip(tmp_path):
    f = TopoField.synthetic(bbox=(-10, 30, 10, 40), resolution=2.0)
    p = str(tmp_path / "topo.json")
    f.to_json(p)
    g = TopoField.from_json(p)
    np.testing.assert_array_equal(g.values, f.values)
    assert g.bbox == f.bbox
    assert g.resolution == f.resolution
    assert g.elevation.shape == f.values.shape[1:]


def test_topofield_shape_validation():
    with pytest.raises(ValueError):
        TopoField(times=[bm.GEBCO_STATIC_TIME], lats=[0.0], lons=[0.0],
                  values=np.zeros((2, 1, 1)), resolution=1.0)


# ---------------------------------------------------------------------------
# fetch_gebco (offline: monkeypatched tile supplier)
# ---------------------------------------------------------------------------

def _patch_tiles(monkeypatch, mapping):
    """mapping: tile name -> (path, vals)."""
    def fake_ensure(name, timeout=600):
        path, _ = mapping[name]
        return path, "deadbeef"
    monkeypatch.setattr(bm, "ensure_gebco_tile", fake_ensure)


def _patch_tile_index(monkeypatch, entries):
    """entries: list of (name, (lon_min, lat_min, lon_max, lat_max))."""
    monkeypatch.setattr(
        bm, "gebco_tiles_for_bbox", lambda bbox: list(entries))


def test_fetch_gebco_native_small(tmp_path, monkeypatch):
    path, vals = _make_fake_tile(
        tmp_path, "t.tif", 9.0, 9.0, 13.0, 13.0, px_deg=1 / 240)
    _patch_tiles(monkeypatch, {"t.tif": (path, vals)})
    _patch_tile_index(monkeypatch, [("t.tif", (9.0, 9.0, 13.0, 13.0))])
    f = fetch_gebco((10.0, 10.0, 11.0, 11.0), resolution="15s")
    assert f.source == f"gebco-{GEBCO_VERSION}"
    assert f.values.shape == (1, 240, 240)
    assert f.lats[0] > f.lats[-1]  # north -> south
    assert f.lons[0] < f.lons[-1]
    np.testing.assert_allclose(
        f.lons, 10.0 + (np.arange(240) + 0.5) / 240.0, rtol=1e-12)
    np.testing.assert_allclose(
        f.lats, 11.0 - (np.arange(240) + 0.5) / 240.0, rtol=1e-12)
    np.testing.assert_allclose(f.elevation, vals[480:720, 240:480])
    assert f.provenance["grid"] == f"GEBCO_{GEBCO_VERSION}"
    assert f.provenance["subset_bbox"] == [10.0, 10.0, 11.0, 11.0]
    assert f.provenance["static_compilation"] is True
    assert len(f.provenance["tiles"]) == 1


def test_fetch_gebco_block_averaged(tmp_path, monkeypatch):
    path, vals = _make_fake_tile(
        tmp_path, "t.tif", 9.0, 9.0, 13.0, 13.0, px_deg=1 / 240)
    _patch_tiles(monkeypatch, {"t.tif": (path, vals)})
    _patch_tile_index(monkeypatch, [("t.tif", (9.0, 9.0, 13.0, 13.0))])
    f = fetch_gebco((10.0, 10.0, 12.0, 12.0), resolution=0.25)
    assert f.values.shape == (1, 8, 8)
    assert f.resolution == pytest.approx(0.25)
    # Block means equal the reader-level block average of the window.
    tif = bm.TiledGeoTIFF(path)
    col0, row0, ncol, nrow = bm._native_window_for_bbox(
        tif, (10.0, 10.0, 12.0, 12.0))
    win = tif.read_window(col0, row0, ncol, nrow).astype(float)
    np.testing.assert_allclose(f.elevation, _block_average(win, 60),
                               rtol=1e-9)


def test_fetch_gebco_bad_resolutions(tmp_path, monkeypatch):
    path, vals = _make_fake_tile(tmp_path, "t.tif", 9.0, 9.0, 13.0, 13.0)
    _patch_tiles(monkeypatch, {"t.tif": (path, vals)})
    _patch_tile_index(monkeypatch, [("t.tif", (9.0, 9.0, 13.0, 13.0))])
    with pytest.raises(ValueError):
        fetch_gebco((10.0, 10.0, 11.0, 11.0), resolution=0.31)
    with pytest.raises(ValueError):
        fetch_gebco((10.0, 10.0, 11.0, 11.0), resolution="ultra-hd")
    with pytest.raises(ValueError):
        fetch_gebco((10.0, 10.0, 11.0, 11.0), resolution=0)
    # 7 * 15s is a multiple of native but does not divide the 90° tiles.
    with pytest.raises(ValueError):
        fetch_gebco((10.0, 10.0, 11.0, 11.0), resolution=7.0 / 240.0)


def test_fetch_gebco_native_cap(tmp_path, monkeypatch):
    path, vals = _make_fake_tile(tmp_path, "t.tif", 9.0, 9.0, 13.0, 13.0)
    _patch_tiles(monkeypatch, {"t.tif": (path, vals)})
    _patch_tile_index(monkeypatch, [("t.tif", (9.0, 9.0, 13.0, 13.0))])
    monkeypatch.setattr(bm, "GEBCO_MAX_NATIVE_WINDOW", 10)
    with pytest.raises(ValueError, match="safety cap"):
        fetch_gebco((10.0, 10.0, 11.0, 11.0), resolution="15s")


def test_fetch_gebco_antimeridian_mosaic(tmp_path, monkeypatch):
    p1, v1 = _make_fake_tile(tmp_path, "t1.tif", 178.0, 0.0, 180.0, 3.0,
                              seed=11)
    p2, v2 = _make_fake_tile(tmp_path, "t2.tif", -180.0, 0.0, -178.0, 3.0,
                              seed=12)
    _patch_tiles(monkeypatch, {"t1.tif": (p1, v1), "t2.tif": (p2, v2)})
    _patch_tile_index(monkeypatch, [
        ("t1.tif", (178.0, 0.0, 180.0, 3.0)),
        ("t2.tif", (-180.0, 0.0, -178.0, 3.0)),
    ])
    f = fetch_gebco((179.0, 1.0, -179.0, 2.0), resolution=0.25)
    assert f.values.shape == (1, 4, 8)
    assert len(f.provenance["tiles"]) == 2
    # West strip (179..180) then wrapped east strip (-180..-179).
    np.testing.assert_allclose(
        f.lats, [1.875, 1.625, 1.375, 1.125], rtol=1e-12)
    np.testing.assert_allclose(
        f.lons,
        [179.125, 179.375, 179.625, 179.875,
         -179.875, -179.625, -179.375, -179.125], rtol=1e-12)
    np.testing.assert_allclose(
        f.elevation[:, :4],
        _block_average(v1[240:480, 240:480].astype(float), 60), rtol=1e-9)
    np.testing.assert_allclose(
        f.elevation[:, 4:],
        _block_average(v2[240:480, 0:240].astype(float), 60), rtol=1e-9)


# ---------------------------------------------------------------------------
# Natural Earth
# ---------------------------------------------------------------------------

def test_ne_url_construction():
    assert bm._ne_zip_url("110m", "coastline") == (
        "https://naturalearth.s3.amazonaws.com/110m_physical/"
        "ne_110m_coastline.zip")
    assert bm._ne_zip_url("50m", "countries") == (
        "https://naturalearth.s3.amazonaws.com/50m_cultural/"
        "ne_50m_admin_0_countries.zip")


def test_ne_invalid_args():
    with pytest.raises(ValueError):
        ensure_naturalearth_zip("30m", "coastline")
    with pytest.raises(ValueError):
        fetch_naturalearth((-10, 30, 10, 40), layers=("roads",))


def test_read_shp_polyline():
    shp = _build_shp([[[(0.0, 0.0), (1.0, 1.0), (2.0, 0.5)]]], 3)
    recs = _read_shp_parts(shp, allowed=(3,))
    assert len(recs) == 1
    st, rings = recs[0]
    assert st == 3
    assert rings[0] == [(0.0, 0.0), (1.0, 1.0), (2.0, 0.5)]


def test_read_shp_bad_header():
    with pytest.raises(ValueError):
        _read_shp_parts(b"\x00" * 200, allowed=(3,))


def test_rings_to_polygons_hole_assignment():
    exterior = [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0),
                (0.0, 0.0)]  # CW in x-right/y-up -> negative area
    assert bm._signed_area(exterior) < 0
    hole = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0),
            (3.0, 3.0)]  # CCW -> positive area
    assert bm._signed_area(hole) > 0
    polys = _rings_to_polygons([exterior, hole])
    assert len(polys) == 1
    assert len(polys[0]) == 2


def test_clip_polyline_bbox():
    line = [(-5.0, 0.0), (5.0, 0.0), (15.0, 0.0)]
    segs = _clip_polyline_bbox(line, (0.0, -1.0, 10.0, 1.0))
    assert len(segs) == 1
    assert segs[0][0][0] == pytest.approx(0.0)
    assert segs[0][-1][0] == pytest.approx(10.0)
    # Fully outside -> no segments.
    assert _clip_polyline_bbox([(20.0, 0.0), (30.0, 0.0)],
                               (0.0, -1.0, 10.0, 1.0)) == []


def test_fetch_naturalearth_offline(cache_env, monkeypatch):
    line = [[(-20.0, 35.0), (0.0, 36.0), (20.0, 35.0)]]
    shp = _build_shp([line], 3)
    payload = _build_ne_zip(shp)

    def fake_ensure(scale, layer, timeout=120):
        path = bm._ne_cache_path(scale, layer)
        path = str(path)
        with open(path, "wb") as fh:
            fh.write(payload)
        digest = hashlib.sha256(payload).hexdigest()
        with open(path + ".sha256", "w") as fh:
            fh.write(digest + "\n")
        return path, digest

    monkeypatch.setattr(bm, "ensure_naturalearth_zip", fake_ensure)
    out = fetch_naturalearth((-10.0, 30.0, 10.0, 40.0), scale="110m",
                             layers=("coastline",))
    fc = out["coastline"]
    assert fc["type"] == "FeatureCollection"
    assert len(fc["features"]) >= 1
    assert all(f["geometry"]["type"] == "LineString"
               for f in fc["features"])
    xs = [x for f in fc["features"]
          for x, _ in f["geometry"]["coordinates"]]
    assert min(xs) >= -10.0 and max(xs) <= 10.0  # clipped to bbox
    assert fc["properties"]["scale"] == "110m"


def test_ensure_naturalearth_zip_cache_hit(cache_env, monkeypatch):
    path = bm._ne_cache_path("110m", "coastline")
    payload = b"PK fake"
    with open(path, "wb") as fh:
        fh.write(payload)
    digest = hashlib.sha256(payload).hexdigest()
    with open(path + ".sha256", "w") as fh:
        fh.write(digest + "\n")

    def _boom(*a, **k):
        raise AssertionError("network must not be touched on cache hit")

    monkeypatch.setattr(bm, "_http_get", _boom)
    got_path, got_digest = ensure_naturalearth_zip("110m", "coastline")
    assert got_path == path and got_digest == digest


def test_synthetic_naturalearth():
    out = synthetic_naturalearth()
    assert set(out) == {"coastline", "countries"}
    assert out["coastline"]["features"][0]["geometry"]["type"] == \
        "LineString"
    assert out["countries"]["features"][0]["geometry"]["type"] == "Polygon"
    out2 = synthetic_naturalearth()
    assert out == out2
    with pytest.raises(ValueError):
        synthetic_naturalearth(layers=("lakes",))


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def test_basemap_cache_dir_env_override(cache_env):
    d = basemap_cache_dir()
    assert d == str(cache_env / "cache" / "basemap")
    assert os.path.isdir(d)


def test_sidecar_ok_missing_and_mismatch(tmp_path):
    p = tmp_path / "x.bin"
    assert not bm._sidecar_ok(str(p))
    p.write_bytes(b"data")
    (tmp_path / "x.bin.sha256").write_text("0" * 64 + "\n")
    assert not bm._sidecar_ok(str(p))
    bm._write_sidecar(str(p))
    assert bm._sidecar_ok(str(p))


def test_gebco_version_pinned():
    assert GEBCO_VERSION == "2024"


# ---------------------------------------------------------------------------
# Official GEBCO 2024 tile (recorded 2026-09-27; offline, cache-only)
# ---------------------------------------------------------------------------

_REAL_TILE = os.path.join(
    bm.basemap_cache_dir(), "gebco", "gebco_2024",
    "gebco_2024_n90.0_s0.0_w-180.0_e-90.0.tif")
_needs_real_tile = pytest.mark.skipif(
    not os.path.exists(_REAL_TILE),
    reason="official GEBCO tile not in the local cache (offline-safe skip)")


@_needs_real_tile
class TestRealGebcoTile:
    """Regression tests against the official 2024 tile entry.

    The tile was range-extracted from the CEDA-hosted
    ``gebco_2024_geotiff.zip`` on 2026-09-27. These tests never touch
    the network: they read the cached tile in place.
    """

    def test_dimensions(self):
        with bm.TiledGeoTIFF(_REAL_TILE) as t:
            assert (t.width, t.height) == (21600, 21600)

    def test_stripped_uncompressed_layout(self):
        with bm.TiledGeoTIFF(_REAL_TILE) as t:
            assert t._tiled is False
            assert t._compression == 1
            assert t.tile_h == 1          # one row per strip
            assert t._nby == 21600

    def test_georeferencing_convention(self):
        with bm.TiledGeoTIFF(_REAL_TILE) as t:
            # Positive Y pixel scale per the GeoTIFF spec; rows run
            # north -> south from the top-left tiepoint.
            assert t._sx == pytest.approx(1 / 240)
            assert t._sy == pytest.approx(-1 / 240)
            assert t.pixel_to_geo(0, 0) == pytest.approx((-180.0, 90.0))
            lon, lat = t.pixel_to_geo(t.width, t.height)
            assert (lon, lat) == pytest.approx((-90.0, 0.0))

    def test_window_read_plausible_depths(self):
        with bm.TiledGeoTIFF(_REAL_TILE) as t:
            w = t.read_window(0, 0, 8, 4)
            assert w.shape == (4, 8) and w.dtype == np.int16
            # NE Pacific abyssal plain: a few km deep, negative.
            assert -6000 < w.mean() < -2000

    def test_file_backed_memory(self):
        # Reading an 8x4 window must not materialize the ~933 MB tile.
        import tracemalloc
        tracemalloc.start()
        with bm.TiledGeoTIFF(_REAL_TILE) as t:
            before = tracemalloc.get_traced_memory()[0]
            t.read_window(10800, 10800, 64, 64)
            after = tracemalloc.get_traced_memory()[0]
        tracemalloc.stop()
        assert after - before < 50 * 1024 * 1024

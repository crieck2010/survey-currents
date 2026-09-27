"""Tests for the NSIDC G02135 v4.0 sea-ice adapter (currents.sea_ice).

Fully offline: HTTP is mocked, and concentration GeoTIFFs are built
in-memory with a minimal stdlib TIFF writer (the real files are
uncompressed single-band 16-bit TIFFs, so this matches the wire
format byte-for-byte).
"""

import datetime as dt
import hashlib
import json
import math
import struct

import numpy as np
import pytest

from currents.sea_ice import (
    HUGHES_A, HUGHES_E2, NSIDC_BASE, NSIDC_CELL_M, NSIDC_GRIDS, NSIDC_START,
    IceField, decode_concentration, fetch_nsidc_sic, nsidc_conc_url,
    pick_hemisphere, ps_forward, read_concentration_geotiff,
    reproject_to_latlon, resolve_native_grid, target_grid, validate_ice_bbox)


# ---------------------------------------------------------------------------
# Minimal uncompressed 16-bit single-band TIFF writer (test fixture)
# ---------------------------------------------------------------------------

def make_tiff(raw: np.ndarray, x0: float | None = None,
              y0_top: float | None = None, dx: float = 25000.0) -> bytes:
    """Encode ``raw`` (uint16, (ny, nx)) as a minimal little-endian TIFF.

    Includes the strip tags the reader needs; georeferencing tags only
    when ``x0`` is given.
    """
    raw = np.asarray(raw, dtype="<u2")
    ny, nx = raw.shape
    strips: list[tuple[bytes, int]] = []
    rows_per_strip = 13
    pix = raw.tobytes()
    for r0 in range(0, ny, rows_per_strip):
        strips.append((pix[r0 * nx * 2:(r0 + rows_per_strip) * nx * 2], r0))
    buf = bytearray(b"II*\x00")
    buf += struct.pack("<I", 8)

    tags: list[tuple[int, int, int, object]] = [
        (256, 3, 1, nx), (257, 3, 1, ny), (258, 3, 1, 16),
        (259, 3, 1, 1), (277, 3, 1, 1), (278, 3, 1, rows_per_strip),
    ]
    # strip payloads come right after the IFD; IFD size is known
    n_tags = len(tags) + 2 + (2 if x0 is not None else 0)
    ifd_size = 2 + n_tags * 12 + 4
    data_off = 8 + ifd_size
    offs = []
    for (sbytes, _r0) in strips:
        offs.append(data_off)
        data_off += len(sbytes)
    geo_off = None
    if x0 is not None:
        geo_off = data_off
        data_off += 24 + 48
    tags += [(273, 4, len(offs), None), (279, 3, len(offs), None)]
    off_arr_off = data_off
    data_off += 4 * len(offs)
    cnt_arr_off = data_off
    data_off += 2 * len(offs)
    if x0 is not None:
        tags += [(33550, 12, 3, None), (33922, 12, 6, None)]

    # lay out
    tags = sorted(tags, key=lambda t: t[0])
    buf += struct.pack("<H", len(tags))
    extra = bytearray()
    for tag, typ, cnt, val in tags:
        if tag == 273:
            buf += struct.pack("<HHI4s", tag, typ, cnt, struct.pack("<I", off_arr_off))
        elif tag == 279:
            buf += struct.pack("<HHI4s", tag, typ, cnt, struct.pack("<I", cnt_arr_off))
        elif tag == 33550:
            buf += struct.pack("<HHI4s", tag, typ, cnt, struct.pack("<I", geo_off))
        elif tag == 33922:
            buf += struct.pack("<HHI4s", tag, typ, cnt, struct.pack("<I", geo_off + 24))
        else:
            buf += struct.pack("<HHI4s", tag, typ, cnt, struct.pack("<H", int(val)) + b"\x00\x00")
    buf += struct.pack("<I", 0)
    for (sbytes, _r0) in strips:
        buf += sbytes
    if x0 is not None:
        buf += struct.pack("<3d", dx, dx, 0.0)
        buf += struct.pack("<6d", 0.0, 0.0, 0.0, x0, y0_top, 0.0)
    buf += struct.pack("<%dI" % len(offs), *offs)
    buf += struct.pack("<%dH" % len(offs), *[len(s[0]) for s in strips])
    return bytes(buf)


def native_grid(hemisphere: str = "north") -> dict:
    g = NSIDC_GRIDS[hemisphere]
    return {"nx": g["nx"], "ny": g["ny"], "x0": g["x0"],
            "y0_top": g["y0_top"], "dx": g["dx"], "dy": g["dx"]}


# ---------------------------------------------------------------------------
# URL builder
# ---------------------------------------------------------------------------

def test_nsidc_conc_url_north():
    url = nsidc_conc_url(dt.date(2026, 2, 1), "north")
    assert url == ("https://noaadata.apps.nsidc.org/NOAA/G02135/north/daily/"
                   "geotiff/2026/02_Feb/N_20260201_concentration_v4.0.tif")


def test_nsidc_conc_url_south():
    url = nsidc_conc_url(dt.date(2026, 9, 15), "south")
    assert url == ("https://noaadata.apps.nsidc.org/NOAA/G02135/south/daily/"
                   "geotiff/2026/09_Sep/S_20260915_concentration_v4.0.tif")


def test_nsidc_conc_url_smmr_era():
    url = nsidc_conc_url(dt.date(1978, 11, 2), "north")
    assert "1978/11_Nov/N_19781102_concentration_v4.0.tif" in url


def test_nsidc_conc_url_bad_hemisphere():
    with pytest.raises(ValueError):
        nsidc_conc_url(dt.date(2026, 1, 1), "equator")


# ---------------------------------------------------------------------------
# pick_hemisphere / validate_ice_bbox
# ---------------------------------------------------------------------------

def test_pick_hemisphere_auto_north():
    assert pick_hemisphere((-180.0, 66.0, 180.0, 90.0)) == "north"


def test_pick_hemisphere_auto_south():
    assert pick_hemisphere((-180.0, -75.0, 180.0, -55.0)) == "south"


def test_pick_hemisphere_equator_crossing_raises():
    with pytest.raises(ValueError, match="equator"):
        pick_hemisphere((-10.0, -5.0, 10.0, 5.0))


def test_pick_hemisphere_explicit():
    assert pick_hemisphere((-10.0, -5.0, 10.0, 5.0), "south") == "south"
    assert pick_hemisphere((-180.0, 66.0, 180.0, 90.0), "NORTH") == "north"


def test_pick_hemisphere_bad_value():
    with pytest.raises(ValueError):
        pick_hemisphere((-180.0, 66.0, 180.0, 90.0), "sideways")


def test_validate_ice_bbox_antimeridian_allowed():
    assert validate_ice_bbox((170.0, 66.0, -170.0, 80.0)) == (170.0, 66.0, -170.0, 80.0)


def test_validate_ice_bbox_bad_lat():
    with pytest.raises(ValueError):
        validate_ice_bbox((-180.0, 80.0, 180.0, 66.0))


# ---------------------------------------------------------------------------
# decode_concentration
# ---------------------------------------------------------------------------

def test_decode_concentration_scale():
    out = decode_concentration(np.array([[0, 500, 803, 1000]], dtype=np.uint16))
    assert out.tolist() == [[0.0, 50.0, 80.3, 100.0]]


def test_decode_concentration_flags_nan():
    out = decode_concentration(np.array([[2510, 2530, 2540, 2550]], dtype=np.uint16))
    assert np.isnan(out).all()


def test_decode_concentration_shape_and_dtype():
    raw = np.zeros((448, 304), dtype=np.uint16)
    out = decode_concentration(raw)
    assert out.shape == (448, 304) and out.dtype == float


# ---------------------------------------------------------------------------
# read_concentration_geotiff
# ---------------------------------------------------------------------------

def test_read_tiff_roundtrip():
    raw = np.arange(12, dtype=np.uint16).reshape(3, 4)
    payload = make_tiff(raw, x0=-3850000.0, y0_top=5850000.0)
    back, grid = read_concentration_geotiff(payload)
    assert (back == raw).all()
    assert grid["x0"] == -3850000.0 and grid["y0_top"] == 5850000.0
    assert grid["dx"] == 25000.0 and grid["nx"] == 4 and grid["ny"] == 3


def test_read_tiff_no_geotags():
    raw = np.zeros((5, 6), dtype=np.uint16)
    _back, grid = read_concentration_geotiff(make_tiff(raw))
    assert grid["x0"] is None  # resolve_native_grid fills these in


def test_read_tiff_rejects_non_tiff():
    with pytest.raises(ValueError, match="TIFF"):
        read_concentration_geotiff(b"<html>error</html>")


def test_read_tiff_rejects_wrong_dims():
    # File with no georeferencing tags whose dims don't match the north
    # grid constants -> resolve_native_grid must refuse.
    raw = np.zeros((5, 6), dtype=np.uint16)
    _back, grid = read_concentration_geotiff(make_tiff(raw))
    with pytest.raises(ValueError, match="no georeferencing"):
        resolve_native_grid(grid, "north")


# ---------------------------------------------------------------------------
# resolve_native_grid
# ---------------------------------------------------------------------------

def test_resolve_native_grid_fills_from_constants():
    grid = {"nx": 304.0, "ny": 448.0, "x0": None, "y0_top": None,
            "dx": None, "dy": None}
    out = resolve_native_grid(grid, "north")
    assert out["x0"] == NSIDC_GRIDS["north"]["x0"]
    assert out["dx"] == NSIDC_CELL_M


def test_resolve_native_grid_keeps_file_tags():
    grid = {"nx": 4.0, "ny": 3.0, "x0": 1.0, "y0_top": 2.0,
            "dx": 5.0, "dy": 5.0}
    out = resolve_native_grid(grid, "north")
    assert out["x0"] == 1.0  # file tags win


# ---------------------------------------------------------------------------
# ps_forward
# ---------------------------------------------------------------------------

def test_ps_forward_poles():
    x, y = ps_forward(90.0, 0.0, "north")
    assert (x, y) == pytest.approx((0.0, 0.0))
    x, y = ps_forward(-90.0, 0.0, "south")
    assert (x, y) == pytest.approx((0.0, 0.0))


def test_ps_forward_north_70_at_center():
    # At (70N, -45) the point lies due south of the north pole on the
    # central meridian -> x = 0, y = -a*m1.
    x, y = ps_forward(70.0, -45.0, "north")
    phi1 = math.radians(70.0)
    e = math.sqrt(HUGHES_E2)
    m1 = math.cos(phi1) / math.sqrt(1 - HUGHES_E2 * math.sin(phi1) ** 2)
    assert x == pytest.approx(0.0, abs=1e-6)
    assert y == pytest.approx(-HUGHES_A * m1, rel=1e-9)


def test_ps_forward_south_sign():
    x, y = ps_forward(-70.0, 0.0, "south")
    phi1 = math.radians(70.0)
    e = math.sqrt(HUGHES_E2)
    m1 = math.cos(phi1) / math.sqrt(1 - HUGHES_E2 * math.sin(phi1) ** 2)
    assert x == pytest.approx(0.0, abs=1e-6)
    assert y == pytest.approx(HUGHES_A * m1, rel=1e-9)  # north-up map


def test_ps_forward_bad_hemisphere():
    with pytest.raises(ValueError):
        ps_forward(80.0, 0.0, "east")


# ---------------------------------------------------------------------------
# reproject_to_latlon
# ---------------------------------------------------------------------------

def _ice_blob_raw(nx: int = 304, ny: int = 448) -> np.ndarray:
    """Full-resolution north grid: 100% ice everywhere except open ocean
    (0) in rows 0..99, land (2540) at pixel (0, 0)."""
    raw = np.full((ny, nx), 1000, dtype=np.uint16)
    raw[:100, :] = 0
    raw[0, 0] = 2540
    return raw


def test_reproject_shape():
    lats, lons = target_grid((-180.0, 66.0, 180.0, 90.0), 2.0)
    out = reproject_to_latlon(_ice_blob_raw(), native_grid("north"),
                              lats, lons, "north")
    assert out.shape == (len(lats), len(lons))


def test_reproject_nearest_neighbor_edge_preserved():
    lats, lons = target_grid((-180.0, 66.0, 180.0, 90.0), 1.0)
    out = reproject_to_latlon(_ice_blob_raw(), native_grid("north"),
                              lats, lons, "north")
    # synthetic boundary between row 99 and row 100 maps near 67.7N;
    # every cell must be exactly 0 or 100 (nearest-neighbor, no blur)
    unique = np.unique(out[~np.isnan(out)])
    assert set(unique.tolist()) <= {0.0, 100.0}


def test_reproject_flags_become_nan():
    raw = _ice_blob_raw()
    # Flag a block straddling the pole pixel (154, 234) -> inside the
    # arctic target domain, so its target cells must come out NaN.
    raw[214:254, 134:174] = 2540
    lats, lons = target_grid((-180.0, 66.0, 180.0, 90.0), 4.0)
    out = reproject_to_latlon(raw, native_grid("north"), lats, lons, "north")
    assert np.isnan(out).any()
    iy = int(np.argmin(np.abs(lats - 88.0)))
    ix = int(np.argmin(np.abs(lons - (-45.0))))
    assert np.isnan(out[iy, ix])


def test_reproject_outside_grid_is_nan():
    lats, lons = target_grid((-180.0, -90.0, 180.0, -30.0), 4.0)
    out = reproject_to_latlon(_ice_blob_raw(), native_grid("north"),
                              lats, lons, "north")
    assert np.isnan(out).all()  # southern bbox on the north grid


# ---------------------------------------------------------------------------
# target_grid
# ---------------------------------------------------------------------------

def test_target_grid_centers():
    lats, lons = target_grid((-10.0, 60.0, 10.0, 70.0), 5.0)
    assert lats.tolist() == [62.5, 67.5]
    assert lons.tolist() == [-7.5, -2.5, 2.5, 7.5]


def test_target_grid_antimeridian_wrap():
    lats, lons = target_grid((170.0, 66.0, -170.0, 80.0), 10.0)
    assert lons[0] == pytest.approx(175.0)
    assert lons[-1] == pytest.approx(185.0)


def test_target_grid_bad_resolution():
    with pytest.raises(ValueError):
        target_grid((-180.0, 66.0, 180.0, 90.0), 0.0)


# ---------------------------------------------------------------------------
# IceField
# ---------------------------------------------------------------------------

def test_icefield_synthetic_shape_and_determinism():
    a = IceField.synthetic()
    b = IceField.synthetic()
    assert a.shape == (5, 24, 360)
    assert np.array_equal(a.values, b.values, equal_nan=True)
    assert a.hemisphere == "north"
    c = IceField.synthetic(bbox=(-180.0, -75.0, 180.0, -55.0))
    assert c.hemisphere == "south"


def test_icefield_synthetic_values_in_range():
    f = IceField.synthetic()
    valid = f.values[~np.isnan(f.values)]
    assert valid.min() >= 0.0 and valid.max() <= 100.0


def test_icefield_select_time():
    f = IceField.synthetic(start="2024-02-01", end="2024-02-10")
    sub = f.select_time("2024-02-03", "2024-02-05")
    assert len(sub) == 3
    assert sub.times[0].date() == dt.date(2024, 2, 3)


def test_icefield_select_bbox():
    f = IceField.synthetic()
    sub = f.select_bbox((-100.0, 70.0, -80.0, 80.0))
    assert sub.lons.min() >= -100.0 and sub.lons.max() <= -80.0
    assert sub.lats.min() >= 70.0 and sub.lats.max() <= 80.0
    with pytest.raises(ValueError):
        f.select_bbox((-10.0, -80.0, 10.0, -70.0))


def test_icefield_json_roundtrip(tmp_path):
    f = IceField.synthetic()
    p = str(tmp_path / "ice.json")
    assert f.to_json(p) == p
    back = IceField.from_json(p)
    assert back.hemisphere == f.hemisphere
    assert np.array_equal(back.values, f.values, equal_nan=True)
    assert back.times == f.times


def test_icefield_time_range_and_bounds():
    f = IceField.synthetic()
    t0, t1 = f.time_range
    assert t0.date() == dt.date(2024, 2, 1) and t1.date() == dt.date(2024, 2, 5)
    assert f.bounds == (-180.0, 66.0, 180.0, 90.0)


# ---------------------------------------------------------------------------
# fetch_nsidc_sic (HTTP mocked)
# ---------------------------------------------------------------------------

class _FakeDownload:
    def __init__(self, payload: bytes):
        self.payload = payload
        self.urls: list[str] = []

    def __call__(self, url: str, timeout: int = 300) -> bytes:
        self.urls.append(url)
        if "missing" in url:
            raise RuntimeError("HTTP 404")
        return self.payload


def test_fetch_urls_and_provenance(monkeypatch):
    import currents.sea_ice as sea_ice
    raw = _ice_blob_raw()
    payload = make_tiff(raw, x0=-3850000.0, y0_top=5850000.0)
    fake = _FakeDownload(payload)
    monkeypatch.setattr(sea_ice, "_download_bytes", fake)
    field = fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0),
                            "2026-02-01", "2026-02-03", resolution=4.0)
    assert len(field) == 3
    assert all(url.startswith(NSIDC_BASE) for url in fake.urls)
    assert "N_20260201_concentration_v4.0.tif" in fake.urls[0]
    assert "N_20260203_concentration_v4.0.tif" in fake.urls[2]
    prov = field.provenance
    assert prov["n_files"] == 3
    assert prov["hemisphere"] == "north"
    assert prov["authentication"] == "none (public keyless HTTPS)"
    assert prov["files"][0]["sha256"] == hashlib.sha256(payload).hexdigest()
    assert "nearest-neighbor" in prov["reprojection"]


def test_fetch_skips_missing_days(monkeypatch):
    import currents.sea_ice as sea_ice
    raw = _ice_blob_raw()
    payload = make_tiff(raw, x0=-3850000.0, y0_top=5850000.0)
    monkeypatch.setattr(sea_ice, "_download_bytes",
                        lambda url, timeout=300: (_ for _ in ()).throw(
                            RuntimeError("HTTP 404")) if "0202" in url else payload)
    field = fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0),
                            "2026-02-01", "2026-02-03", resolution=4.0)
    assert len(field) == 2
    assert field.times[0].date() == dt.date(2026, 2, 1)
    assert field.times[1].date() == dt.date(2026, 2, 3)
    assert len(field.provenance["skipped_days"]) == 1


def test_fetch_stride_days(monkeypatch):
    import currents.sea_ice as sea_ice
    payload = make_tiff(_ice_blob_raw(), x0=-3850000.0, y0_top=5850000.0)
    fake = _FakeDownload(payload)
    monkeypatch.setattr(sea_ice, "_download_bytes", fake)
    field = fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0),
                            "2026-02-01", "2026-02-07",
                            stride_days=3, resolution=4.0)
    assert [u[-31:-23] for u in fake.urls] == ["20260201", "20260204", "20260207"]
    assert len(field) == 3


def test_fetch_southern_hemisphere_url(monkeypatch):
    import currents.sea_ice as sea_ice
    raw = np.full((332, 316), 1000, dtype=np.uint16)
    payload = make_tiff(raw, x0=-3950000.0, y0_top=4350000.0)
    fake = _FakeDownload(payload)
    monkeypatch.setattr(sea_ice, "_download_bytes", fake)
    field = fetch_nsidc_sic((-180.0, -75.0, 180.0, -55.0),
                            "2026-02-01", "2026-02-01", resolution=4.0)
    assert field.hemisphere == "south"
    assert "/south/daily/geotiff/" in fake.urls[0]
    assert "S_20260201_concentration_v4.0.tif" in fake.urls[0]
    assert field.values.shape == (1, 5, 90)


def test_fetch_all_missing_raises(monkeypatch):
    import currents.sea_ice as sea_ice
    monkeypatch.setattr(
        sea_ice, "_download_bytes",
        lambda url, timeout=300: (_ for _ in ()).throw(RuntimeError("HTTP 404")))
    with pytest.raises(ValueError, match="no NSIDC sea-ice files"):
        fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0),
                        "2026-02-01", "2026-02-02")


def test_fetch_bad_dates_and_range():
    with pytest.raises(ValueError, match="starts 1978-11-01"):
        fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0), "1970-01-01", "1970-01-02")
    with pytest.raises(ValueError, match="before start"):
        fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0), "2026-02-02", "2026-02-01")


def test_fetch_accepts_datetime_and_str_forms(monkeypatch):
    import currents.sea_ice as sea_ice
    payload = make_tiff(_ice_blob_raw(), x0=-3850000.0, y0_top=5850000.0)
    fake = _FakeDownload(payload)
    monkeypatch.setattr(sea_ice, "_download_bytes", fake)
    field = fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0),
                            dt.datetime(2026, 2, 1, 12, 0),
                            dt.date(2026, 2, 2), resolution=4.0)
    assert len(field) == 2


def test_fetch_values_are_decoded_percent(monkeypatch):
    import currents.sea_ice as sea_ice
    raw = np.full((448, 304), 803, dtype=np.uint16)
    payload = make_tiff(raw, x0=-3850000.0, y0_top=5850000.0)
    monkeypatch.setattr(sea_ice, "_download_bytes", _FakeDownload(payload))
    field = fetch_nsidc_sic((-180.0, 66.0, 180.0, 90.0),
                            "2026-02-01", "2026-02-01", resolution=4.0)
    valid = field.values[~np.isnan(field.values)]
    assert np.allclose(valid, 80.3)

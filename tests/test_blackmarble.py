"""Tests for the NASA Black Marble night-lights adapter (currents.blackmarble).

Fully offline: CMR discovery and downloads are mocked, HDF5 fixtures
are built with h5py in-memory, and LightsField.synthetic() covers the
field-model tests.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import io
import json

import numpy as np
import pytest

import currents.blackmarble as bm
from currents.blackmarble import (
    CredentialsMissing,
    LightsField,
    _block_average_to_grid,
    _parse_blackmarble_bytes,
    bm_sample_days,
    cmr_granule_query,
    discover_granule_url,
    fetch_blackmarble,
    normalize_bm_product,
    tile_bounds,
    tile_for_lon_lat,
    tile_token,
    tiles_for_bbox,
)


def _h5py():
    return pytest.importorskip("h5py")


# ---------------------------------------------------------------------------
# product normalization
# ---------------------------------------------------------------------------


def test_normalize_product_ok():
    assert normalize_bm_product("daily") == "daily"
    assert normalize_bm_product(" Daily ") == "daily"
    assert normalize_bm_product("DAILY") == "daily"


def test_normalize_product_bad():
    with pytest.raises(ValueError, match="only 'daily'"):
        normalize_bm_product("monthly")
    with pytest.raises(ValueError, match="only 'daily'"):
        normalize_bm_product("nrt")


# ---------------------------------------------------------------------------
# tile math (verified against live CMR granule footprints 2026-09-26:
# h07v10 = lon -110..-100, lat -20..-10)
# ---------------------------------------------------------------------------


def test_tile_for_lon_lat_known():
    assert tile_for_lon_lat(-105.0, -15.0) == (7, 10)
    assert tile_for_lon_lat(-110.0, -10.0) == (7, 10)


def test_tile_for_lon_lat_corners():
    assert tile_for_lon_lat(-180.0, 90.0) == (0, 0)
    assert tile_for_lon_lat(179.9, -89.9) == (35, 17)
    assert tile_for_lon_lat(0.0, 0.0) == (18, 9)


def test_tile_bounds_roundtrip():
    lon_min, lat_min, lon_max, lat_max = tile_bounds(7, 10)
    assert (lon_min, lat_min, lon_max, lat_max) == (-110.0, -20.0, -100.0, -10.0)
    assert tile_for_lon_lat(lon_min + 0.1, lat_max - 0.1) == (7, 10)


def test_tile_bounds_bad():
    with pytest.raises(ValueError, match="out of range"):
        tile_bounds(36, 0)
    with pytest.raises(ValueError, match="out of range"):
        tile_bounds(0, 18)


def test_tile_token():
    assert tile_token(7, 10) == "h07v10"
    assert tile_token(0, 0) == "h00v00"
    assert tile_token(35, 17) == "h35v17"


def test_tiles_for_bbox_single_tile():
    assert tiles_for_bbox((-110.0, -20.0, -100.0, -10.0)) == [(7, 10)]


def test_tiles_for_bbox_us():
    tiles = tiles_for_bbox((-125.0, 25.0, -66.0, 49.0))
    assert len(tiles) == 21
    assert tiles == sorted(tiles)


def test_tiles_for_bbox_antimeridian():
    tiles = tiles_for_bbox((170.0, 20.0, -170.0, 40.0))
    hs = {h for h, _ in tiles}
    assert 35 in hs and 0 in hs  # wraps the date line


def test_tiles_for_bbox_global():
    assert len(tiles_for_bbox((-180.0, -90.0, 180.0, 90.0))) == 36 * 18


# ---------------------------------------------------------------------------
# CMR discovery (keyless) — URL shape and token matching
# ---------------------------------------------------------------------------


def test_cmr_granule_query():
    day = dt.date(2024, 6, 1)
    url = cmr_granule_query(day, 7, 10)
    assert url.startswith(bm.CMR_GRANULES_URL)
    assert "short_name=VNP46A2" in url
    assert "version=2" in url
    assert "2024-06-01T00%3A00%3A00Z" in url or "2024-06-01" in url
    # tile footprint passed as the search bounding box
    assert "-110" in url and "-20" in url


def _fake_cmr(entries):
    return {"feed": {"entry": entries}}


def _cmr_entry(href):
    return {"links": [{"rel": "http://esipfed.org/ns/fedsearch/1.1/data#",
                       "href": href},
                      {"rel": "via",
                       "href": "https://example.invalid/browse"}]}


def test_discover_granule_url_match(monkeypatch):
    day = dt.date(2024, 6, 1)
    good = ("https://data.laadsdaac.earthdatacloud.nasa.gov/prod-lads/"
            "VNP46A2/VNP46A2.A2024153.h07v10.002.2024153120000.h5")
    other = ("https://data.laadsdaac.earthdatacloud.nasa.gov/prod-lads/"
             "VNP46A2/VNP46A2.A2024153.h07v11.002.2024153120000.h5")
    monkeypatch.setattr(
        bm, "_cmr_get_json",
        lambda url, timeout=60: _fake_cmr([_cmr_entry(other), _cmr_entry(good)]))
    assert discover_granule_url(day, 7, 10) == good


def test_discover_granule_url_no_match(monkeypatch):
    day = dt.date(2024, 6, 1)
    other = ("https://data.laadsdaac.earthdatacloud.nasa.gov/prod-lads/"
             "VNP46A2/VNP46A2.A2024153.h07v11.002.2024153120000.h5")
    monkeypatch.setattr(
        bm, "_cmr_get_json",
        lambda url, timeout=60: _fake_cmr([_cmr_entry(other)]))
    assert discover_granule_url(day, 7, 10) is None


def test_discover_granule_url_empty(monkeypatch):
    monkeypatch.setattr(
        bm, "_cmr_get_json", lambda url, timeout=60: _fake_cmr([]))
    assert discover_granule_url(dt.date(2024, 6, 1), 7, 10) is None


# ---------------------------------------------------------------------------
# date validation
# ---------------------------------------------------------------------------


def test_bm_sample_days():
    days = bm_sample_days(dt.date(2024, 1, 1), dt.date(2024, 1, 5))
    assert len(days) == 5
    assert bm_sample_days(dt.date(2024, 1, 1), dt.date(2024, 1, 5),
                          stride_days=2) == [dt.date(2024, 1, 1),
                                             dt.date(2024, 1, 3),
                                             dt.date(2024, 1, 5)]
    with pytest.raises(ValueError, match="stride_days"):
        bm_sample_days(dt.date(2024, 1, 1), dt.date(2024, 1, 2),
                       stride_days=0)


def _fake_field(monkeypatch, payload, url="https://example.invalid/x.h5"):
    monkeypatch.setattr(bm, "discover_granule_url",
                        lambda day, h, v, timeout=60: url)
    monkeypatch.setattr(bm, "_bm_download",
                        lambda u, opener, timeout=300: payload)
    monkeypatch.setattr(bm, "_blackmarble_opener", lambda: object())


def _tile_payload(h5py, h=7, v=10, ny=20, nx=20, value=10.0,
                  sds=bm.BM_SDS, fill=None, scale=1.0, offset=0.0):
    buf = io.BytesIO()
    with h5py.File(buf, "w") as fh:
        data = np.full((ny, nx), (value - offset) / scale)
        if fill is not None:
            data[0, 0] = fill
        ds = fh.create_dataset(sds, data=data)
        ds.attrs["_FillValue"] = fill if fill is not None else -1.0
        ds.attrs["scale_factor"] = scale
        ds.attrs["add_offset"] = offset
        ds.attrs["units"] = bm.BM_UNITS
    return buf.getvalue()


def test_fetch_rejects_pre_record_dates(monkeypatch):
    _h5py()
    monkeypatch.setattr(bm, "_blackmarble_opener", lambda: object())
    with pytest.raises(ValueError, match="2012-01-19"):
        fetch_blackmarble((-110.0, -20.0, -100.0, -10.0),
                          "2000-01-01", "2000-01-05")


def test_fetch_rejects_bad_product():
    with pytest.raises(ValueError, match="only 'daily'"):
        fetch_blackmarble((-110.0, -20.0, -100.0, -10.0),
                          "2024-01-01", "2024-01-02", product="monthly")


def test_fetch_rejects_bad_resolution():
    with pytest.raises(ValueError, match="resolution"):
        fetch_blackmarble((-110.0, -20.0, -100.0, -10.0),
                          "2024-01-01", "2024-01-02", resolution=0.0)


# ---------------------------------------------------------------------------
# HDF5 parsing
# ---------------------------------------------------------------------------


def test_parse_tile_basic():
    h5py = _h5py()
    payload = _tile_payload(h5py, value=10.0)
    lats, lons, radiance, sds = _parse_blackmarble_bytes(payload, 7, 10)
    assert sds == bm.BM_SDS
    assert lats.shape == (20,) and lons.shape == (20,)
    assert radiance.shape == (20, 20)
    # lats descend (row 0 = north), lons ascend
    assert lats[0] > lats[-1] and lons[0] < lons[-1]
    assert np.nanmax(radiance) == pytest.approx(10.0)
    # tile spans exactly 10x10
    assert lons[-1] - lons[0] < 10.0 and lats[0] - lats[-1] < 10.0
    assert lons[0] >= -110.0 and lons[-1] <= -100.0
    assert lats[0] <= -10.0 and lats[-1] >= -20.0


def test_parse_tile_fill_and_scale():
    h5py = _h5py()
    payload = _tile_payload(h5py, value=10.0, fill=65535.0,
                            scale=0.1, offset=1.0)
    _, _, radiance, _ = _parse_blackmarble_bytes(payload, 7, 10)
    assert np.isnan(radiance[0, 0])  # fill masked
    assert radiance[0, 1] == pytest.approx(10.0)  # scale/offset applied


def test_parse_tile_fallback_sds():
    h5py = _h5py()
    payload = _tile_payload(h5py, sds=bm.BM_SDS_FALLBACK, value=5.0)
    _, _, radiance, sds = _parse_blackmarble_bytes(payload, 7, 10)
    assert sds == bm.BM_SDS_FALLBACK
    assert np.nanmax(radiance) == pytest.approx(5.0)


def test_parse_tile_no_sds_raises():
    h5py = _h5py()
    payload = _tile_payload(h5py, sds="Something_Else", value=5.0)
    with pytest.raises(ValueError, match="no NTL science dataset"):
        _parse_blackmarble_bytes(payload, 7, 10)


# ---------------------------------------------------------------------------
# block averaging
# ---------------------------------------------------------------------------


def test_block_average_single_tile():
    h5py = _h5py()
    payload = _tile_payload(h5py, h=7, v=10, ny=20, nx=20, value=10.0)
    t_lats, t_lons, radiance, _ = _parse_blackmarble_bytes(payload, 7, 10)
    # output grid: the exact tile bbox at 5° resolution -> 2x2 cells
    acc_sum, acc_count = _block_average_to_grid(
        t_lats, t_lons, radiance, -110.0, -10.0, 5.0, 2, 2)
    assert acc_count.shape == (2, 2)
    assert np.all(acc_count == 100)  # 10x10 native pixels per cell
    assert np.allclose(acc_sum / acc_count, 10.0)


def test_block_average_nan_aware():
    h5py = _h5py()
    payload = _tile_payload(h5py, h=7, v=10, ny=20, nx=20, value=10.0)
    t_lats, t_lons, radiance, _ = _parse_blackmarble_bytes(payload, 7, 10)
    radiance[:, :] = np.nan
    radiance[0, 0] = 4.0
    acc_sum, acc_count = _block_average_to_grid(
        t_lats, t_lons, radiance, -110.0, -10.0, 5.0, 2, 2)
    assert acc_count[0, 0] == 1 and acc_count[1, 1] == 0
    assert acc_sum[0, 0] == pytest.approx(4.0)


# ---------------------------------------------------------------------------
# LightsField model
# ---------------------------------------------------------------------------


def test_lights_field_synthetic():
    field = LightsField.synthetic()
    assert field.shape[0] == 5
    assert field.units == bm.BM_UNITS
    assert field.source == "synthetic"
    t0, t1 = field.time_range
    assert (t1 - t0).days == 4
    # city blobs are bright, background is NaN
    assert np.nanmax(field.values) > 10.0
    assert np.isnan(field.values).any()


def test_lights_field_synthetic_deterministic():
    a = LightsField.synthetic(seed=3)
    b = LightsField.synthetic(seed=3)
    assert np.allclose(np.nan_to_num(a.values), np.nan_to_num(b.values))


def test_lights_field_shape_mismatch():
    with pytest.raises(ValueError):
        LightsField(times=[dt.datetime(2024, 1, 1, tzinfo=dt.timezone.utc)],
                    lats=np.array([0.0]), lons=np.array([0.0, 1.0]),
                    values=np.zeros((1, 1, 1)))


def test_lights_field_select_time():
    field = LightsField.synthetic(start="2024-01-01", end="2024-01-05")
    sub = field.select_time("2024-01-02", "2024-01-03")
    assert len(sub) == 2
    assert sub.times[0].date() == dt.date(2024, 1, 2)


def test_lights_field_select_bbox():
    field = LightsField.synthetic(bbox=(-125.0, 25.0, -66.0, 49.0),
                                  resolution=1.0)
    sub = field.select_bbox((-120.0, 30.0, -110.0, 40.0))
    assert sub.bounds == (-120.0, 30.0, -110.0, 40.0)
    assert sub.shape[1] == 10 and sub.shape[2] == 10


def test_lights_field_json_roundtrip(tmp_path):
    field = LightsField.synthetic()
    path = str(tmp_path / "lights.json")
    field.to_json(path)
    back = LightsField.from_json(path)
    assert back.shape == field.shape
    assert back.units == field.units
    assert np.allclose(np.nan_to_num(back.values),
                       np.nan_to_num(field.values))
    assert back.times == field.times


def test_lights_field_total_radiance():
    field = LightsField.synthetic()
    totals = field.total_radiance
    assert totals.shape == (len(field),)
    assert np.all(totals > 0)


# ---------------------------------------------------------------------------
# fetch_blackmarble end-to-end (mocked network)
# ---------------------------------------------------------------------------


def test_fetch_blackmarble_mosaic(monkeypatch):
    h5py = _h5py()
    payload = _tile_payload(h5py, h=7, v=10, ny=20, nx=20, value=10.0)
    _fake_field(monkeypatch, payload)
    field = fetch_blackmarble((-110.0, -20.0, -100.0, -10.0),
                              "2024-06-01", "2024-06-02",
                              resolution=5.0)
    assert field.shape == (2, 2, 2)
    assert np.allclose(field.values, 10.0)
    assert field.units == bm.BM_UNITS
    prov = field.provenance
    assert prov["product"] == "VNP46A2"
    assert prov["n_files"] == 2
    assert len(prov["files"]) == 2
    rec = prov["files"][0]
    assert rec["url"] == "https://example.invalid/x.h5"
    assert rec["tile"] == "h07v10"
    assert rec["sha256"] == hashlib.sha256(payload).hexdigest()
    assert prov["day_tiles"]["2024-06-01"] == {"expected": 1, "retrieved": 1}
    assert prov["retrieved_at"]


def test_fetch_blackmarble_skips_missing_tiles(monkeypatch):
    h5py = _h5py()
    payload = _tile_payload(h5py, h=7, v=10, ny=20, nx=20, value=10.0)
    # two tiles in bbox, one has no CMR granule
    monkeypatch.setattr(
        bm, "discover_granule_url",
        lambda day, h, v, timeout=60: ("https://example.invalid/x.h5"
                                      if (h, v) == (7, 10) else None))
    monkeypatch.setattr(bm, "_bm_download",
                        lambda u, opener, timeout=300: payload)
    monkeypatch.setattr(bm, "_blackmarble_opener", lambda: object())
    field = fetch_blackmarble((-110.0, -20.0, -90.0, -10.0),
                              "2024-06-01", "2024-06-01",
                              resolution=5.0)
    assert field.shape[0] == 1
    assert any("no CMR granule" in s for s in field.provenance["skipped"])
    assert field.provenance["day_tiles"]["2024-06-01"]["retrieved"] == 1


def test_fetch_blackmarble_no_granules_raises(monkeypatch):
    monkeypatch.setattr(bm, "discover_granule_url",
                        lambda day, h, v, timeout=60: None)
    monkeypatch.setattr(bm, "_blackmarble_opener", lambda: object())
    with pytest.raises(ValueError, match="no Black Marble granules"):
        fetch_blackmarble((-110.0, -20.0, -100.0, -10.0),
                          "2024-06-01", "2024-06-01")


def test_fetch_blackmarble_credentials_missing(monkeypatch):
    _h5py()

    def _no_creds():
        raise CredentialsMissing("nope")

    monkeypatch.setattr(bm, "_blackmarble_opener", _no_creds)
    with pytest.raises(CredentialsMissing):
        fetch_blackmarble((-110.0, -20.0, -100.0, -10.0),
                          "2024-06-01", "2024-06-02")


def test_fetch_blackmarble_stride(monkeypatch):
    h5py = _h5py()
    payload = _tile_payload(h5py, h=7, v=10, ny=20, nx=20, value=10.0)
    _fake_field(monkeypatch, payload)
    field = fetch_blackmarble((-110.0, -20.0, -100.0, -10.0),
                              "2024-06-01", "2024-06-05",
                              resolution=5.0, stride_days=2)
    assert [t.date().isoformat() for t in field.times] == [
        "2024-06-01", "2024-06-03", "2024-06-05"]


def test_cmr_http_error_becomes_runtimeerror(monkeypatch):
    import urllib.error

    def _boom(url, timeout=60):
        raise urllib.error.HTTPError(url, 500, "boom", {}, None)

    monkeypatch.setattr(urllib.request, "urlopen", _boom)
    with pytest.raises(RuntimeError, match="CMR granule search failed"):
        discover_granule_url(dt.date(2024, 6, 1), 7, 10)


def test_main_demo(capsys):
    bm.main_demo()
    out = capsys.readouterr().out
    assert "synthetic LightsField" in out

"""Tests for currents.currents_global (OSCAR v2.0 + CMEMS currents).

Fully offline: OSCAR NetCDF payloads are fabricated in-memory with
netCDF4; the CMR search, the authenticated download, and the CMEMS
toolbox are all monkeypatched. No live network in this suite.
"""

import datetime as _dt
import json
import urllib.error

import numpy as np
import pytest

import currents.currents_global as cg
from currents.models import CurrentField


# ---------------------------------------------------------------------------
# fabricated OSCAR NetCDF payloads (on-wire dim order: time, longitude, latitude)
# ---------------------------------------------------------------------------

def _make_oscar_payload(path, nx, ny, lon0_idx=0, lat0_idx=0):
    """Write a minimal OSCAR-shaped NetCDF file; return its bytes."""
    import netCDF4 as nc
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("longitude", nx)
        ds.createDimension("latitude", ny)
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units = "days since 1990-01-01 00:00:00"
        tv[:] = [0.0]
        ds.createVariable("longitude", "f4", ("longitude",))[:] = \
            (lon0_idx + np.arange(nx, dtype=float)) * 0.25
        ds.createVariable("latitude", "f4", ("latitude",))[:] = \
            -89.75 + (lat0_idx + np.arange(ny, dtype=float)) * 0.25
        u = ds.createVariable("u", "f4", ("time", "longitude", "latitude"),
                              fill_value=-999.0)
        v = ds.createVariable("v", "f4", ("time", "longitude", "latitude"),
                              fill_value=-999.0)
        uu = np.arange(nx * ny, dtype=np.float32).reshape(1, nx, ny) * 0.01
        vv = (np.arange(nx * ny, dtype=np.float32).reshape(1, nx, ny)
              * 0.02 + 0.5)
        uu[0, 0, 0] = -999.0  # one masked cell
        u[:] = uu
        v[:] = vv
    with open(path, "rb") as fh:
        return fh.read()


def _make_oscar_payload_unsorted_lon(path):
    """Payload whose longitude axis crosses the 0/360 seam unsorted."""
    import netCDF4 as nc
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("longitude", 3)
        ds.createDimension("latitude", 2)
        ds.createVariable("time", "f8", ("time",))[:] = [0.0]
        ds.createVariable("longitude", "f4", ("longitude",))[:] = \
            [359.5, 0.0, 0.5]
        ds.createVariable("latitude", "f4", ("latitude",))[:] = [10.0, 10.25]
        u = ds.createVariable("u", "f4", ("time", "longitude", "latitude"),
                              fill_value=-999.0)
        v = ds.createVariable("v", "f4", ("time", "longitude", "latitude"),
                              fill_value=-999.0)
        # u column j carries the value j (tracks the permutation)
        u[:] = np.array([[[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]],
                        dtype=np.float32)
        v[:] = np.zeros((1, 3, 2), dtype=np.float32)
    with open(path, "rb") as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# collection picking
# ---------------------------------------------------------------------------

_TODAY = _dt.date(2026, 9, 26)


def test_oscar_collection_for_final():
    assert cg.oscar_collection_for(_dt.date(2020, 1, 1), _TODAY) == "final"
    assert cg.oscar_collection_for(_dt.date(1993, 1, 1), _TODAY) == "final"


def test_oscar_collection_for_interim():
    assert cg.oscar_collection_for(_dt.date(2026, 6, 1), _TODAY) == "interim"
    assert cg.oscar_collection_for(_dt.date(2025, 5, 1), _TODAY) == "interim"


def test_oscar_collection_for_nrt():
    assert cg.oscar_collection_for(_dt.date(2026, 9, 20), _TODAY) == "nrt"
    assert cg.oscar_collection_for(_TODAY, _TODAY) == "nrt"


def test_oscar_collection_for_boundary_final():
    # age == 540 days is NOT older than 540 -> interim
    assert cg.oscar_collection_for(
        _TODAY - _dt.timedelta(days=540), _TODAY) == "interim"
    assert cg.oscar_collection_for(
        _TODAY - _dt.timedelta(days=541), _TODAY) == "final"


def test_oscar_collection_for_boundary_interim():
    # age == 45 days is NOT older than 45 -> nrt
    assert cg.oscar_collection_for(
        _TODAY - _dt.timedelta(days=45), _TODAY) == "nrt"
    assert cg.oscar_collection_for(
        _TODAY - _dt.timedelta(days=46), _TODAY) == "interim"


def test_oscar_collection_for_before_record_raises():
    with pytest.raises(ValueError, match="1993-01-01"):
        cg.oscar_collection_for(_dt.date(1992, 12, 31), _TODAY)


def test_oscar_collection_for_future_raises():
    with pytest.raises(ValueError, match="after today"):
        cg.oscar_collection_for(_TODAY + _dt.timedelta(days=1), _TODAY)


# ---------------------------------------------------------------------------
# granule titles / service URLs
# ---------------------------------------------------------------------------

def test_oscar_granule_title_all_tiers():
    day = _dt.date(2024, 1, 15)
    assert cg.oscar_granule_title("final", day) == \
        "oscar_currents_final_20240115.nc"
    assert cg.oscar_granule_title("interim", day) == \
        "oscar_currents_interim_20240115.nc"
    assert cg.oscar_granule_title("nrt", day) == \
        "oscar_currents_nrt_20240115.nc"


def test_oscar_granule_title_unknown_collection_raises():
    with pytest.raises(ValueError, match="unknown OSCAR collection"):
        cg.oscar_granule_title("bogus", _dt.date(2024, 1, 15))


def test_oscar_service_url():
    url = cg.oscar_service_url("final", _dt.date(2024, 1, 15))
    assert url == ("https://opendap.earthdata.nasa.gov/collections/"
                   "C2098858642-POCLOUD/granules/"
                   "oscar_currents_final_20240115.nc")


def test_oscar_subset_urls_shape():
    day = _dt.date(2024, 1, 15)
    urls = cg.oscar_subset_urls(day, "final", (-81.0, 25.0, -55.0, 43.0))
    assert len(urls) == 1
    url = urls[0]
    assert ".nc.nc" not in url  # title already ends in .nc; no double suffix
    assert ".nc?u[0:1:0]" in url
    assert url.startswith(
        "https://opendap.earthdata.nasa.gov/collections/C2098858642-POCLOUD/")
    assert "u[0:1:0]" in url and "v[0:1:0]" in url
    # on-wire dim order is (time, longitude, latitude): lon index first
    assert "u[0:1:0][1116:1:1220][459:1:531]" in url


def test_oscar_subset_urls_base_override():
    day = _dt.date(2024, 1, 15)
    urls = cg.oscar_subset_urls(day, "final", (-81.0, 25.0, -55.0, 43.0),
                                base_url="https://cmr.example/opendap/x")
    assert urls[0].startswith("https://cmr.example/opendap/x?")
    assert "u[0:1:0]" in urls[0]


# ---------------------------------------------------------------------------
# longitude windows
# ---------------------------------------------------------------------------

def test_oscar_lon_windows_simple():
    assert cg.oscar_lon_windows((-81.0, 25.0, -55.0, 43.0)) == \
        [(279.0, 305.0)]


def test_oscar_lon_windows_antimeridian_single():
    # (170, ., -170, .) wraps to one 0-360 window
    assert cg.oscar_lon_windows((170.0, -10.0, -170.0, 10.0)) == \
        [(170.0, 190.0)]


def test_oscar_lon_windows_wide_splits():
    # (-170, ., 170, .) crosses 360 in 0-360 space -> two windows
    assert cg.oscar_lon_windows((-170.0, -10.0, 170.0, 10.0)) == \
        [(190.0, 359.75), (0.0, 170.0)]


def test_oscar_lon_windows_full_globe():
    assert cg.oscar_lon_windows((-180.0, -60.0, 180.0, 60.0)) == \
        [(0.0, 359.75)]


def test_oscar_index_windows_math():
    wins = cg.oscar_index_windows((-81.0, 25.0, -55.0, 43.0))
    assert wins == [(1116, 1220, 459, 531)]
    for i0, i1, j0, j1 in wins:
        assert 0 <= i0 <= i1 < cg.OSCAR_NLON
        assert 0 <= j0 <= j1 < cg.OSCAR_NLAT


def test_oscar_index_windows_wrap_count():
    wins = cg.oscar_index_windows((-170.0, -10.0, 170.0, 10.0))
    assert len(wins) == 2


# ---------------------------------------------------------------------------
# NetCDF parsing
# ---------------------------------------------------------------------------

def test_parse_oscar_bytes_transposes(tmp_path):
    payload = _make_oscar_payload(str(tmp_path / "t.nc"), nx=4, ny=3,
                                  lon0_idx=1116, lat0_idx=459)
    lats, lons, u, v = cg._parse_oscar_bytes(payload)
    # on-wire (1, nx_lon=4, ny_lat=3) -> canonical (1, ny=3, nx=4)
    assert u.shape == (1, 3, 4)
    assert v.shape == (1, 3, 4)
    # masked fill cell: on-wire u[0, 0, 0] -> canonical u[0, 0, 0]
    assert bool(u.mask[0, 0, 0])
    # value check: on-wire uu[0, ix, iy] = (ix*ny+iy)*0.01
    # canonical u[0, iy, ix]
    assert float(u[0, 2, 3]) == pytest.approx((3 * 3 + 2) * 0.01)
    assert float(v[0, 1, 1]) == pytest.approx((1 * 3 + 1) * 0.02 + 0.5)
    assert lats.shape == (3,) and lons.shape == (4,)


def test_parse_oscar_bytes_lon_normalization(tmp_path):
    payload = _make_oscar_payload_unsorted_lon(str(tmp_path / "t.nc"))
    lats, lons, u, v = cg._parse_oscar_bytes(payload)
    # 359.5 -> -0.5, 0.0 -> 0.0, 0.5 -> 0.5, sorted increasing
    assert list(lons) == pytest.approx([-0.5, 0.0, 0.5])
    # u column j carried value j; after the sort the -0.5 column (was
    # index 0) comes first
    assert [float(u[0, 0, j]) for j in range(3)] == \
        pytest.approx([0.0, 1.0, 2.0])


def test_parse_oscar_bytes_wrong_dim_order_raises(tmp_path):
    import netCDF4 as nc
    path = str(tmp_path / "bad.nc")
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("lat", 2)
        ds.createDimension("lon", 3)
        ds.createVariable("time", "f8", ("time",))[:] = [0.0]
        ds.createVariable("lat", "f4", ("lat",))[:] = [10.0, 10.25]
        ds.createVariable("lon", "f4", ("lon",))[:] = [0.0, 0.25, 0.5]
        ds.createVariable("u", "f4", ("time", "lat", "lon"))[:] = \
            np.zeros((1, 2, 3), dtype=np.float32)
        ds.createVariable("v", "f4", ("time", "lat", "lon"))[:] = \
            np.zeros((1, 2, 3), dtype=np.float32)
    with open(path, "rb") as fh:
        payload = fh.read()
    with pytest.raises(ValueError, match="dimensions"):
        cg._parse_oscar_bytes(payload)


def test_parse_oscar_bytes_missing_var_raises(tmp_path):
    import netCDF4 as nc
    path = str(tmp_path / "nov.nc")
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("longitude", 2)
        ds.createDimension("latitude", 2)
        ds.createVariable("time", "f8", ("time",))[:] = [0.0]
        ds.createVariable("longitude", "f4", ("longitude",))[:] = [0.0, 0.25]
        ds.createVariable("latitude", "f4", ("latitude",))[:] = [10.0, 10.25]
        ds.createVariable("v", "f4",
                          ("time", "longitude", "latitude"))[:] = \
            np.zeros((1, 2, 2), dtype=np.float32)
    with open(path, "rb") as fh:
        payload = fh.read()
    with pytest.raises(ValueError, match="missing variable 'u'"):
        cg._parse_oscar_bytes(payload)


def test_concat_oscar_parts_sorts_lon():
    lats = np.array([10.0, 10.25])
    u1 = np.ma.array(np.ones((1, 2, 2)))
    v1 = np.ma.array(np.zeros((1, 2, 2)))
    u2 = np.ma.array(np.ones((1, 2, 2)) * 2.0)
    v2 = np.ma.array(np.zeros((1, 2, 2)))
    lats_o, lons_o, u_o, v_o = cg._concat_oscar_parts([
        (lats, np.array([5.0, 10.0]), u1, v1),
        (lats, np.array([-170.0, -160.0]), u2, v2),
    ])
    assert list(lons_o) == [-170.0, -160.0, 5.0, 10.0]
    assert u_o.shape == (1, 2, 4)
    assert float(u_o[0, 0, 0]) == 2.0  # the -170 column came from part 2
    assert float(u_o[0, 0, 3]) == 1.0


# ---------------------------------------------------------------------------
# CMR discovery (mocked download)
# ---------------------------------------------------------------------------

def _fake_cmr_payload():
    return json.dumps({
        "feed": {"entry": [
            {"title": "oscar_currents_final_20240115.nc",
             "time_start": "2024-01-15T00:00:00Z",
             "time_end": "2024-01-15T23:59:59Z",
             "links": [
                 {"href": "https://opendap.earthdata.nasa.gov/x",
                  "rel": "OPENDAP DATA"},
                 {"href": "https://example/data.nc"}]},
        ]}
    }).encode()


def test_cmr_search_oscar_granules_url_and_parse(monkeypatch):
    seen = {}

    def fake_download(url, timeout=300):
        seen["url"] = url
        return _fake_cmr_payload()

    monkeypatch.setattr(cg, "_download_bytes", fake_download)
    granules = cg.cmr_search_oscar_granules(
        "final", "2024-01-15", "2024-01-16")
    assert "collection_concept_id=C2098858642-POCLOUD" in seen["url"].replace(
        "%3A", ":").replace("%2F", "/") or "C2098858642" in seen["url"]
    assert len(granules) == 1
    assert granules[0]["title"] == "oscar_currents_final_20240115.nc"
    assert granules[0]["opendap_url"] == \
        "https://opendap.earthdata.nasa.gov/x"


def test_cmr_search_oscar_granules_unknown_collection_raises():
    with pytest.raises(ValueError, match="unknown OSCAR collection"):
        cg.cmr_search_oscar_granules("bogus", "2024-01-15", "2024-01-16")


def test_oscar_match_granule_by_title():
    granules = [{"title": "oscar_currents_final_20240115.nc",
                 "time_start": "2024-01-15T00:00:00Z",
                 "opendap_url": None, "links": []}]
    hit = cg.oscar_match_granule(_dt.date(2024, 1, 15), granules)
    assert hit["title"] == "oscar_currents_final_20240115.nc"


def test_oscar_match_granule_missing_raises():
    with pytest.raises(RuntimeError, match="no OSCAR granule"):
        cg.oscar_match_granule(_dt.date(2024, 1, 16), [])


# ---------------------------------------------------------------------------
# fetch_oscar validation + credentials
# ---------------------------------------------------------------------------

def test_fetch_oscar_start_before_record_raises():
    with pytest.raises(ValueError, match="1993-01-01"):
        cg.fetch_oscar((-80, 20, -60, 40), "1990-01-01", "2024-01-15")


def test_fetch_oscar_stride_zero_raises():
    with pytest.raises(ValueError, match="stride_days"):
        cg.fetch_oscar((-80, 20, -60, 40), "2024-01-15", "2024-01-20",
                       stride_days=0)


def test_fetch_oscar_no_credentials_raises(monkeypatch):
    monkeypatch.setattr(cg, "earthdata_credentials", lambda: None)
    with pytest.raises(cg.CredentialsMissing, match="Earthdata Login"):
        cg.fetch_oscar((-80, 20, -60, 40), "2024-01-15", "2024-01-20")


def test_fetch_oscar_full_assembly_mocked(monkeypatch, tmp_path):
    monkeypatch.setattr(cg, "earthdata_credentials",
                        lambda: ("user", "pass"))
    granules = [{"title": "oscar_currents_final_20240115.nc",
                 "time_start": "2024-01-15T00:00:00Z",
                 "opendap_url": None, "links": []}]
    monkeypatch.setattr(cg, "cmr_search_oscar_granules",
                        lambda coll, s, e, page_size=200: granules)
    bbox = (-81.0, 25.0, -55.0, 43.0)
    i0, i1, j0, j1 = cg.oscar_index_windows(bbox)[0]
    payload = _make_oscar_payload(str(tmp_path / "g.nc"), nx=i1 - i0 + 1,
                                  ny=j1 - j0 + 1, lon0_idx=i0, lat0_idx=j0)
    monkeypatch.setattr(cg, "_oscar_get", lambda opener, url, ctx: payload)

    field = cg.fetch_oscar(bbox, "2024-01-15", "2024-01-15")
    assert isinstance(field, CurrentField)
    assert field.times == ["2024-01-15T00:00:00+00:00"]
    assert field.temperature is None
    ny, nx = j1 - j0 + 1, i1 - i0 + 1
    assert field.u.shape == (1, ny, nx)
    assert field.source == "podaac/OSCAR_L4_OC_V2.0"
    prov = field.provenance
    assert prov["collections_used"] == ["final"]
    assert prov["dates"][0]["collection"] == "final"
    assert prov["dates"][0]["granule_title"] == \
        "oscar_currents_final_20240115.nc"
    assert len(prov["service_urls"]) == 1
    assert "C2098858642-POCLOUD" in prov["service_urls"][0]
    assert len(prov["combined_sha256"]) == 64
    assert "retrieved_at" in prov
    assert prov["units"].startswith("m/s")
    # corrections from the verified-source audit are on the record
    assert any("ERDDAP" in c for c in prov["corrections"])


# ---------------------------------------------------------------------------
# _oscar_get error mapping
# ---------------------------------------------------------------------------

class _FakeOpener:
    def __init__(self, exc):
        self._exc = exc

    def open(self, req, timeout=None):
        raise self._exc


def test_oscar_get_401_raises_credentials_missing():
    exc = urllib.error.HTTPError("https://x", 401, "Unauthorized", {}, None)
    with pytest.raises(cg.CredentialsMissing, match="rejected the credentials"):
        cg._oscar_get(_FakeOpener(exc), "https://x", "ctx")


def test_oscar_get_other_http_raises_runtime():
    exc = urllib.error.HTTPError("https://x", 500, "Server Error", {}, None)
    with pytest.raises(RuntimeError, match="HTTP 500"):
        cg._oscar_get(_FakeOpener(exc), "https://x", "ctx")


# ---------------------------------------------------------------------------
# fetch_cmems_currents
# ---------------------------------------------------------------------------

def _fake_cmems_field(nt=5):
    return CurrentField.synthetic(nt=nt, ny=4, nx=5, seed=3,
                                  source="cmems:test")


def test_fetch_cmems_currents_stride_and_provenance(monkeypatch, tmp_path):
    import currents.cmems as cmems
    calls = {}

    def fake_subset(preset, bbox, start, end, work_dir, variables=None,
                    **kw):
        calls.update(preset=preset, bbox=bbox, start=start, end=end,
                     work_dir=work_dir, variables=variables)
        return str(tmp_path / "subset.nc")

    monkeypatch.setattr(cmems, "subset_cmems", fake_subset)
    monkeypatch.setattr(cmems, "parse_cmems_netcdf",
                        lambda path: _fake_cmems_field())

    field = cg.fetch_cmems_currents((-80, 20, -60, 40), "2024-01-01",
                                    "2024-01-05", stride_days=2)
    assert calls["preset"] == cg.CMEMS_CURRENTS_PRESET
    assert calls["variables"] == ["uo", "vo", "thetao"]
    assert len(field.times) == 3  # 5 daily steps, every 2nd kept
    assert field.u.shape[0] == 3
    prov = field.provenance
    assert prov["n_timesteps_total"] == 5
    assert prov["n_timesteps_kept"] == 3
    assert prov["stride_days"] == 2
    assert "potential temperature" in prov["temperature_note"]
    assert prov["netcdf_path"] == str(tmp_path / "subset.nc")


def test_fetch_cmems_currents_temperature_carried(monkeypatch, tmp_path):
    import currents.cmems as cmems
    monkeypatch.setattr(cmems, "subset_cmems",
                        lambda *a, **k: str(tmp_path / "s.nc"))
    fake = _fake_cmems_field()
    monkeypatch.setattr(cmems, "parse_cmems_netcdf", lambda path: fake)
    field = cg.fetch_cmems_currents((-80, 20, -60, 40), "2024-01-01",
                                    "2024-01-02")
    assert field.temperature is not None
    assert field.temperature.shape == field.u.shape


def test_fetch_cmems_currents_toolbox_missing(monkeypatch):
    import currents.cmems as cmems

    def boom():
        raise ImportError("CMEMS access needs the copernicusmarine toolbox")

    monkeypatch.setattr(cmems, "require_toolbox", boom)
    with pytest.raises(ImportError, match="copernicusmarine"):
        cg.fetch_cmems_currents((-80, 20, -60, 40), "2024-01-01",
                                "2024-01-02")


def test_fetch_cmems_currents_bad_stride_raises():
    with pytest.raises(ValueError, match="stride_days"):
        cg.fetch_cmems_currents((-80, 20, -60, 40), "2024-01-01",
                                "2024-01-02", stride_days=0)


def test_fetch_cmems_currents_start_after_end_raises():
    with pytest.raises(ValueError, match="after end"):
        cg.fetch_cmems_currents((-80, 20, -60, 40), "2024-01-05",
                                "2024-01-01")


# ---------------------------------------------------------------------------
# synthetic fixtures + demo
# ---------------------------------------------------------------------------

def test_oscar_synthetic_no_temperature():
    f = cg.oscar_synthetic(nt=2)
    assert f.temperature is None
    assert f.u.shape == (2, 6, 8)
    assert "OSCAR" in f.source
    assert f.provenance["synthetic"] is True


def test_cmems_currents_synthetic_has_temperature():
    f = cg.cmems_currents_synthetic(nt=2)
    assert f.temperature is not None
    assert f.temperature.shape == f.u.shape
    assert "cmems" in f.source


def test_main_demo_runs(capsys):
    cg.main_demo()
    out = capsys.readouterr().out
    assert "[oscar]" in out and "[cmems-currents]" in out

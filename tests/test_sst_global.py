"""Tests for the global-SST module (NOAA OISST v2.1 + NASA JPL MUR v4.1).

Fully offline: NetCDF payloads are fabricated in-memory with netCDF4,
HTTP is mocked, and credentials are faked via monkeypatched env vars.
"""

import datetime as dt
import hashlib
import io
import json
import os
import sys

import numpy as np
import pytest

from currents import sst_global
from currents.sst_global import (
    CredentialsMissing,
    SstField,
    _chunk_date_range,
    _parse_mur_bytes,
    _parse_mur_grid,
    _parse_oisst_bytes,
    earthdata_credentials,
    fetch_mur,
    fetch_oisst,
    mur_granule_title,
    mur_opendap_url,
    mur_subset_urls,
    oisst_lon_windows,
    oisst_sst_urls,
    validate_sst_bbox,
)


# ---------------------------------------------------------------------------
# fabricated NetCDF payloads
# ---------------------------------------------------------------------------

def _make_oisst_payload(path, nt=3, ny=5, nx=6, lon0=280.0):
    """Write a minimal OISST-shaped NetCDF file; return its bytes."""
    import netCDF4 as nc
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("time", nt)
        ds.createDimension("zlev", 1)
        ds.createDimension("latitude", ny)
        ds.createDimension("longitude", nx)
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units = "days since 1981-09-01 12:00:00"
        tv[:] = np.arange(nt, dtype=float) * 30.0
        ds.createVariable("zlev", "f4", ("zlev",))[:] = [0.0]
        ds.createVariable("latitude", "f4", ("latitude",))[:] = \
            np.linspace(20.0, 40.0, ny)
        ds.createVariable("longitude", "f4", ("longitude",))[:] = \
            lon0 + np.arange(nx, dtype=float) * 0.25
        sv = ds.createVariable("sst", "f4", ("time", "zlev", "latitude", "longitude"),
                               fill_value=sst_global.OISST_FILL)
        data = 25.0 + np.arange(nt * ny * nx, dtype=np.float32).reshape(nt, 1, ny, nx) * 0.01
        data[0, 0, 0, 0] = sst_global.OISST_FILL  # one masked cell
        sv[:] = data
    with open(path, "rb") as fh:
        return fh.read()


def _make_mur_payload(path, ny=5, nx=7, kelvin=True):
    """Write a minimal MUR-shaped NetCDF file; return its bytes."""
    import netCDF4 as nc
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("lat", ny)
        ds.createDimension("lon", nx)
        ds.createVariable("time", "f8", ("time",))[:] = [0.0]
        ds.createVariable("lat", "f4", ("lat",))[:] = np.linspace(20.0, 24.0, ny)
        ds.createVariable("lon", "f4", ("lon",))[:] = np.linspace(-80.0, -74.0, nx)
        sv = ds.createVariable("analysed_sst", "f4", ("time", "lat", "lon"),
                               fill_value=-32768.0)
        base = 290.0 if kelvin else 17.0
        sv[:] = base + np.arange(ny * nx, dtype=np.float32).reshape(1, ny, nx) * 0.01
        if kelvin:
            sv.units = "kelvin"
    with open(path, "rb") as fh:
        return fh.read()


MUR_DAS = """Attributes {
 lat {
  String long_name "latitude";
  Float64 actual_range 20.0, 24.0;
 }
 lon {
  Float64 actual_range -80.0, -74.0;
 }
}"""

MUR_DDS = """Dataset {
 Float64 lat[lat = 5];
 Float64 lon[lon = 7];
 Float32 analysed_sst[time = 1][lat = 5][lon = 7];
} fake;"""


# ---------------------------------------------------------------------------
# bbox validation
# ---------------------------------------------------------------------------

def test_validate_sst_bbox_ok():
    assert validate_sst_bbox((-80, 20, -60, 40)) == (-80, 20, -60, 40)


def test_validate_sst_bbox_antimeridian_ok():
    # dateline-crossing boxes are valid; they wrap internally
    assert validate_sst_bbox((170, -10, -170, 10)) == (170, -10, -170, 10)


def test_validate_sst_bbox_bad_latitude():
    with pytest.raises(ValueError, match="latitudes"):
        validate_sst_bbox((-80, -95, -60, 40))


def test_validate_sst_bbox_bad_longitude():
    with pytest.raises(ValueError, match="longitudes"):
        validate_sst_bbox((-190, 20, -60, 40))


def test_validate_sst_bbox_descending_lat():
    with pytest.raises(ValueError, match="ascend"):
        validate_sst_bbox((-80, 40, -60, 20))


# ---------------------------------------------------------------------------
# OISST longitude windows + URL construction
# ---------------------------------------------------------------------------

def test_oisst_lon_windows_normal():
    assert oisst_lon_windows((-80, 20, -60, 40)) == [(280.0, 300.0)]


def test_oisst_lon_windows_antimeridian_wraps():
    assert oisst_lon_windows((170, -10, -170, 10)) == [(170.0, 190.0)]


def test_oisst_lon_windows_cross_360_splits():
    # (-170..170) becomes (190..530) in 0-360 space: two requests, no silent loss
    assert oisst_lon_windows((-170, -10, 170, 10)) == [
        (190.0, 359.875), (0.125, 170.0)]


def test_oisst_lon_windows_full_globe():
    assert oisst_lon_windows((-180, -90, 180, 90)) == [(0.125, 359.875)]


def test_oisst_sst_urls_include_zlev():
    urls = oisst_sst_urls((-80, 20, -60, 40), "2020-01-01", "2020-01-31")
    assert len(urls) == 1
    assert "[(0.0)]" in urls[0]  # singleton depth axis is required by ERDDAP
    assert "ncdcOisst21Agg" in urls[0]
    assert "[(280.0):(300.0)]" in urls[0]


def test_oisst_sst_urls_stride_and_time():
    urls = oisst_sst_urls((-80, 20, -60, 40), "2020-01-01", "2020-01-31",
                          stride_days=7)
    assert ":7:(" in urls[0]
    assert "(2020-01-01T12:00:00Z)" in urls[0]
    assert "(2020-01-31T12:00:00Z)" in urls[0]


def test_oisst_sst_urls_two_for_cross_360():
    urls = oisst_sst_urls((-170, -10, 170, 10), "2020-01-01", "2020-01-31")
    assert len(urls) == 2


def test_oisst_sst_urls_reject_bad_stride():
    with pytest.raises(ValueError, match="stride_days"):
        oisst_sst_urls((-80, 20, -60, 40), "2020-01-01", "2020-01-31",
                       stride_days=0)


def test_oisst_sst_urls_reject_reversed_dates():
    with pytest.raises(ValueError, match="after end"):
        oisst_sst_urls((-80, 20, -60, 40), "2020-02-01", "2020-01-01")


def test_oisst_sst_urls_reject_before_record():
    with pytest.raises(ValueError, match="1981"):
        oisst_sst_urls((-80, 20, -60, 40), "1970-01-01", "1970-02-01")


def test_chunk_date_range_short_is_single():
    d0 = dt.datetime(2020, 1, 1, 12)
    d1 = dt.datetime(2021, 6, 1, 12)
    assert _chunk_date_range(d0, d1, 5) == [(d0, d1)]


def test_chunk_date_range_splits_long():
    d0 = dt.datetime(2000, 1, 1, 12)
    d1 = dt.datetime(2019, 12, 31, 12)
    chunks = _chunk_date_range(d0, d1, 5)
    assert len(chunks) == 4
    assert chunks[0][0] == d0 and chunks[-1][1] == d1
    # contiguous, no gaps
    for (_a, b), (c, _d) in zip(chunks, chunks[1:]):
        assert c == b + dt.timedelta(seconds=1)


# ---------------------------------------------------------------------------
# OISST parsing
# ---------------------------------------------------------------------------

def test_parse_oisst_bytes_shape_and_lons(tmp_path):
    payload = _make_oisst_payload(str(tmp_path / "oisst.nc"))
    lats, lons, times, sst = _parse_oisst_bytes(payload)
    assert sst.shape == (3, 5, 6)
    assert lons[0] < -70 and lons[-1] < -60  # 280..281.25 -> -80..-78.75
    assert np.all(np.diff(lons) > 0)         # increasing after normalization
    assert len(times) == 3
    assert times[0].startswith("1981-09-01")


def test_parse_oisst_bytes_masks_fill(tmp_path):
    payload = _make_oisst_payload(str(tmp_path / "oisst.nc"))
    _, _, _, sst = _parse_oisst_bytes(payload)
    assert bool(sst.mask[0, 0, 0])           # the injected fill cell
    assert not bool(sst.mask[0, 0, 1])


def test_parse_oisst_bytes_wrong_dims(tmp_path):
    import netCDF4 as nc
    p = str(tmp_path / "bad.nc")
    with nc.Dataset(p, "w") as ds:
        ds.createDimension("time", 2)
        ds.createVariable("sst", "f4", ("time",))[:] = [1.0, 2.0]
    with open(p, "rb") as fh:
        payload = fh.read()
    with pytest.raises(ValueError, match="dimensions"):
        _parse_oisst_bytes(payload)


# ---------------------------------------------------------------------------
# fetch_oisst with mocked HTTP
# ---------------------------------------------------------------------------

def test_fetch_oisst_mocked(tmp_path, monkeypatch):
    payload = _make_oisst_payload(str(tmp_path / "oisst.nc"))
    calls = []

    def fake_download(url):
        calls.append(url)
        return payload

    monkeypatch.setattr(sst_global, "_download_bytes", fake_download)
    f = fetch_oisst((-80, 20, -60, 40), "2020-01-01", "2020-01-31",
                    stride_days=30)
    assert isinstance(f, SstField)
    assert f.sst.shape == (3, 5, 6)
    assert f.source == "noaa-oisst/ncdcOisst21Agg"
    assert len(calls) == 1
    assert f.provenance["urls"] == calls
    assert f.provenance["sha256"] == hashlib.sha256(payload).hexdigest()
    assert f.provenance["units"] == "degree_C"


def test_fetch_oisst_antimeridian_two_requests(tmp_path, monkeypatch):
    # each lon window returns its own longitude range, like the real server
    payload_west = _make_oisst_payload(str(tmp_path / "w.nc"), lon0=190.0)
    payload_east = _make_oisst_payload(str(tmp_path / "e.nc"), lon0=0.125)
    calls = []

    def fake_download(url):
        calls.append(url)
        return payload_west if "(190.0)" in url else payload_east

    monkeypatch.setattr(sst_global, "_download_bytes", fake_download)
    f = fetch_oisst((-170, -10, 170, 10), "2020-01-01", "2020-01-31")
    assert len(calls) == 2                      # split lon windows
    assert f.sst.shape == (3, 5, 12)
    assert np.all(np.diff(f.lons) > 0)          # concatenated + sorted
    assert f.lons[0] == pytest.approx(-170.0)
    assert f.provenance["lon_windows"] == [[190.0, 359.875], [0.125, 170.0]]


def test_fetch_oisst_network_error_wrapped(tmp_path, monkeypatch):
    def boom(url):
        raise RuntimeError("download failed after 3 attempts")
    monkeypatch.setattr(sst_global, "_download_bytes", boom)
    with pytest.raises(RuntimeError, match="OISST request failed"):
        fetch_oisst((-80, 20, -60, 40), "2020-01-01", "2020-01-31")


# ---------------------------------------------------------------------------
# SstField model
# ---------------------------------------------------------------------------

def test_sstfield_shape_mismatch_raises():
    with pytest.raises(ValueError, match="does not match"):
        SstField(sst=np.ma.zeros((3, 3, 4)), times=["t1", "t2"],
                 lats=np.zeros(3), lons=np.zeros(4))


def test_sstfield_spatial_mean_ignores_masked():
    sst = np.ma.array([[[1.0, 2.0], [3.0, 4.0]]],
                      mask=[[[False, False], [False, True]]])
    f = SstField(sst=sst, times=["2020-01-01T12:00:00+00:00"],
                 lats=np.array([0.0, 1.0]), lons=np.array([0.0, 1.0]))
    assert f.spatial_mean(0) == pytest.approx(2.0)


def test_sstfield_select_time():
    f = SstField.synthetic(nt=4)
    one = f.select_time(2)
    assert one.sst.shape == (1, 6, 8)
    assert one.times == [f.times[2]]


def test_sstfield_select_bbox():
    f = SstField.synthetic(nt=2, ny=6, nx=8,
                           lats=(-30, 30), lons=(-60, 60))
    sub = f.select_bbox((-60, -30, 0, 0))
    assert sub.sst.shape == (2, 3, 4)
    assert sub.bounds[0] == pytest.approx(-60.0)


def test_sstfield_select_bbox_no_overlap_raises():
    f = SstField.synthetic()
    with pytest.raises(ValueError, match="does not overlap"):
        f.select_bbox((100, 80, 120, 85))


def test_sstfield_json_roundtrip(tmp_path):
    f = SstField.synthetic(nt=3, seed=5)
    p = str(tmp_path / "f.json")
    f.to_json(p)
    g = SstField.from_json(p)
    assert g.sst.shape == f.sst.shape
    assert g.times == f.times
    assert np.allclose(np.ma.filled(g.sst, np.nan),
                       np.ma.filled(f.sst, np.nan), equal_nan=True)
    assert g.provenance == f.provenance


def test_sstfield_from_dict_missing_keys():
    with pytest.raises(ValueError, match="missing keys"):
        SstField.from_dict({"sst": []})


def test_sstfield_synthetic_deterministic():
    a = SstField.synthetic(nt=3, seed=9)
    b = SstField.synthetic(nt=3, seed=9)
    assert np.ma.allequal(a.sst, b.sst)
    assert a.times == b.times


def test_main_demo_runs(capsys):
    sst_global.main_demo()
    out = capsys.readouterr().out
    assert "[oisst]" in out and "[mur]" in out


# ---------------------------------------------------------------------------
# MUR: granule naming, grid parsing, subset URLs
# ---------------------------------------------------------------------------

def test_mur_granule_title_format():
    assert mur_granule_title("2026-09-25") == \
        "20260925090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1"


def test_mur_opendap_url_uses_collection():
    url = mur_opendap_url("2026-09-25")
    assert sst_global.MUR_COLLECTION_ID in url
    assert url.startswith("https://opendap.earthdata.nasa.gov/collections/")


def test_parse_mur_grid_ok():
    grid = _parse_mur_grid(MUR_DAS, MUR_DDS)
    assert grid == {"lat": (5, 20.0, 24.0), "lon": (7, -80.0, -74.0)}


def test_parse_mur_grid_missing_axis_raises():
    with pytest.raises(ValueError, match="missing axis"):
        _parse_mur_grid(MUR_DAS, "Dataset { Float64 lat[lat = 5]; } x;")


def test_mur_subset_urls_constraint():
    grid = _parse_mur_grid(MUR_DAS, MUR_DDS)
    urls = mur_subset_urls("2026-09-25", (-80, 20, -74, 24), grid)
    assert len(urls) == 1
    assert ".nc?analysed_sst[0:1:0][" in urls[0]
    assert "[0:1:4][0:1:6]" in urls[0]  # [time][lat 0..4][lon 0..6]


def test_mur_subset_urls_antimeridian_two():
    das = MUR_DAS.replace("20.0, 24.0", "-10.0, 10.0").replace(
        "-80.0, -74.0", "-179.99, 179.99")
    dds = MUR_DDS.replace("lat = 5", "lat = 21").replace("lon = 7", "lon = 36000")
    grid = _parse_mur_grid(das, dds)
    urls = mur_subset_urls("2026-09-25", (170, -10, -170, 10), grid)
    assert len(urls) == 2


def test_parse_mur_bytes_kelvin_to_celsius(tmp_path):
    payload = _make_mur_payload(str(tmp_path / "mur.nc"), kelvin=True)
    lats, lons, sst = _parse_mur_bytes(payload)
    assert sst.shape == (1, 5, 7)
    assert float(sst[0, 0, 0]) == pytest.approx(16.85, abs=0.01)
    assert lats[0] == pytest.approx(20.0)


# ---------------------------------------------------------------------------
# MUR credentials
# ---------------------------------------------------------------------------

def test_earthdata_credentials_from_env(monkeypatch):
    monkeypatch.setenv("EARTHDATA_USERNAME", "user1")
    monkeypatch.setenv("EARTHDATA_PASSWORD", "pw1")
    assert earthdata_credentials() == ("user1", "pw1")


def test_earthdata_credentials_none(monkeypatch):
    monkeypatch.delenv("EARTHDATA_USERNAME", raising=False)
    monkeypatch.delenv("EARTHDATA_PASSWORD", raising=False)
    monkeypatch.setenv("HOME", "/nonexistent-home-for-test")
    monkeypatch.setenv("NETRC", "/nonexistent-home-for-test/.netrc")
    assert earthdata_credentials() is None


def test_fetch_mur_requires_credentials(monkeypatch):
    monkeypatch.setattr(sst_global, "earthdata_credentials", lambda: None)
    with pytest.raises(CredentialsMissing, match="Earthdata Login"):
        fetch_mur((-80, 20, -74, 24), "2024-01-01", "2024-01-31")


def _fake_cmr_granules(*stamps):
    """Build fake CMR discovery results for date stamps like '20240101'."""
    return [{
        "title": f"{s}090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1",
        "time_start": None,
        "time_end": None,
        "opendap_url": None,
        "links": [],
    } for s in stamps]


def test_fetch_mur_mocked(tmp_path, monkeypatch):
    payload = _make_mur_payload(str(tmp_path / "mur.nc"), kelvin=True)
    urls_seen = []

    def fake_get(opener, url, context):
        urls_seen.append(url)
        if url.endswith(".das"):
            return MUR_DAS.encode()
        if url.endswith(".dds"):
            return MUR_DDS.encode()
        return payload

    monkeypatch.setattr(sst_global, "_earthdata_opener", lambda: object())
    monkeypatch.setattr(
        sst_global, "cmr_search_mur_granules",
        lambda start, end: _fake_cmr_granules("20240101", "20240103",
                                             "20240105"))
    monkeypatch.setattr(sst_global, "_opendap_get", fake_get)
    f = fetch_mur((-80, 20, -74, 24), "2024-01-01", "2024-01-05",
                  stride_days=2)
    assert isinstance(f, SstField)
    assert f.sst.shape == (3, 5, 7)          # 3 granules, degC after conversion
    assert float(f.sst[0, 0, 0]) == pytest.approx(16.85, abs=0.01)
    assert f.source == "nasa-jpl-mur/MUR-JPL-L4-GLOB-v4.1"
    assert f.provenance["n_granules"] == 3
    assert len(f.provenance["granule_titles"]) == 3
    assert f.provenance["granule_titles"][0].startswith("20240101090000-JPL")
    assert any(u.endswith(".das") for u in urls_seen)
    assert f.times[0].endswith("09:00:00+00:00")


def test_fetch_mur_http401_becomes_credentials_missing(monkeypatch):
    import urllib.error

    def fake_opener():
        return object()

    def fake_get(opener, url, context):
        raise sst_global.CredentialsMissing("rejected (HTTP 401)")

    monkeypatch.setattr(sst_global, "_earthdata_opener", fake_opener)
    monkeypatch.setattr(
        sst_global, "cmr_search_mur_granules",
        lambda start, end: _fake_cmr_granules("20240101"))
    monkeypatch.setattr(sst_global, "_opendap_get", fake_get)
    with pytest.raises(CredentialsMissing):
        fetch_mur((-80, 20, -74, 24), "2024-01-01", "2024-01-05")


def test_opendap_get_401_maps_to_credentials_missing(monkeypatch):
    import urllib.error

    class FakeOpener:
        def open(self, req, timeout=None):
            raise urllib.error.HTTPError(req.full_url, 401, "Unauthorized",
                                         {}, io.BytesIO(b""))

    with pytest.raises(CredentialsMissing, match="HTTP 401"):
        sst_global._opendap_get(FakeOpener(), "https://x.example/y", "test")


def test_fetch_mur_rejects_before_record(monkeypatch):
    monkeypatch.setattr(sst_global, "earthdata_credentials",
                        lambda: ("u", "p"))
    with pytest.raises(ValueError, match="2002"):
        fetch_mur((-80, 20, -74, 24), "1990-01-01", "1990-02-01")


# ---------------------------------------------------------------------------
# CMR granule discovery (mocked HTTP; verified live 2026-09-26)
# ---------------------------------------------------------------------------

CMR_FEED = {
    "feed": {"entry": [
        {
            "title": "20260924090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1",
            "time_start": "2026-09-23T21:00:00.000Z",
            "time_end": "2026-09-24T21:00:00.000Z",
            "links": [
                {"rel": "http://esipfed.org/ns/fedsearch/1.1/data#",
                 "href": "https://opendap.earthdata.nasa.gov/collections/"
                         "C1996881146-POCLOUD/granules/20260924090000-JPL-"
                         "L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1"},
                {"rel": "via",
                 "href": "https://cmr.earthdata.nasa.gov/meta"},
            ],
        },
        {
            "title": "20260925090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1",
            "time_start": "2026-09-24T21:00:00.000Z",
            "time_end": "2026-09-25T21:00:00.000Z",
            "links": [{"rel": "data",
                       "href": "https://example.invalid/plain-data-file"}],
        },
    ]}
}


def _mock_cmr(monkeypatch, feed=CMR_FEED):
    seen = {}

    def fake_download(url, timeout=300):
        seen["url"] = url
        return json.dumps(feed).encode()

    monkeypatch.setattr(sst_global, "_download_bytes", fake_download)
    return seen


def test_cmr_search_parses_granules(monkeypatch):
    _mock_cmr(monkeypatch)
    granules = sst_global.cmr_search_mur_granules("2026-09-24", "2026-09-26")
    assert len(granules) == 2
    first, second = granules
    assert first["title"].startswith("20260924090000-JPL")
    assert first["time_start"] == "2026-09-23T21:00:00.000Z"
    assert first["opendap_url"].startswith(
        "https://opendap.earthdata.nasa.gov/collections/C1996881146-POCLOUD")
    assert second["opendap_url"] is None  # no OPeNDAP link advertised


def test_cmr_search_targets_mur_collection(monkeypatch):
    seen = _mock_cmr(monkeypatch)
    sst_global.cmr_search_mur_granules("2026-09-24", "2026-09-26")
    assert "collection_concept_id=C1996881146-POCLOUD" in seen["url"]
    assert "granules.json" in seen["url"]
    assert "temporal=" in seen["url"]


def test_cmr_search_non_json_raises(monkeypatch):
    monkeypatch.setattr(sst_global, "_download_bytes",
                        lambda url, timeout=300: b"<html>not json</html>")
    with pytest.raises(RuntimeError, match="non-JSON"):
        sst_global.cmr_search_mur_granules("2026-09-24", "2026-09-26")


def test_mur_match_granules_by_title_stamp(monkeypatch):
    _mock_cmr(monkeypatch)
    granules = sst_global.cmr_search_mur_granules("2026-09-24", "2026-09-26")
    matched = sst_global.mur_match_granules(
        [dt.date(2026, 9, 24), dt.date(2026, 9, 25)], granules)
    assert [m["title"][:8] for m in matched] == ["20260924", "20260925"]


def test_mur_match_granules_gap_raises(monkeypatch):
    _mock_cmr(monkeypatch)
    granules = sst_global.cmr_search_mur_granules("2026-09-24", "2026-09-26")
    with pytest.raises(RuntimeError, match="no MUR granule for 2026-09-26"):
        sst_global.mur_match_granules([dt.date(2026, 9, 26)], granules)


def test_mur_match_granules_time_start_fallback():
    granules = [{
        "title": "WEIRD-TITLE-NO-STAMP",
        "time_start": "2026-09-24T21:00:00.000Z",
        "time_end": None, "opendap_url": None, "links": [],
    }]
    matched = sst_global.mur_match_granules([dt.date(2026, 9, 24)], granules)
    assert matched[0]["title"] == "WEIRD-TITLE-NO-STAMP"


def test_mur_service_url_prefers_cmr_link():
    url, source = sst_global._mur_service_url({
        "title": "T", "opendap_url": "https://cmr.example.invalid/opendap/x"})
    assert url == "https://cmr.example.invalid/opendap/x"
    assert source == "cmr-link"


def test_mur_service_url_constructed_fallback():
    url, source = sst_global._mur_service_url({
        "title": "20260924090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1",
        "opendap_url": None})
    assert url == sst_global.mur_opendap_url("2026-09-24")
    assert source == "constructed"


def test_mur_subset_urls_accepts_discovered_base_url():
    grid = _parse_mur_grid(MUR_DAS, MUR_DDS)
    urls = sst_global.mur_subset_urls(
        "2026-09-24", (-80, 20, -74, 24), grid,
        base_url="https://cmr.example.invalid/opendap/granule-XYZ")
    assert len(urls) == 1
    assert urls[0].startswith(
        "https://cmr.example.invalid/opendap/granule-XYZ.nc?analysed_sst")


def test_fetch_mur_uses_discovered_cmr_links(tmp_path, monkeypatch):
    """fetch_mur fetches through the CMR-advertised service URLs, not
    constructed names: the discovered link below could never be guessed."""
    payload = _make_mur_payload(str(tmp_path / "mur.nc"), kelvin=True)
    urls_seen = []
    cmr_link = "https://cmr.example.invalid/opendap/granule-XYZ"

    def fake_get(opener, url, context):
        urls_seen.append(url)
        if url.endswith(".das"):
            return MUR_DAS.encode()
        if url.endswith(".dds"):
            return MUR_DDS.encode()
        return payload

    monkeypatch.setattr(sst_global, "_earthdata_opener", lambda: object())
    monkeypatch.setattr(
        sst_global, "cmr_search_mur_granules",
        lambda start, end: [{
            "title": "20240101090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1",
            "time_start": "2023-12-31T21:00:00.000Z",
            "time_end": None,
            "opendap_url": cmr_link,
            "links": [],
        }])
    monkeypatch.setattr(sst_global, "_opendap_get", fake_get)
    f = fetch_mur((-80, 20, -74, 24), "2024-01-01", "2024-01-01")
    subset_urls = [u for u in urls_seen if ".nc?" in u]
    assert subset_urls, "expected subset requests"
    assert all(u.startswith(cmr_link + ".nc?") for u in subset_urls)
    assert f.provenance["opendap_url_source"] == ["cmr-link"]
    assert f.provenance["cmr_search"]["n_results"] == 1
    assert f.times[0].endswith("09:00:00+00:00")  # from title stamp


def test_fetch_mur_discovery_gap_is_honest(monkeypatch):
    monkeypatch.setattr(sst_global, "_earthdata_opener", lambda: object())
    monkeypatch.setattr(sst_global, "cmr_search_mur_granules",
                        lambda start, end: [])
    with pytest.raises(RuntimeError, match="no MUR granule"):
        fetch_mur((-80, 20, -74, 24), "2024-01-01", "2024-01-05")


# ---------------------------------------------------------------------------
# CLI wiring (offline)
# ---------------------------------------------------------------------------

def test_cli_has_new_commands():
    import argparse
    from currents.cli import build_parser
    p = build_parser()
    subs = [a for a in p._actions
            if isinstance(a, argparse._SubParsersAction)]
    choices = subs[0].choices
    for cmd in ("fetch-oisst", "fetch-mur", "sst-synthetic"):
        assert cmd in choices
        assert callable(choices[cmd].get_default("func"))


def test_cli_sst_synthetic_offline(tmp_path, monkeypatch, capsys):
    from currents.cli import main
    monkeypatch.chdir(tmp_path)
    rc = main(["sst-synthetic", "--nt", "2", "--seed", "3",
               "--out", "demo_sst"])
    assert rc == 0
    with open("demo_sst.json") as fh:
        data = json.load(fh)
    assert len(data["times"]) == 2
    assert data["source"] == "synthetic"
    assert "wrote synthetic global SST field" in capsys.readouterr().out

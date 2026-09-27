"""Tests for currents.streamgages (USGS Water Services / NWIS streamgages).

The module is fully offline-testable: the URL builders, RDB/JSON parsers,
and GageField model are pure, and the download seam (``_cached_get_text``)
is mocked for the fetch tests.

Live-network tests are skipped unless ``SURVEY_CURRENTS_LIVE=1``.
"""

import datetime as dt
import json as _json
import os

import numpy as np
import pytest

import currents.streamgages as sg
from currents.streamgages import (GageField, GageRecord, fetch_usgs,
                                  usgs_cache_dir)

BBOX = (-83.5, 42.0, -82.0, 43.5)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

RDB_FIXTURE = """#
# USGS Water Services site service
#
agency_cd\tsite_no\tstation_nm\tsite_tp_cd\tdec_lat_va\tdec_long_va\tdrain_area_va\thuc_cd\ttz_cd\talt_va
5s\t15s\t50s\t7s\t16s\t16s\t8s\t16s\t6s\t8s
USGS\t04159130\tHURON RIVER AT FLAT ROCK MI\tST\t42.10812397\t-83.28017320\t\t04100009\tEST\t
USGS\t04160050\tLONG LAKE OUTLET AT LONG LAKE NY\tST\t43.97600000\t-74.42000000\t51.2\t02020003\tEST\t1490
USGS\t04159500\tFLAT ROCK CREEK NEAR OAKWOOD MI\tST\t42.00000000\t-83.40000000\t\t04100009\tEDT\t
"""


def _dv_fixture() -> str:
    """dv JSON: station 04159130 with 5 days of 00060 + 00065,
    station 04159500 with a missing day (sentinel) on 00060, and no
    timeSeries for 04160050 (no data in window)."""
    def _ts(site_no, name, param, unit, values):
        return {
            "sourceInfo": {
                "siteName": name,
                "siteCode": [{"value": site_no, "agencyCode": "USGS"}],
            },
            "variable": {
                "variableCode": [{"value": param}],
                "unit": {"unitCode": unit},
            },
            "values": [{
                "value": [
                    {"value": v, "qualifiers": [q],
                     "dateTime": f"2026-09-{day:02d}T00:00:00.000"}
                    for (day, v, q) in values
                ]
            }],
        }
    body = {"value": {"timeSeries": [
        _ts("04159130", "HURON RIVER AT FLAT ROCK MI", "00060", "ft3/s",
            [(20, "312", "P"), (21, "298", "P"), (22, "285", "P"),
             (23, "277", "P"), (24, "271", "P")]),
        _ts("04159130", "HURON RIVER AT FLAT ROCK MI", "00065", "ft",
            [(20, "4.12", "P"), (21, "4.05", "P"), (22, "3.99", "P"),
             (23, "3.95", "P"), (24, "3.92", "P")]),
        _ts("04159500", "FLAT ROCK CREEK NEAR OAKWOOD MI", "00060", "ft3/s",
            [(20, "41", "A"), (21, "-999999", "P"), (22, "39", "A"),
             (23, "38", "A"), (24, "37", "A")]),
    ]}}
    return _json.dumps(body)


def _mock_usgs(monkeypatch, tmp_path, rdb=RDB_FIXTURE, dv=None,
               missing=False):
    """Mock the download seam for fetch_usgs.

    ``rdb``: site-service response text (default fixture).
    ``dv``: dv-service response text (default built fixture).
    ``missing``: when True, dv returns no timeSeries at all.
    """
    dv_text = (_json.dumps({"value": {"timeSeries": []}})
               if missing else (dv if dv is not None else _dv_fixture()))
    calls = []

    def fake_cached_get(url, cache_dir, max_age_days, refresh=False):
        calls.append(url)
        if "nwis/site/" in url:
            return rdb, {"url": url, "sha256": "abc", "cache_hit": False,
                         "downloaded": True}
        return dv_text, {"url": url, "sha256": "def", "cache_hit": False,
                         "downloaded": True}

    monkeypatch.setattr(sg, "_cached_get_text", fake_cached_get)
    monkeypatch.setattr(sg, "usgs_cache_dir", lambda: str(tmp_path))
    return calls


# ---------------------------------------------------------------------------
# URL builders
# ---------------------------------------------------------------------------

def test_site_inventory_url():
    url = sg.site_inventory_url(BBOX)
    assert url.startswith(sg.USGS_SITE_BASE)
    assert "format=rdb" in url
    assert "hasDataTypeCd=dv" in url
    assert "parameterCd=00060" in url
    assert "siteOutput=expanded" in url
    assert "-83.5,42.0,-82.0,43.5" in url


def test_site_inventory_url_multi_parameter():
    url = sg.site_inventory_url(BBOX, ("00060", "00065"))
    assert "parameterCd=00060,00065" in url


def test_dv_request_url():
    url = sg.dv_request_url(["04159130", "04160050"], "2026-09-20",
                            "2026-09-24")
    assert url.startswith(sg.USGS_DV_BASE)
    assert "format=json" in url
    assert "sites=04159130,04160050" in url
    assert "startDT=2026-09-20" in url
    assert "endDT=2026-09-24" in url
    assert "parameterCd=00060" in url
    # Daily MEAN statistic must be explicit (dv default, stated for
    # determinism).
    assert "statCd=00003" in url


# ---------------------------------------------------------------------------
# Parsers (pure, offline)
# ---------------------------------------------------------------------------

def test_parse_site_rdb():
    sites = sg._parse_site_rdb(RDB_FIXTURE)
    assert len(sites) == 3
    first = sites[0]
    assert first["site_no"] == "04159130"
    assert first["site_name"] == "HURON RIVER AT FLAT ROCK MI"
    assert first["lat"] == pytest.approx(42.10812397)
    assert first["lon"] == pytest.approx(-83.28017320)
    assert first["drain_area_sqmi"] is None  # empty cell -> None
    assert first["huc"] == "04100009"
    second = sites[1]
    assert second["drain_area_sqmi"] == pytest.approx(51.2)
    assert second["alt_ft"] == pytest.approx(1490.0)


def test_parse_site_rdb_empty():
    assert sg._parse_site_rdb("#\n# no sites\n") == []
    assert sg._parse_site_rdb("") == []


def test_parse_dv_json():
    parsed = sg._parse_dv_json(_json.loads(_dv_fixture()))
    q60 = parsed["04159130"]["00060"]
    assert q60["unit"] == "ft3/s"
    assert len(q60["dates"]) == 5
    assert q60["dates"][0] == dt.date(2026, 9, 20)
    assert list(q60["values"]) == [312.0, 298.0, 285.0, 277.0, 271.0]
    assert q60["qualifiers"][0] == ["P"]
    gh = parsed["04159130"]["00065"]
    assert gh["unit"] == "ft"
    # Missing sentinel -> NaN, recorded with its qualifier.
    gap = parsed["04159500"]["00060"]
    assert np.isnan(gap["values"][1])
    assert gap["qualifiers"][1] == ["P"]
    assert gap["values"][0] == pytest.approx(41.0)


def test_parse_dv_json_empty_timeseries_skipped():
    parsed = sg._parse_dv_json({"value": {"timeSeries": []}})
    assert parsed == {}


def test_parse_dv_json_bad_date_skipped():
    payload = {"value": {"timeSeries": [{
        "sourceInfo": {"siteCode": [{"value": "04159130"}]},
        "variable": {"variableCode": [{"value": "00060"}],
                     "unit": {"unitCode": "ft3/s"}},
        "values": [{"value": [
            {"value": "100", "qualifiers": ["P"],
             "dateTime": "2026-09-20T00:00:00.000"},
            {"value": "200", "qualifiers": ["P"], "dateTime": "not-a-date"},
        ]}],
    }]}}
    parsed = sg._parse_dv_json(payload)
    assert len(parsed["04159130"]["00060"]["dates"]) == 1


# ---------------------------------------------------------------------------
# GageRecord
# ---------------------------------------------------------------------------

def _record():
    dates = [dt.date(2026, 9, 20) + dt.timedelta(days=k) for k in range(5)]
    return GageRecord(
        site_no="04159130", site_name="HURON RIVER AT FLAT ROCK MI",
        lat=42.1, lon=-83.28, huc="04100009", drain_area_sqmi=None,
        tz="EST",
        series={"00060": {
            "dates": dates,
            "values": np.array([312., 298., np.nan, 277., 271.]),
            "unit": "ft3/s",
            "qualifiers": [["P"], ["P"], ["P"], ["P"], ["P"]]}})


def test_gagerecord_latest():
    rec = _record()
    d, v = rec.latest()
    assert d == dt.date(2026, 9, 24)
    assert v == pytest.approx(271.0)
    assert rec.n_valid() == 4


def test_gagerecord_latest_all_nan():
    rec = _record()
    rec.series["00060"]["values"] = np.full(5, np.nan)
    assert rec.latest() == (None, None)
    assert rec.percentile_of_record() is None


def test_gagerecord_missing_days():
    rec = _record()
    missing = rec.missing_days()
    assert missing == [dt.date(2026, 9, 22)]
    # Bounded range honours gaps outside the series range too.
    missing = rec.missing_days(start="2026-09-20", end="2026-09-26")
    assert missing == [dt.date(2026, 9, 22), dt.date(2026, 9, 25),
                       dt.date(2026, 9, 26)]


def test_gagerecord_percentile_of_record():
    rec = _record()
    # Latest value 271 is the lowest of the 4 valid values -> 0.0.
    assert rec.percentile_of_record() == pytest.approx(0.0)
    # Single valid value -> 50 by convention.
    rec.series["00060"]["values"] = np.array(
        [np.nan, np.nan, np.nan, np.nan, 271.0])
    assert rec.percentile_of_record() == pytest.approx(50.0)


def test_gagerecord_to_si():
    rec = _record()
    si = rec.to_si()
    assert si.series["00060"]["unit"] == "m3/s"
    assert si.series["00060"]["values"][0] == pytest.approx(
        312.0 * sg.FT3S_TO_M3S)
    assert np.isnan(si.series["00060"]["values"][2])
    assert rec.series["00060"]["unit"] == "ft3/s"  # original untouched


def test_gagerecord_to_si_gage_height():
    dates = [dt.date(2026, 9, 20)]
    rec = GageRecord(site_no="x", site_name="x", lat=0.0, lon=0.0,
                     series={"00065": {"dates": dates,
                                       "values": np.array([4.0]),
                                       "unit": "ft",
                                       "qualifiers": [[]]}})
    assert rec.to_si().series["00065"]["unit"] == "m"
    assert rec.to_si().series["00065"]["values"][0] == pytest.approx(1.2192)


def test_gagerecord_roundtrip_nan():
    rec = _record()
    rt = GageRecord.from_dict(rec.to_dict())
    assert rt.site_no == "04159130"
    assert np.isnan(rt.series["00060"]["values"][2])
    assert rt.series["00060"]["values"][0] == pytest.approx(312.0)
    assert rt.drain_area_sqmi is None


# ---------------------------------------------------------------------------
# GageField
# ---------------------------------------------------------------------------

def _field():
    return GageField.synthetic(n_sites=3, parameters=("00060", "00065"))


def test_synthetic_deterministic():
    a = GageField.synthetic(seed=11)
    b = GageField.synthetic(seed=11)
    assert a.site_numbers() == b.site_numbers()
    va = a.records[0].series["00060"]["values"]
    vb = b.records[0].series["00060"]["values"]
    np.testing.assert_allclose(va, vb, equal_nan=True)


def test_synthetic_named_site_and_gap():
    f = GageField.synthetic(n_sites=4, parameters=("00060",))
    assert f.site_numbers()[0] == "01111110"
    assert f.records[0].site_name == "SYNTHETIC RIVER AT TESTVILLE"
    gap = f.records[1].missing_days()
    assert len(gap) == 10  # the engineered gap block
    assert f.source == "synthetic"


def test_synthetic_00065_included():
    f = GageField.synthetic(n_sites=2, parameters=("00060", "00065"))
    assert "00065" in f.records[0].series
    assert f.records[0].series["00065"]["unit"] == "ft"


def test_gagefield_regional_median():
    f = GageField.synthetic(n_sites=2, parameters=("00060",), seed=11)
    dates, med = f.regional_median()
    assert len(dates) == (f.end - f.start).days + 1
    va = f.records[0].series["00060"]["values"]
    vb = f.records[1].series["00060"]["values"]
    expected = np.nanmedian(np.vstack([va, vb]), axis=0)
    np.testing.assert_allclose(med, expected, equal_nan=True)
    # The engineered gap (site 2, days 20..29) does not zero the median:
    assert np.isfinite(med[25])


def test_gagefield_select_site():
    f = _field()
    one = f.select_site(f.site_numbers()[0])
    assert len(one) == 1
    assert one.records[0].site_no == f.site_numbers()[0]
    assert len(f.select_site("00000000")) == 0


def test_gagefield_to_si():
    f = _field()
    si = f.to_si()
    assert si.units == "si"
    assert si.unit_label == "m³/s"
    v0 = f.records[0].series["00060"]["values"][0]
    assert si.records[0].series["00060"]["values"][0] == pytest.approx(
        v0 * sg.FT3S_TO_M3S)
    assert f.unit_label == "ft³/s"  # original untouched


def test_gagefield_json_roundtrip():
    f = _field()
    tmp = "/tmp/_test_gage_roundtrip.json"
    f.to_json(tmp)
    rt = GageField.from_json(tmp)
    assert rt.site_numbers() == f.site_numbers()
    va = f.records[1].series["00060"]["values"]
    vb = rt.records[1].series["00060"]["values"]
    np.testing.assert_allclose(va, vb, equal_nan=True)  # NaN survives
    assert rt.units == f.units
    os.remove(tmp)


def test_gagefield_bad_dates():
    with pytest.raises(ValueError):
        GageField(records=[], bbox=BBOX, start="2024-08-31", end="2024-06-01")


# ---------------------------------------------------------------------------
# fetch_usgs (mocked download seam)
# ---------------------------------------------------------------------------

def test_fetch_usgs_mocked(monkeypatch, tmp_path):
    calls = _mock_usgs(monkeypatch, tmp_path)
    field = fetch_usgs(BBOX, "2026-09-20", "2026-09-24",
                       parameters=("00060", "00065"), min_record_days=3)
    # 04159130 (5 days x 2 params) and 04159500 (5 days x 1 param) kept;
    # 04160050 has no timeSeries -> counted as no-data.
    assert sorted(field.site_numbers()) == ["04159130", "04159500"]
    rec = field.records[0]
    assert rec.site_name == "HURON RIVER AT FLAT ROCK MI"
    assert rec.lat == pytest.approx(42.10812397)
    assert set(rec.series) == {"00060", "00065"}
    assert rec.n_valid("00060") == 5
    assert len(rec.missing_days("00060")) == 0
    # Sentinel day on 04159500 recorded as a missing day, never filled.
    gap_rec = [r for r in field.records
               if r.site_no == "04159500"][0]
    assert gap_rec.n_valid("00060") == 4
    assert gap_rec.missing_days("00060") == [dt.date(2026, 9, 21)]

    prov = field.provenance
    assert prov["source"] == "usgs"
    assert prov["n_sites_discovered"] == 3
    assert prov["n_sites_requested"] == 3
    assert prov["n_sites_with_data"] == 2
    assert prov["n_sites_no_data"] == 1
    assert prov["n_sites_excluded_short_record"] == 0
    assert prov["missing_days_total"] == 1
    assert "00003 (daily mean" in prov["statistic"]
    assert any("nwis/site/" in u for u in prov["request_urls"])
    assert any("nwis/dv/" in u for u in prov["request_urls"])
    assert prov["units"] == "native"
    assert len(calls) == 2  # one inventory + one dv batch


def test_fetch_usgs_min_record_days_excludes(monkeypatch, tmp_path):
    _mock_usgs(monkeypatch, tmp_path)
    field = fetch_usgs(BBOX, "2026-09-20", "2026-09-24",
                       min_record_days=5)
    assert field.site_numbers() == ["04159130"]  # 4-day site excluded
    assert field.provenance["n_sites_excluded_short_record"] == 1


def test_fetch_usgs_site_limit(monkeypatch, tmp_path):
    _mock_usgs(monkeypatch, tmp_path)
    field = fetch_usgs(BBOX, "2026-09-20", "2026-09-24",
                       min_record_days=3, site_limit=1)
    assert field.site_numbers() == ["04159130"]  # sorted, first only


def test_fetch_usgs_empty_field(monkeypatch, tmp_path):
    _mock_usgs(monkeypatch, tmp_path, missing=True)
    field = fetch_usgs(BBOX, "2026-09-20", "2026-09-24")
    assert len(field) == 0
    assert "empty_reason" in field.provenance
    assert "US-only" in field.provenance["empty_reason"]


def test_fetch_usgs_empty_inventory(monkeypatch, tmp_path):
    _mock_usgs(monkeypatch, tmp_path, rdb="#\n")
    field = fetch_usgs(BBOX, "2026-09-20", "2026-09-24")
    assert len(field) == 0
    assert field.provenance["n_sites_discovered"] == 0


def test_fetch_usgs_si_units(monkeypatch, tmp_path):
    _mock_usgs(monkeypatch, tmp_path)
    field = fetch_usgs(BBOX, "2026-09-20", "2026-09-24",
                       min_record_days=3, units="si")
    assert field.units == "si"
    assert field.unit_label == "m³/s"
    assert field.records[0].series["00060"]["values"][0] == pytest.approx(
        312.0 * sg.FT3S_TO_M3S)
    assert field.provenance["units"] == "si"


def test_fetch_usgs_validation(monkeypatch, tmp_path):
    _mock_usgs(monkeypatch, tmp_path)
    with pytest.raises(ValueError, match="unknown USGS parameter"):
        fetch_usgs(BBOX, "2026-09-20", "2026-09-24",
                   parameters=("99999",))
    with pytest.raises(ValueError, match="units must be"):
        fetch_usgs(BBOX, "2026-09-20", "2026-09-24", units="furlongs")
    with pytest.raises(ValueError, match="before start"):
        fetch_usgs(BBOX, "2026-09-24", "2026-09-20")
    with pytest.raises(ValueError, match="in the future"):
        fetch_usgs(BBOX, "2026-09-20", "2099-01-01")
    with pytest.raises(ValueError, match="predates the NWIS"):
        fetch_usgs(BBOX, "1800-01-01", "1800-12-31")
    with pytest.raises(ValueError, match="site_limit"):
        fetch_usgs(BBOX, "2026-09-20", "2026-09-24", site_limit=0)


# ---------------------------------------------------------------------------
# Cache behaviour (mocked download)
# ---------------------------------------------------------------------------

def test_usgs_cache_dir_env(monkeypatch, tmp_path):
    monkeypatch.setenv("SURVEY_CURRENTS_CACHE", str(tmp_path))
    d = usgs_cache_dir()
    assert d.endswith(os.path.join("streamgages"))
    assert os.path.isdir(d)


def test_cached_get_text_caches(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(
        sg, "_download_bytes",
        lambda url, timeout=120: calls.append(url) or b"site-data")
    text, prov = sg._cached_get_text("https://example.test/inv", str(tmp_path),
                                     max_age_days=7)
    assert text == "site-data"
    assert prov["downloaded"] is True
    assert len(prov["sha256"]) == 64
    # Second call: cache hit, no download.
    text2, prov2 = sg._cached_get_text("https://example.test/inv",
                                       str(tmp_path), max_age_days=7)
    assert text2 == "site-data"
    assert prov2["cache_hit"] is True
    assert len(calls) == 1


def test_cached_get_text_corrupt_redownload(monkeypatch, tmp_path):
    monkeypatch.setattr(sg, "_download_bytes",
                        lambda url, timeout=120: b"v2")
    path = str(tmp_path)
    t1, _ = sg._cached_get_text("https://example.test/inv", path,
                                max_age_days=7)
    assert t1 == "v2"
    key = sg._cache_key("https://example.test/inv")
    with open(f"{path}/{key}.txt", "w", encoding="utf-8") as fh:
        fh.write("corrupted")
    t2, prov2 = sg._cached_get_text("https://example.test/inv", path,
                                    max_age_days=7)
    assert t2 == "v2"
    assert prov2["downloaded"] is True


def test_cached_get_text_max_age_forces_refetch(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(
        sg, "_download_bytes",
        lambda url, timeout=120: calls.append(url) or b"data")
    sg._cached_get_text("https://example.test/inv", str(tmp_path),
                        max_age_days=7)
    sg._cached_get_text("https://example.test/inv", str(tmp_path),
                        max_age_days=0)
    assert len(calls) == 2


# ---------------------------------------------------------------------------
# CLI (offline)
# ---------------------------------------------------------------------------

def test_cli_usgs_synthetic(tmp_path):
    from currents import cli
    ns = type("NS", (), {"bbox": "-83.5,42.0,-82.0,43.5",
                         "start": "2024-06-01", "end": "2024-06-10",
                         "n_sites": 2, "parameters": "00060,00065",
                         "seed": 11, "source": "synthetic",
                         "out": str(tmp_path / "syn")})()
    assert cli.cmd_usgs_synthetic(ns) == 0
    f = GageField.from_json(str(tmp_path / "syn.json"))
    assert len(f) == 2
    assert "00065" in f.records[0].series


def test_cli_fetch_usgs_mocked(tmp_path, monkeypatch):
    from currents import cli
    _mock_usgs(monkeypatch, tmp_path)
    ns = type("NS", (), {"bbox": "-83.5,42.0,-82.0,43.5",
                         "start": "2026-09-20", "end": "2026-09-24",
                         "parameters": "00060", "min_record_days": 3,
                         "site_limit": 200, "units": "native",
                         "out": str(tmp_path / "usgs")})()
    assert cli.cmd_fetch_usgs(ns) == 0
    f = GageField.from_json(str(tmp_path / "usgs.json"))
    assert len(f) == 2


# ---------------------------------------------------------------------------
# Live (skipped unless SURVEY_CURRENTS_LIVE=1)
# ---------------------------------------------------------------------------

_needs_live = pytest.mark.skipif(
    os.environ.get("SURVEY_CURRENTS_LIVE") != "1",
    reason="live USGS download (set SURVEY_CURRENTS_LIVE=1)")


@_needs_live
class TestLiveUsgs:
    def test_fetch_huron_bbox(self, tmp_path):
        field = fetch_usgs((-83.5, 42.0, -83.0, 42.5), "2026-09-20",
                           "2026-09-24", min_record_days=3,
                           cache_dir=str(tmp_path))
        assert len(field) >= 1
        assert field.provenance["n_sites_discovered"] >= 1
        rec = field.records[0]
        assert rec.site_no and rec.site_name
        assert np.isfinite(rec.lat) and np.isfinite(rec.lon)

    def test_fetch_outside_us_empty(self, tmp_path):
        # NWIS is US-only: an empty mid-Atlantic bbox yields an honest
        # empty field, never fabricated data.
        field = fetch_usgs((-50.0, 35.0, -49.0, 36.0), "2026-09-20",
                           "2026-09-24", cache_dir=str(tmp_path))
        assert len(field) == 0
        assert "empty_reason" in field.provenance

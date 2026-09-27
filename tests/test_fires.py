"""Offline tests for currents.fires (NASA FIRMS active-fire adapter).

All tests run with no network and no FIRMS MAP_KEY: the HTTP layer is
faked by monkeypatching ``currents.fires._download_bytes``, and the
NRT/SP product rule takes an explicit ``today``.
"""

from __future__ import annotations

import datetime as _dt
import json

import numpy as np
import pytest

from currents import fires
from currents.fires import (
    FIRMS_INSTRUMENTS,
    CredentialsMissing,
    FireField,
    fetch_firms,
    firms_area_url,
    firms_area_windows,
    firms_map_key,
    firms_product_for,
    firms_window_chunks,
    redact_map_key,
    validate_firms_bbox,
)


# ---------------------------------------------------------------------------
# sample payloads
# ---------------------------------------------------------------------------

VIIRS_CSV = """latitude,longitude,bright_ti4,bright_ti5,frp,scan,track,acq_date,acq_time,satellite,instrument,confidence,version,daynight
34.0522,-118.2437,320.5,290.1,12.4,0.4,0.4,2024-08-01,1930,Suomi-NPP,VIIRS,h,2.0NRT,D
34.0600,-118.2500,335.2,295.0,45.8,0.4,0.4,2024-08-01,1930,Suomi-NPP,VIIRS,h,2.0NRT,D
35.1000,-119.0000,305.0,285.2,3.1,0.5,0.5,2024-08-02,0215,NOAA-20,VIIRS,n,2.0NRT,N
"""

MODIS_CSV = """latitude,longitude,brightness,scan,track,acq_date,acq_time,satellite,instrument,confidence,version,bright_t31,frp,daynight
-3.1000,-60.0200,312.4,1.0,1.0,2024-08-01,1430,Terra,MODIS,85,6.1NRT,298.1,22.5,D
"""

ERROR_PAYLOAD = "Invalid MAP_KEY."
EMPTY_CSV = "latitude,longitude,bright_ti4,bright_ti5,frp,scan,track,acq_date,acq_time,satellite,instrument,confidence,version,daynight\n"

BBOX = (-125.0, 32.0, -114.0, 42.0)


def _fake_download_factory(payloads):
    """Return a fake _download_bytes serving payloads in order."""
    calls = []

    def fake(url, timeout=300):
        calls.append(url)
        idx = min(len(calls) - 1, len(payloads) - 1)
        return payloads[idx].encode("utf-8")

    fake.calls = calls
    return fake


@pytest.fixture
def viirs(monkeypatch):
    fake = _fake_download_factory([VIIRS_CSV])
    monkeypatch.setattr(fires, "_download_bytes", fake)
    monkeypatch.setenv("FIRMS_MAP_KEY", "TESTKEY123")
    return fake


# ---------------------------------------------------------------------------
# pure helpers
# ---------------------------------------------------------------------------

class TestProductTiering:
    def test_recent_date_uses_nrt(self):
        today = _dt.date(2024, 8, 15)
        assert firms_product_for("VIIRS_SNPP", _dt.date(2024, 8, 1),
                                 today=today) == "VIIRS_SNPP_NRT"

    def test_old_date_uses_sp(self):
        today = _dt.date(2024, 8, 15)
        assert firms_product_for("VIIRS_SNPP", _dt.date(2024, 1, 10),
                                 today=today) == "VIIRS_SNPP_SP"

    def test_modis_products(self):
        today = _dt.date(2024, 8, 15)
        assert firms_product_for("MODIS", _dt.date(2024, 8, 10),
                                 today=today) == "MODIS_NRT"
        assert firms_product_for("modis", _dt.date(2020, 5, 5),
                                 today=today) == "MODIS_SP"

    def test_unknown_instrument(self):
        with pytest.raises(ValueError, match="unknown FIRMS instrument"):
            firms_product_for("GOES", _dt.date(2024, 1, 1),
                              today=_dt.date(2024, 8, 15))

    def test_before_record_start(self):
        with pytest.raises(ValueError, match="no data before"):
            firms_product_for("VIIRS_SNPP", _dt.date(2010, 1, 1),
                              today=_dt.date(2024, 8, 15))

    def test_instruments_table_sane(self):
        for inst, fam in FIRMS_INSTRUMENTS.items():
            assert fam["nrt"].endswith("_NRT")
            assert fam["sp"].endswith("_SP")


class TestWindows:
    def test_single_window(self):
        chunks = firms_window_chunks(_dt.date(2024, 8, 1), _dt.date(2024, 8, 3))
        assert chunks == [(_dt.date(2024, 8, 1), _dt.date(2024, 8, 3))]

    def test_five_day_boundary(self):
        chunks = firms_window_chunks(_dt.date(2024, 8, 1), _dt.date(2024, 8, 5))
        assert chunks == [(_dt.date(2024, 8, 1), _dt.date(2024, 8, 5))]

    def test_six_days_splits(self):
        chunks = firms_window_chunks(_dt.date(2024, 8, 1), _dt.date(2024, 8, 6))
        assert chunks == [(_dt.date(2024, 8, 1), _dt.date(2024, 8, 5)),
                          (_dt.date(2024, 8, 6), _dt.date(2024, 8, 6))]

    def test_long_range(self):
        chunks = firms_window_chunks(_dt.date(2024, 8, 1), _dt.date(2024, 8, 31))
        assert all((b - a).days + 1 <= 5 for a, b in chunks)
        assert chunks[0][0] == _dt.date(2024, 8, 1)
        assert chunks[-1][1] == _dt.date(2024, 8, 31)

    def test_reversed_dates(self):
        with pytest.raises(ValueError, match="before start"):
            firms_window_chunks(_dt.date(2024, 8, 5), _dt.date(2024, 8, 1))


class TestBbox:
    def test_valid(self):
        assert validate_firms_bbox(BBOX) == BBOX

    def test_antimeridian_allowed(self):
        assert validate_firms_bbox((170.0, 40.0, -170.0, 50.0)) == (
            170.0, 40.0, -170.0, 50.0)

    def test_bad_lat(self):
        with pytest.raises(ValueError):
            validate_firms_bbox((-125.0, 42.0, -114.0, 32.0))

    def test_area_windows_split(self):
        wins = firms_area_windows((170.0, 40.0, -170.0, 50.0))
        assert wins == [(170.0, 40.0, 180.0, 50.0),
                        (-180.0, 40.0, -170.0, 50.0)]

    def test_area_windows_passthrough(self):
        assert firms_area_windows(BBOX) == [BBOX]


class TestUrl:
    def test_url_shape(self):
        url = firms_area_url("KEY", "VIIRS_SNPP_NRT", BBOX, 5,
                             _dt.date(2024, 8, 1))
        assert url == ("https://firms.modaps.eosdis.nasa.gov/api/area/csv/"
                       "KEY/VIIRS_SNPP_NRT/"
                       "-125.0000,32.0000,-114.0000,42.0000/5/2024-08-01")

    def test_redaction(self):
        url = firms_area_url("SECRETKEY", "VIIRS_SNPP_NRT", BBOX, 5,
                             _dt.date(2024, 8, 1))
        red = redact_map_key(url, "SECRETKEY")
        assert "SECRETKEY" not in red
        assert "<redacted>" in red
        assert "VIIRS_SNPP_NRT" in red


class TestMapKey:
    def test_env_key(self, monkeypatch):
        monkeypatch.setenv("FIRMS_MAP_KEY", "  abc123  ")
        assert firms_map_key() == "abc123"

    def test_explicit_wins(self, monkeypatch):
        monkeypatch.setenv("FIRMS_MAP_KEY", "envkey")
        assert firms_map_key("explicit") == "explicit"

    def test_missing_raises(self, monkeypatch):
        monkeypatch.delenv("FIRMS_MAP_KEY", raising=False)
        with pytest.raises(CredentialsMissing, match="FIRMS_MAP_KEY"):
            firms_map_key()


# ---------------------------------------------------------------------------
# fetch_firms (faked HTTP)
# ---------------------------------------------------------------------------

class TestFetchFirms:
    def test_basic(self, viirs):
        today = _dt.date.today()
        d0, d1 = today - _dt.timedelta(days=5), today - _dt.timedelta(days=3)
        field = fetch_firms(BBOX, d0.isoformat(), d1.isoformat())
        assert len(field) == 3
        assert field.source == "firms"
        assert field.instruments == ("VIIRS_SNPP",)
        assert field.products == (firms_product_for("VIIRS_SNPP", d0),)
        # oldest-first
        assert field.times == sorted(field.times)
        assert field.frp[1] == pytest.approx(45.8)
        assert field.confidence[0] == "h"
        assert field.satellite[2] == "NOAA-20"
        assert field.daynight[2] == "N"
        # brightness comes from bright_ti4 for VIIRS
        assert field.brightness[0] == pytest.approx(320.5)

    def test_modis_schema(self, monkeypatch):
        monkeypatch.setattr(fires, "_download_bytes",
                            _fake_download_factory([MODIS_CSV]))
        monkeypatch.setenv("FIRMS_MAP_KEY", "TESTKEY123")
        today = _dt.date.today()
        d0 = (today - _dt.timedelta(days=4)).isoformat()
        field = fetch_firms((-65.0, -8.0, -55.0, 2.0), d0,
                            d0, instruments=("MODIS",))
        assert len(field) == 1
        assert field.brightness[0] == pytest.approx(312.4)
        assert field.confidence[0] == "85"
        assert field.products == (
            firms_product_for("MODIS", _dt.date.fromisoformat(d0)),)

    def test_multi_window_request_count(self, monkeypatch):
        fake = _fake_download_factory([VIIRS_CSV, VIIRS_CSV])
        monkeypatch.setattr(fires, "_download_bytes", fake)
        monkeypatch.setenv("FIRMS_MAP_KEY", "TESTKEY123")
        field = fetch_firms(BBOX, "2024-08-01", "2024-08-07")
        # 7 days -> two 5-day windows, one request each
        assert len(fake.calls) == 2
        assert "/5/2024-08-01" in fake.calls[0]
        assert "/2/2024-08-06" in fake.calls[1]
        assert field.provenance["n_requests"] == 2

    def test_nrt_sp_tiering_across_chunks(self, monkeypatch):
        fake = _fake_download_factory([VIIRS_CSV] * 3)
        monkeypatch.setattr(fires, "_download_bytes", fake)
        monkeypatch.setenv("FIRMS_MAP_KEY", "TESTKEY123")
        # Span the NRT/SP boundary relative to "today" is hard without
        # time travel; instead force products explicitly.
        field = fetch_firms(BBOX, "2024-08-01", "2024-08-07",
                            products=("VIIRS_SNPP_SP",))
        assert field.products == ("VIIRS_SNPP_SP",)
        assert "VIIRS_SNPP_SP" in fake.calls[0]

    def test_empty_result(self, monkeypatch):
        monkeypatch.setattr(fires, "_download_bytes",
                            _fake_download_factory([EMPTY_CSV]))
        monkeypatch.setenv("FIRMS_MAP_KEY", "TESTKEY123")
        field = fetch_firms(BBOX, "2024-08-01", "2024-08-01")
        assert len(field) == 0
        assert field.time_range == (None, None)

    def test_invalid_key_error(self, monkeypatch):
        monkeypatch.setattr(fires, "_download_bytes",
                            _fake_download_factory([ERROR_PAYLOAD]))
        monkeypatch.setenv("FIRMS_MAP_KEY", "BADKEY")
        with pytest.raises(ValueError, match="did not return CSV"):
            fetch_firms(BBOX, "2024-08-01", "2024-08-01")

    def test_missing_key(self, monkeypatch):
        monkeypatch.delenv("FIRMS_MAP_KEY", raising=False)
        with pytest.raises(CredentialsMissing):
            fetch_firms(BBOX, "2024-08-01", "2024-08-01")

    def test_provenance_redacted(self, viirs):
        field = fetch_firms(BBOX, "2024-08-01", "2024-08-02")
        prov = field.provenance
        assert prov["map_key"] == "<redacted>"
        for req in prov["requests"]:
            assert "TESTKEY123" not in req["url"]
            assert "<redacted>" in req["url"]
            assert len(req["sha256"]) == 64
        assert prov["bbox"] == list(BBOX)
        assert prov["time_window"] == ["2024-08-01", "2024-08-02"]

    def test_before_coverage(self, monkeypatch):
        monkeypatch.setenv("FIRMS_MAP_KEY", "TESTKEY123")
        with pytest.raises(ValueError, match="coverage starts"):
            fetch_firms(BBOX, "1999-01-01", "1999-02-01")

    def test_unknown_instrument(self, monkeypatch):
        monkeypatch.setenv("FIRMS_MAP_KEY", "TESTKEY123")
        with pytest.raises(ValueError, match="unknown FIRMS instrument"):
            fetch_firms(BBOX, "2024-08-01", "2024-08-02",
                        instruments=("GOES16",))


# ---------------------------------------------------------------------------
# FireField model
# ---------------------------------------------------------------------------

class TestFireField:
    def _field(self):
        return FireField.synthetic(bbox=BBOX, start="2024-08-01",
                                   end="2024-08-03", n=40, seed=7)

    def test_synthetic_deterministic(self):
        a = FireField.synthetic(seed=7)
        b = FireField.synthetic(seed=7)
        assert np.array_equal(a.lats, b.lats)
        assert a.times == b.times

    def test_len_and_bounds(self):
        f = self._field()
        assert len(f) == 40
        assert f.bounds == BBOX
        lo, hi = f.time_range
        assert lo.date() >= _dt.date(2024, 8, 1)
        assert hi.date() <= _dt.date(2024, 8, 3)

    def test_select_time(self):
        f = self._field()
        sub = f.select_time("2024-08-02", "2024-08-02")
        assert all(t.date() == _dt.date(2024, 8, 2) for t in sub.times)
        assert 0 < len(sub) <= len(f)

    def test_select_bbox(self):
        f = self._field()
        sub = f.select_bbox((-125.0, 32.0, -120.0, 42.0))
        assert all(v <= -120.0 for v in sub.lons)

    def test_density_grid_counts(self):
        f = FireField(
            times=[_dt.datetime(2024, 8, 1, 12, tzinfo=_dt.timezone.utc)] * 3,
            lats=np.array([34.0, 34.05, 36.0]),
            lons=np.array([-120.0, -119.95, -118.0]),
            brightness=np.array([320.0, 330.0, 310.0]),
            frp=np.array([10.0, 20.0, 5.0]),
            confidence=["h", "h", "n"],
            satellite=["Suomi-NPP"] * 3,
            instrument=["VIIRS"] * 3,
            daynight=["D"] * 3,
            bbox=(-121.0, 33.0, -117.0, 37.0),
        )
        grid = f.to_density_grid(resolution=1.0)
        assert grid["values"].shape == (1, 4, 4)
        assert grid["values"].sum() == 3.0
        assert grid["weighting"] == "count"
        # the two close-together detections share a cell -> count 2
        assert grid["values"].max() == 2.0
        assert grid["lats"].shape == (4,)
        assert grid["lons"].shape == (4,)

    def test_density_grid_frp_weighted(self):
        f = self._field()
        grid = f.to_density_grid(resolution=0.5, frp_weighted=True)
        assert grid["weighting"] == "frp_MW"
        assert grid["values"].sum() == pytest.approx(
            float(np.nansum(f.frp)), rel=1e-6)

    def test_density_grid_per_day(self):
        f = self._field()
        grid = f.to_density_grid(resolution=1.0)
        assert len(grid["times"]) == 3  # Aug 1..3
        assert grid["values"].shape[0] == 3

    def test_json_roundtrip(self, tmp_path):
        f = self._field()
        path = str(tmp_path / "fires.json")
        f.to_json(path)
        g = FireField.from_json(path)
        assert g.times == f.times
        assert np.array_equal(g.lats, f.lats)
        assert np.array_equal(g.frp, f.frp)
        assert g.confidence == f.confidence
        assert g.bbox == f.bbox
        assert g.provenance == f.provenance

    def test_dict_roundtrip(self):
        f = self._field()
        g = FireField.from_dict(json.loads(json.dumps(f.to_dict())))
        assert len(g) == len(f)
        assert g.instruments == f.instruments

"""Tests for currents.gfs_wind — fully offline.

Live-network tests are never in this suite: the recorded GRIB2 fixtures
(``tests/fixtures/gfs_wind_20261001_5x5.grib2`` — a real 2026-10-01 00z
f000 NOMADS subregion response, 2t/10u/10v over -100..-99, 40..41 — and
``tests/fixtures/gfs_wind_20261001_f001_5x5.grib2``, the same window's
f001 hourly step, recorded live 2026-10-01)
cover parsing, and the HTTP layer is mocked for fetch tests. Tests that
need cfgrib skip cleanly when it is not installed.
"""

import datetime as dt
import hashlib
import json
import os
import urllib.error

import numpy as np
import pytest

from currents import gfs_wind
from currents.gfs_wind import (
    GFS_CYCLES,
    GFS_RETENTION_DAYS,
    NOMADS_FILTER_BASE,
    GfsWindField,
    UnavailableRangeError,
    fetch_gfs_wind,
    gfs_filter_url,
    gfs_lon_windows,
    gfs_retention_window,
    gfs_sample_plan,
    validate_gfs_bbox,
)

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures",
                       "gfs_wind_20261001_5x5.grib2")
F001_FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures",
                            "gfs_wind_20261001_f001_5x5.grib2")
FIXTURE_DATE = dt.date(2026, 10, 1)
FIXTURE_BBOX = (-100.0, 40.0, -99.0, 41.0)


def _fixture_bytes() -> bytes:
    with open(FIXTURE, "rb") as fh:
        return fh.read()


def _f001_fixture_bytes() -> bytes:
    with open(F001_FIXTURE, "rb") as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# URL construction (no network)
# ---------------------------------------------------------------------------


def test_filter_url_exact():
    url = gfs_filter_url(FIXTURE_DATE, "00", FIXTURE_BBOX)
    assert url == (
        NOMADS_FILTER_BASE
        + "?dir=/gfs.20261001/00/atmos"
        + "&file=gfs.t00z.pgrb2.0p25.f000"
        + "&subregion=on"
        + "&leftlon=-100.0&rightlon=-99.0&toplat=41.0&bottomlat=40.0"
        + "&lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on"
        + "&lev_2_m_above_ground=on&var_TMP=on"
    )


def test_filter_url_exact_f001():
    # The recorded f001 fixture was fetched with exactly this URL.
    url = gfs_filter_url(FIXTURE_DATE, "00", FIXTURE_BBOX, forecast_hour=1)
    assert url == (
        NOMADS_FILTER_BASE
        + "?dir=/gfs.20261001/00/atmos"
        + "&file=gfs.t00z.pgrb2.0p25.f001"
        + "&subregion=on"
        + "&leftlon=-100.0&rightlon=-99.0&toplat=41.0&bottomlat=40.0"
        + "&lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on"
        + "&lev_2_m_above_ground=on&var_TMP=on"
    )


def test_filter_url_hourly_filename_segments():
    # Zero-padded 3-digit forecast steps (verified live 2026-10-01).
    url1 = gfs_filter_url(FIXTURE_DATE, "00", FIXTURE_BBOX, forecast_hour=1)
    assert "file=gfs.t00z.pgrb2.0p25.f001" in url1
    assert "dir=/gfs.20261001/00/atmos" in url1
    url23 = gfs_filter_url(FIXTURE_DATE, "12", FIXTURE_BBOX, forecast_hour=23)
    assert "file=gfs.t12z.pgrb2.0p25.f023" in url23
    url120 = gfs_filter_url(FIXTURE_DATE, "18", FIXTURE_BBOX, forecast_hour=120)
    assert "file=gfs.t18z.pgrb2.0p25.f120" in url120


def test_filter_url_rejects_out_of_range_hour():
    with pytest.raises(ValueError):
        gfs_filter_url(FIXTURE_DATE, "00", FIXTURE_BBOX, forecast_hour=121)
    with pytest.raises(ValueError):
        gfs_filter_url(FIXTURE_DATE, "00", FIXTURE_BBOX, forecast_hour=-1)


def test_filter_url_uses_subregion_and_f000():
    # subregion=on is what makes the filter honor the box (verified
    # 2026-10-01); f000 is the analysis, not a forecast hour.
    url = gfs_filter_url(dt.date(2026, 9, 25), "18", (-130.0, 25.0, -65.0, 50.0))
    assert "subregion=on" in url
    assert "gfs.t18z.pgrb2.0p25.f000" in url
    assert "dir=/gfs.20260925/18/atmos" in url
    assert "&leftlon=-130.0&rightlon=-65.0&toplat=50.0&bottomlat=25.0" in url


def test_filter_url_rejects_bad_cycle():
    with pytest.raises(ValueError):
        gfs_filter_url(FIXTURE_DATE, "03", FIXTURE_BBOX)


def test_gfs_cycles_are_the_four_analysis_cycles():
    assert GFS_CYCLES == ("00", "06", "12", "18")


# ---------------------------------------------------------------------------
# bbox / windows
# ---------------------------------------------------------------------------


def test_validate_gfs_bbox_roundtrip():
    assert validate_gfs_bbox(FIXTURE_BBOX) == FIXTURE_BBOX


def test_validate_gfs_bbox_rejects_out_of_range():
    with pytest.raises(ValueError):
        validate_gfs_bbox((-200.0, 40.0, -99.0, 41.0))


def test_lon_windows_single_for_normal_bbox():
    assert gfs_lon_windows(FIXTURE_BBOX) == [FIXTURE_BBOX]


def test_lon_windows_split_antimeridian():
    windows = gfs_lon_windows((170.0, 30.0, -170.0, 50.0))
    assert windows == [(170.0, 30.0, 180.0, 50.0),
                       (-180.0, 30.0, -170.0, 50.0)]
    # Each split window gets its own filter URL (no leftlon > rightlon).
    for w in windows:
        url = gfs_filter_url(FIXTURE_DATE, "00", w)
        assert "leftlon=170.0&rightlon=180.0" in url or \
            "leftlon=-180.0&rightlon=-170.0" in url


# ---------------------------------------------------------------------------
# retention + sampling (no network)
# ---------------------------------------------------------------------------


def test_retention_window_is_ten_days_ending_today():
    oldest, newest = gfs_retention_window(dt.date(2026, 10, 1))
    assert oldest == dt.date(2026, 9, 22)
    assert newest == dt.date(2026, 10, 1)
    assert (newest - oldest).days == GFS_RETENTION_DAYS - 1


def test_dates_outside_retention_raise_honestly():
    # One day older than the verified 2026-09-22 floor.
    with pytest.raises(UnavailableRangeError) as ei:
        fetch_gfs_wind(FIXTURE_BBOX, "2026-09-21", "2026-09-21")
    assert "retention" in str(ei.value).lower()
    with pytest.raises(UnavailableRangeError):
        fetch_gfs_wind(FIXTURE_BBOX, "2026-10-02", "2026-10-02")  # future


def test_start_after_end_rejected():
    with pytest.raises(ValueError):
        fetch_gfs_wind(FIXTURE_BBOX, "2026-09-30", "2026-09-29")


def test_stride_days_must_be_positive_int():
    with pytest.raises(ValueError):
        fetch_gfs_wind(FIXTURE_BBOX, "2026-09-30", "2026-10-01", stride_days=0)
    with pytest.raises(ValueError):
        fetch_gfs_wind(FIXTURE_BBOX, "2026-09-30", "2026-10-01",
                       stride_days=1.5)


def test_sample_plan_daily_stride():
    d0 = dt.datetime(2026, 9, 28, tzinfo=dt.timezone.utc)
    d1 = dt.datetime(2026, 10, 1, tzinfo=dt.timezone.utc)
    plan = gfs_sample_plan(d0, d1, 1)
    assert plan == [dt.date(2026, 9, 28), dt.date(2026, 9, 29),
                    dt.date(2026, 9, 30), dt.date(2026, 10, 1)]
    plan2 = gfs_sample_plan(d0, d1, 2)
    assert plan2 == [dt.date(2026, 9, 28), dt.date(2026, 9, 30)]


# ---------------------------------------------------------------------------
# GRIB parsing against the recorded fixture (needs cfgrib)
# ---------------------------------------------------------------------------


def test_parse_fixture_grids():
    pytest.importorskip("cfgrib")
    parsed = gfs_wind._parse_grib_payload(
        _fixture_bytes(), "fixture-url", FIXTURE_DATE, "00")
    assert sorted(parsed) == ["t2m", "u10", "v10"]
    for key, part in parsed.items():
        assert part["values"].shape == (5, 5)
        assert part["data_date"] == "20261001"
        assert part["data_time"] == "00"
    # Lats increasing, lons normalized to -180..180.
    assert parsed["u10"]["lats"][0] == pytest.approx(40.0)
    assert parsed["u10"]["lats"][-1] == pytest.approx(41.0)
    assert parsed["u10"]["lons"][0] == pytest.approx(-100.0)
    assert parsed["u10"]["lons"][-1] == pytest.approx(-99.0)


def test_parse_fixture_values_are_real_gfs():
    # The fixture is a real NOMADS response: 2t ~291 K (October Iowa),
    # wind components within physical range.
    pytest.importorskip("cfgrib")
    parsed = gfs_wind._parse_grib_payload(
        _fixture_bytes(), "fixture-url", FIXTURE_DATE, "00")
    t2m = parsed["t2m"]["values"]
    assert 280.0 < t2m.mean() < 300.0  # Kelvin on the wire
    assert np.abs(parsed["u10"]["values"]).max() < 60.0
    assert np.abs(parsed["v10"]["values"]).max() < 60.0


def test_parse_rejects_wrong_day():
    pytest.importorskip("cfgrib")
    with pytest.raises(RuntimeError) as ei:
        gfs_wind._parse_grib_payload(
            _fixture_bytes(), "fixture-url", dt.date(2026, 9, 30), "00")
    assert "wrong day" in str(ei.value)


def test_verify_coverage_rejects_wrong_region():
    pytest.importorskip("cfgrib")
    parsed = gfs_wind._parse_grib_payload(
        _fixture_bytes(), "fixture-url", FIXTURE_DATE, "00")
    with pytest.raises(RuntimeError) as ei:
        gfs_wind._verify_coverage(parsed, (10.0, 40.0, 11.0, 41.0),
                                  "fixture-url")
    assert "does not cover" in str(ei.value)


def test_crop_to_window_snaps_to_grid():
    pytest.importorskip("cfgrib")
    parsed = gfs_wind._parse_grib_payload(
        _fixture_bytes(), "fixture-url", FIXTURE_DATE, "00")
    lats, lons, out = gfs_wind._crop_to_window(parsed, FIXTURE_BBOX)
    assert lats.shape == (5,) and lons.shape == (5,)
    assert out["t2m"].shape == (5, 5)


# ---------------------------------------------------------------------------
# fetch with a mocked HTTP layer (needs cfgrib for the real parse)
# ---------------------------------------------------------------------------


def _mock_cached_get(monkeypatch, tmp_path, payload=None, calls=None):
    """Replace the HTTP + parse layers: serve the fixtures, record the URLs.

    ``_cached_get_bytes`` serves the recorded fixture bytes for every
    URL (f000 bytes, or the recorded f001 bytes when the URL names
    ``f001``), so URL construction, provenance hashes, and the cache
    are exercised for real. ``_parse_grib_payload`` is stubbed to build
    grids consistent with the *requested* window/date/cycle/hour out of
    the fixture's values: for the exact recorded (window, day, cycle,
    hour) triples the real parser runs (proving the fixtures parse
    end-to-end); for other windows/hours the stub mirrors the real
    parser's contract (lats increasing, lons -180..180, ``valid_time``
    = cycle + forecast hour) without needing a distinct recorded
    payload per combination.
    """
    payload = _fixture_bytes() if payload is None else payload
    f001_payload = _f001_fixture_bytes()
    seen = [] if calls is None else calls
    real_parse = gfs_wind._parse_grib_payload

    def fake_get(url, cache_dir):
        seen.append(url)
        if "f001" in url:
            return f001_payload, False
        return payload, False

    def fake_parse(payload_bytes, url, day, cycle, forecast_hour=0):
        if (f"dir=/gfs.{FIXTURE_DATE.strftime('%Y%m%d')}/{cycle}/atmos" in url
                and "leftlon=-100.0&rightlon=-99.0" in url
                and day == FIXTURE_DATE
                and forecast_hour in (0, 1)):
            # Recorded triples: f000 (2026-10-01 00z) and f001
            # (2026-10-01 00z) fixtures, parsed for real.
            want = f001_payload if forecast_hour == 1 else payload
            return real_parse(want, url, day, cycle,
                              forecast_hour=forecast_hour)
        # Stub: derive the window from the URL and reuse the fixture's
        # values (the real parser's wrong-day/wrong-step refusals are
        # covered against real bytes below; the stub only exercises
        # fetch plumbing).
        import urllib.parse as _up
        q = _up.parse_qs(_up.urlparse(url).query)
        minx, maxx = float(q["leftlon"][0]), float(q["rightlon"][0])
        miny, maxy = float(q["bottomlat"][0]), float(q["toplat"][0])
        ref = real_parse(payload, url, FIXTURE_DATE, "00", forecast_hour=0)
        lats = np.linspace(miny, maxy, 5)
        lons = np.linspace(minx, maxx, 5)
        cycle_dt = dt.datetime(day.year, day.month, day.day, int(cycle),
                               tzinfo=dt.timezone.utc)
        valid = (cycle_dt + dt.timedelta(hours=forecast_hour)).isoformat()
        return {k: {"values": np.asarray(ref[k]["values"]),
                    "lats": lats, "lons": lons,
                    "data_date": day.strftime("%Y%m%d"),
                    "data_time": cycle,
                    "valid_time": valid}
                for k in ("u10", "v10", "t2m")}

    monkeypatch.setattr(gfs_wind, "_cached_get_bytes", fake_get)
    monkeypatch.setattr(gfs_wind, "_parse_grib_payload", fake_parse)
    monkeypatch.setattr(gfs_wind, "gfs_cache_dir",
                        lambda work_dir=None: str(tmp_path))
    return seen


def test_fetch_field_shape_units_and_times(monkeypatch, tmp_path):
    pytest.importorskip("cfgrib")
    _mock_cached_get(monkeypatch, tmp_path)
    f = fetch_gfs_wind(FIXTURE_BBOX, "2026-09-30", "2026-10-01",
                       stride_days=1, cycle="00")
    assert isinstance(f, GfsWindField)
    assert f.grids["u10"].shape == (2, 5, 5)
    assert f.grids["v10"].shape == (2, 5, 5)
    assert f.grids["t2m"].shape == (2, 5, 5)
    assert f.units == {"u10": "m/s", "v10": "m/s", "t2m": "\u00b0C"}
    assert f.times == ["2026-09-30T00:00:00+00:00",
                       "2026-10-01T00:00:00+00:00"]
    assert f.cycle == "00"
    # 2t converted Kelvin -> Celsius on ingest.
    assert 0.0 < f.grids["t2m"].mean() < 30.0


def test_fetch_default_forecast_hours_is_f000_only():
    # The default (0,) path is backwards compatible: one f000 analysis
    # per sampled day, timestamps at the cycle hour.
    import inspect
    sig = inspect.signature(fetch_gfs_wind)
    assert sig.parameters["forecast_hours"].default == (0,)


def test_fetch_multi_hour_chronological_times(monkeypatch, tmp_path):
    pytest.importorskip("cfgrib")
    seen = _mock_cached_get(monkeypatch, tmp_path)
    f = fetch_gfs_wind(FIXTURE_BBOX, "2026-10-01", "2026-10-01",
                       cycle="00", forecast_hours=(0, 1, 6))
    assert f.grids["u10"].shape == (3, 5, 5)
    # Timestamps carry the forecast hour: valid time, not the cycle.
    assert f.times == ["2026-10-01T00:00:00+00:00",
                       "2026-10-01T01:00:00+00:00",
                       "2026-10-01T06:00:00+00:00"]
    # One NOMADS request per hour, exact filename segments.
    assert len(seen) == 3
    assert "gfs.t00z.pgrb2.0p25.f000" in seen[0]
    assert "gfs.t00z.pgrb2.0p25.f001" in seen[1]
    assert "gfs.t00z.pgrb2.0p25.f006" in seen[2]
    # Provenance is per timestep: URL + SHA-256 + byte count each.
    reqs = f.provenance["requests"]
    assert [r["forecast_hour"] for r in reqs] == [0, 1, 6]
    assert [r["url"] for r in reqs] == seen
    assert all(r["sha256"] and r["n_bytes"] > 0 for r in reqs)
    assert f.provenance["forecast_hours"] == [0, 1, 6]
    assert f.provenance["n_requests"] == 3
    assert f.provenance["n_bytes"] == sum(r["n_bytes"] for r in reqs)


def test_fetch_forecast_hours_normalized_to_sorted_order(monkeypatch,
                                                        tmp_path):
    # Scrambled and duplicated input assembles chronologically.
    pytest.importorskip("cfgrib")
    _mock_cached_get(monkeypatch, tmp_path)
    f = fetch_gfs_wind(FIXTURE_BBOX, "2026-10-01", "2026-10-01",
                       forecast_hours=(6, 0, 1, 0))
    assert f.times == ["2026-10-01T00:00:00+00:00",
                       "2026-10-01T01:00:00+00:00",
                       "2026-10-01T06:00:00+00:00"]
    assert f.provenance["forecast_hours"] == [0, 1, 6]


def test_fetch_multi_day_multi_hour_is_chronological(monkeypatch, tmp_path):
    pytest.importorskip("cfgrib")
    _mock_cached_get(monkeypatch, tmp_path)
    f = fetch_gfs_wind(FIXTURE_BBOX, "2026-09-30", "2026-10-01",
                       forecast_hours=(0, 6))
    assert f.grids["u10"].shape == (4, 5, 5)
    assert f.times == ["2026-09-30T00:00:00+00:00",
                       "2026-09-30T06:00:00+00:00",
                       "2026-10-01T00:00:00+00:00",
                       "2026-10-01T06:00:00+00:00"]


def test_forecast_hours_validation_rejects_bad_hours():
    bad_hours = [(121,), (-1,), (1.5,), (True,), ("1",), ()]
    for hours in bad_hours:
        with pytest.raises(ValueError):
            fetch_gfs_wind(FIXTURE_BBOX, "2026-09-30", "2026-10-01",
                           forecast_hours=hours)
    with pytest.raises(ValueError):
        fetch_gfs_wind(FIXTURE_BBOX, "2026-09-30", "2026-10-01",
                       forecast_hours=True)


def test_parse_f001_fixture_carries_valid_time():
    # The recorded f001 GRIB: dataDate/dataTime stay at the cycle, the
    # valid time (cycle + 1h) is what lands on the times axis.
    pytest.importorskip("cfgrib")
    parsed = gfs_wind._parse_grib_payload(
        _f001_fixture_bytes(), "fixture-url", FIXTURE_DATE, "00",
        forecast_hour=1)
    for key, part in parsed.items():
        assert part["valid_time"] == "2026-10-01T01:00:00+00:00"
        assert part["data_date"] == "20261001"
        assert part["data_time"] == "00"
    # The same bytes refused when parsed as the wrong step.
    with pytest.raises(RuntimeError) as ei:
        gfs_wind._parse_grib_payload(
            _f001_fixture_bytes(), "fixture-url", FIXTURE_DATE, "00",
            forecast_hour=0)
    assert "stepRange" in str(ei.value)


def test_fetch_kelvin_to_celsius_and_fahrenheit_conventions():
    pytest.importorskip("cfgrib")
    parsed = gfs_wind._parse_grib_payload(
        _fixture_bytes(), "fixture-url", FIXTURE_DATE, "00")
    kelvin = float(parsed["t2m"]["values"][2, 2])
    f = GfsWindField.synthetic(nt=1)
    # Synthetic t2m is °C already; check the conversion identity on a
    # hand-built field instead.
    t2m_c = np.ma.array([[[kelvin - 273.15]]])
    u = np.ma.array([[[3.0]]])
    v = np.ma.array([[[4.0]]])
    hand = GfsWindField(
        grids={"u10": u, "v10": v, "t2m": t2m_c},
        times=["2026-10-01T00:00:00+00:00"],
        lats=np.array([40.0]), lons=np.array([-100.0]))
    assert float(hand.grids["t2m"][0, 0, 0]) == pytest.approx(kelvin - 273.15)
    # air_temperature is °F (the dark_strands strand-color convention).
    assert hand.temperature_unit == "\u00b0F"
    assert float(hand.air_temperature[0, 0, 0]) == pytest.approx(
        (kelvin - 273.15) * 9.0 / 5.0 + 32.0)
    # wind_speed is the vector magnitude.
    assert float(hand.wind_speed[0, 0, 0]) == pytest.approx(5.0)


def test_fetch_provenance_records_urls_hashes_bytes(monkeypatch, tmp_path):
    pytest.importorskip("cfgrib")
    seen = _mock_cached_get(monkeypatch, tmp_path)
    f = fetch_gfs_wind(FIXTURE_BBOX, "2026-10-01", "2026-10-01")
    prov = f.provenance
    assert prov["n_requests"] == 1
    assert len(seen) == 1
    req = prov["requests"][0]
    assert req["url"] == seen[0]
    assert req["url"].startswith(NOMADS_FILTER_BASE)
    assert req["sha256"] == hashlib.sha256(_fixture_bytes()).hexdigest()
    assert req["n_bytes"] == len(_fixture_bytes())
    assert prov["sha256"] == hashlib.sha256(
        hashlib.sha256(_fixture_bytes()).digest()).hexdigest()
    assert prov["n_bytes"] == len(_fixture_bytes())
    assert "retrieved_at" in prov
    assert prov["access"] == "keyless"
    assert prov["retention_days"] == GFS_RETENTION_DAYS


def test_fetch_antimeridian_uses_two_windows(monkeypatch, tmp_path):
    pytest.importorskip("cfgrib")
    seen = _mock_cached_get(monkeypatch, tmp_path)
    f = fetch_gfs_wind((170.0, 40.0, -170.0, 41.0), "2026-10-01", "2026-10-01")
    assert len(seen) == 2  # one request per split window
    assert "leftlon=170.0" in seen[0]
    assert "leftlon=-180.0" in seen[1]
    # Grids concatenated along longitude (5 + 5 columns here since the
    # fixture mock returns the same 5x5 payload for both windows).
    assert f.grids["u10"].shape == (1, 5, 10)


def test_fetch_http_404_becomes_unavailable_range_error(monkeypatch, tmp_path):
    # Bypass the cache helper and fail at urlopen: the adapter must turn
    # a NOMADS 404 into UnavailableRangeError, never a silent skip.
    import urllib.request

    def boom(req, timeout=None):
        raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(urllib.request, "urlopen", boom)
    monkeypatch.setattr(gfs_wind, "gfs_cache_dir",
                        lambda work_dir=None: str(tmp_path))
    with pytest.raises(UnavailableRangeError) as ei:
        fetch_gfs_wind(FIXTURE_BBOX, "2026-10-01", "2026-10-01")
    assert "404" in str(ei.value)


def test_fetch_hour_404_becomes_unavailable_range_error(monkeypatch,
                                                       tmp_path):
    # A 404 on ONE (day, hour) must fail loudly with the exact URL —
    # never skip the hour and never pad it silently.
    pytest.importorskip("cfgrib")
    import urllib.request

    class Resp:
        def __init__(self, payload):
            self.payload = payload

        def read(self):
            return self.payload

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def selective(req, timeout=None):
        if "f006" in req.full_url:
            raise urllib.error.HTTPError(req.full_url, 404, "Not Found",
                                         {}, None)
        return Resp(_fixture_bytes())

    monkeypatch.setattr(urllib.request, "urlopen", selective)
    monkeypatch.setattr(gfs_wind, "gfs_cache_dir",
                        lambda work_dir=None: str(tmp_path))
    with pytest.raises(UnavailableRangeError) as ei:
        fetch_gfs_wind(FIXTURE_BBOX, "2026-10-01", "2026-10-01",
                       forecast_hours=(0, 6))
    assert "404" in str(ei.value)
    assert "f006" in str(ei.value)


def test_fetch_uses_cache_on_second_call(monkeypatch, tmp_path):
    pytest.importorskip("cfgrib")
    calls = []
    real = gfs_wind._cached_get_bytes

    def counting(url, cache_dir):
        calls.append(url)
        return real(url, cache_dir)

    monkeypatch.setattr(gfs_wind, "_cached_get_bytes", counting)
    # Pre-seed the cache by writing through the real helper once, with
    # urlopen mocked to serve the fixture.
    import urllib.request

    class Resp:
        def __init__(self, payload):
            self.payload = payload

        def read(self):
            return self.payload

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        return Resp(_fixture_bytes())

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    cache = str(tmp_path)
    monkeypatch.setattr(gfs_wind, "gfs_cache_dir",
                        lambda work_dir=None: cache)
    fetch_gfs_wind(FIXTURE_BBOX, "2026-10-01", "2026-10-01")
    assert len(calls) == 1
    # Second fetch: served from cache, urlopen not hit again. Replace
    # urlopen with a bomb to prove it.
    def bomb(req, timeout=None):
        raise AssertionError("network hit on cached fetch")

    monkeypatch.setattr(urllib.request, "urlopen", bomb)
    f = fetch_gfs_wind(FIXTURE_BBOX, "2026-10-01", "2026-10-01")
    assert f.provenance["n_cached"] == 1
    assert f.grids["u10"].shape == (1, 5, 5)


def test_require_cfgrib_error_is_actionable(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "cfgrib.messages" or name.startswith("cfgrib"):
            raise ImportError("no cfgrib")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError) as ei:
        gfs_wind._require_cfgrib()
    assert "survey-currents[gfs]" in str(ei.value)


# ---------------------------------------------------------------------------
# GfsWindField model
# ---------------------------------------------------------------------------


def test_field_rejects_missing_grids():
    with pytest.raises(ValueError):
        GfsWindField(
            grids={"u10": np.ma.zeros((1, 2, 2))},
            times=["2026-10-01T00:00:00+00:00"],
            lats=np.array([40.0, 41.0]), lons=np.array([-100.0, -99.0]))


def test_field_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        GfsWindField(
            grids={k: np.ma.zeros((2, 2, 2)) for k in ("u10", "v10", "t2m")},
            times=["2026-10-01T00:00:00+00:00"],
            lats=np.array([40.0, 41.0]), lons=np.array([-100.0, -99.0]))


def test_overlay_grids_and_spatial_mean():
    f = GfsWindField.synthetic(nt=2, seed=3)
    assert sorted(f.overlay_grids) == ["t2m"]
    assert f.overlay_grid("t2m", 0).shape == f.grids["t2m"][0].shape
    with pytest.raises(KeyError):
        f.overlay_grid("msl", 0)
    assert f.spatial_mean("wind", 0) == pytest.approx(
        float(np.ma.masked_invalid(f.wind_speed[0]).mean()))
    assert f.spatial_mean("t2m", 1) < f.spatial_mean("air_temperature", 1)
    with pytest.raises(KeyError):
        f.spatial_mean("bogus", 0)


def test_select_time_and_bbox():
    f = GfsWindField.synthetic(nt=3, ny=6, nx=8,
                               lats=(-30.0, 30.0), lons=(-60.0, 60.0))
    one = f.select_time(1)
    assert one.grids["u10"].shape == (1, 6, 8)
    assert one.times == [f.times[1]]
    sub = f.select_bbox((-30.0, -10.0, 30.0, 10.0))
    assert sub.lats[0] >= -10.0 and sub.lats[-1] <= 10.0
    assert sub.lons[0] >= -30.0 and sub.lons[-1] <= 30.0
    with pytest.raises(ValueError):
        f.select_bbox((100.0, 40.0, 110.0, 50.0))


def test_dict_json_roundtrip():
    f = GfsWindField.synthetic(nt=2, seed=5)
    d = f.to_dict()
    assert "air_temperature" in d and "temperature_unit" in d
    assert "grids" in d and "values" in d
    g = GfsWindField.from_dict(json.loads(json.dumps(d)))
    assert g.times == f.times
    assert g.cycle == f.cycle
    np.testing.assert_allclose(np.ma.filled(g.grids["u10"], np.nan),
                               np.ma.filled(f.grids["u10"], np.nan))
    np.testing.assert_allclose(np.ma.filled(g.air_temperature, np.nan),
                               np.ma.filled(f.air_temperature, np.nan))
    with pytest.raises(ValueError):
        GfsWindField.from_dict({"times": []})


def test_json_file_roundtrip(tmp_path):
    f = GfsWindField.synthetic(nt=2, seed=6)
    path = str(tmp_path / "gfs.json")
    f.to_json(path)
    g = GfsWindField.from_json(path)
    assert g.times == f.times
    assert g.grids["v10"].shape == f.grids["v10"].shape


def test_synthetic_is_deterministic():
    a = GfsWindField.synthetic(seed=9)
    b = GfsWindField.synthetic(seed=9)
    np.testing.assert_array_equal(np.ma.getdata(a.grids["u10"]),
                                  np.ma.getdata(b.grids["u10"]))


def test_main_demo_runs(capsys):
    gfs_wind.main_demo()
    out = capsys.readouterr().out
    assert "[gfs-wind]" in out and "m/s" in out

"""Tests for the IBTrACS tropical-cyclone adapter (fully offline).

Parsing is exercised through :func:`storms._parse_ibtracs_dataset` with
duck-typed fake datasets — no netCDF4 needed. Downloads are mocked.
Live-network tests are skipped unless ``SURVEY_CURRENTS_LIVE=1``.
"""

import datetime as dt
import json
import os

import numpy as np
import pytest

from currents import storms as st
from currents.storms import (StormField, StormTrack, fetch_ibtracs,
                             ibtracs_filename, ibtracs_url, sshs_category,
                             storm_cache_dir)


# ---------------------------------------------------------------------------
# Fake NetCDF dataset
# ---------------------------------------------------------------------------

class _FakeVar:
    def __init__(self, data, fill=None):
        self._data = np.asarray(data)
        if fill is not None:
            self._FillValue = fill

    def __getitem__(self, key):
        return self._data[key]


class _FakeDS:
    def __init__(self, variables):
        self.variables = variables


def _chars(strings, width):
    arr = np.zeros((len(strings), width), dtype="S1")
    for i, s in enumerate(strings):
        b = s.encode("ascii")[:width]
        arr[i, :len(b)] = np.frombuffer(b, dtype="S1")
    return arr


def _chars3(strings, n2, width):
    """(n, n2, width) char array from a list of per-storm value lists."""
    arr = np.zeros((len(strings), n2, width), dtype="S1")
    for i, vals in enumerate(strings):
        for j, s in enumerate(vals):
            b = s.encode("ascii")[:width]
            arr[i, j, :len(b)] = np.frombuffer(b, dtype="S1")
    return arr


def _iso(times):
    return _chars3([[t.strftime("%Y-%m-%d %H:%M:%S")] for t in times], 1, 19)[:, 0, :]


def _make_ds():
    """3 fake storms x 8 obs slots."""
    n, m = 3, 8
    numobs = np.array([8, 5, 5], dtype=np.int16)
    sid = _chars(["2005236N12345", "2005240N23456", "2006001S34567"], 13)
    name = _chars(["TESTKATRINA", "UNNAMED", "TESTDATELINE"], 128)
    season = np.array([2005, 2005, 2006], dtype=np.int16)
    basin = _chars3([["NA"] * m, ["NA"] * m, ["SP"] * m], m, 2)

    t0 = dt.datetime(2005, 8, 24, tzinfo=dt.timezone.utc)
    iso0 = [t0 + dt.timedelta(hours=6 * k) for k in range(8)]
    t1 = dt.datetime(2005, 9, 1, tzinfo=dt.timezone.utc)
    iso1 = [t1 + dt.timedelta(hours=6 * k) for k in range(5)]
    t2 = dt.datetime(2006, 1, 10, tzinfo=dt.timezone.utc)
    iso2 = [t2 + dt.timedelta(hours=6 * k) for k in range(5)]
    iso = np.zeros((n, m, 19), dtype="S1")
    for i, times in enumerate((iso0, iso1, iso2)):
        for j, t in enumerate(times):
            b = t.strftime("%Y-%m-%d %H:%M:%S").encode("ascii")
            iso[i, j, :len(b)] = np.frombuffer(b, dtype="S1")

    lat = np.full((n, m), np.nan)
    lon = np.full((n, m), np.nan)
    lat[0] = [23.0, 24.1, 25.2, 26.0, 27.1, 28.2, 29.0, 30.1]
    lon[0] = [-75.0, -76.5, -78.0, -79.5, -81.0, -82.5, -84.0, -85.5]
    lat[1, :5] = [15.0, 16.0, 17.0, 18.0, 19.0]
    lon[1, :5] = [-65.0, -66.0, -67.0, -68.0, -69.0]
    # Dateline crosser in 0..360 longitudes (as the archive mixes).
    lat[2, :5] = [-15.0, -16.0, -17.0, -18.0, -19.0]
    lon[2, :5] = [170.0, 175.0, 185.0, 190.0, 195.0]

    usa_wind = np.full((n, m), -9999, dtype=np.int16)
    wmo_wind = np.full((n, m), -9999, dtype=np.int16)
    usa_wind[0] = [35, 45, 65, 90, 110, 100, 80, 50]
    usa_wind[1, :5] = [30, 35, 40, 35, 30]
    # Storm 2: usa missing on obs 1 -> wmo fallback.
    usa_wind[2, :5] = [40, -9999, 70, 85, 60]
    wmo_wind[2, :5] = [38, 55, 68, 80, 58]
    usa_pres = np.full((n, m), -9999, dtype=np.int16)
    wmo_pres = np.full((n, m), -9999, dtype=np.int16)
    usa_pres[0] = [1002, 998, 990, 980, 970, 975, 985, 995]
    usa_pres[1, :5] = [1008, 1006, 1004, 1006, 1008]
    usa_pres[2, :5] = [1000, -9999, 990, 985, 992]
    wmo_pres[2, :5] = [999, 996, 989, 984, 991]

    variables = {
        "numobs": _FakeVar(numobs),
        "sid": _FakeVar(sid),
        "name": _FakeVar(name),
        "season": _FakeVar(season),
        "basin": _FakeVar(basin),
        "iso_time": _FakeVar(iso),
        "lat": _FakeVar(lat),
        "lon": _FakeVar(lon),
        "usa_wind": _FakeVar(usa_wind, fill=-9999),
        "wmo_wind": _FakeVar(wmo_wind, fill=-9999),
        "usa_pres": _FakeVar(usa_pres, fill=-9999),
        "wmo_pres": _FakeVar(wmo_pres, fill=-9999),
    }
    return _FakeDS(variables)


_NA_BBOX = (-90.0, 10.0, -60.0, 35.0)


def _parse(**kw):
    ds = _make_ds()
    args = dict(bbox=_NA_BBOX, d0=dt.date(2005, 8, 1),
                d1=dt.date(2005, 9, 30))
    args.update(kw)
    return st._parse_ibtracs_dataset(ds, **args)


# ---------------------------------------------------------------------------
# Category mapping
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("wind,cat", [
    (0, "TD"), (33.9, "TD"), (34, "TS"), (63.9, "TS"),
    (64, "C1"), (82.9, "C1"), (83, "C2"), (95.9, "C2"),
    (96, "C3"), (112.9, "C3"), (113, "C4"), (136.9, "C4"),
    (137, "C5"), (200, "C5"),
])
def test_sshs_category_boundaries(wind, cat):
    assert sshs_category(wind) == cat


def test_sshs_category_unknown():
    assert sshs_category(float("nan")) == "unknown"
    assert sshs_category(-5) == "unknown"
    for code in ("TD", "TS", "C1", "C2", "C3", "C4", "C5", "unknown"):
        assert code in st.SSHS_COLORS and code in st.SSHS_LABELS


# ---------------------------------------------------------------------------
# Constants / cache helpers
# ---------------------------------------------------------------------------

def test_urls_and_filenames():
    assert ibtracs_filename(False) == "IBTrACS.since1980.v04r01.nc"
    assert ibtracs_filename(True) == "IBTrACS.ALL.v04r01.nc"
    assert ibtracs_url(False).startswith("https://www.ncei.noaa.gov/")
    assert ibtracs_url(False).endswith(ibtracs_filename(False))
    assert st.IBTRACS_VERSION == "v04r01"


def test_storm_cache_dir_env(monkeypatch, tmp_path):
    monkeypatch.setenv("SURVEY_CURRENTS_CACHE", str(tmp_path))
    assert storm_cache_dir() == os.path.join(str(tmp_path), "storms")
    assert os.path.isdir(storm_cache_dir())


def test_norm_lon():
    assert st._norm_lon(185.0) == pytest.approx(-175.0)
    assert st._norm_lon(-185.0) == pytest.approx(175.0)
    assert st._norm_lon(170.0) == pytest.approx(170.0)
    assert st._norm_lon(360.0) == pytest.approx(0.0)


def test_chars_to_str():
    arr = _chars(["KATRINA", "UNNAMED"], 13)
    assert st._chars_to_str(arr[0]) == "KATRINA"
    assert st._chars_to_str(arr[1]) == "UNNAMED"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def test_parse_basic():
    tracks = _parse()
    # Storm 0 (NA, in bbox) and storm 1 (in bbox) kept; storm 2 is in the
    # South Pacific, outside the NA bbox.
    assert [t.name for t in tracks] == ["TESTKATRINA", "UNNAMED"]
    t = tracks[0]
    assert t.sid == "2005236N12345"
    assert t.season == 2005
    assert t.basin == "NA"
    assert len(t) == 8
    assert list(t.winds) == [35, 45, 65, 90, 110, 100, 80, 50]
    assert t.max_wind == 110.0
    assert t.max_category == "C3"
    assert t.min_pres == 970.0


def test_parse_wind_priority_fallback():
    tracks = _parse(bbox=(-180.0, -25.0, 180.0, -5.0),
                    d0=dt.date(2006, 1, 1), d1=dt.date(2006, 2, 1))
    assert len(tracks) == 1
    t = tracks[0]
    # Obs 1 had usa_wind fill -> wmo_wind 55 used.
    assert list(t.winds) == [40, 55, 70, 85, 60]
    assert list(t.press) == [1000, 996, 990, 985, 992]


def test_parse_dateline_normalization():
    tracks = _parse(bbox=(-180.0, -25.0, 180.0, -5.0),
                    d0=dt.date(2006, 1, 1), d1=dt.date(2006, 2, 1))
    t = tracks[0]
    assert list(t.lons) == pytest.approx([170.0, 175.0, -175.0, -170.0, -165.0])
    assert all(-180.0 <= lo < 180.0 for lo in t.lons)


def test_parse_time_clipping():
    tracks = _parse(d0=dt.date(2005, 8, 25), d1=dt.date(2005, 8, 26))
    t = tracks[0]
    # Storm kept (obs in bbox); fixes clipped to the window (4 of 8).
    assert len(t) == 4
    assert all(dt.date(2005, 8, 25) <= tm.date() <= dt.date(2005, 8, 26)
               for tm in t.times)


def test_parse_bbox_excludes():
    tracks = _parse(bbox=(0.0, 40.0, 10.0, 50.0))
    assert tracks == []


def test_parse_storm_name_match():
    tracks = _parse(storm_name="testkatrina")
    assert [t.name for t in tracks] == ["TESTKATRINA"]
    assert _parse(storm_name="no-such-storm") == []


def test_parse_min_wind():
    tracks = _parse(min_wind=100.0)
    assert [t.name for t in tracks] == ["TESTKATRINA"]
    assert _parse(min_wind=200.0) == []


def test_parse_deterministic_order():
    tracks = _parse()
    seasons = [(t.season, t.times[0]) for t in tracks]
    assert seasons == sorted(seasons)


# ---------------------------------------------------------------------------
# StormField model
# ---------------------------------------------------------------------------

def _field():
    return StormField.synthetic(seed=7)


def test_synthetic():
    f = _field()
    assert len(f) == 3
    assert "TESTALPHA" in [t.name for t in f.tracks]
    assert f.n_obs > 0
    assert f.source == "synthetic"
    lo, hi = f.time_range
    assert lo <= hi


def test_track_round_trip():
    f = _field()
    f2 = StormField.from_dict(f.to_dict())
    assert len(f2) == len(f)
    assert [t.name for t in f2.tracks] == [t.name for t in f.tracks]
    assert f2.tracks[0].times[0] == f.tracks[0].times[0]
    assert f2.tracks[0].max_wind == pytest.approx(f.tracks[0].max_wind)


def test_field_json_round_trip(tmp_path):
    f = _field()
    path = str(tmp_path / "storms.json")
    f.to_json(path)
    f2 = StormField.from_json(path)
    assert len(f2) == 3
    assert f2.provenance["seed"] == 7


def test_rank_by_intensity():
    f = _field()
    ranked = f.rank_by_intensity()
    winds = [t.max_wind for t in ranked]
    assert winds == sorted(winds, reverse=True)


def test_select_time():
    f = _field()
    d0 = f.start + dt.timedelta(days=20)
    sub = f.select_time(d0, f.end)
    assert all(all(d0 <= tm.date() <= f.end for tm in t.times)
               for t in sub.tracks)
    assert sub.start == d0


def test_select_bbox():
    f = _field()
    sub = f.select_bbox((-100.0, 10.0, -60.0, 40.0))
    assert len(sub) == len(f)
    empty = f.select_bbox((0.0, 40.0, 10.0, 50.0))
    assert len(empty) == 0


def test_field_validation():
    with pytest.raises(ValueError):
        StormField(tracks=[], bbox=(-100, 10, -60, 40),
                   start="2024-09-01", end="2024-08-01")


# ---------------------------------------------------------------------------
# fetch_ibtracs validation (offline — fails before any download)
# ---------------------------------------------------------------------------

def test_fetch_validates_dates():
    with pytest.raises(ValueError):
        fetch_ibtracs(_NA_BBOX, "2024-09-01", "2024-08-01")
    with pytest.raises(ValueError):
        fetch_ibtracs(_NA_BBOX, "1970-01-01", "1970-02-01")
    with pytest.raises(ValueError):
        fetch_ibtracs(_NA_BBOX, "2024-08-01", "2024-09-01", min_wind=-1)
    # Full archive accepts 1842+ but still validates order.
    with pytest.raises(ValueError):
        fetch_ibtracs(_NA_BBOX, "1900-01-01", "1899-01-01",
                      full_archive=True)


def test_fetch_requires_netcdf4(monkeypatch, tmp_path):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "netCDF4":
            raise ImportError("no netCDF4")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="netCDF4"):
        st._require_netcdf4()


# ---------------------------------------------------------------------------
# Cache behaviour (mocked download)
# ---------------------------------------------------------------------------

def test_ensure_ibtracs_file_caches(monkeypatch, tmp_path):
    calls = []

    def fake_download(url, timeout=600):
        calls.append(url)
        return b"fake-ibtracs-bytes"

    monkeypatch.setattr(st, "_download_bytes", fake_download)
    monkeypatch.setattr(st, "_http_head_last_modified",
                        lambda url, timeout=60: "Thu, 24 Sep 2026 09:28:46 GMT")
    path, prov = st.ensure_ibtracs_file(cache_dir=str(tmp_path))
    assert os.path.isfile(path)
    assert prov["downloaded"] is True
    assert len(prov["sha256"]) == 64
    assert os.path.isfile(path + ".sha256")
    # Second call: cache hit, no download.
    path2, prov2 = st.ensure_ibtracs_file(cache_dir=str(tmp_path))
    assert path2 == path
    assert prov2["cache_hit"] is True
    assert len(calls) == 1


def test_ensure_ibtracs_file_corrupt_redownload(monkeypatch, tmp_path):
    monkeypatch.setattr(st, "_download_bytes", lambda url, timeout=600: b"v2")
    monkeypatch.setattr(st, "_http_head_last_modified",
                        lambda url, timeout=60: None)
    path, _ = st.ensure_ibtracs_file(cache_dir=str(tmp_path))
    with open(path, "wb") as fh:
        fh.write(b"corrupted")
    path2, prov2 = st.ensure_ibtracs_file(cache_dir=str(tmp_path))
    assert path2 == path
    assert prov2["downloaded"] is True  # corrupt entry re-downloaded
    assert open(path, "rb").read() == b"v2"


def test_ensure_ibtracs_file_refresh(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(st, "_download_bytes",
                        lambda url, timeout=600: calls.append(url) or b"data")
    monkeypatch.setattr(st, "_http_head_last_modified",
                        lambda url, timeout=60: "lm")
    st.ensure_ibtracs_file(cache_dir=str(tmp_path))
    st.ensure_ibtracs_file(cache_dir=str(tmp_path), refresh=True)
    assert len(calls) == 2


# ---------------------------------------------------------------------------
# CLI (offline)
# ---------------------------------------------------------------------------

def test_cli_storms_synthetic(tmp_path, capsys):
    from currents import cli
    ns = type("NS", (), {"bbox": "-100,10,-60,40", "start": "2024-08-01",
                         "end": "2024-09-15", "n_storms": 2, "seed": 7,
                         "source": "synthetic",
                         "out": str(tmp_path / "syn")})()
    assert cli.cmd_storms_synthetic(ns) == 0
    path = str(tmp_path / "syn.json")
    assert os.path.exists(path)
    f = StormField.from_json(path)
    assert len(f) == 2


def test_cli_fetch_ibtracs_mocked(tmp_path, capsys, monkeypatch):
    from currents import cli
    import currents.storms as sm

    def fake_fetch(bbox, start, end, min_wind=None, storm_name=None,
                   full_archive=False):
        assert storm_name == "testkatrina"
        return StormField.synthetic()

    monkeypatch.setattr(sm, "fetch_ibtracs", fake_fetch)
    ns = type("NS", (), {"bbox": "-100,10,-60,40", "start": "2024-08-01",
                         "end": "2024-09-15", "min_wind": None,
                         "storm_name": "testkatrina", "full_archive": False,
                         "out": str(tmp_path / "fetch")})()
    assert cli.cmd_fetch_ibtracs(ns) == 0
    out = capsys.readouterr().out
    assert "IBTrACS" in out
    assert os.path.exists(str(tmp_path / "fetch.json"))


# ---------------------------------------------------------------------------
# Live (skipped unless SURVEY_CURRENTS_LIVE=1)
# ---------------------------------------------------------------------------

_needs_live = pytest.mark.skipif(
    os.environ.get("SURVEY_CURRENTS_LIVE") != "1",
    reason="live NCEI download (set SURVEY_CURRENTS_LIVE=1)")


@_needs_live
class TestLiveIbtracs:
    def test_fetch_katrina_2005(self):
        pytest.importorskip("netCDF4")
        field = fetch_ibtracs((-100.0, 20.0, -60.0, 40.0),
                              "2005-08-23", "2005-08-31",
                              storm_name="katrina")
        assert len(field) >= 1
        assert field.tracks[0].name.upper() == "KATRINA"
        assert field.tracks[0].max_wind >= 64
        assert field.provenance["sha256"]
        assert field.provenance["ibtracs_version"] == "v04r01"

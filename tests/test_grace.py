"""Tests for the CSR GRACE/GRACE-FO RL06.3 adapter (fully offline).

Parsing is exercised through :func:`grace._parse_grace_dataset` with
duck-typed fake datasets — no netCDF4 needed. Downloads are mocked.
Live-network tests are skipped unless ``SURVEY_CURRENTS_LIVE=1``.
"""

import datetime as dt
import json
import os

import numpy as np
import pytest

from currents import grace as gr
from currents.grace import (GRACE_BASELINE, GRACE_EPOCH, GRACE_MASK_FILE,
                            GRACE_MASK_URL, GRACE_RECORD_START,
                            GRACE_SOLUTION_FILE, GRACE_SOLUTION_URL,
                            GRACE_UNITS, WaterField, fetch_grace,
                            grace_cache_dir, grace_file_url)


# ---------------------------------------------------------------------------
# Fake NetCDF datasets
# ---------------------------------------------------------------------------

class _FakeVar:
    def __init__(self, data):
        # asanyarray preserves np.ma.masked_array semantics, like a real
        # netCDF4 variable slice.
        self._data = np.asanyarray(data)
        self.dimensions = tuple(f"dim{i}" for i in range(self._data.ndim))

    def __getitem__(self, key):
        return self._data[key]


class _FakeDS:
    """Context-manager fake with a plain-dict ``variables``."""

    def __init__(self, variables):
        self.variables = dict(variables)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _days(y, m, d=1):
    return (dt.date(y, m, d) - GRACE_EPOCH).days


def _fake_solution(months, lats, lons, value=1.0):
    """Fake CSR solution dataset for the given (year, month) labels.

    ``months`` is a list of ``(year, month)`` tuples in time order.
    The frame for month k holds the constant ``value + k``.
    """
    times, bounds = [], []
    for (y, m) in months:
        first = _days(y, m, 1)
        last = _days(y, m, 1) + 31  # bounds are only midpoint-labeled
        times.append(first + 15)
        bounds.append([first, last])
    frames = [np.full((len(lats), len(lons)), value + k, dtype=float)
              for k in range(len(months))]
    return _FakeDS({
        "time": _FakeVar(times),
        "time_bounds": _FakeVar(bounds),
        "lat": _FakeVar(lats),
        "lon": _FakeVar(lons),
        "lwe_thickness": _FakeVar(np.stack(frames)),
    })


def _fake_mask(lats, lons, land_where=None):
    """Fake land mask: 1 = land, 0 = ocean."""
    arr = np.zeros((len(lats), len(lons)))
    if land_where is None:
        arr[:] = 1.0
    else:
        arr[np.ix_(*land_where)] = 1.0
    return np.asarray(arr)


# ---------------------------------------------------------------------------
# Constants / cache plumbing
# ---------------------------------------------------------------------------

def test_solution_url_is_keyless_csr():
    assert GRACE_SOLUTION_URL == (
        "https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/"
        "CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc")
    assert GRACE_MASK_URL == (
        "https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/"
        "CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc")


def test_grace_file_url_which():
    assert grace_file_url("solution") == GRACE_SOLUTION_URL
    assert grace_file_url("mask") == GRACE_MASK_URL
    with pytest.raises(ValueError):
        grace_file_url("bogus")


def test_grace_cache_dir_env(monkeypatch, tmp_path):
    monkeypatch.setenv("SURVEY_CURRENTS_CACHE", str(tmp_path))
    assert grace_cache_dir() == os.path.join(str(tmp_path), "grace")
    assert os.path.isdir(grace_cache_dir())


def test_grace_units_and_baseline_constants():
    assert GRACE_UNITS == "cm"
    assert "2004" in GRACE_BASELINE and "2009" in GRACE_BASELINE
    assert GRACE_RECORD_START == dt.date(2002, 4, 1)


# ---------------------------------------------------------------------------
# Month arithmetic
# ---------------------------------------------------------------------------

def test_month_label_from_time_bounds():
    # Bounds [94, 120] days since 2002-01-01 midpoint to 2002-04-16.
    assert gr._month_label(107, 94, 120) == (2002, 4)
    assert gr._month_label(140, 120, 150) == (2002, 5)


def test_window_months_inclusive():
    months = gr._window_months(dt.date(2019, 11, 3), dt.date(2020, 2, 27))
    assert months == [(2019, 11), (2019, 12), (2020, 1), (2020, 2)]


def test_norm_lon_360():
    assert gr._norm_lon_360(-125.0) == 235.0
    assert gr._norm_lon_360(10.0) == 10.0
    assert gr._norm_lon_360(-180.0) == 180.0


def test_lon_windows_normal():
    assert gr._lon_windows_360((-125.0, 30.0, -110.0, 45.0)) == [(235.0, 250.0)]


def test_lon_windows_antimeridian():
    # 170E..170W is contiguous on the 0..360 axis (170..190): single
    # window, no wrap needed.
    assert gr._lon_windows_360((170.0, -50.0, -170.0, -30.0)) == [(170.0, 190.0)]


def test_lon_windows_prime_meridian_wrap():
    # -10..10 straddles the prime meridian: wraps on 0..360.
    assert gr._lon_windows_360((-10.0, -50.0, 10.0, -30.0)) == [
        (350.0, 360.0), (0.0, 10.0)]


def test_lon_windows_full_globe():
    assert gr._lon_windows_360((-180.0, -90.0, 180.0, 90.0)) == [(0.0, 360.0)]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def test_parse_months_labeled_from_bounds():
    ds = _fake_solution([(2020, 1), (2020, 2), (2020, 4)],
                        [-1.0, 1.0], [0.0, 90.0])
    parsed = gr._parse_grace_dataset(ds, _fake_mask([-1.0, 1.0], [0.0, 90.0]),
                                     (-1.0, -2.0, 91.0, 2.0),
                                     dt.date(2020, 1, 1), dt.date(2020, 4, 30))
    # Requested window 2020-01..2020-04: 4 frames, 2020-03 a gap.
    assert [t.isoformat() for t in parsed["times"]] == [
        "2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01"]
    assert [t.isoformat() for t in parsed["gap_months"]] == ["2020-03-01"]
    assert parsed["values"].shape == (4, 2, 2)
    # Frames hold value + k from the fake: Jan=1, Feb=2, Mar=NaN, Apr=3.
    assert parsed["values"][0, 0, 0] == 1.0
    assert parsed["values"][1, 0, 0] == 2.0
    assert np.all(np.isnan(parsed["values"][2]))
    assert parsed["values"][3, 0, 0] == 3.0


def test_parse_gap_frames_are_never_interpolated():
    ds = _fake_solution([(2020, 1)], [-1.0], [0.0])
    parsed = gr._parse_grace_dataset(ds, _fake_mask([-1.0], [0.0]),
                                     (-1.0, -2.0, 1.0, 2.0),
                                     dt.date(2020, 1, 1), dt.date(2020, 3, 31))
    assert len(parsed["times"]) == 3  # window stays calendar-complete
    assert len(parsed["gap_months"]) == 2
    assert np.all(np.isnan(parsed["values"][1]))
    assert np.all(np.isnan(parsed["values"][2]))


def test_parse_applies_land_mask():
    lats = [-1.0, 1.0]
    lons = [0.0, 90.0]
    ds = _fake_solution([(2020, 1)], lats, lons, value=5.0)
    # Only the north-east cell is land.
    mask = _fake_mask(lats, lons, land_where=([1], [1]))
    parsed = gr._parse_grace_dataset(ds, mask, (-1.0, -2.0, 91.0, 2.0),
                                     dt.date(2020, 1, 1), dt.date(2020, 1, 31))
    frame = parsed["values"][0]
    assert frame[1, 1] == 5.0
    assert np.isnan(frame[0, 0]) and np.isnan(frame[0, 1]) and np.isnan(frame[1, 0])


def test_parse_normalizes_lons_to_minus180_180():
    lats = [-1.0]
    lons = [0.0, 90.0, 180.0, 270.0]  # file axis 0..360
    ds = _fake_solution([(2020, 1)], lats, lons)
    parsed = gr._parse_grace_dataset(ds, _fake_mask(lats, lons),
                                     (-180.0, -2.0, 180.0, 2.0),
                                     dt.date(2020, 1, 1), dt.date(2020, 1, 31))
    # 180 stays 180 (strict > 180); nothing maps to -180 here.
    assert list(parsed["lons"]) == [-90.0, 0.0, 90.0, 180.0]
    assert (np.diff(parsed["lons"]) > 0).all()


def test_parse_antimeridian_bbox():
    lats = [-1.0]
    lons = [0.0, 90.0, 170.0, 180.0, 190.0, 270.0]
    ds = _fake_solution([(2020, 1)], lats, lons)
    parsed = gr._parse_grace_dataset(ds, _fake_mask(lats, lons),
                                     (160.0, -2.0, -160.0, 2.0),
                                     dt.date(2020, 1, 1), dt.date(2020, 1, 31))
    # 160E..160W(=200E) -> single 0..360 window [160, 200]: 170, 180, 190.
    assert list(parsed["lons"]) == [-170.0, 170.0, 180.0]


def test_parse_prime_meridian_wrap_bbox():
    lats = [-1.0]
    lons = [0.0, 90.0, 180.0, 270.0, 355.0]
    ds = _fake_solution([(2020, 1)], lats, lons)
    parsed = gr._parse_grace_dataset(ds, _fake_mask(lats, lons),
                                     (-10.0, -2.0, 10.0, 2.0),
                                     dt.date(2020, 1, 1), dt.date(2020, 1, 31))
    # -10..10 -> windows [350, 360) and [0, 10]: 0 and 355 match.
    assert list(parsed["lons"]) == [-5.0, 0.0]


def test_parse_masked_array_data():
    lats = [-1.0]
    lons = [0.0]
    data = np.ma.masked_array([[[7.0]]], mask=[[[True]]])
    ds = _FakeDS({
        "time": _FakeVar([_days(2020, 1, 15)]),
        "time_bounds": _FakeVar([[_days(2020, 1, 1), _days(2020, 2, 1)]]),
        "lat": _FakeVar(lats), "lon": _FakeVar(lons),
        "lwe_thickness": data,
    })
    parsed = gr._parse_grace_dataset(ds, _fake_mask(lats, lons),
                                     (-1.0, -2.0, 1.0, 2.0),
                                     dt.date(2020, 1, 1), dt.date(2020, 1, 31))
    assert np.all(np.isnan(parsed["values"][0]))


def test_parse_empty_lat_selection_raises():
    ds = _fake_solution([(2020, 1)], [-1.0], [0.0])
    with pytest.raises(ValueError):
        gr._parse_grace_dataset(ds, _fake_mask([-1.0], [0.0]),
                                (0.0, 50.0, 1.0, 60.0),
                                dt.date(2020, 1, 1), dt.date(2020, 1, 31))


def test_parse_mask_shape_mismatch_raises():
    ds = _fake_solution([(2020, 1)], [-1.0], [0.0])
    with pytest.raises(ValueError):
        gr._parse_grace_dataset(ds, np.zeros((2, 2)),
                                (-1.0, -2.0, 1.0, 2.0),
                                dt.date(2020, 1, 1), dt.date(2020, 1, 31))


# ---------------------------------------------------------------------------
# WaterField
# ---------------------------------------------------------------------------

def _waterfield(**kw):
    return WaterField.synthetic(bbox=(-125.0, 30.0, -110.0, 45.0),
                                start="2020-01-01", end="2020-12-01",
                                resolution=2.5, **kw)


def test_synthetic_deterministic():
    a = _waterfield()
    b = _waterfield()
    assert np.allclose(a.values, b.values, equal_nan=True)
    assert a.gap_months == [dt.date(2020, 6, 1)]
    assert a.units == "cm"
    assert "2004" in a.anomaly_baseline


def test_synthetic_shape_and_gaps():
    f = _waterfield()
    assert len(f) == 12
    assert f.n_gap_months == 1
    assert f.shape == (12, 7, 7)  # 30..45 and -125..-110 at 2.5°
    assert f.is_gap(5)
    assert not f.is_gap(0)


def test_synthetic_gap_frame_all_nan():
    f = _waterfield()
    assert np.all(np.isnan(f.values[5]))
    assert np.isnan(f.spatial_mean(5))


def test_spatial_mean_ignores_nan():
    f = _waterfield()
    m = f.spatial_mean(0)
    assert np.isfinite(m)
    assert m == pytest.approx(float(np.nanmean(f.values[0])))


def test_len_bounds_shape_time_range():
    f = _waterfield()
    assert len(f) == 12
    assert f.bounds == (-125.0, 30.0, -110.0, 45.0)
    assert f.shape == (12, 7, 7)
    assert f.time_range == (dt.date(2020, 1, 1), dt.date(2020, 12, 1))


def test_select_time_keeps_gaps_in_range():
    f = _waterfield()
    sub = f.select_time("2020-05-01", "2020-07-31")
    assert len(sub) == 3
    assert sub.gap_months == [dt.date(2020, 6, 1)]
    out = f.select_time("2020-01-01", "2020-03-31")
    assert out.gap_months == []


def test_select_bbox():
    f = _waterfield()
    sub = f.select_bbox((-125.0, 30.0, -117.5, 37.5))
    assert sub.shape == (12, 4, 4)
    assert sub.bounds == (-125.0, 30.0, -117.5, 37.5)
    # Gap months survive spatial subsetting.
    assert sub.gap_months == [dt.date(2020, 6, 1)]


def test_round_trip_dict():
    f = _waterfield()
    g = WaterField.from_dict(f.to_dict())
    assert np.allclose(g.values, f.values, equal_nan=True)
    assert g.times == f.times
    assert g.gap_months == f.gap_months
    assert g.units == "cm"
    assert g.anomaly_baseline == f.anomaly_baseline
    assert g.provenance == f.provenance


def test_round_trip_json(tmp_path):
    f = _waterfield()
    path = str(tmp_path / "w.json")
    assert f.to_json(path) == path
    g = WaterField.from_json(path)
    assert np.allclose(g.values, f.values, equal_nan=True)
    assert g.times == f.times
    assert g.gap_months == f.gap_months


def test_dict_round_trip_starts_from_empty_gaps():
    f = _waterfield(gap=())
    assert f.gap_months == []
    g = WaterField.from_dict(f.to_dict())
    assert g.gap_months == []


# ---------------------------------------------------------------------------
# fetch_grace validation + end-to-end with fake datasets
# ---------------------------------------------------------------------------

def _fake_fetch_env(monkeypatch, months, lats, lons, mask=None, value=1.0):
    sol_ds = _fake_solution(months, lats, lons, value=value)
    mask_ds = _FakeDS({"LO_val": _FakeVar(
        mask if mask is not None else _fake_mask(lats, lons))})

    class _FakeNC:
        @staticmethod
        def Dataset(path, mode):  # noqa: N802 - matches netCDF4 API
            if path.endswith(GRACE_SOLUTION_FILE):
                return sol_ds
            return mask_ds

    monkeypatch.setattr(gr, "_require_netcdf4", lambda: _FakeNC)
    monkeypatch.setattr(gr, "ensure_grace_files",
                        lambda **kw: ("/fake/" + GRACE_SOLUTION_FILE,
                                      "/fake/" + GRACE_MASK_FILE,
                                      {"solution": {"file_url": GRACE_SOLUTION_URL,
                                                    "sha256": "s" * 64},
                                       "mask": {"file_url": GRACE_MASK_URL,
                                                "sha256": "m" * 64}}))


def test_fetch_rejects_pre_record_start():
    with pytest.raises(ValueError, match="2002-04-01"):
        fetch_grace((-10.0, 30.0, 10.0, 50.0), "2001-01-01", "2002-05-01",
                    cache_dir="/nonexistent-cache")


def test_fetch_rejects_future_end():
    with pytest.raises(ValueError, match="in the future"):
        fetch_grace((-10.0, 30.0, 10.0, 50.0), "2020-01-01", "2100-01-01",
                    cache_dir="/nonexistent-cache")


def test_fetch_rejects_inverted_range():
    with pytest.raises(ValueError, match="before start"):
        fetch_grace((-10.0, 30.0, 10.0, 50.0), "2020-02-01", "2020-01-01",
                    cache_dir="/nonexistent-cache")


def test_fetch_end_to_end_fake(monkeypatch):
    lats = [30.0, 35.0]
    lons = [235.0, 240.0]  # 0..360 axis -> -125, -120
    _fake_fetch_env(monkeypatch, [(2020, 1), (2020, 2), (2020, 4)],
                   lats, lons)
    f = fetch_grace((-125.5, 29.5, -119.5, 35.5),
                    "2020-01-01", "2020-04-30", cache_dir="/fake")
    assert isinstance(f, WaterField)
    assert [t.isoformat() for t in f.times] == [
        "2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01"]
    assert [t.isoformat() for t in f.gap_months] == ["2020-03-01"]
    assert list(f.lons) == [-125.0, -120.0]
    assert f.units == "cm"
    prov = f.provenance
    assert prov["product_version"] == "RL06.3"
    assert "2004" in prov["anomaly_baseline"]
    assert prov["solution_url"] == GRACE_SOLUTION_URL
    assert prov["solution_sha256"] == "s" * 64
    assert prov["mask_url"] == GRACE_MASK_URL
    assert prov["gap_months"] == ["2020-03-01"]
    assert prov["ocean_masked"] is True
    assert prov["cache_hit"] is False


def test_fetch_applies_land_mask_fake(monkeypatch):
    lats = [30.0, 35.0]
    lons = [235.0, 240.0]
    mask = np.array([[1.0, 0.0], [0.0, 0.0]])  # only one land cell
    _fake_fetch_env(monkeypatch, [(2020, 1)], lats, lons, mask=mask, value=9.0)
    f = fetch_grace((-125.5, 29.5, -119.5, 35.5),
                    "2020-01-01", "2020-01-31", cache_dir="/fake")
    frame = f.values[0]
    assert frame[0, 0] == 9.0
    assert np.isnan(frame[0, 1]) and np.isnan(frame[1, 0]) and np.isnan(frame[1, 1])


def test_fetch_mask_variable_fallback(monkeypatch):
    """Mask file without 'LO_val' falls back to its first 2D variable."""
    lats = [30.0]
    lons = [235.0]
    sol_ds = _fake_solution([(2020, 1)], lats, lons)
    mask_ds = _FakeDS({"other_mask": _FakeVar(np.ones((1, 1)))})

    class _FakeNC:
        @staticmethod
        def Dataset(path, mode):  # noqa: N802
            return sol_ds if path.endswith(GRACE_SOLUTION_FILE) else mask_ds

    monkeypatch.setattr(gr, "_require_netcdf4", lambda: _FakeNC)
    monkeypatch.setattr(gr, "ensure_grace_files",
                        lambda **kw: ("/fake/" + GRACE_SOLUTION_FILE,
                                      "/fake/" + GRACE_MASK_FILE,
                                      {"solution": {"file_url": GRACE_SOLUTION_URL,
                                                    "sha256": "s"},
                                       "mask": {"file_url": GRACE_MASK_URL,
                                                "sha256": "m"}}))
    f = fetch_grace((-125.5, 29.5, -119.5, 35.5),
                    "2020-01-01", "2020-01-31", cache_dir="/fake")
    assert np.isfinite(f.values[0, 0, 0])


def test_fetch_mask_without_2d_variable_raises(monkeypatch):
    lats = [30.0]
    lons = [235.0]
    sol_ds = _fake_solution([(2020, 1)], lats, lons)
    mask_ds = _FakeDS({"scalar": _FakeVar(1.0)})

    class _FakeNC:
        @staticmethod
        def Dataset(path, mode):  # noqa: N802
            return sol_ds if path.endswith(GRACE_SOLUTION_FILE) else mask_ds

    monkeypatch.setattr(gr, "_require_netcdf4", lambda: _FakeNC)
    monkeypatch.setattr(gr, "ensure_grace_files",
                        lambda **kw: ("/fake/" + GRACE_SOLUTION_FILE,
                                      "/fake/" + GRACE_MASK_FILE,
                                      {"solution": {"file_url": GRACE_SOLUTION_URL,
                                                    "sha256": "s"},
                                       "mask": {"file_url": GRACE_MASK_URL,
                                                "sha256": "m"}}))
    with pytest.raises(RuntimeError, match="no 2D variable"):
        fetch_grace((-125.5, 29.5, -119.5, 35.5),
                    "2020-01-01", "2020-01-31", cache_dir="/fake")


# ---------------------------------------------------------------------------
# Cache / download plumbing (mocked)
# ---------------------------------------------------------------------------

def _fake_download(payload=b"GRACE-DATA"):
    def _dl(url, timeout=900):
        return payload
    return _dl


def test_ensure_grace_files_downloads_once(monkeypatch, tmp_path):
    monkeypatch.setattr(gr, "_download_bytes", _fake_download())
    monkeypatch.setattr(gr, "_http_head_last_modified", lambda url, timeout=60: None)
    sol, mask, prov = gr.ensure_grace_files(cache_dir=str(tmp_path))
    assert os.path.isfile(sol) and os.path.isfile(mask)
    assert os.path.isfile(sol + ".sha256") and os.path.isfile(mask + ".sha256")
    assert prov["solution"]["downloaded"] is True
    # Second call reuses the verified cache without downloading.
    calls = []

    def _boom(url, timeout=900):
        calls.append(url)
        raise AssertionError("should not download")

    monkeypatch.setattr(gr, "_download_bytes", _boom)
    monkeypatch.setattr(gr, "_http_head_last_modified", lambda url, timeout=60: None)
    sol2, _, prov2 = gr.ensure_grace_files(cache_dir=str(tmp_path))
    assert sol2 == sol
    assert prov2["solution"]["cache_hit"] is True
    assert calls == []


def test_ensure_grace_files_corrupt_cache_redownloads(monkeypatch, tmp_path):
    monkeypatch.setattr(gr, "_download_bytes", _fake_download(b"v1"))
    monkeypatch.setattr(gr, "_http_head_last_modified", lambda url, timeout=60: None)
    sol, _, _ = gr.ensure_grace_files(cache_dir=str(tmp_path))
    # Corrupt the solution file in place.
    with open(sol, "wb") as fh:
        fh.write(b"corrupted")
    monkeypatch.setattr(gr, "_download_bytes", _fake_download(b"v2"))
    sol2, _, prov2 = gr.ensure_grace_files(cache_dir=str(tmp_path))
    assert prov2["solution"]["downloaded"] is True
    with open(sol2, "rb") as fh:
        assert fh.read() == b"v2"


def test_ensure_grace_files_refresh_forces_redownload(monkeypatch, tmp_path):
    monkeypatch.setattr(gr, "_download_bytes", _fake_download(b"v1"))
    monkeypatch.setattr(gr, "_http_head_last_modified", lambda url, timeout=60: None)
    gr.ensure_grace_files(cache_dir=str(tmp_path))
    monkeypatch.setattr(gr, "_download_bytes", _fake_download(b"v2"))
    _, _, prov = gr.ensure_grace_files(cache_dir=str(tmp_path), refresh=True)
    assert prov["solution"]["downloaded"] is True


def test_ensure_stale_cache_revalidates(monkeypatch, tmp_path):
    monkeypatch.setattr(gr, "_download_bytes", _fake_download(b"v1"))
    monkeypatch.setattr(gr, "_http_head_last_modified",
                        lambda url, timeout=60: "Wed, 24 Aug 2026 00:00:00 GMT")
    gr.ensure_grace_files(cache_dir=str(tmp_path), max_cache_age_days=0)
    # Same Last-Modified: no re-download.
    calls = []
    monkeypatch.setattr(gr, "_download_bytes", lambda url, timeout=900:
                        calls.append(url) or b"v2")
    _, _, prov = gr.ensure_grace_files(cache_dir=str(tmp_path),
                                       max_cache_age_days=0)
    assert calls == []
    assert prov["solution"]["cache_hit"] is True


def test_fetch_requires_netcdf4(monkeypatch):
    """The optional-dependency error is raised before any download."""
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "netCDF4":
            raise ImportError("no netCDF4")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    monkeypatch.setattr(gr, "ensure_grace_files",
                        lambda **kw: ("/x", "/y",
                                      {"solution": {"file_url": "", "sha256": ""},
                                       "mask": {"file_url": "", "sha256": ""}}))
    with pytest.raises(ImportError, match="netCDF4"):
        fetch_grace((-10.0, 30.0, 10.0, 50.0), "2020-01-01", "2020-02-01",
                    cache_dir="/fake")


# ---------------------------------------------------------------------------
# Live test (opt-in)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(os.environ.get("SURVEY_CURRENTS_LIVE") != "1",
                    reason="live CSR download; set SURVEY_CURRENTS_LIVE=1")
def test_live_fetch_small_bbox():
    f = fetch_grace((-125.0, 30.0, -110.0, 45.0), "2019-01-01", "2019-12-31")
    assert len(f) == 12
    assert f.units == "cm"
    assert (np.diff(f.lons) > 0).all() and f.lons.min() >= -180.0
    # The 2017-07..2018-05 mission gap must be present when requested.
    g = fetch_grace((-125.0, 30.0, -110.0, 45.0), "2017-01-01", "2018-12-31")
    assert [t.isoformat() for t in g.gap_months][:3] == [
        "2017-02-01", "2017-07-01", "2017-08-01"]
    assert "2017-10-01" in [t.isoformat() for t in g.gap_months]
    assert "2018-05-01" in [t.isoformat() for t in g.gap_months]
    assert np.all(np.isnan(g.values[g.times.index(dt.date(2017, 10, 1))]))

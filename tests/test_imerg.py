"""Tests for the GPM IMERG precipitation adapter (currents.imerg).

Fully offline: HTTP/auth are mocked, HDF5 fixtures are built with h5py
in-memory (skipped when h5py is missing), and RainField.synthetic()
covers the field-model tests.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import io
import json

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

import currents.imerg as imerg
from currents.imerg import (
    CredentialsMissing,
    RainField,
    accumulate_daily,
    fetch_imerg,
    imerg_file_url,
    imerg_granule_name,
    imerg_index_windows,
    imerg_sample_days,
    normalize_imerg_accumulate,
    normalize_imerg_run,
)


# ---------------------------------------------------------------------------
# normalization helpers
# ---------------------------------------------------------------------------


def test_normalize_run_ok():
    assert normalize_imerg_run("late") == "late"
    assert normalize_imerg_run(" Late ") == "late"
    assert normalize_imerg_run("EARLY") == "early"
    assert normalize_imerg_run("Final") == "final"


def test_normalize_run_bad():
    with pytest.raises(ValueError, match="run must be one of"):
        normalize_imerg_run("nrt")


def test_normalize_accumulate_ok():
    assert normalize_imerg_accumulate("daily") == "daily"
    assert normalize_imerg_accumulate("NATIVE") == "native"


def test_normalize_accumulate_bad():
    with pytest.raises(ValueError, match="accumulate must be one of"):
        normalize_imerg_accumulate("monthly")


# ---------------------------------------------------------------------------
# granule naming (verified against the live GES DISC catalog 2026-09-26)
# ---------------------------------------------------------------------------


def test_granule_name_final():
    day = dt.date(2024, 1, 1)
    assert imerg_granule_name(day, 0, "final") == (
        "3B-HHR.MS.MRG.3IMERG.20240101-S000000-E002959.0000.V07B.HDF5")
    assert imerg_granule_name(day, 1, "final") == (
        "3B-HHR.MS.MRG.3IMERG.20240101-S003000-E005959.0030.V07B.HDF5")
    assert imerg_granule_name(day, 47, "final") == (
        "3B-HHR.MS.MRG.3IMERG.20240101-S233000-E235959.2330.V07B.HDF5")


def test_granule_name_early_late_infix():
    day = dt.date(2024, 1, 1)
    assert "-E.MS.MRG.3IMERG" in imerg_granule_name(day, 0, "early")
    assert "-L.MS.MRG.3IMERG" in imerg_granule_name(day, 0, "late")
    assert imerg_granule_name(day, 5, "early").endswith(".0230.V07B.HDF5")


def test_granule_name_bad_slot():
    with pytest.raises(ValueError, match="half_hour must be 0..47"):
        imerg_granule_name(dt.date(2024, 1, 1), 48, "late")
    with pytest.raises(ValueError, match="half_hour must be 0..47"):
        imerg_granule_name(dt.date(2024, 1, 1), -1, "late")


def test_file_url_layout():
    url = imerg_file_url(dt.date(2024, 1, 1), 0, "late")
    assert url == ("https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/"
                   "GPM_3IMERGHHL.07/2024/001/"
                   "3B-HHR-L.MS.MRG.3IMERG.20240101-S000000-E002959.0000.V07B.HDF5")
    url = imerg_file_url(dt.date(2024, 3, 5), 12, "early")
    assert "/GPM_3IMERGHHE.07/2024/065/" in url
    url = imerg_file_url(dt.date(2024, 3, 5), 12, "final")
    assert "/GPM_3IMERGHH.07/2024/065/" in url


def test_sample_days():
    days = imerg_sample_days(dt.date(2024, 1, 1), dt.date(2024, 1, 5), 2)
    assert days == [dt.date(2024, 1, 1), dt.date(2024, 1, 3), dt.date(2024, 1, 5)]
    with pytest.raises(ValueError, match="stride_days"):
        imerg_sample_days(dt.date(2024, 1, 1), dt.date(2024, 1, 2), 0)


# ---------------------------------------------------------------------------
# fake HDF5 fixture (mimics the GES DISC granule layout)
# ---------------------------------------------------------------------------


def _fake_imerg_payload(seed=3, ny=18, nx=36, order="time-lat-lon"):
    """Build a minimal IMERG-like HDF5 granule in memory.

    ``order`` is the precipitation dataspace order: the live .das
    reports (time, lon, lat), so both orders are exercised.
    """
    rng = np.random.default_rng(seed)
    lats = np.linspace(-89.95, 89.95, ny)
    lons = np.linspace(-179.95, 179.95, nx)
    rates = np.abs(rng.normal(0.5, 1.5, (ny, nx))).astype(np.float32)
    rates[0, 0] = -9999.9  # a fill cell
    buf = io.BytesIO()
    with h5py.File(buf, "w") as f:
        grid = f.create_group("Grid")
        grid.create_dataset("lat", data=lats)
        grid.create_dataset("lon", data=lons)
        if order == "time-lat-lon":
            grid.create_dataset("precipitation", data=rates[None, :, :])
        else:  # the live .das order: (time, lon, lat)
            grid.create_dataset("precipitation",
                                data=np.ascontiguousarray(rates[None, :, :].transpose(0, 2, 1)))
        grid["precipitation"].attrs["_FillValue"] = np.float32(-9999.9)
    return buf.getvalue(), lats, lons, rates


def test_parse_imerg_bytes_lat_lon_order():
    payload, lats, lons, rates = _fake_imerg_payload(order="time-lat-lon")
    pl, pn, values = imerg._parse_imerg_bytes(payload)
    assert pl.shape == (18,) and pn.shape == (36,)
    assert values.shape == (18, 36)
    assert np.isnan(values[0, 0])  # fill -> NaN
    np.testing.assert_allclose(values[1:, 1:], rates[1:, 1:], rtol=1e-5)


def test_parse_imerg_bytes_das_order():
    # The live .das reports dims (time, lon, lat): the parser must still
    # land on (nlat, nlon) via the size-based axis resolution.
    payload, lats, lons, rates = _fake_imerg_payload(order="time-lon-lat")
    pl, pn, values = imerg._parse_imerg_bytes(payload)
    assert values.shape == (18, 36)
    assert np.isnan(values[0, 0])
    np.testing.assert_allclose(values[1:, 1:], rates[1:, 1:], rtol=1e-5)


def test_parse_imerg_bytes_bad_magic():
    with pytest.raises(ValueError, match="not an HDF5 payload"):
        imerg._parse_imerg_bytes(b"<html>not hdf5</html>")


def test_parse_imerg_bytes_missing_dataset():
    buf = io.BytesIO()
    with h5py.File(buf, "w") as f:
        f.create_dataset("Grid/lat", data=np.linspace(-89.95, 89.95, 18))
    with pytest.raises(ValueError, match="missing dataset"):
        imerg._parse_imerg_bytes(buf.getvalue())


def test_axis_permutation_ambiguous():
    with pytest.raises(ValueError, match="ambiguous"):
        imerg._imerg_axis_permutation((1, 1800, 1800), 1800, 3600)
    with pytest.raises(ValueError, match="3-D dataspace"):
        imerg._imerg_axis_permutation((1800, 3600), 1800, 3600)


def test_index_windows_basic():
    lats = np.linspace(-89.95, 89.95, 1800)
    lons = np.linspace(-179.95, 179.95, 3600)
    wins = imerg_index_windows((-125.0, 25.0, -66.0, 49.0), lats, lons)
    assert len(wins) == 1
    i0, i1, j0, j1 = wins[0]
    # i1/j1 are exclusive slice ends: the included cells lie within the bbox
    assert lons[i0] >= -125.0 and lons[i1 - 1] <= -66.0
    assert lats[j0] >= 25.0 and lats[j1 - 1] <= 49.0


def test_index_windows_antimeridian():
    lats = np.linspace(-89.95, 89.95, 1800)
    lons = np.linspace(-179.95, 179.95, 3600)
    wins = imerg_index_windows((170.0, -10.0, -170.0, 10.0), lats, lons)
    assert len(wins) == 2


# ---------------------------------------------------------------------------
# daily accumulation
# ---------------------------------------------------------------------------


def _frames(day, nslots, ny=4, nx=5, seed=1):
    rng = np.random.default_rng(seed)
    times = [dt.datetime(day.year, day.month, day.day, tzinfo=dt.timezone.utc)
             + dt.timedelta(minutes=30 * i) for i in range(nslots)]
    frames = [np.abs(rng.normal(1.0, 0.5, (ny, nx))) for _ in range(nslots)]
    return times, frames


def test_accumulate_daily_sums_rates():
    day = dt.date(2024, 6, 1)
    times, frames = _frames(day, 48)
    frames[0][:] = 2.0  # 2 mm/hr for one slot -> 1 mm
    frames[1][:] = np.nan  # missing slot contributes 0
    dtimes, daily = accumulate_daily(times, frames)
    assert len(dtimes) == 1
    assert dtimes[0] == dt.datetime(2024, 6, 1, tzinfo=dt.timezone.utc)
    expected = sum(np.where(np.isnan(f), 0.0, f) * 0.5 for f in frames)
    np.testing.assert_allclose(daily[0], expected)


def test_accumulate_daily_all_nan_stays_nan():
    day = dt.date(2024, 6, 1)
    times, frames = _frames(day, 4)
    for f in frames:
        f[:, :] = np.nan
    _dt2, daily = accumulate_daily(times, frames)
    assert np.isnan(daily[0]).all()


def test_accumulate_daily_multi_day():
    t1, f1 = _frames(dt.date(2024, 6, 1), 48, seed=1)
    t2, f2 = _frames(dt.date(2024, 6, 2), 48, seed=2)
    dtimes, daily = accumulate_daily(t1 + t2, f1 + f2)
    assert len(dtimes) == 2
    assert daily.shape == (2, 4, 5)


def test_accumulate_daily_length_mismatch():
    with pytest.raises(ValueError, match="differ in length"):
        accumulate_daily([dt.datetime(2024, 6, 1)], [])


# ---------------------------------------------------------------------------
# RainField model
# ---------------------------------------------------------------------------


def test_rainfield_synthetic_deterministic():
    a = RainField.synthetic(seed=7)
    b = RainField.synthetic(seed=7)
    np.testing.assert_array_equal(a.values, b.values)
    assert a.units == "mm/day"
    assert a.accumulate == "daily"
    assert np.isnan(a.values).any()  # the NaN corner


def test_rainfield_synthetic_native():
    f = RainField.synthetic(accumulate="native",
                            start="2024-01-01", end="2024-01-02")
    assert f.units == "mm/hr"
    assert len(f.times) == 96  # 2 days x 48 slots


def test_rainfield_shape_validation():
    with pytest.raises(ValueError):
        RainField(times=[dt.datetime(2024, 1, 1)],
                  lats=np.array([0.0, 1.0]), lons=np.array([0.0]),
                  values=np.zeros((1, 3, 1)))


def test_rainfield_select_time():
    f = RainField.synthetic(start="2024-01-01", end="2024-01-05")
    sub = f.select_time("2024-01-02", "2024-01-03")
    assert len(sub.times) == 2
    assert sub.times[0].date() == dt.date(2024, 1, 2)


def test_rainfield_select_bbox():
    f = RainField.synthetic()
    sub = f.select_bbox((-120.0, 30.0, -100.0, 40.0))
    assert sub.lons.min() >= -120.0 and sub.lons.max() <= -100.0
    assert sub.lats.min() >= 30.0 and sub.lats.max() <= 40.0
    assert sub.values.shape[1:] == (len(sub.lats), len(sub.lons))


def test_rainfield_select_bbox_no_overlap():
    f = RainField.synthetic()
    with pytest.raises(ValueError, match="no overlap"):
        f.select_bbox((150.0, 30.0, 160.0, 40.0))


def test_rainfield_json_roundtrip(tmp_path):
    f = RainField.synthetic(start="2024-01-01", end="2024-01-02")
    f.provenance["note"] = "roundtrip"
    path = str(tmp_path / "rain.json")
    f.to_json(path)
    g = RainField.from_json(path)
    assert g.times == f.times
    assert g.units == f.units and g.run == f.run
    np.testing.assert_allclose(g.values, f.values, equal_nan=True)
    assert g.provenance["note"] == "roundtrip"


def test_rainfield_from_dict_missing_keys():
    with pytest.raises(ValueError, match="missing keys"):
        RainField.from_dict({"times": []})


def test_rainfield_total():
    f = RainField.synthetic(start="2024-01-01", end="2024-01-01")
    assert f.total(0) == pytest.approx(float(np.nansum(f.values[0])))


# ---------------------------------------------------------------------------
# fetch_imerg (mocked network + auth)
# ---------------------------------------------------------------------------

_BBOX = (-125.0, 25.0, -66.0, 49.0)


def _install_fake_download(monkeypatch, payload, calls):
    def fake_download(url, opener, timeout=300):
        calls.append(url)
        return payload
    monkeypatch.setattr(imerg, "_imerg_download", fake_download)
    monkeypatch.setattr(imerg, "_imerg_opener", lambda: object())


def test_fetch_imerg_daily(monkeypatch):
    payload, _lats, _lons, _rates = _fake_imerg_payload()
    calls = []
    _install_fake_download(monkeypatch, payload, calls)
    f = fetch_imerg(_BBOX, "2024-01-01", "2024-01-02", run="late")
    assert len(calls) == 96  # 2 days x 48 slots
    assert "GPM_3IMERGHHL.07" in calls[0]
    assert f.units == "mm/day"
    assert f.accumulate == "daily"
    assert f.run == "late"
    assert len(f.times) == 2
    assert f.values.shape[0] == 2
    prov = f.provenance
    assert prov["n_files"] == 96
    assert prov["run"] == "late"
    assert prov["day_coverage"]["2024-01-01"] == {"expected": 48, "retrieved": 48}
    assert len(prov["files"]) == 96
    assert all(len(r["sha256"]) == 64 for r in prov["files"])
    assert prov["files"][0]["sha256"] == hashlib.sha256(payload).hexdigest()
    assert "retrieved_utc" in prov


def test_fetch_imerg_native(monkeypatch):
    payload, _lats, _lons, _rates = _fake_imerg_payload()
    calls = []
    _install_fake_download(monkeypatch, payload, calls)
    f = fetch_imerg(_BBOX, "2024-01-01", "2024-01-01",
                    accumulate="native", run="final")
    assert f.units == "mm/hr"
    assert len(f.times) == 48
    assert (f.times[1] - f.times[0]) == dt.timedelta(minutes=30)


def test_fetch_imerg_skips_404_slots(monkeypatch):
    payload, _lats, _lons, _rates = _fake_imerg_payload()
    calls = []

    def fake_download(url, opener, timeout=300):
        calls.append(url)
        if "S003000" in url:  # one slot 404s
            raise FileNotFoundError("404")
        return payload

    monkeypatch.setattr(imerg, "_imerg_download", fake_download)
    monkeypatch.setattr(imerg, "_imerg_opener", lambda: object())
    f = fetch_imerg(_BBOX, "2024-01-01", "2024-01-01")
    prov = f.provenance
    assert prov["n_files"] == 47
    assert prov["day_coverage"]["2024-01-01"] == {"expected": 48, "retrieved": 47}
    assert len(prov["skipped_slots"]) == 1
    assert len(f.times) == 1  # the day still accumulates


def test_fetch_imerg_no_granules_raises(monkeypatch):
    def fake_download(url, opener, timeout=300):
        raise FileNotFoundError("404")
    monkeypatch.setattr(imerg, "_imerg_download", fake_download)
    monkeypatch.setattr(imerg, "_imerg_opener", lambda: object())
    with pytest.raises(ValueError, match="no IMERG granules retrieved"):
        fetch_imerg(_BBOX, "2024-01-01", "2024-01-02")


def test_fetch_imerg_credentials_missing(monkeypatch):
    def no_creds():
        raise CredentialsMissing("no creds")
    monkeypatch.setattr(imerg, "_imerg_opener", no_creds)
    with pytest.raises(CredentialsMissing):
        fetch_imerg(_BBOX, "2024-01-01", "2024-01-02")


def test_fetch_imerg_bad_run():
    with pytest.raises(ValueError, match="run must be one of"):
        fetch_imerg(_BBOX, "2024-01-01", "2024-01-02", run="nrt")


def test_fetch_imerg_before_record():
    with pytest.raises(ValueError, match="1999-01-01"):
        fetch_imerg(_BBOX, "1999-01-01", "1999-02-01", run="final")
    with pytest.raises(ValueError, match="1999-06-01"):
        fetch_imerg(_BBOX, "1999-06-01", "1999-07-01", run="early")


def test_fetch_imerg_end_before_start():
    with pytest.raises(ValueError, match="end .* is before start"):
        fetch_imerg(_BBOX, "2024-02-01", "2024-01-01")


def test_fetch_imerg_future_end():
    future = (dt.date.today() + dt.timedelta(days=10)).isoformat()
    with pytest.raises(ValueError, match="in the future"):
        fetch_imerg(_BBOX, "2024-01-01", future)


def test_credentials_missing_message():
    msg = str(imerg.CredentialsMissing("x"))
    assert "urs.earthdata.nasa.gov/users/new" in imerg.CredentialsMissing.__doc__
    assert "EARTHDATA_USERNAME" in imerg.CredentialsMissing.__doc__


def test_opener_raises_without_creds(monkeypatch):
    monkeypatch.setattr("currents.sst_global.earthdata_credentials",
                        lambda: None)
    with pytest.raises(CredentialsMissing, match="Earthdata Login"):
        imerg._imerg_opener()


def test_download_401_becomes_credentials_missing(monkeypatch):
    import urllib.error
    import urllib.request

    class FakeOpener:
        def open(self, req, timeout=None):
            raise urllib.error.HTTPError(req.full_url, 401, "Unauthorized",
                                         {}, None)

    with pytest.raises(CredentialsMissing, match="rejected the credentials"):
        imerg._imerg_download("https://example.com/x", FakeOpener())


def test_download_404_becomes_file_not_found(monkeypatch):
    import urllib.error

    class FakeOpener:
        def open(self, req, timeout=None):
            raise urllib.error.HTTPError(req.full_url, 404, "Not Found",
                                         {}, None)

    with pytest.raises(FileNotFoundError, match="granule not found"):
        imerg._imerg_download("https://example.com/x", FakeOpener())


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_has_imerg_commands():
    from currents.cli import build_parser
    p = build_parser()
    cmds = {a.dest for a in p._subparsers._group_actions[0]._choices_actions}
    # argparse internals: simpler to check the choices dict
    choices = p._subparsers._group_actions[0].choices
    assert "fetch-imerg" in choices
    assert "rain-synthetic" in choices


def test_cli_rain_synthetic(tmp_path, monkeypatch, capsys):
    from currents import cli
    monkeypatch.chdir(tmp_path)
    rc = cli.main(["rain-synthetic", "--bbox=-125,25,-66,49",
                   "--start", "2024-01-01", "--end", "2024-01-02",
                   "--out", "rain_test"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "wrote synthetic IMERG rain field" in out
    data = json.loads(open("rain_test.json").read())
    assert data["units"] == "mm/day"
    assert data["accumulate"] == "daily"


def test_main_demo(capsys):
    imerg.main_demo()
    assert "[imerg] synthetic field" in capsys.readouterr().out

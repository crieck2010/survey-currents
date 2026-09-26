"""Tests for currents.era5 — fully offline (fake cdsapi via sys.modules)."""

import datetime as dt
import json
import sys
import types

import numpy as np
import pytest

from currents import era5
from currents.era5 import (ERA5_DATASET, ERA5_VARIABLES, Era5Field,
                           era5_area_windows, era5_cds_names,
                           era5_month_chunks, era5_request, era5_sample_plan,
                           fetch_era5, normalize_era5_variables,
                           validate_era5_bbox)


# ---------------------------------------------------------------------------
# fake cdsapi: writes a deterministic NetCDF matching the request
# ---------------------------------------------------------------------------

_LONG_TO_SHORT = {
    "10m_u_component_of_wind": "u10",
    "10m_v_component_of_wind": "v10",
    "mean_sea_level_pressure": "msl",
    "2m_temperature": "t2m",
    "total_precipitation": "tp",
}
_ONWIRE = {"u10": 5.0, "v10": 2.0, "msl": 101300.0, "t2m": 290.0, "tp": 0.001}

RECORDED_REQUESTS = []


def _build_fake_era5_nc(path, request, time_var="valid_time"):
    import netCDF4
    n, w, s, e = (float(x) for x in request["area"])
    year, month = int(request["year"]), int(request["month"])
    lats = np.arange(n, s - 1e-9, -0.25)          # descending, like the CDS
    w360, e360 = w % 360.0, e % 360.0
    if e360 <= w360:
        e360 += 360.0
    lons = np.arange(w360, e360 + 1e-9, 0.25)     # 0-360, like the CDS
    stamps = []
    for day in request["day"]:
        for hhmm in request["time"]:
            stamps.append(dt.datetime(year, month, int(day),
                                      int(hhmm[:2]), tzinfo=dt.timezone.utc))
    epoch = dt.datetime(1900, 1, 1, tzinfo=dt.timezone.utc)
    hours = [(t - epoch).total_seconds() / 3600.0 for t in stamps]
    ds = netCDF4.Dataset(path, "w")
    ds.createDimension("time", len(stamps))
    ds.createDimension("latitude", len(lats))
    ds.createDimension("longitude", len(lons))
    tv = ds.createVariable(time_var, "f8", ("time",))
    tv.units = "hours since 1900-01-01 00:00:00"
    tv[:] = hours
    ds.createVariable("latitude", "f8", ("latitude",))[:] = lats
    ds.createVariable("longitude", "f8", ("longitude",))[:] = lons
    for long_name in request["variable"]:
        short = _LONG_TO_SHORT[long_name]
        var = ds.createVariable(short, "f4",
                                ("time", "latitude", "longitude"))
        var[:] = np.full((len(stamps), len(lats), len(lons)),
                         _ONWIRE[short], dtype="f4")
    ds.close()


class _FakeResult:
    def __init__(self, request):
        self.request = request

    def download(self, target):
        RECORDED_REQUESTS.append(self.request)
        _build_fake_era5_nc(target, self.request)


class _FakeClient:
    def __init__(self, *a, **k):
        pass

    def retrieve(self, dataset, request):
        assert dataset == ERA5_DATASET, dataset
        return _FakeResult(request)


@pytest.fixture
def fake_cdsapi(monkeypatch):
    RECORDED_REQUESTS.clear()
    module = types.ModuleType("cdsapi")
    module.Client = _FakeClient
    monkeypatch.setitem(sys.modules, "cdsapi", module)
    return module


@pytest.fixture
def failing_cdsapi(monkeypatch):
    class _BadClient:
        def __init__(self, *a, **k):
            raise RuntimeError("Missing/incomplete configuration")
    module = types.ModuleType("cdsapi")
    module.Client = _BadClient
    monkeypatch.setitem(sys.modules, "cdsapi", module)
    return module


# ---------------------------------------------------------------------------
# validation / request building (no cdsapi needed)
# ---------------------------------------------------------------------------


def test_normalize_single_string():
    assert normalize_era5_variables("wind") == ["wind"]


def test_normalize_list_dedupes():
    assert normalize_era5_variables(["msl", "wind", "msl"]) == ["msl", "wind"]


def test_normalize_unknown_raises():
    with pytest.raises(ValueError, match="unknown ERA5 variable"):
        normalize_era5_variables(["wind", "clouds"])


def test_normalize_empty_raises():
    with pytest.raises(ValueError, match="at least one variable"):
        normalize_era5_variables([])


def test_validate_bbox_ok():
    assert validate_era5_bbox((-98.0, 29.0, -97.0, 31.0)) == (-98.0, 29.0, -97.0, 31.0)


def test_validate_bbox_bad_length():
    with pytest.raises(ValueError, match="4 numbers"):
        validate_era5_bbox((-98.0, 29.0, -97.0))


def test_validate_bbox_lat_out_of_range():
    with pytest.raises(ValueError, match="within \\[-90, 90\\]"):
        validate_era5_bbox((-98.0, -95.0, -97.0, 31.0))


def test_validate_bbox_lat_not_ascending():
    with pytest.raises(ValueError, match="ascend"):
        validate_era5_bbox((-98.0, 31.0, -97.0, 29.0))


def test_sample_plan_subdaily():
    d0 = dt.datetime(2024, 1, 6, tzinfo=dt.timezone.utc)
    plan = era5_sample_plan(d0, d0, 6)
    assert plan == [(dt.date(2024, 1, 6),
                     ["00:00", "06:00", "12:00", "18:00"])]


def test_sample_plan_daily():
    d0 = dt.datetime(2024, 1, 6, tzinfo=dt.timezone.utc)
    d1 = dt.datetime(2024, 1, 8, tzinfo=dt.timezone.utc)
    plan = era5_sample_plan(d0, d1, 24)
    assert [p[0] for p in plan] == [dt.date(2024, 1, 6), dt.date(2024, 1, 7),
                                    dt.date(2024, 1, 8)]
    assert all(p[1] == ["12:00"] for p in plan)


def test_sample_plan_weekly():
    d0 = dt.datetime(2024, 1, 1, tzinfo=dt.timezone.utc)
    d1 = dt.datetime(2024, 1, 22, tzinfo=dt.timezone.utc)
    plan = era5_sample_plan(d0, d1, 168)
    assert [p[0].day for p in plan] == [1, 8, 15, 22]


def test_sample_plan_bad_stride():
    d0 = dt.datetime(2024, 1, 6, tzinfo=dt.timezone.utc)
    with pytest.raises(ValueError, match="stride_hours must be >= 1"):
        era5_sample_plan(d0, d0, 0)


def test_month_chunks_split():
    d0 = dt.datetime(2024, 1, 30, tzinfo=dt.timezone.utc)
    d1 = dt.datetime(2024, 2, 2, tzinfo=dt.timezone.utc)
    chunks = era5_month_chunks(era5_sample_plan(d0, d1, 24))
    assert [(c[0], c[1], c[2]) for c in chunks] == [
        (2024, 1, ["30", "31"]), (2024, 2, ["01", "02"])]
    assert chunks[0][3] == ["12:00"]


def test_area_windows_single():
    assert era5_area_windows((-98.0, 29.0, -97.0, 31.0)) == [(31.0, -98.0, 29.0, -97.0)]


def test_area_windows_antimeridian():
    assert era5_area_windows((170.0, 20.0, -170.0, 30.0)) == [
        (30.0, 170.0, 20.0, 180.0), (30.0, -180.0, 20.0, -170.0)]


def test_area_windows_global():
    assert era5_area_windows((-180.0, -90.0, 180.0, 90.0)) == [(90.0, -180.0, -90.0, 180.0)]


def test_cds_names_wind_pair():
    assert era5_cds_names(["wind"]) == ["10m_u_component_of_wind",
                                        "10m_v_component_of_wind"]
    assert era5_cds_names(["wind", "msl"]) == [
        "10m_u_component_of_wind", "10m_v_component_of_wind",
        "mean_sea_level_pressure"]


def test_request_shape():
    req = era5_request(["wind", "t2m"], 2024, 1, ["06"], ["00:00", "12:00"],
                       (31.0, -98.0, 29.0, -97.0))
    assert req["product_type"] == "reanalysis"
    assert req["variable"] == ["10m_u_component_of_wind",
                               "10m_v_component_of_wind", "2m_temperature"]
    assert req["area"] == [31.0, -98.0, 29.0, -97.0]  # N, W, S, E
    assert req["grid"] == "0.25/0.25"
    assert req["data_format"] == "netcdf"
    assert req["download_format"] == "unarchived"


def test_require_cdsapi_missing(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "cdsapi":
            raise ImportError("no module")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="pip install cdsapi"):
        era5._require_cdsapi()


def test_client_failure_raises_credentials_missing(failing_cdsapi):
    with pytest.raises(era5.CredentialsMissing, match="Copernicus CDS"):
        fetch_era5("wind", (-98.0, 29.0, -97.0, 31.0), "2024-01-06",
                   "2024-01-06", stride_hours=24)


# ---------------------------------------------------------------------------
# fetch_era5 against the fake cdsapi
# ---------------------------------------------------------------------------


def test_fetch_wind_msl_shapes_and_units(fake_cdsapi):
    f = fetch_era5(["wind", "msl"], (-98.0, 29.0, -97.0, 31.0),
                   "2024-01-06", "2024-01-06", stride_hours=6)
    assert isinstance(f, Era5Field)
    assert f.variables == ["wind", "msl"]
    assert f.base_variable == "wind"
    assert f.grids["u10"].shape == (4, 9, 5)  # 4 steps, lats 29..31, lons -98..-97
    assert f.units == {"u10": "m/s", "v10": "m/s", "msl": "hPa"}
    # unit conversions: Pa -> hPa, m/s untouched
    assert float(f.grids["msl"].mean()) == pytest.approx(1013.0)
    assert float(f.grids["u10"].mean()) == pytest.approx(5.0)
    # wind speed = sqrt(5^2 + 2^2)
    assert float(f.values.mean()) == pytest.approx((25 + 4) ** 0.5)


def test_fetch_lon_lat_normalization(fake_cdsapi):
    f = fetch_era5("t2m", (-98.0, 29.0, -97.0, 31.0),
                   "2024-01-06", "2024-01-06", stride_hours=24)
    assert f.lats[0] < f.lats[-1]                      # ascending, not N->S
    assert f.lats[0] == pytest.approx(29.0)
    assert f.lons[0] == pytest.approx(-98.0)           # -180..180, not 0..360
    assert f.lons[-1] == pytest.approx(-97.0)
    assert float(f.grids["t2m"].mean()) == pytest.approx(290.0 - 273.15)
    assert [t[:10] for t in f.times] == ["2024-01-06"]
    assert f.times[0].endswith("12:00:00+00:00")


def test_fetch_tp_units(fake_cdsapi):
    f = fetch_era5("tp", (-98.0, 29.0, -97.0, 31.0),
                   "2024-01-06", "2024-01-06", stride_hours=24)
    assert float(f.grids["tp"].mean()) == pytest.approx(1.0)  # m -> mm


def test_fetch_provenance(fake_cdsapi):
    f = fetch_era5("msl", (-98.0, 29.0, -97.0, 31.0),
                   "2024-01-06", "2024-01-06", stride_hours=24)
    prov = f.provenance
    assert prov["cds_dataset"] == ERA5_DATASET
    assert prov["n_requests"] == 1
    assert len(prov["sha256"]) == 64
    assert prov["n_bytes"] > 0
    assert "retrieved_at" in prov
    assert prov["variables_requested"] == ["msl"]
    assert prov["units"] == {"msl": "hPa"}
    assert len(RECORDED_REQUESTS) == 1
    assert RECORDED_REQUESTS[0]["variable"] == ["mean_sea_level_pressure"]


def test_fetch_month_chunking_two_requests(fake_cdsapi):
    f = fetch_era5("t2m", (-98.0, 29.0, -97.0, 31.0),
                   "2024-01-31", "2024-02-01", stride_hours=24)
    assert len(f.times) == 2
    assert len(RECORDED_REQUESTS) == 2
    assert {r["month"] for r in RECORDED_REQUESTS} == {"01", "02"}
    assert f.provenance["time_chunks"] == 2


def test_fetch_antimeridian_two_windows(fake_cdsapi):
    f = fetch_era5("msl", (170.0, 20.0, -170.0, 30.0),
                   "2024-01-06", "2024-01-06", stride_hours=24)
    assert len(RECORDED_REQUESTS) == 2
    assert f.lons[0] == pytest.approx(-180.0)
    assert f.lons[-1] == pytest.approx(179.75)  # 180.0 aliases to -180.0
    assert np.all(np.diff(f.lons) > 0)          # strictly increasing, no seam dup
    assert f.grids["msl"].shape == (1, 41, 81)  # 41 + 41 - 1 seam dup


def test_fetch_future_end_raises(fake_cdsapi):
    future = (dt.date.today() + dt.timedelta(days=30)).isoformat()
    with pytest.raises(ValueError, match="in the future"):
        fetch_era5("wind", (-98.0, 29.0, -97.0, 31.0), "2024-01-06", future,
                   stride_hours=24)


def test_fetch_before_1940_raises(fake_cdsapi):
    with pytest.raises(ValueError, match="ERA5 starts 1940"):
        fetch_era5("wind", (-98.0, 29.0, -97.0, 31.0), "1939-12-31",
                   "1940-01-02", stride_hours=24)


def test_fetch_start_after_end_raises(fake_cdsapi):
    with pytest.raises(ValueError, match="after end"):
        fetch_era5("wind", (-98.0, 29.0, -97.0, 31.0), "2024-01-07",
                   "2024-01-06", stride_hours=24)


def test_parse_time_fallback(tmp_path):
    path = str(tmp_path / "era5.nc")
    _build_fake_era5_nc(path, {
        "area": [31.0, -98.0, 29.0, -97.0], "year": "2024", "month": "01",
        "day": ["06"], "time": ["12:00"], "variable": ["2m_temperature"],
    }, time_var="time")
    with open(path, "rb") as fh:
        payload = fh.read()
    lats, lons, times, grids = era5._parse_era5_bytes(payload, ["t2m"])
    assert times[0].startswith("2024-01-06T12:00:00")
    assert set(grids) == {"t2m"}


def test_parse_missing_variable_raises(tmp_path):
    path = str(tmp_path / "era5.nc")
    _build_fake_era5_nc(path, {
        "area": [31.0, -98.0, 29.0, -97.0], "year": "2024", "month": "01",
        "day": ["06"], "time": ["12:00"], "variable": ["2m_temperature"],
    })
    with open(path, "rb") as fh:
        payload = fh.read()
    with pytest.raises(ValueError, match="missing variable 'msl'"):
        era5._parse_era5_bytes(payload, ["msl"])


# ---------------------------------------------------------------------------
# Era5Field model
# ---------------------------------------------------------------------------


def test_field_synthetic_shapes():
    f = Era5Field.synthetic(variables=("wind", "msl"), nt=3)
    assert set(f.grids) == {"u10", "v10", "msl"}
    assert f.grids["msl"].shape == (3, 6, 8)
    assert f.base_variable == "wind"
    assert f.units["msl"] == "hPa"


def test_field_synthetic_unknown_variable():
    with pytest.raises(ValueError, match="unknown ERA5 variable"):
        Era5Field.synthetic(variables=("clouds",))


def test_field_values_wind_is_speed():
    f = Era5Field.synthetic(variables=("wind",), nt=2, seed=3)
    expected = np.ma.sqrt(f.grids["u10"] ** 2 + f.grids["v10"] ** 2)
    assert np.allclose(np.ma.filled(f.values, np.nan),
                       np.ma.filled(expected, np.nan))


def test_field_values_single_variable():
    f = Era5Field.synthetic(variables=("t2m",), nt=2)
    assert f.values is f.grids["t2m"]


def test_field_overlay_grids_excludes_base():
    f = Era5Field.synthetic(variables=("wind", "msl"), nt=2)
    assert sorted(f.overlay_grids) == ["msl"]
    g = f.overlay_grid("msl", 0)
    assert g.shape == (6, 8)


def test_field_overlay_grid_missing_raises():
    f = Era5Field.synthetic(variables=("wind",), nt=2)
    with pytest.raises(KeyError, match="no overlay 'msl'"):
        f.overlay_grid("msl", 0)


def test_field_spatial_mean():
    f = Era5Field.synthetic(variables=("wind", "t2m"), nt=2, seed=5)
    assert f.spatial_mean("wind", 0) == pytest.approx(float(f.wind_speed[0].mean()))
    assert f.spatial_mean("t2m", 1) == pytest.approx(float(f.grids["t2m"][1].mean()))


def test_field_select_time_and_bbox():
    f = Era5Field.synthetic(variables=("wind", "msl"), nt=4)
    one = f.select_time(2)
    assert one.grids["msl"].shape == (1, 6, 8)
    assert one.times == [f.times[2]]
    sub = f.select_bbox((-60.0, -30.0, 0.0, 0.0))
    assert sub.grids["u10"].shape == (4, 3, 4)
    with pytest.raises(ValueError, match="does not overlap"):
        f.select_bbox((100.0, 50.0, 110.0, 60.0))


def test_field_roundtrip_json(tmp_path):
    f = Era5Field.synthetic(variables=("wind", "msl", "tp"), nt=2)
    path = str(tmp_path / "era5.json")
    f.to_json(path)
    f2 = Era5Field.from_json(path)
    assert f2.variables == ["wind", "msl", "tp"]
    assert f2.base_variable == "wind"
    assert f2.grids["tp"].shape == (2, 6, 8)
    assert f2.times == f.times


def test_field_shape_mismatch_raises():
    f = Era5Field.synthetic(nt=2)
    bad = dict(f.grids)
    bad["msl"] = bad["msl"][:, :, :4]
    with pytest.raises(ValueError, match="does not match"):
        Era5Field(grids=bad, times=f.times, lats=f.lats, lons=f.lons,
                  variables=["wind", "msl"], base_variable="wind")


def test_field_base_variable_must_be_requested():
    f = Era5Field.synthetic(nt=2)
    with pytest.raises(ValueError, match="not in variables"):
        Era5Field(grids=f.grids, times=f.times, lats=f.lats, lons=f.lons,
                  variables=["wind", "msl"], base_variable="t2m")


def test_main_demo(capsys):
    era5.main_demo()
    out = capsys.readouterr().out
    assert "[era5]" in out and "mean wind=" in out


def test_cli_synthetic(tmp_path, monkeypatch):
    from currents import cli
    out = tmp_path / "demo"
    monkeypatch.setattr(sys, "argv",
                        ["survey-currents", "era5-synthetic",
                         "--variables", "wind,t2m", "--nt", "2",
                         "--out", str(out)])
    assert cli.main() == 0
    with open(str(out) + ".json") as fh:
        data = json.load(fh)
    assert data["variables"] == ["wind", "t2m"]
    assert len(data["times"]) == 2


def test_cli_fetch(tmp_path, monkeypatch, fake_cdsapi):
    from currents import cli
    out = tmp_path / "fetch"
    monkeypatch.setattr(sys, "argv",
                        ["survey-currents", "fetch-era5",
                         "--variables", "msl",
                         "--bbox=-98,29,-97,31",  # = form: raw - looks like a flag
                         "--start", "2024-01-06", "--end", "2024-01-06",
                         "--stride-hours", "24",
                         "--out", str(out)])
    assert cli.main() == 0
    with open(str(out) + ".json") as fh:
        data = json.load(fh)
    assert data["variables"] == ["msl"]
    assert data["provenance"]["n_requests"] == 1

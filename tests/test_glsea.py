"""Tests for the NOAA GLSEA module. Fully offline except the two
explicitly network-guarded tests at the bottom (skipped without network)."""

import datetime as dt
import hashlib
import io
import json
import socket
import sys
import types
import urllib.request

import numpy as np
import numpy.ma as ma
import pytest

from currents import glsea
from currents.glsea import (LAKE_COLUMNS, GlseaField, LakeSeries,
                            fetch_glsea_lake_averages, fetch_glsea_sst,
                            glsea_averages_url, glsea_sst_url,
                            parse_glsea_averages_csv, validate_glsea_bbox)


# ---------------------------------------------------------------------------
# date coercion
# ---------------------------------------------------------------------------

def test_coerce_date_from_date():
    assert glsea._coerce_date(dt.date(2026, 1, 5)) == dt.date(2026, 1, 5)


def test_coerce_date_from_datetime():
    assert glsea._coerce_date(dt.datetime(2026, 1, 5, 9, 30)) == dt.date(2026, 1, 5)


def test_coerce_date_from_iso_string():
    assert glsea._coerce_date("2026-01-05") == dt.date(2026, 1, 5)
    assert glsea._coerce_date("2026-01-05T12:00:00Z") == dt.date(2026, 1, 5)


def test_coerce_date_bad():
    with pytest.raises(ValueError, match="cannot parse"):
        glsea._coerce_date("not-a-date")


def test_erddap_iso_date_defaults_to_noon():
    # GLSEA daily timesteps are stamped 12:00 UTC
    assert glsea._erddap_iso("2016-01-01") == "2016-01-01T12:00:00Z"


def test_erddap_iso_naive_datetime_assumed_utc():
    assert glsea._erddap_iso(dt.datetime(2016, 1, 2, 3, 4, 5)) == \
        "2016-01-02T03:04:05Z"


# ---------------------------------------------------------------------------
# bbox validation
# ---------------------------------------------------------------------------

def test_validate_bbox_ok():
    assert validate_glsea_bbox((-88.0, 41.5, -86.0, 44.5)) == \
        (-88.0, 41.5, -86.0, 44.5)


def test_validate_bbox_west_of_floor():
    with pytest.raises(ValueError, match="west of the GLSEA grid's longitude floor"):
        validate_glsea_bbox((-93.0, 41.5, -86.0, 44.5))


def test_validate_bbox_east_of_edge():
    with pytest.raises(ValueError, match="east of the GLSEA grid's eastern edge"):
        validate_glsea_bbox((-88.0, 41.5, -75.0, 44.5))


def test_validate_bbox_lat_outside():
    with pytest.raises(ValueError, match="outside the GLSEA grid's latitude range"):
        validate_glsea_bbox((-88.0, 35.0, -86.0, 44.5))


def test_validate_bbox_lon_descending():
    with pytest.raises(ValueError, match="longitude constraint must be ascending"):
        validate_glsea_bbox((-86.0, 41.5, -88.0, 44.5))


def test_validate_bbox_lat_descending():
    with pytest.raises(ValueError, match="latitude constraint must be ascending"):
        validate_glsea_bbox((-88.0, 44.5, -86.0, 41.5))


def test_validate_bbox_wrong_arity():
    with pytest.raises(ValueError, match="needs \\(lon_min, lat_min, lon_max, lat_max\\)"):
        validate_glsea_bbox((-88.0, 41.5, -86.0))


def test_validate_bbox_exact_floor_is_allowed():
    # the floor itself is on-grid
    assert validate_glsea_bbox((glsea.GLSEA_LON_MIN, 41.5, -86.0, 44.5))[0] == \
        pytest.approx(glsea.GLSEA_LON_MIN)


# ---------------------------------------------------------------------------
# URL construction
# ---------------------------------------------------------------------------

def test_glsea_sst_url():
    url = glsea_sst_url((-88.0, 41.5, -86.0, 44.5), "2016-01-01", "2016-03-01")
    assert url == (
        "https://apps.glerl.noaa.gov/erddap/griddap/GLSEA_ACSPO_GCS.nc?"
        "sst[(2016-01-01T12:00:00Z):30:(2016-03-01T12:00:00Z)]"
        "[(41.5):1:(44.5)][(-88.0):1:(-86.0)]")


def test_glsea_sst_url_stride():
    url = glsea_sst_url((-88.0, 41.5, -86.0, 44.5), "2016-01-01", "2016-01-10",
                       stride_days=1)
    assert "):1:(" in url.split("?sst")[1].split("][")[0]
    assert url.endswith("][(-88.0):1:(-86.0)]")


def test_glsea_sst_url_start_after_end():
    with pytest.raises(ValueError, match="is after end"):
        glsea_sst_url((-88.0, 41.5, -86.0, 44.5), "2016-03-01", "2016-01-01")


def test_glsea_sst_url_bad_stride():
    with pytest.raises(ValueError, match="stride_days must be >= 1"):
        glsea_sst_url((-88.0, 41.5, -86.0, 44.5), "2016-01-01", "2016-01-10",
                      stride_days=0)


def test_glsea_sst_url_rejects_bad_bbox():
    with pytest.raises(ValueError, match="longitude floor"):
        glsea_sst_url((-95.0, 41.5, -86.0, 44.5), "2016-01-01", "2016-01-10")


def test_glsea_averages_url():
    assert glsea_averages_url("superior", 2020) == (
        "https://apps.glerl.noaa.gov/erddap/tabledap/glsea_avgtemps_3.csv"
        "?Year,Day,Sup&Year%3E=2020")


def test_glsea_averages_url_lake_mapping():
    assert "Day,Mich&" in glsea_averages_url("MICHIGAN", 2021)
    assert "Day,Ont&" in glsea_averages_url(" ontario ", 2021)


def test_glsea_averages_url_unknown_lake():
    with pytest.raises(ValueError, match="unknown lake"):
        glsea_averages_url("baikal", 2020)


def test_lake_columns_mapping():
    assert LAKE_COLUMNS == {"superior": "Sup", "michigan": "Mich",
                            "huron": "Huron", "erie": "Erie", "ontario": "Ont"}


# ---------------------------------------------------------------------------
# avgtemps CSV parsing (real fixture lines captured live 2026-09-26)
# ---------------------------------------------------------------------------

_AVG_FIXTURE = """Year,Day,Sup,Mich,Huron,Erie,Ont
,,,,,,
2026,001,3.13,3.91,2.93,1.95,3.78
2026,002,3.03,3.75,2.76,1.68,3.55
2026,003,2.93,3.61,2.64,1.5,3.5
2026,004,2.79,,2.45,1.25,3.44
2026,005,2.69,3.53,2.39,1.25,3.37
"""


def _only_col(csv_text, col):
    """Reduce a multi-column fixture to Year,Day,<col> like the engine's URL."""
    out = ["Year,Day,Value", ",,"]
    for row in csv_text.splitlines()[2:]:
        cells = row.split(",")
        idx = {"Sup": 2, "Mich": 3, "Huron": 4, "Erie": 5, "Ont": 6}[col]
        out.append(f"{cells[0]},{cells[1]},{cells[idx]}")
    return "\n".join(out) + "\n"


def test_parse_averages_csv_basic():
    s = parse_glsea_averages_csv(_only_col(_AVG_FIXTURE, "Sup"), "superior",
                                 dt.date(2026, 1, 1), dt.date(2026, 12, 31))
    assert isinstance(s, LakeSeries)
    assert s.lake == "superior"
    assert s.dates == ["2026-01-01", "2026-01-02", "2026-01-03",
                       "2026-01-04", "2026-01-05"]
    assert s.temps == [3.13, 3.03, 2.93, 2.79, 2.69]


def test_parse_averages_csv_skips_missing_temp():
    # row 004 has an empty Mich value -> skipped for michigan
    s = parse_glsea_averages_csv(_only_col(_AVG_FIXTURE, "Mich"), "michigan",
                                 dt.date(2026, 1, 1), dt.date(2026, 12, 31))
    assert s.temps == [3.91, 3.75, 3.61, 3.53]
    assert s.dates == ["2026-01-01", "2026-01-02", "2026-01-03", "2026-01-05"]


def test_parse_averages_csv_window_filter():
    s = parse_glsea_averages_csv(_only_col(_AVG_FIXTURE, "Erie"), "erie",
                                 dt.date(2026, 1, 2), dt.date(2026, 1, 3))
    assert s.dates == ["2026-01-02", "2026-01-03"]
    assert s.temps == [1.68, 1.5]


def test_parse_averages_csv_unknown_lake():
    with pytest.raises(ValueError, match="unknown lake"):
        parse_glsea_averages_csv(_only_col(_AVG_FIXTURE, "Sup"), "baikal",
                                 dt.date(2026, 1, 1), dt.date(2026, 1, 2))


def test_parse_averages_csv_empty_window():
    s = parse_glsea_averages_csv(_only_col(_AVG_FIXTURE, "Sup"), "superior",
                                 dt.date(2027, 1, 1), dt.date(2027, 1, 2))
    assert s.n == 0
    assert s.dates == [] and s.temps == []


def test_lakeseries_length_mismatch():
    with pytest.raises(ValueError, match="lengths differ"):
        LakeSeries(lake="superior", dates=["2026-01-01"],
                   temps=[3.1, 3.0])


def test_lakeseries_unknown_lake():
    with pytest.raises(ValueError, match="unknown lake"):
        LakeSeries(lake="baikal", dates=[], temps=[])


def test_lakeseries_json_roundtrip(tmp_path):
    s = parse_glsea_averages_csv(_only_col(_AVG_FIXTURE, "Ont"), "ontario",
                                 dt.date(2026, 1, 1), dt.date(2026, 1, 31))
    p = str(tmp_path / "s.json")
    s.to_json(p)
    s2 = LakeSeries.from_json(p)
    assert s2.lake == "ontario"
    assert s2.dates == s.dates and s2.temps == s.temps
    assert abs(s2.mean() - s.mean()) < 1e-12


def test_lakeseries_mean_empty():
    with pytest.raises(ValueError, match="empty series"):
        LakeSeries(lake="erie", dates=[], temps=[]).mean()


# ---------------------------------------------------------------------------
# GlseaField model
# ---------------------------------------------------------------------------

def test_synthetic_shape_and_mask():
    f = GlseaField.synthetic(nt=3, ny=5, nx=6, seed=3)
    assert f.sst.shape == (3, 5, 6)
    assert isinstance(f.sst, ma.MaskedArray)
    assert f.sst.mask.sum() > 0          # corner "land" cells masked
    assert len(f.times) == 3
    assert f.bounds == (f.lons[0], f.lats[0], f.lons[-1], f.lats[-1])
    assert f.provenance["synthetic"] is True
    assert list(f.lons) == sorted(f.lons) and list(f.lats) == sorted(f.lats)


def test_synthetic_deterministic():
    a = GlseaField.synthetic(seed=11)
    b = GlseaField.synthetic(seed=11)
    assert ma.allequal(a.sst, b.sst)


def test_glseafield_shape_mismatch():
    f = GlseaField.synthetic(nt=2)
    with pytest.raises(ValueError, match="does not match"):
        GlseaField(sst=f.sst[:1], times=f.times, lats=f.lats, lons=f.lons)


def test_spatial_mean_ignores_mask():
    f = GlseaField.synthetic(nt=1, seed=5)
    mean = f.spatial_mean(0)
    assert mean == pytest.approx(ma.mean(ma.masked_invalid(f.sst[0])))


def test_select_time():
    f = GlseaField.synthetic(nt=4, seed=5)
    one = f.select_time(2)
    assert one.sst.shape == (1, 6, 8)
    assert one.times == [f.times[2]]
    assert ma.allequal(one.sst[0], f.sst[2])


def test_select_bbox():
    f = GlseaField.synthetic(nt=2, lats=(41.5, 46.5), lons=(-92.4, -84.5),
                             ny=6, nx=8, seed=5)
    sub = f.select_bbox((-90.0, 42.0, -87.0, 45.0))
    assert sub.sst.shape[0] == 2
    assert sub.lons[0] >= -90.0 and sub.lons[-1] <= -87.0
    assert sub.lats[0] >= 42.0 and sub.lats[-1] <= 45.0


def test_select_bbox_no_overlap():
    f = GlseaField.synthetic(nt=1, seed=5)
    with pytest.raises(ValueError, match="does not overlap"):
        f.select_bbox((-80.0, 30.0, -79.0, 31.0))


def test_glseafield_json_roundtrip(tmp_path):
    f = GlseaField.synthetic(nt=2, seed=9)
    p = str(tmp_path / "g.json")
    f.to_json(p)
    g = GlseaField.from_json(p)
    assert g.times == f.times
    assert g.source == f.source
    assert ma.allequal(ma.masked_invalid(g.sst), f.sst)
    assert g.bounds == pytest.approx(f.bounds)


def test_glseafield_from_dict_missing_keys():
    with pytest.raises(ValueError, match="missing keys"):
        GlseaField.from_dict({"sst": [], "times": []})


# ---------------------------------------------------------------------------
# NetCDF parsing with a fake netCDF4 (no binary dependency)
# ---------------------------------------------------------------------------

class _FakeVar:
    def __init__(self, data, dims=(), units=None, fill=None):
        self._data = np.asarray(data)
        self.dimensions = tuple(dims)
        self._units = units
        self._fill = fill

    def ncattrs(self):
        attrs = []
        if self._fill is not None:
            attrs.append("_FillValue")
        if self._units is not None:
            attrs.append("units")
        return attrs

    def getncattr(self, name):
        if name == "_FillValue":
            return self._fill
        if name == "units":
            return self._units
        raise AttributeError(name)

    def __getitem__(self, key):
        return self._data[key]


class _FakeDS:
    def __init__(self, variables):
        self.variables = variables

    def close(self):
        pass


def _fake_glsea_ds(nt=2, ny=3, nx=4):
    lat = np.linspace(42.0, 44.0, ny)
    lon = np.linspace(-90.0, -88.0, nx)
    sst = np.full((nt, ny, nx), 17.5)
    sst[:, 0, 0] = -99999.0  # fill sentinel -> masked
    return _FakeDS({
        "sst": _FakeVar(sst, ("time", "latitude", "longitude"), fill=-99999.0),
        "latitude": _FakeVar(lat, ("latitude",)),
        "longitude": _FakeVar(lon, ("longitude",)),
        "time": _FakeVar(np.array([1451606400, 1451692800]),
                         ("time",), units="seconds since 1970-01-01T00:00:00Z"),
    })


def test_parse_glsea_dataset_masks_fill(monkeypatch):
    # fake netCDF4.num2date so no binary dependency is needed
    fake_nc = types.ModuleType("netCDF4")
    fake_nc.num2date = lambda v, u, calendar="standard": \
        dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc) + \
        dt.timedelta(seconds=float(v))
    monkeypatch.setitem(sys.modules, "netCDF4", fake_nc)
    lats, lons, times, sst = glsea._parse_glsea_dataset(_fake_glsea_ds())
    assert list(lats) == pytest.approx([42.0, 43.0, 44.0])
    assert list(lons) == pytest.approx([-90.0, -89.33333333333333, -88.66666666666667,
                                        -88.0])
    assert times == ["2016-01-01T00:00:00+00:00", "2016-01-02T00:00:00+00:00"]
    assert sst.shape == (2, 3, 4)
    assert bool(sst.mask[0, 0, 0]) is True     # fill sentinel masked
    assert bool(sst.mask[0, 1, 1]) is False
    assert sst[0, 1, 1] == pytest.approx(17.5)


def test_parse_glsea_dataset_bad_rank():
    ds = _FakeDS({"sst": _FakeVar(np.full((3, 4), 1.0), ("lat", "lon")),
                  "latitude": _FakeVar([1, 2, 3], ("latitude",)),
                  "longitude": _FakeVar([1, 2, 3, 4], ("longitude",)),
                  "time": _FakeVar([0], ("time",))})
    with pytest.raises(ValueError, match="expected sst\\(time, lat, lon\\)"):
        glsea._parse_glsea_dataset(ds)


def test_require_netcdf4_error_names_pip(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "netCDF4":
            raise ImportError("No module named 'netCDF4'")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="pip install netCDF4"):
        glsea._require_netcdf4()


# ---------------------------------------------------------------------------
# fetch_* with mocked network (urlopen stubbed; parse stubbed)
# ---------------------------------------------------------------------------

class _FakeResp:
    def __init__(self, payload: bytes):
        self._buf = io.BytesIO(payload)

    def read(self, n=-1):
        return self._buf.read(n)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_fetch_glsea_sst_mocked(monkeypatch):
    seen = {}
    payload = b"fake-netcdf-bytes"

    def fake_urlopen(url, timeout=300):
        seen["url"] = url
        return _FakeResp(payload)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    synth = GlseaField.synthetic(nt=2, seed=1)
    monkeypatch.setattr(glsea, "parse_glsea_bytes",
                        lambda b: (synth.lats, synth.lons, synth.times, synth.sst))
    f = fetch_glsea_sst((-88.0, 41.5, -86.0, 44.5), "2016-01-01", "2016-03-01")
    assert f.source == "noaa-glsea/GLSEA_ACSPO_GCS"
    assert len(f.times) == 2
    assert seen["url"] == glsea_sst_url((-88.0, 41.5, -86.0, 44.5),
                                       "2016-01-01", "2016-03-01")
    prov = f.provenance
    assert prov["url"] == seen["url"]
    assert prov["sha256"] == hashlib.sha256(payload).hexdigest()
    assert prov["retrieved_at"]  # ISO timestamp present
    assert prov["units"] == "degree_C"
    assert prov["stride_days"] == 30
    # bbox validation still applies before any download
    with pytest.raises(ValueError, match="longitude floor"):
        fetch_glsea_sst((-95.0, 41.5, -86.0, 44.5), "2016-01-01", "2016-03-01")


def test_fetch_glsea_averages_mocked(monkeypatch):
    seen = {}
    payload = _only_col(_AVG_FIXTURE, "Huron").encode()

    def fake_urlopen(url, timeout=300):
        seen["url"] = url
        return _FakeResp(payload)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    s = fetch_glsea_lake_averages("huron", "2026-01-01", "2026-01-31")
    assert s.lake == "huron"
    assert s.n == 5
    assert s.temps[0] == pytest.approx(2.93)
    assert seen["url"] == glsea_averages_url("huron", 2026)
    prov = s.provenance
    assert prov["url"] == seen["url"]
    assert prov["sha256"] == hashlib.sha256(payload).hexdigest()
    assert prov["retrieved_at"]
    assert prov["lake_column"] == "Huron"
    assert prov["rows_returned"] == 5


def test_fetch_glsea_averages_unknown_lake_no_network(monkeypatch):
    # lake validation happens before any download
    def boom(url, timeout=300):
        raise AssertionError("must not hit the network")

    monkeypatch.setattr(urllib.request, "urlopen", boom)
    with pytest.raises(ValueError, match="unknown lake"):
        fetch_glsea_lake_averages("baikal", "2026-01-01", "2026-01-31")


def test_fetch_glsea_averages_start_after_end():
    with pytest.raises(ValueError, match="is after end"):
        fetch_glsea_lake_averages("erie", "2026-02-01", "2026-01-01")


# ---------------------------------------------------------------------------
# package exports + CLI wiring
# ---------------------------------------------------------------------------

def test_package_exports():
    import currents
    for name in ("GlseaField", "LakeSeries", "fetch_glsea_sst",
                 "fetch_glsea_lake_averages", "LAKE_COLUMNS",
                 "validate_glsea_bbox", "glsea_sst_url", "glsea_averages_url"):
        assert name in currents.__all__
        assert hasattr(currents, name)
    assert currents.__version__ == "0.6.0"


def test_cli_has_glsea_subcommands():
    from currents.cli import build_parser
    p = build_parser()
    subs = {a.dest for a in p._subparsers._group_actions[0]._choices_actions}
    assert {"fetch-glsea-sst", "fetch-glsea-averages", "glsea-synthetic"} <= subs


def test_cli_glsea_synthetic_offline(tmp_path, monkeypatch, capsys):
    from currents.cli import main
    monkeypatch.chdir(tmp_path)
    rc = main(["glsea-synthetic", "--nt", "3", "--out", "g"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "synthetic GLSEA SST field (3 steps)" in out
    data = json.load(open("g.json"))
    assert len(data["times"]) == 3
    assert data["source"] == "synthetic"


# ---------------------------------------------------------------------------
# live-network tests (skipped when offline)
# ---------------------------------------------------------------------------

def _has_network() -> bool:
    try:
        socket.create_connection(("apps.glerl.noaa.gov", 443), timeout=5).close()
        return True
    except OSError:
        return False


_NEEDS_NET = pytest.mark.skipif(not _has_network(),
                                 reason="no network to apps.glerl.noaa.gov")

try:
    import netCDF4  # noqa: F401
    _HAS_NETCDF4 = True
except ImportError:
    _HAS_NETCDF4 = False

_NEEDS_NC = pytest.mark.skipif(not _HAS_NETCDF4,
                               reason="netCDF4 not installed")


@_NEEDS_NET
@_NEEDS_NC
def test_live_fetch_glsea_sst(tmp_path):
    # tiny subset: 1-day window, small bbox over western Lake Superior
    f = fetch_glsea_sst((-91.0, 46.5, -90.0, 47.0), "2026-01-01", "2026-01-02",
                        stride_days=1)
    assert f.sst.shape[0] >= 1
    assert f.sst.shape[1:] == (f.lats.shape[0], f.lons.shape[0])
    assert f.lons[0] >= glsea.GLSEA_LON_MIN
    assert "sha256" in f.provenance


@_NEEDS_NET
def test_live_fetch_glsea_averages():
    s = fetch_glsea_lake_averages("superior", "2026-01-01", "2026-01-07")
    assert s.n >= 5  # allow a couple of missing days
    assert all(-5.0 < t < 35.0 for t in s.temps)

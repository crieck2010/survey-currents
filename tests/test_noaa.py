"""Tests for the NOAA OFS module. No network: S3 and xarray are mocked."""

import io
import json
import sys
import types
import urllib.request

import numpy as np
import pytest

from currents import noaa_ofs
from currents.models import CurrentField


# ---------------------------------------------------------------------------
# path / key helpers (pure, offline)
# ---------------------------------------------------------------------------

def test_archive_prefix():
    assert noaa_ofs.archive_prefix("2026-09-16") == "OFS/netcdf/202609/"


def test_recent_prefix():
    assert noaa_ofs.recent_prefix("2026-09-16") == "OFS.20260916"


def test_filter_keys():
    keys = [
        "OFS/netcdf/202609/nos.lmhofs.fields.f006.20260916.t00z.nc",
        "OFS/netcdf/202609/nos.lmhofs.fields.f012.20260916.t00z.nc",
        "OFS/netcdf/202609/nos.glofs.fields.f006.20260916.t00z.nc",
        "OFS/netcdf/202609/nos.lmhofs.fields.f006.20260916.t06z.nc",
        "OFS/netcdf/202609/readme.txt",
    ]
    out = noaa_ofs.filter_keys(keys, "LMHOFS", "2026-09-16", "00", 6)
    assert out == ["OFS/netcdf/202609/nos.lmhofs.fields.f006.20260916.t00z.nc"]
    # without hour: both matching hours
    out2 = noaa_ofs.filter_keys(keys, "LMHOFS", "2026-09-16", "00")
    assert len(out2) == 2


def test_select_files_hour_outside_horizon():
    with pytest.raises(ValueError, match="outside.*horizon"):
        noaa_ofs.select_files(["k.nc"], "GLOFS", "2026-09-16", "00", [61])


def test_select_files_missing():
    with pytest.raises(FileNotFoundError, match="no NetCDF"):
        noaa_ofs.select_files(["other.nc"], "GLOFS", "2026-09-16", "00", [6])


def test_select_files_ok():
    keys = ["OFS/netcdf/202609/nos.glofs.fields.f006.20260916.t00z.nc"]
    out = noaa_ofs.select_files(keys, "GLOFS", "2026-09-16", "00", [6])
    assert out == {6: keys[0]}


def test_list_field_files_bad_cycle():
    with pytest.raises(ValueError, match="not in"):
        noaa_ofs.list_field_files("GLOFS", "2026-09-16", cycle="03")


def test_list_field_files_unknown_ofs():
    with pytest.raises(KeyError):
        noaa_ofs.list_field_files("NOPE", "2026-09-16")


# ---------------------------------------------------------------------------
# anonymous S3 over stdlib (urlopen mocked)
# ---------------------------------------------------------------------------

_LIST_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <IsTruncated>false</IsTruncated>
  <Contents><Key>OFS/netcdf/202609/a.nc</Key></Contents>
  <Contents><Key>OFS/netcdf/202609/b.nc</Key></Contents>
</ListBucketResult>"""


class _FakeResp:
    def __init__(self, payload: bytes):
        self._buf = io.BytesIO(payload)

    def read(self, n=-1):
        return self._buf.read(n)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _no_boto3(monkeypatch):
    monkeypatch.setattr(noaa_ofs, "_boto3_client", lambda: None)


def test_s3_list_stdlib(monkeypatch):
    _no_boto3(monkeypatch)
    monkeypatch.setattr(urllib.request, "urlopen",
                        lambda url, timeout=60: _FakeResp(_LIST_XML.encode()))
    keys = noaa_ofs.s3_list("noaa-nos-ofs-pds", "OFS/netcdf/202609/")
    assert keys == ["OFS/netcdf/202609/a.nc", "OFS/netcdf/202609/b.nc"]


def test_s3_download_stdlib(monkeypatch, tmp_path):
    _no_boto3(monkeypatch)
    payload = b"\x89HDF fake-bytes"
    monkeypatch.setattr(urllib.request, "urlopen",
                        lambda url, timeout=300: _FakeResp(payload))
    dest = str(tmp_path / "sub" / "f.nc")
    out = noaa_ofs.s3_download("bkt", "some/key.nc", dest)
    assert out == dest
    assert open(dest, "rb").read() == payload


def test_s3_uses_boto3_when_available(monkeypatch, tmp_path):
    class FakeClient:
        def download_file(self, bucket, key, dest):
            open(dest, "wb").write(b"via-boto3")
    monkeypatch.setattr(noaa_ofs, "_boto3_client", lambda: FakeClient())
    dest = str(tmp_path / "f.nc")
    noaa_ofs.s3_download("bkt", "k", dest)
    assert open(dest, "rb").read() == b"via-boto3"


# ---------------------------------------------------------------------------
# NetCDF parsing with a fake xarray (no netCDF4 needed)
# ---------------------------------------------------------------------------

class _FakeDA:
    def __init__(self, values, dims, name=""):
        self._values = np.asarray(values)  # keep dtype (incl. datetime64)
        self.dims = tuple(dims)
        self.name = name

    @property
    def values(self):
        return self._values

    def astype(self, dt):
        return _FakeDA(self._values.astype(dt), self.dims, self.name)

    def isel(self, sel):
        arr = self._values
        dims = list(self.dims)
        for dim, idx in sel.items():
            ax = dims.index(dim)
            arr = np.take(arr, [idx], axis=ax)
            # keep dim (xarray isel drops it, but our surface() squeezes anyway)
        return _FakeDA(arr, self.dims, self.name)


class _FakeDS:
    def __init__(self, variables):
        self.variables = variables
        self.coords = {}
        self.attrs = {}

    def __getitem__(self, name):
        return self.variables[name]

    def close(self):
        pass


def _fake_xr(variables):
    mod = types.ModuleType("xarray")
    mod.open_dataset = lambda path: _FakeDS(variables)
    return mod


def _regular_vars(descending_lat=False):
    lat = np.array([46.0, 44.0, 42.0]) if descending_lat else np.array([42.0, 44.0, 46.0])
    lon = np.array([-92.0, -90.0, -88.0])
    shape = (1, 3, 3)
    return {
        "u": _FakeDA(np.full(shape, 0.3), ("time", "lat", "lon"), "u"),
        "v": _FakeDA(np.full(shape, -0.1), ("time", "lat", "lon"), "v"),
        "temp": _FakeDA(np.full(shape, 17.5), ("time", "lat", "lon"), "temp"),
        "lat": _FakeDA(lat, ("lat",), "lat"),
        "lon": _FakeDA(lon, ("lon",), "lon"),
        "time": _FakeDA(np.array(["2026-09-16T00:00:00"],
                                 dtype="datetime64[s]"), ("time",), "time"),
    }


def test_parse_regular_grid(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "xarray", _fake_xr(_regular_vars()))
    f = noaa_ofs.parse_ofs_netcdf(str(tmp_path / "x.nc"))
    assert isinstance(f, CurrentField)
    assert f.u.shape == (1, 3, 3)
    assert f.u[0, 0, 0] == pytest.approx(0.3)
    assert f.temperature[0, 0, 0] == pytest.approx(17.5)
    # lats stored increasing even though file order may vary
    assert list(f.lats) == sorted(f.lats)


def test_parse_descending_lat_reordered(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "xarray", _fake_xr(_regular_vars(descending_lat=True)))
    f = noaa_ofs.parse_ofs_netcdf(str(tmp_path / "x.nc"))
    assert list(f.lats) == [42.0, 44.0, 46.0]


def test_parse_unstructured_raises(monkeypatch, tmp_path):
    # FVCOM-style: u/v on a node dimension, not on lat/lon
    node = np.arange(5)
    variables = {
        "u": _FakeDA(np.full((1, 5), 0.2), ("time", "node"), "u"),
        "v": _FakeDA(np.full((1, 5), 0.1), ("time", "node"), "v"),
        "lon": _FakeDA(np.linspace(-92, -88, 5), ("node",), "lon"),
        "lat": _FakeDA(np.linspace(42, 46, 5), ("node",), "lat"),
    }
    monkeypatch.setitem(sys.modules, "xarray", _fake_xr(variables))
    with pytest.raises(ValueError, match="unstructured"):
        noaa_ofs.parse_ofs_netcdf(str(tmp_path / "x.nc"))


def test_parse_missing_variables(monkeypatch, tmp_path):
    variables = {
        "lon": _FakeDA(np.array([-92.0]), ("lon",), "lon"),
        "lat": _FakeDA(np.array([42.0]), ("lat",), "lat"),
    }
    monkeypatch.setitem(sys.modules, "xarray", _fake_xr(variables))
    with pytest.raises(KeyError):
        noaa_ofs.parse_ofs_netcdf(str(tmp_path / "x.nc"))


def test_parse_needs_xarray(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "xarray", None)
    # force the import machinery to fail
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "xarray":
            raise ImportError("no xarray")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="survey-currents\\[noaa\\]"):
        noaa_ofs.parse_ofs_netcdf(str(tmp_path / "x.nc"))


# ---------------------------------------------------------------------------
# stacking + end-to-end fetch (mocked)
# ---------------------------------------------------------------------------

def test_stack_steps():
    a = CurrentField.synthetic(nt=1, seed=1)
    b = CurrentField.synthetic(nt=1, seed=2)
    s = noaa_ofs.stack_steps([a, b])
    assert s.u.shape == (2, 6, 8)
    assert len(s.times) == 2
    assert len(s.provenance["steps"]) == 2


def test_stack_steps_grid_mismatch():
    a = CurrentField.synthetic(nt=1, seed=1)
    b = CurrentField.synthetic(nt=1, ny=9, seed=2)
    with pytest.raises(ValueError, match="Grids differ|grids differ"):
        noaa_ofs.stack_steps([a, b])


def test_stack_steps_empty():
    with pytest.raises(ValueError, match="no steps"):
        noaa_ofs.stack_steps([])


def test_fetch_currents_mocked(monkeypatch, tmp_path):
    keys = ["OFS/netcdf/202609/nos.glofs.fields.f000.20260916.t00z.nc",
            "OFS/netcdf/202609/nos.glofs.fields.f001.20260916.t00z.nc"]
    monkeypatch.setattr(noaa_ofs, "list_field_files", lambda *a, **k: keys)
    dl_paths = []
    real_select = noaa_ofs.select_files

    def fake_download(bucket, key, dest):
        dl_paths.append(dest)
        open(dest, "wb").write(b"nc")
        return dest

    monkeypatch.setattr(noaa_ofs, "s3_download", fake_download)
    monkeypatch.setitem(sys.modules, "xarray", _fake_xr(_regular_vars()))
    field = noaa_ofs.fetch_currents("GLOFS", "2026-09-16", cycle="00",
                                    hours=[0, 1], work_dir=str(tmp_path))
    assert field.u.shape == (2, 3, 3)
    assert field.forecast_hours == [0, 1]
    assert field.source == "noaa-ofs/GLOFS"
    assert len(dl_paths) == 2
    # provenance sidecars were written for both downloads
    for p in dl_paths:
        assert (p + ".provenance.json")
        import os
        assert os.path.exists(p + ".provenance.json")

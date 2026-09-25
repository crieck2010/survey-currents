"""Tests for the CMEMS module. No network: the toolbox is faked."""

import os
import sys
import types

import numpy as np
import pytest

from currents import cmems
from currents.models import CurrentField


def test_presets():
    p = cmems.get_preset("global-physics-daily")
    assert p.dataset_id == "cmems_mod_glo_phy_anfc_0.083deg_P1D-m"
    assert p.variables == ("uo", "vo", "thetao")
    with pytest.raises(KeyError, match="unknown CMEMS preset"):
        cmems.get_preset("nope")


def test_require_toolbox_missing(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "copernicusmarine":
            raise ImportError("no toolbox")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="survey-currents\\[cmems\\]"):
        cmems.require_toolbox()


def test_require_toolbox_no_credentials(monkeypatch, tmp_path):
    fake = types.ModuleType("copernicusmarine")
    monkeypatch.setitem(sys.modules, "copernicusmarine", fake)
    monkeypatch.setenv("HOME", str(tmp_path))  # no cred file there
    monkeypatch.delenv("COPERNICUSMARINE_SERVICE_USERNAME", raising=False)
    with pytest.raises(RuntimeError, match="CMEMS credentials not found"):
        cmems.require_toolbox()


def test_require_toolbox_env_credentials(monkeypatch, tmp_path):
    fake = types.ModuleType("copernicusmarine")
    monkeypatch.setitem(sys.modules, "copernicusmarine", fake)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("COPERNICUSMARINE_SERVICE_USERNAME", "user")
    assert cmems.require_toolbox() is fake


def test_subset_cmems_calls_toolbox(monkeypatch, tmp_path):
    calls = {}

    class FakeCM:
        @staticmethod
        def subset(**kwargs):
            calls.update(kwargs)
            out = os.path.join(kwargs["output_directory"],
                               kwargs["output_filename"])
            open(out, "wb").write(b"nc")

    monkeypatch.setattr(cmems, "require_toolbox", lambda: FakeCM)
    out = cmems.subset_cmems(
        "global-physics-daily",
        bbox=(-92.0, 41.0, -84.0, 47.0),
        start="2026-09-16T00:00:00", end="2026-09-17T00:00:00",
        work_dir=str(tmp_path))
    assert out.endswith(".nc")
    assert os.path.exists(out)
    assert calls["dataset_id"] == "cmems_mod_glo_phy_anfc_0.083deg_P1D-m"
    assert calls["variables"] == ["uo", "vo", "thetao"]
    assert calls["minimum_longitude"] == -92.0
    assert calls["maximum_latitude"] == 47.0
    # provenance sidecar written
    assert os.path.exists(out + ".provenance.json")


class _FakeDA:
    def __init__(self, values, dims, name=""):
        self._values = np.asarray(values, dtype=float)
        self.dims = tuple(dims)
        self.name = name

    @property
    def values(self):
        return self._values

    def astype(self, dt):
        return _FakeDA(self._values.astype(dt), self.dims, self.name)

    def isel(self, sel):
        arr = self._values
        for dim, idx in sel.items():
            arr = np.take(arr, [idx], axis=self.dims.index(dim))
        return _FakeDA(arr, self.dims, self.name)


class _FakeDS:
    def __init__(self, variables):
        self.variables = variables
        self.attrs = {"dataset_id": "cmems_mod_glo_phy_anfc_0.083deg_P1D-m"}

    def __getitem__(self, name):
        return self.variables[name]

    def __contains__(self, name):
        return name in self.variables

    def close(self):
        pass


def _fake_xr():
    mod = types.ModuleType("xarray")
    lat = np.array([41.0, 42.0])
    lon = np.array([-92.0, -91.0])
    mod.open_dataset = lambda path: _FakeDS({
        "uo": _FakeDA(np.full((2, 1, 2, 2), 0.2), ("time", "depth", "latitude", "longitude"), "uo"),
        "vo": _FakeDA(np.full((2, 1, 2, 2), 0.1), ("time", "depth", "latitude", "longitude"), "vo"),
        "thetao": _FakeDA(np.full((2, 1, 2, 2), 19.0), ("time", "depth", "latitude", "longitude"), "thetao"),
        "latitude": _FakeDA(lat, ("latitude",), "latitude"),
        "longitude": _FakeDA(lon, ("longitude",), "longitude"),
        "time": _FakeDA(np.array(["2026-09-16", "2026-09-17"], dtype="datetime64[D]"), ("time",), "time"),
        "depth": _FakeDA(np.array([0.5]), ("depth",), "depth"),
    })
    return mod


def test_parse_cmems_netcdf(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "xarray", _fake_xr())
    f = cmems.parse_cmems_netcdf(str(tmp_path / "x.nc"))
    assert isinstance(f, CurrentField)
    assert f.u.shape == (2, 2, 2)
    assert f.temperature[0, 0, 0] == pytest.approx(19.0)
    assert f.source.startswith("cmems:")
    assert len(f.times) == 2

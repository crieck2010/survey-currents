"""Tests for COG export and the interop helpers (offline)."""

import sys

import numpy as np
import pytest

from currents import convert
from currents.interop import (CurrentsPassProvider, _bbox_of_geometry,
                              align_to_thermal_zone)
from currents.models import CurrentField


# ---------------------------------------------------------------------------
# convert
# ---------------------------------------------------------------------------

def test_band_names_contract():
    assert convert.BAND_NAMES == ("speed", "u", "v", "temperature")


def test_timestep_stamp():
    assert convert._timestep_stamp("2026-09-16T06:00:00+00:00") == "20260916T0600"
    assert convert._timestep_stamp("not-a-date")  # falls back, never raises


def test_export_cogs_needs_rasterio(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "rasterio":
            raise ImportError("no rasterio")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    f = CurrentField.synthetic(nt=1, seed=1)
    with pytest.raises(ImportError, match="survey-currents\\[raster\\]"):
        convert.export_cogs(f, "/tmp/nowhere")


rasterio = None  # imported lazily inside rasterio-only tests


def _require_rasterio():
    return pytest.importorskip("rasterio", reason="rasterio not installed")


def test_export_cogs_round_trip(tmp_path):
    rio = _require_rasterio()
    f = CurrentField.synthetic(nt=2, seed=1)
    paths = convert.export_cogs(f, str(tmp_path), prefix="t")
    assert len(paths) == 2
    with rio.open(paths[0]) as ds:
        assert ds.count == 4
        assert ds.descriptions == ("speed", "u", "v", "temperature")
        assert ds.crs.to_string() == "EPSG:4326"
        speed = ds.read(1)
    # latitude axis flipped to north-up: compare against flipped field
    np.testing.assert_allclose(speed, f.speed()[0][::-1, :], rtol=1e-5)
    with rio.open(paths[1]) as ds:
        temp = ds.read(4)
    np.testing.assert_allclose(temp, f.temperature[1][::-1, :], rtol=1e-5)


def test_export_cogs_no_temperature(tmp_path):
    rio = _require_rasterio()
    f = CurrentField.synthetic(nt=1, seed=1)
    f.temperature = None
    paths = convert.export_cogs(f, str(tmp_path), prefix="t")
    with rio.open(paths[0]) as ds:
        band4 = ds.read(4)
    assert np.all(np.isnan(band4))


# ---------------------------------------------------------------------------
# interop geometry helper
# ---------------------------------------------------------------------------

def test_bbox_of_geometry_dict():
    g = {"bbox": [-92.0, 41.0, -84.0, 47.0]}
    assert _bbox_of_geometry(g) == [-92.0, 41.0, -84.0, 47.0]


def test_bbox_of_geojson_polygon():
    g = {"type": "Polygon", "coordinates": [[[-92, 41], [-84, 41],
                                             [-84, 47], [-92, 47],
                                             [-92, 41]]]}
    assert _bbox_of_geometry(g) == [-92.0, 41.0, -84.0, 47.0]


def test_bbox_of_bad_geometry():
    with pytest.raises(ValueError):
        _bbox_of_geometry({"type": "Point"})


# ---------------------------------------------------------------------------
# CurrentsPassProvider
# ---------------------------------------------------------------------------

def _provider():
    f = CurrentField.synthetic(nt=3, seed=1, step_hours=6)
    return CurrentsPassProvider(f, "lake-michigan",
                                {"bbox": [-92.5, 41.5, -84.5, 46.5]})


def test_provider_list_passes():
    p = _provider()
    passes = p.list_passes()
    assert len(passes) == 3
    assert passes[0].pass_id.startswith("currents-lake-michigan-f000-")
    assert passes[0].date == "2026-09-16"
    assert [x.date for x in passes] == sorted(x.date for x in passes)


def test_provider_list_passes_since():
    p = _provider()
    passes = p.list_passes(since="2026-09-16")
    # step_hours=6 keeps everything on 2026-09-16; since is exclusive
    assert passes == []


def test_provider_metrics():
    p = _provider()
    passes = p.list_passes()
    m = p.metrics(passes[1].pass_id)
    assert set(m) == {"speed_mean", "u_mean", "v_mean", "temp_mean"}
    assert m["speed_mean"] == pytest.approx(p.field.zonal_mean(1, "speed"))


def test_provider_metrics_unknown():
    p = _provider()
    with pytest.raises(KeyError):
        p.metrics("currents-lake-michigan-f999-2026-09-16")


def test_provider_metrics_no_temperature():
    f = CurrentField.synthetic(nt=2, seed=1)
    f.temperature = None
    p = CurrentsPassProvider(f, "s", {"bbox": [-92.5, 41.5, -84.5, 46.5]})
    m = p.metrics(p.list_passes()[0].pass_id)
    assert "temp_mean" not in m


# ---------------------------------------------------------------------------
# thermal zone alignment
# ---------------------------------------------------------------------------

def test_align_to_thermal_zone():
    f = CurrentField.synthetic(nt=4, seed=1)
    zone = {"id": "erie-west",
            "geometry": {"bbox": [-92.0, 42.0, -86.0, 45.0]}}
    out = align_to_thermal_zone(f, zone)
    assert out["zone_id"] == "erie-west"
    assert len(out["timesteps"]) == 4
    assert len(out["temp_mean_series"]) == 4
    assert out["timesteps"][0]["temp_mean_c"] == pytest.approx(
        out["temp_mean_series"][0])


def test_align_to_thermal_zone_window():
    f = CurrentField.synthetic(nt=4, seed=1)
    zone = {"id": "z", "geometry": {"bbox": [-92.0, 42.0, -86.0, 45.0]}}
    out = align_to_thermal_zone(f, zone, start="2026-09-17")
    assert out["timesteps"] == []
    assert out["temp_mean_series"] == []

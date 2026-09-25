"""Tests for the canonical model and the OFS registry (offline)."""

import json

import numpy as np
import pytest

from currents.models import (OFS_REGISTRY, CurrentField, get_ofs_model,
                             list_ofs_models)


def test_registry_has_expected_models():
    assert get_ofs_model("LMHOFS").horizon_hours == 120
    assert get_ofs_model("glofs").horizon_hours == 60  # case-insensitive
    assert get_ofs_model("CBOFS").region == "Chesapeake Bay"
    assert len(OFS_REGISTRY) >= 10


def test_registry_unknown_code():
    with pytest.raises(KeyError, match="unknown OFS code"):
        get_ofs_model("NOPE")


def test_list_ofs_models_sorted():
    codes = [m.code for m in list_ofs_models()]
    assert codes == sorted(codes)


def test_synthetic_shape_and_units():
    f = CurrentField.synthetic(nt=3, ny=4, nx=5, seed=1)
    assert f.u.shape == (3, 4, 5)
    assert f.v.shape == (3, 4, 5)
    assert f.temperature.shape == (3, 4, 5)
    assert len(f.times) == 3
    assert f.lats[0] < f.lats[-1] and f.lons[0] < f.lons[-1]
    assert f.bounds == (pytest.approx(-92.5), pytest.approx(41.5),
                        pytest.approx(-84.5), pytest.approx(46.5))


def test_synthetic_deterministic():
    a = CurrentField.synthetic(seed=3)
    b = CurrentField.synthetic(seed=3)
    assert np.array_equal(a.u, b.u)
    c = CurrentField.synthetic(seed=4)
    assert not np.array_equal(a.u, c.u)


def test_speed_and_direction():
    f = CurrentField.synthetic(nt=1, seed=1)
    speed = f.speed()
    assert speed.shape == f.u.shape
    assert np.all(speed >= 0)
    np.testing.assert_allclose(speed[0], np.hypot(f.u[0], f.v[0]))
    d = f.direction_deg()
    assert np.all((d >= 0) & (d < 360))


def test_shape_mismatch_raises():
    f = CurrentField.synthetic(nt=2, ny=3, nx=4)
    with pytest.raises(ValueError, match="does not match"):
        CurrentField(u=np.zeros((2, 3, 5)), v=f.v, temperature=None,
                     times=f.times, lats=f.lats, lons=f.lons)
    with pytest.raises(ValueError, match="does not match"):
        CurrentField(u=f.u, v=f.v, temperature=np.zeros((2, 3, 3)),
                     times=f.times, lats=f.lats, lons=f.lons)


def test_select_time():
    f = CurrentField.synthetic(nt=3, seed=1)
    one = f.select_time(1)
    assert one.u.shape == (1, 6, 8)
    assert one.times == [f.times[1]]
    np.testing.assert_array_equal(one.u[0], f.u[1])


def test_select_bbox():
    f = CurrentField.synthetic(nt=2, seed=1)
    sub = f.select_bbox((-90.0, 42.0, -86.0, 45.0))
    assert sub.bounds[0] >= -90.0 and sub.bounds[2] <= -86.0
    assert sub.u.shape[0] == 2
    assert len(sub.times) == 2


def test_select_bbox_no_overlap():
    f = CurrentField.synthetic(seed=1)
    with pytest.raises(ValueError, match="does not overlap"):
        f.select_bbox((0.0, 0.0, 1.0, 1.0))


def test_zonal_mean():
    f = CurrentField.synthetic(nt=2, seed=1)
    m = f.zonal_mean(0, "speed")
    assert m == pytest.approx(float(np.nanmean(np.hypot(f.u[0], f.v[0]))))
    t = f.zonal_mean(1, "temperature")
    assert t == pytest.approx(float(np.nanmean(f.temperature[1])))
    with pytest.raises(ValueError, match="unknown stat"):
        f.zonal_mean(0, "salinity")


def test_zonal_mean_no_temperature():
    f = CurrentField.synthetic(seed=1)
    f.temperature = None
    with pytest.raises(ValueError, match="no temperature"):
        f.zonal_mean(0, "temperature")


def test_dict_round_trip():
    f = CurrentField.synthetic(nt=2, seed=2)
    f.provenance = {"source": "unit"}
    g = CurrentField.from_dict(json.loads(json.dumps(f.to_dict())))
    np.testing.assert_array_equal(g.u, f.u)
    np.testing.assert_array_equal(g.temperature, f.temperature)
    assert g.times == f.times
    assert g.provenance == {"source": "unit"}


def test_dict_missing_keys():
    with pytest.raises(ValueError, match="missing keys"):
        CurrentField.from_dict({"u": []})


def test_json_round_trip(tmp_path):
    f = CurrentField.synthetic(nt=2, seed=5)
    p = str(tmp_path / "field.json")
    f.to_json(p)
    g = CurrentField.from_json(p)
    np.testing.assert_array_equal(g.v, f.v)
    assert g.crs == f.crs

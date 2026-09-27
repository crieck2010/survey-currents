"""Offline tests for currents.oceancolor.

Everything here runs with stdlib + numpy only (netCDF4-gated tests
skip when it is absent). Live CoastWatch requests are never made:
``_download_bytes`` is monkeypatched to return a synthetic NetCDF
payload built in-memory.
"""

import datetime as _dt
import json
import os
import subprocess
import sys

import numpy as np
import pytest

from currents import oceancolor as oc
from currents.oceancolor import (
    OceanColorField,
    fetch_oceancolor,
    fetch_oceancolor_coastwatch,
    oceancolor_urls,
)

netcdf4 = pytest.importorskip("netCDF4")  # noqa: F401  (skip whole module without it)


# ---------------------------------------------------------------------------
# synthetic CoastWatch-like NetCDF payload (offline fixture)
# ---------------------------------------------------------------------------


def _make_payload(nt=2, ny=6, nx=8, lons=None, lats=None, seed=12):
    """Build a CoastWatch-shaped NetCDF payload in memory; return bytes."""
    import tempfile
    if lons is None:
        lons = np.linspace(-80.0, -40.0, nx)
    if lats is None:
        lats = np.linspace(20.0, 50.0, ny)
    rng = np.random.default_rng(seed)
    path = os.path.join(tempfile.mkdtemp(), "chl.nc")
    ds = netcdf4.Dataset(path, mode="w")
    ds.createDimension("time", nt)
    ds.createDimension("altitude", 1)
    ds.createDimension("latitude", ny)
    ds.createDimension("longitude", nx)
    tvar = ds.createVariable("time", "f8", ("time",))
    tvar.units = "days since 2012-01-01 00:00:00"
    tvar[:] = np.arange(0, nt * 30, 30)
    ds.createVariable("altitude", "f4", ("altitude",))[:] = [0.0]
    ds.createVariable("latitude", "f4", ("latitude",))[:] = lats
    ds.createVariable("longitude", "f4", ("longitude",))[:] = lons
    cvar = ds.createVariable("chlor_a", "f4",
                             ("time", "altitude", "latitude", "longitude"),
                             fill_value=-999.0)
    cvar.units = "mg m^-3"
    grid = 0.3 * 10.0 ** rng.normal(0, 0.3, size=(nt, 1, ny, nx))
    grid[:, :, :, 0] = -999.0  # fake cloud column
    cvar[:] = grid
    ds.close()
    with open(path, "rb") as fh:
        data = fh.read()
    os.unlink(path)
    return data


@pytest.fixture()
def payload():
    return _make_payload()


@pytest.fixture()
def patched_download(payload, monkeypatch, tmp_path):
    """Monkeypatch _download_bytes -> payload; cache into tmp_path."""
    calls = {"n": 0}

    def fake(url, timeout=300):
        calls["n"] += 1
        assert url.startswith("https://coastwatch.noaa.gov/erddap/griddap/")
        assert "chlor_a" in url
        return payload

    monkeypatch.setattr(oc, "_download_bytes", fake)
    return calls


BBOX = (-80.0, 20.0, -40.0, 50.0)


# ---------------------------------------------------------------------------
# URL builder
# ---------------------------------------------------------------------------


class TestOceancolorUrls:
    def test_monthly_url_shape(self):
        urls = oceancolor_urls(BBOX, "2024-01-01", "2024-03-01")
        assert len(urls) == 1
        url = urls[0]
        # default sensor is now MODIS Aqua
        assert url.startswith(
            "https://coastwatch.noaa.gov/erddap/griddap/"
            "erdMH1chlamday_R2022SQ.nc?chlor_a[")
        assert "(0.0)" in url  # singleton altitude axis pinned

    def test_sensor_cadence_select_datasets(self):
        assert "erdMH1chlamday_R2022SQ" in oceancolor_urls(
            BBOX, "2024-01-01", "2024-02-01")[0]  # default = modis-aqua
        assert "nesdisVHNSQchlaMonthly" in oceancolor_urls(
            BBOX, "2024-01-01", "2024-02-01",
            sensor="viirs-snpp")[0]
        assert "nesdisVHNSQchlaWeekly" in oceancolor_urls(
            BBOX, "2024-01-01", "2024-02-01",
            sensor="viirs-snpp", cadence="weekly")[0]
        assert "pmlEsaCCI60OceanColorMonthly" in oceancolor_urls(
            BBOX, "2000-01-01", "2000-02-01",
            sensor="multi")[0]

    def test_unknown_sensor_cadence_raises(self):
        with pytest.raises(ValueError, match="known combinations"):
            oceancolor_urls(BBOX, "2024-01-01", "2024-02-01",
                            sensor="modis-aqua", cadence="weekly")

    def test_dates_before_record_raise(self):
        # default sensor is MODIS Aqua: record starts 2002-07-01
        with pytest.raises(ValueError, match="starts 2002"):
            oceancolor_urls(BBOX, "2000-01-01", "2000-02-01")
        with pytest.raises(ValueError, match="starts 2012"):
            oceancolor_urls(BBOX, "2000-01-01", "2000-02-01",
                            sensor="viirs-snpp")

    def test_start_after_end_raises(self):
        with pytest.raises(ValueError, match="after end"):
            oceancolor_urls(BBOX, "2024-03-01", "2024-01-01")

    def test_future_end_raises(self):
        future = (_dt.date.today() + _dt.timedelta(days=30)).isoformat()
        with pytest.raises(ValueError, match="in the future"):
            oceancolor_urls(BBOX, "2024-01-01", future)

    def test_bad_bbox_raises(self):
        with pytest.raises(ValueError, match="[Bb]box|longitude"):
            oceancolor_urls((200.0, 20.0, -40.0, 50.0),
                            "2024-01-01", "2024-02-01")

    def test_bad_stride_raises(self):
        with pytest.raises(ValueError, match="stride"):
            oceancolor_urls(BBOX, "2024-01-01", "2024-02-01", stride=0)

    def test_mirror_base_override(self):
        urls = oceancolor_urls(BBOX, "2024-01-01", "2024-02-01",
                               erddap_base="https://mirror.example/erddap")
        assert urls[0].startswith("https://mirror.example/erddap/griddap/")


# ---------------------------------------------------------------------------
# payload parsing
# ---------------------------------------------------------------------------


class TestParsePayload:
    def test_parse_grid_shape_and_coords(self, payload):
        lats, lons, times, chl = oc._parse_coastwatch_bytes(payload)
        assert lats.shape == (6,) and lons.shape == (8,)
        assert len(times) == 2
        assert chl.shape == (2, 6, 8)

    def test_lons_normalized_sorted(self):
        data = _make_payload(lons=np.linspace(300.0, 320.0, 8))  # 0..360 style
        _, lons, _, chl = oc._parse_coastwatch_bytes(data)
        assert np.all(lons <= 180.0) and np.all(lons >= -180.0)
        assert np.all(np.diff(lons) > 0)
        # column sorting is carried into the grid
        assert chl.shape == (2, 6, 8)

    def test_fill_value_is_masked(self, payload):
        _, _, _, chl = oc._parse_coastwatch_bytes(payload)
        assert bool(chl.mask[..., 0].all())  # fake cloud column stays NaN
        assert not bool(chl.mask[..., 1:].any())

    def test_missing_variable_raises(self):
        data = _make_payload()
        # rebuild without chlor_a: strip the variable
        import tempfile
        path = os.path.join(tempfile.mkdtemp(), "no_chl.nc")
        src = netcdf4.Dataset(path, mode="r", memory=data)
        dst_path = os.path.join(tempfile.mkdtemp(), "stripped.nc")
        dst = netcdf4.Dataset(dst_path, mode="w")
        for name, dim in src.dimensions.items():
            dst.createDimension(name, len(dim))
        for name, var in src.variables.items():
            if name == "chlor_a":
                continue
            v = dst.createVariable(name, var.dtype, var.dimensions)
            v[:] = var[:]
            v.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
        dst.close()
        src.close()
        with open(dst_path, "rb") as fh:
            stripped = fh.read()
        with pytest.raises(ValueError, match="missing variable"):
            oc._parse_coastwatch_bytes(stripped)


# ---------------------------------------------------------------------------
# fetch (offline, monkeypatched transport)
# ---------------------------------------------------------------------------


class TestFetchCoastwatch:
    def test_fetch_returns_field(self, patched_download, tmp_path):
        f = fetch_oceancolor_coastwatch(
            BBOX, "2024-01-01", "2024-03-01", work_dir=str(tmp_path))
        assert isinstance(f, OceanColorField)
        assert f.values.shape == (2, 6, 8)
        # default sensor is now MODIS Aqua
        assert f.source == "noaa-coastwatch/erdMH1chlamday_R2022SQ"
        assert f.provenance["units"] == "mg m^-3"
        assert f.provenance["cadence"] == "monthly"
        assert f.provenance["sensor"] == "MODIS (Aqua)"
        assert f.provenance["dataset_id"] == "erdMH1chlamday_R2022SQ"
        assert len(f.provenance["urls"]) == 1
        assert len(f.provenance["sha256"]) == 64
        assert "retrieved_at" in f.provenance
        assert len(f.provenance["gap_fractions"]) == 2
        # fake cloud column: 1/8 of cells masked
        assert f.provenance["gap_fractions"][0] == pytest.approx(0.125)
        assert patched_download["n"] == 1

    def test_cache_reuse_skips_download(self, patched_download, tmp_path):
        fetch_oceancolor_coastwatch(BBOX, "2024-01-01", "2024-03-01",
                                    work_dir=str(tmp_path))
        assert patched_download["n"] == 1
        fetch_oceancolor_coastwatch(BBOX, "2024-01-01", "2024-03-01",
                                    work_dir=str(tmp_path))
        assert patched_download["n"] == 1  # still one download
        # sidecar holds the sha256
        nc_files = [p for p in os.listdir(tmp_path) if p.endswith(".nc")]
        assert len(nc_files) == 1
        side = os.path.join(tmp_path, nc_files[0] + ".sha256")
        assert os.path.exists(side)

    def test_corrupt_cache_recovers(self, patched_download, tmp_path):
        fetch_oceancolor_coastwatch(BBOX, "2024-01-01", "2024-03-01",
                                    work_dir=str(tmp_path))
        nc_files = [p for p in os.listdir(tmp_path) if p.endswith(".nc")]
        with open(os.path.join(tmp_path, nc_files[0]), "wb") as fh:
            fh.write(b"corrupted")
        f = fetch_oceancolor_coastwatch(BBOX, "2024-01-01", "2024-03-01",
                                        work_dir=str(tmp_path))
        assert patched_download["n"] == 2  # re-downloaded
        assert f.values.shape == (2, 6, 8)

    def test_network_failure_is_actionable(self, monkeypatch, tmp_path):
        def boom(url, timeout=300):
            raise ConnectionError("nope")

        monkeypatch.setattr(oc, "_download_bytes", boom)
        with pytest.raises(RuntimeError, match="CoastWatch ocean-color"):
            fetch_oceancolor_coastwatch(BBOX, "2024-01-01", "2024-02-01",
                                        work_dir=str(tmp_path))

    def test_product_and_source_gates(self, tmp_path):
        with pytest.raises(ValueError, match="product"):
            fetch_oceancolor(BBOX, "2024-01-01", "2024-02-01",
                             product="kd490")
        with pytest.raises(ValueError, match="source"):
            fetch_oceancolor(BBOX, "2024-01-01", "2024-02-01",
                             source="mars")

    def test_auto_fallback_without_creds_is_honest(
            self, monkeypatch, tmp_path):
        def boom(url, timeout=300):
            raise ConnectionError("down")

        monkeypatch.setattr(oc, "_download_bytes", boom)
        # deterministic: no Earthdata credentials in the test env
        from currents.sst_global import earthdata_credentials  # noqa: F401
        monkeypatch.setattr("currents.sst_global.earthdata_credentials",
                            lambda: None)
        with pytest.raises(RuntimeError, match="every source"):
            fetch_oceancolor(BBOX, "2024-01-01", "2024-02-01",
                             source="auto", work_dir=str(tmp_path))

    def test_cache_key_differs_by_cadence(self, patched_download, tmp_path):
        fetch_oceancolor_coastwatch(BBOX, "2024-01-01", "2024-02-01",
                                    work_dir=str(tmp_path))
        fetch_oceancolor_coastwatch(BBOX, "2024-01-01", "2024-02-01",
                                    sensor="viirs-snpp", cadence="weekly",
                                    work_dir=str(tmp_path))
        assert patched_download["n"] == 2  # distinct cache entries


# ---------------------------------------------------------------------------
# OceanColorField model
# ---------------------------------------------------------------------------


class TestOceanColorField:
    def test_synthetic_shape_and_bounds(self):
        f = OceanColorField.synthetic()
        assert f.values.shape == (4, 6, 8)
        assert f.lats.shape == (6,) and f.lons.shape == (8,)
        assert len(f.times) == 4
        assert f.bounds == pytest.approx((-60.0, -30.0, 60.0, 30.0))

    def test_synthetic_gap_and_positive_values(self):
        f = OceanColorField.synthetic(seed=12, gap_fraction=0.25)
        assert f.gap_fraction(0) == pytest.approx(0.25, abs=0.05)
        vals = f.values[0].compressed()
        assert np.all(vals > 0)

    def test_spatial_median_ignores_nans(self):
        f = OceanColorField.synthetic()
        med = f.spatial_median(0)
        assert np.isfinite(med) and med > 0

    def test_shape_mismatch_raises(self):
        f = OceanColorField.synthetic()
        with pytest.raises(ValueError, match="does not match"):
            OceanColorField(values=np.zeros((3, 3, 3)), times=f.times[:2],
                            lats=f.lats[:3], lons=f.lons[:3])

    def test_chl_alias(self):
        f = OceanColorField.synthetic()
        assert f.chl is f.values

    def test_select_time(self):
        f = OceanColorField.synthetic(nt=4)
        g = f.select_time(2)
        assert g.values.shape == (1, 6, 8)
        assert g.times == [f.times[2]]

    def test_select_bbox(self):
        f = OceanColorField.synthetic()
        g = f.select_bbox((-60.0, -30.0, 0.0, 0.0))
        assert g.lons[-1] <= 0.0 and g.lats[-1] <= 0.0
        assert g.values.shape[0] == 4

    def test_select_bbox_no_overlap_raises(self):
        f = OceanColorField.synthetic()
        with pytest.raises(ValueError, match="does not overlap"):
            f.select_bbox((100.0, 40.0, 120.0, 50.0))

    def test_dict_roundtrip(self):
        f = OceanColorField.synthetic()
        g = OceanColorField.from_dict(f.to_dict())
        assert g.values.shape == f.values.shape
        assert np.ma.allequal(g.values, f.values)
        assert g.times == f.times

    def test_dict_missing_keys_raise(self):
        with pytest.raises(ValueError, match="missing keys"):
            OceanColorField.from_dict({"values": [[1.0]]})
        with pytest.raises(ValueError, match="missing 'values'"):
            OceanColorField.from_dict({"times": ["x"], "lats": [0.0],
                                       "lons": [0.0]})

    def test_json_roundtrip(self, tmp_path):
        f = OceanColorField.synthetic()
        path = str(tmp_path / "chl.json")
        f.to_json(path)
        g = OceanColorField.from_json(path)
        assert np.ma.allequal(g.values, f.values)
        assert json.load(open(path))["provenance"]


# ---------------------------------------------------------------------------
# CLI smoke tests (offline)
# ---------------------------------------------------------------------------


def _run_cli(*args, env=None):
    cmd = [sys.executable, "-m", "currents.cli", *args]
    e = dict(os.environ)
    e["SURVEY_CURRENTS_CACHE"] = e.get("SURVEY_CURRENTS_CACHE", "")
    if env:
        e.update(env)
    return subprocess.run(cmd, capture_output=True, text=True, env=e)


class TestOceancolorCli:
    def test_synthetic_cli_runs_offline(self, tmp_path):
        out = str(tmp_path / "syn")  # --out is a path prefix; CLI appends .json
        res = _run_cli("oceancolor-synthetic", "--out", out)
        assert res.returncode == 0, res.stderr
        assert os.path.exists(out + ".json")

    def test_fetch_cli_needs_args(self):
        res = _run_cli("fetch-oceancolor")
        assert res.returncode != 0


# ---------------------------------------------------------------------------
# OBPG fallback (offline, mocked transport)
# ---------------------------------------------------------------------------


def _make_obpg_payload(nt=1, ny=36, nx=72, flipped=True, seed=7):
    """Build a SeaDAS L3-mapped-shaped NetCDF payload; return bytes."""
    import tempfile
    rng = np.random.default_rng(seed)
    path = os.path.join(tempfile.mkdtemp(), "obpg.nc")
    ds = netcdf4.Dataset(path, mode="w")
    ds.createDimension("lat", ny)
    ds.createDimension("lon", nx)
    ds.createDimension("time", nt)
    lat = ds.createVariable("lat", "f4", ("lat",))
    lon = ds.createVariable("lon", "f4", ("lon",))
    la = np.linspace(-89.9, 89.9, ny)
    if flipped:  # SeaDAS L3 files run north-to-south
        la = la[::-1]
    lat[:] = la
    lon[:] = np.linspace(-179.9, 179.9, nx)
    tvar = ds.createVariable("time", "f8", ("time",))
    tvar.units = "days since 2024-01-01 00:00:00"
    tvar[:] = np.arange(nt) * 30.0 + 15.0
    v = ds.createVariable("chlor_a", "f4", ("time", "lat", "lon"),
                          fill_value=-999.0)
    v.units = "mg m^-3"
    grid = np.exp(rng.normal(-1, 1, (nt, ny, nx))).astype("f4")
    grid[:, :, 24] = -999.0  # fake cloud column (inside the test BBOX)
    v[:] = grid
    ds.close()
    with open(path, "rb") as fh:
        data = fh.read()
    os.unlink(path)
    return data


@pytest.fixture()
def obpg_payload():
    return _make_obpg_payload()


@pytest.fixture()
def patched_obpg(monkeypatch, obpg_payload):
    """Mock the Earthdata-authenticated opener -> synthetic OBPG payload."""
    from currents.sst_global import earthdata_credentials  # noqa: F401

    class _Resp:
        def __init__(self, data):
            self._data = data

        def read(self):
            return self._data

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    class _Opener:
        def open(self, url, timeout=None):
            assert url.startswith(oc.OBPG_DIRECTACCESS_BASE), url
            assert url.endswith(".nc"), url
            return _Resp(obpg_payload)

    monkeypatch.setattr("currents.sst_global.earthdata_credentials",
                        lambda: ("user", "pass"))
    monkeypatch.setattr(oc, "_obpg_opener", lambda: _Opener())


class TestDefaultSignature:
    def test_public_defaults_are_modis_aqua(self):
        import inspect
        sig = inspect.signature(fetch_oceancolor)
        assert sig.parameters["sensor"].default == "modis-aqua"
        assert sig.parameters["product"].default == "chlorophyll-a"
        assert sig.parameters["cadence"].default == "monthly"
        assert sig.parameters["source"].default == "coastwatch"
        sig2 = inspect.signature(fetch_oceancolor_coastwatch)
        assert sig2.parameters["sensor"].default == "modis-aqua"

    def test_obpg_fallback_is_a_named_source(self):
        with pytest.raises(ValueError, match="source"):
            fetch_oceancolor(BBOX, "2024-01-01", "2024-02-01",
                             source="obpg-mirror")  # not a real source


class TestObpgFileNames:
    def test_monthly_names_follow_documented_convention(self):
        names = oc.obpg_file_names("2024-01-10", "2024-03-20", "monthly")
        assert names == [
            "AQUA_MODIS.20240101_20240131.L3m.MO.CHL.chlor_a.4km.nc",
            "AQUA_MODIS.20240201_20240229.L3m.MO.CHL.chlor_a.4km.nc",
            "AQUA_MODIS.20240301_20240331.L3m.MO.CHL.chlor_a.4km.nc",
        ]

    def test_weekly_names_use_obpg_8day_bins(self):
        names = oc.obpg_file_names("2024-01-10", "2024-01-20", "weekly")
        assert names == [
            "AQUA_MODIS.20240109_20240116.L3m.8D.CHL.chlor_a.4km.nc",
            "AQUA_MODIS.20240117_20240124.L3m.8D.CHL.chlor_a.4km.nc",
        ]

    def test_daily_names_are_single_date(self):
        names = oc.obpg_file_names("2024-01-10", "2024-01-12", "daily")
        assert names == [
            "AQUA_MODIS.20240110.L3m.DAY.CHL.chlor_a.4km.nc",
            "AQUA_MODIS.20240111.L3m.DAY.CHL.chlor_a.4km.nc",
            "AQUA_MODIS.20240112.L3m.DAY.CHL.chlor_a.4km.nc",
        ]

    def test_unknown_cadence_raises(self):
        with pytest.raises(ValueError, match="cadence"):
            oc.obpg_file_names("2024-01-01", "2024-01-31", "hourly")

    def test_pre_record_dates_raise(self, patched_obpg):
        with pytest.raises(ValueError, match="2002-07-04"):
            oc.fetch_oceancolor_obpg(BBOX, "1999-01-01", "1999-02-01",
                                     work_dir="/tmp/nope")


class TestObpgParse:
    def test_parse_3d_masks_fill_and_decodes_time(self, obpg_payload):
        lats, lons, times, chl = oc._parse_obpg_bytes(
            obpg_payload, "http://x/test.nc")
        assert chl.shape == (1, 36, 72)
        assert chl.mask[:, :, 24].all()  # fake cloud column masked
        assert not chl.mask[:, :, 0].any()  # and nothing else is
        assert times and times[0].startswith("2024-01-16")

    def test_parse_2d_gains_time_axis(self):
        import tempfile
        path = os.path.join(tempfile.mkdtemp(), "flat.nc")
        ds = netcdf4.Dataset(path, mode="w")
        ds.createDimension("lat", 4)
        ds.createDimension("lon", 5)
        ds.createVariable("lat", "f4", ("lat",))[:] = np.linspace(-1, 1, 4)
        ds.createVariable("lon", "f4", ("lon",))[:] = np.linspace(0, 4, 5)
        v = ds.createVariable("chlor_a", "f4", ("lat", "lon"))
        v[:] = np.full((4, 5), 0.5, dtype="f4")
        ds.close()
        with open(path, "rb") as fh:
            payload = fh.read()
        os.unlink(path)
        _lats, _lons, times, chl = oc._parse_obpg_bytes(
            payload, "http://x/flat.nc")
        assert chl.shape == (1, 4, 5)
        assert times == []

    def test_missing_variable_raises(self):
        import tempfile
        path = os.path.join(tempfile.mkdtemp(), "bad.nc")
        ds = netcdf4.Dataset(path, mode="w")
        ds.createDimension("x", 2)
        ds.close()
        with open(path, "rb") as fh:
            payload = fh.read()
        os.unlink(path)
        with pytest.raises(ValueError, match="chlor_a"):
            oc._parse_obpg_bytes(payload, "http://x/bad.nc")


class TestFetchObpg:
    def test_fetch_subsets_bbox_and_records_provenance(
            self, patched_obpg, tmp_path):
        f = oc.fetch_oceancolor_obpg(BBOX, "2024-01-01", "2024-03-01",
                                     work_dir=str(tmp_path))
        assert isinstance(f, OceanColorField)
        assert f.values.shape == (3, 6, 8)
        assert f.values.shape[0] == 3  # Jan, Feb, Mar monthly files
        assert (np.diff(f.lats) > 0).all()  # normalized south-to-north
        assert f.source.startswith("nasa-obpg/")
        assert f.provenance["cmr_collection"] == (
            "MODISA_L3m_CHL (C3380709133-OB_CLOUD)")
        assert all(u.startswith(oc.OBPG_DIRECTACCESS_BASE)
                   for u in f.provenance["urls"])
        assert all(n.endswith(".nc") for n in f.provenance["urls"])
        assert f.provenance["live_verified"] is False
        assert f.gap_fraction(0) == pytest.approx(1 / 8, rel=0.05)

    def test_fetch_via_source_kwarg(self, patched_obpg, tmp_path):
        f = fetch_oceancolor(BBOX, "2024-01-01", "2024-02-01",
                             source="obpg", work_dir=str(tmp_path))
        assert isinstance(f, OceanColorField)
        assert f.values.shape[0] == 2  # Jan + Feb monthly files

    def test_no_credentials_raises_credentials_missing(self, monkeypatch,
                                                       tmp_path):
        from currents.sst_global import CredentialsMissing
        monkeypatch.setattr("currents.sst_global.earthdata_credentials",
                            lambda: None)
        with pytest.raises(CredentialsMissing, match="Earthdata Login"):
            fetch_oceancolor(BBOX, "2024-01-01", "2024-02-01",
                             source="obpg", work_dir=str(tmp_path))

    def test_rejected_credentials_raise_credentials_missing(
            self, monkeypatch, tmp_path):
        import urllib.error
        from currents.sst_global import CredentialsMissing
        monkeypatch.setattr("currents.sst_global.earthdata_credentials",
                            lambda: ("user", "wrong"))

        class _BadOpener:
            def open(self, url, timeout=None):
                raise urllib.error.HTTPError(url, 401, "Unauthorized",
                                             {}, None)

        monkeypatch.setattr(oc, "_obpg_opener", lambda: _BadOpener())
        with pytest.raises(CredentialsMissing, match="401"):
            oc.fetch_oceancolor_obpg(BBOX, "2024-01-01", "2024-02-01",
                                     work_dir=str(tmp_path))

    def test_auto_names_all_sources_when_everything_fails(
            self, monkeypatch, tmp_path):
        from currents.sst_global import CredentialsMissing

        def boom(url, timeout=300):
            raise ConnectionError("down")

        monkeypatch.setattr(oc, "_download_bytes", boom)
        monkeypatch.setattr("currents.sst_global.earthdata_credentials",
                            lambda: None)
        with pytest.raises(RuntimeError) as excinfo:
            fetch_oceancolor(BBOX, "2024-01-01", "2024-02-01",
                             source="auto", work_dir=str(tmp_path))
        msg = str(excinfo.value)
        assert "CoastWatch" in msg and "OBPG" in msg and "CMEMS" in msg

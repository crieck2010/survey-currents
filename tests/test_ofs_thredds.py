"""Tests for currents.ofs_thredds — fully offline.

Live-network tests are never in this suite. The fixtures are real
responses recorded 2026-10-01 from the CO-OPS THREDDS server:

* ``ofs_thredds_catalog_20261001.xml`` / ``ofs_thredds_catalog_20260930.xml``
  — real SSCOFS day catalogs, trimmed to the handful of regulargrid files
  the tests address.
* ``ofs_thredds_sscofs.das`` — real DAS for
  ``sscofs.t15z.20261001.regulargrid.n003.nc`` (time units, _FillValue).
* ``ofs_thredds_subset.dds`` — real unconstrained DDS for the same file.
* ``ofs_thredds_lat_col.ascii`` / ``lon_row`` / ``lat_row0`` / ``lon_col0``
  — real 1-D grid-vector probes (1553 / 1519 points).
* ``ofs_thredds_subset_n003.dods`` / ``.ascii`` — real 3x4 surface subset
  (y 600..602, x 1300..1303) of the n003 file, binary and ASCII forms.
* ``ofs_thredds_subset_f003.dods`` / ``.ascii`` — same window of the f003
  file (valid 2026-10-01T18:00Z).
* ``ofs_thredds_time_calib.ascii`` — real numeric ``time`` of
  ``sscofs.t09z.20260930.regulargrid.n003.nc`` (valid 2026-09-30T06:00Z).

The HTTP layer (``ofs_thredds._http_get_bytes``) is monkeypatched to serve
these fixtures by URL pattern; nothing here touches the network.
"""

import datetime as dt
import os
import re
import struct

import numpy as np
import pytest

from currents import ofs_thredds
from currents.ofs_thredds import (
    OFS_THREDDS_MODELS,
    THREDDS_BASE,
    UnavailableRangeError,
    _calibrate_nowcast_span,
    _pick_file,
    day_catalog_url,
    dods_base_url,
    fetch_ascii_values,
    fetch_dds,
    fetch_grid_vectors,
    fetch_ofs_thredds,
    index_range,
    list_regulargrid_files,
    nominal_valid_time,
    parse_dods_response,
    parse_fill_value,
    parse_ofs_filename,
    parse_time_units,
    subset_constraint,
)

FIX = os.path.join(os.path.dirname(__file__), "fixtures")

UP_N003 = ("NOAA/SSCOFS/MODELS/2026/10/01/"
           "sscofs.t15z.20261001.regulargrid.n003.nc")
UP_F003 = ("NOAA/SSCOFS/MODELS/2026/10/01/"
           "sscofs.t15z.20261001.regulargrid.f003.nc")

# Test window recorded in the subset fixtures: y 600..602, x 1300..1303,
# i.e. lats 47.37..47.38, lons -123.03..-123.015. The bbox below maps to
# exactly those indices (with margin away from exact grid nodes).
BBOX = (-123.034, 47.366, -123.011, 47.384)


def _read(name: str) -> bytes:
    with open(os.path.join(FIX, name), "rb") as fh:
        return fh.read()


def _fake_get(url: str, timeout: float = 120.0) -> bytes:
    """Serve recorded fixtures by URL pattern (no network)."""
    if "MODELS/2026/10/01/catalog.xml" in url:
        return _read("ofs_thredds_catalog_20261001.xml")
    if "MODELS/2026/09/30/catalog.xml" in url:
        return _read("ofs_thredds_catalog_20260930.xml")
    if url.endswith(".das"):
        return _read("ofs_thredds_sscofs.das")
    if url.endswith(".dds"):
        return _read("ofs_thredds_subset.dds")
    if ".dods?" in url:
        if "regulargrid.n003.nc" in url:
            return _read("ofs_thredds_subset_n003.dods")
        if "regulargrid.f003.nc" in url:
            return _read("ofs_thredds_subset_f003.dods")
        raise AssertionError(f"unexpected .dods URL in test: {url}")
    if ".ascii?" in url:
        if "Latitude%5B0:1:1552%5D%5B0:1:0%5D" in url:
            return _read("ofs_thredds_lat_col.ascii")
        if "Latitude%5B0:1:0%5D%5B0:1:1518%5D" in url:
            return _read("ofs_thredds_lat_row0.ascii")
        if "Longitude%5B0:1:0%5D%5B0:1:1518%5D" in url:
            return _read("ofs_thredds_lon_row.ascii")
        if "Longitude%5B0:1:1552%5D%5B0:1:0%5D" in url:
            return _read("ofs_thredds_lon_col0.ascii")
        if "time%5B0:1:0%5D" in url:
            return _read("ofs_thredds_time_calib.ascii")
        raise AssertionError(f"unexpected .ascii URL in test: {url}")
    raise AssertionError(f"unexpected URL in test: {url}")


@pytest.fixture
def offline(monkeypatch):
    monkeypatch.setattr(ofs_thredds, "_http_get_bytes", _fake_get)


# ---------------------------------------------------------------------------
# URL construction (no network)
# ---------------------------------------------------------------------------


def test_day_catalog_url():
    url = day_catalog_url("SSCOFS", dt.date(2026, 10, 1))
    assert url == (f"{THREDDS_BASE}/catalog/NOAA/SSCOFS/MODELS/"
                   "2026/10/01/catalog.xml")


def test_dods_base_url():
    assert dods_base_url(UP_N003) == f"{THREDDS_BASE}/dodsC/{UP_N003}"


def test_subset_constraint_encodes_brackets():
    ce = subset_constraint(600, 602, 1300, 1303)
    assert "[" not in ce and "]" not in ce
    assert "u_eastward%5B0:1:0%5D%5B0:1:0%5D%5B600:1:602%5D%5B1300:1:1303%5D" in ce
    assert "v_northward" in ce and "temp" in ce and "time%5B0:1:0%5D" in ce


def test_known_models():
    assert "SSCOFS" in OFS_THREDDS_MODELS and "CBOFS" in OFS_THREDDS_MODELS
    assert len(OFS_THREDDS_MODELS) == 12


# ---------------------------------------------------------------------------
# Catalog parsing + filename -> valid time
# ---------------------------------------------------------------------------


def test_list_regulargrid_files(offline):
    files = list_regulargrid_files("SSCOFS", dt.date(2026, 10, 1))
    names = [f["name"] for f in files]
    assert "sscofs.t15z.20261001.regulargrid.n003.nc" in names
    assert all(f["url_path"].startswith("NOAA/SSCOFS/MODELS/2026/10/01/")
               for f in files)


def test_list_regulargrid_files_missing_day(monkeypatch):
    def boom(url, timeout=120.0):
        raise UnavailableRangeError(f"THREDDS request failed for {url}: 404")
    monkeypatch.setattr(ofs_thredds, "_http_get_bytes", boom)
    with pytest.raises(UnavailableRangeError, match="31 days"):
        list_regulargrid_files("SSCOFS", dt.date(2020, 1, 1))


def test_parse_ofs_filename():
    p = parse_ofs_filename("sscofs.t15z.20261001.regulargrid.n006.nc")
    assert (p["code"], p["kind"], p["hour"]) == ("sscofs", "n", 6)
    assert p["cycle"] == dt.datetime(2026, 10, 1, 15,
                                     tzinfo=dt.timezone.utc)
    with pytest.raises(ValueError):
        parse_ofs_filename("not_a_file.nc")


def test_nominal_valid_time_forecast():
    p = parse_ofs_filename("sscofs.t15z.20261001.regulargrid.f003.nc")
    assert nominal_valid_time(p, 6) == dt.datetime(
        2026, 10, 1, 18, tzinfo=dt.timezone.utc)


def test_nominal_valid_time_nowcast_span6():
    p = parse_ofs_filename("sscofs.t15z.20261001.regulargrid.n003.nc")
    assert nominal_valid_time(p, 6) == dt.datetime(
        2026, 10, 1, 12, tzinfo=dt.timezone.utc)


def test_nominal_valid_time_nowcast_span24():
    # WCOFS nowcast window is 24 h (verified live 2026-10-01: n024@t03z
    # is valid at the cycle time).
    p = parse_ofs_filename("wcofs.t03z.20261001.regulargrid.n024.nc")
    assert nominal_valid_time(p, 24) == dt.datetime(
        2026, 10, 1, 3, tzinfo=dt.timezone.utc)


def test_calibrate_nowcast_span(offline):
    index = ofs_thredds._day_file_index(
        "SSCOFS", [dt.date(2026, 9, 30), dt.date(2026, 10, 1)], 120.0)
    span, epoch = _calibrate_nowcast_span("SSCOFS", index, 120.0)
    assert span == 6
    assert epoch == dt.datetime(2018, 1, 1, tzinfo=dt.timezone.utc)


def test_pick_file_prefers_nowcast(offline):
    index = ofs_thredds._day_file_index(
        "SSCOFS", [dt.date(2026, 10, 1)], 120.0)
    # 15:00 is covered by both n006 (nowcast) and... only n006 here.
    e = _pick_file(dt.datetime(2026, 10, 1, 15, tzinfo=dt.timezone.utc),
                   index, 6, "nowcast")
    assert e["name"] == "sscofs.t15z.20261001.regulargrid.n006.nc"
    e = _pick_file(dt.datetime(2026, 10, 1, 12, tzinfo=dt.timezone.utc),
                   index, 6, "nowcast")
    assert e["name"] == "sscofs.t15z.20261001.regulargrid.n003.nc"


def test_pick_file_prefers_latest_cycle(offline):
    # Two forecast files cover 18:00 (t09z f009 and t15z f003): the later
    # cycle (shorter lead time) must win.
    index = ofs_thredds._day_file_index(
        "SSCOFS", [dt.date(2026, 10, 1)], 120.0)
    extra = {"name": "sscofs.t09z.20261001.regulargrid.f009.nc",
             "url_path": ("NOAA/SSCOFS/MODELS/2026/10/01/"
                          "sscofs.t09z.20261001.regulargrid.f009.nc"),
             **parse_ofs_filename("sscofs.t09z.20261001.regulargrid.f009.nc")}
    index[dt.date(2026, 10, 1)].append(extra)
    e = _pick_file(dt.datetime(2026, 10, 1, 18, tzinfo=dt.timezone.utc),
                   index, 6, "nowcast")
    assert e["name"] == "sscofs.t15z.20261001.regulargrid.f003.nc"


def test_pick_file_missing_hour(offline):
    index = ofs_thredds._day_file_index(
        "SSCOFS", [dt.date(2026, 10, 1)], 120.0)
    with pytest.raises(UnavailableRangeError, match="no nowcast file"):
        _pick_file(dt.datetime(2026, 10, 1, 4, tzinfo=dt.timezone.utc),
                   index, 6, "nowcast")


# ---------------------------------------------------------------------------
# DAS / DDS / ASCII parsing
# ---------------------------------------------------------------------------


def test_parse_time_units(offline):
    das = _read("ofs_thredds_sscofs.das").decode()
    assert parse_time_units(das) == dt.datetime(2018, 1, 1,
                                                tzinfo=dt.timezone.utc)


def test_parse_fill_value(offline):
    das = _read("ofs_thredds_sscofs.das").decode()
    assert parse_fill_value(das) == -99999.0


def test_fetch_dds(offline):
    specs = fetch_dds(UP_N003)
    by_name = {n: (t, dict(d)) for n, t, d in specs}
    assert by_name["u_eastward"][0] == "Float32"
    assert by_name["u_eastward"][1] == {"time": 1, "Depth": 37,
                                       "ny": 1553, "nx": 1519}
    assert by_name["time"][1] == {"time": 1}


def test_fetch_ascii_values(offline):
    vals = fetch_ascii_values(UP_N003, "Latitude", ["0:1:1552", "0:1:0"])
    assert len(vals) == 1553
    assert abs(vals[0] - 44.37) < 1e-9
    assert all(b > a for a, b in zip(vals, vals[1:]))


def test_fetch_grid_vectors(offline):
    lats, lons = fetch_grid_vectors(UP_N003, 1553, 1519)
    assert lats.shape == (1553,) and lons.shape == (1519,)
    assert abs(lats[600] - 47.37) < 1e-9
    assert abs(lons[1300] - (-123.03)) < 1e-9


def test_fetch_grid_vectors_rejects_irregular(monkeypatch):
    def fake_ascii(url_path, var, slices, timeout=120.0):
        if var == "Latitude" and slices == ["0:1:1552", "0:1:0"]:
            return [44.0 + 0.005 * i for i in range(1553)]
        if var == "Longitude" and slices == ["0:1:0", "0:1:1518"]:
            return [-129.0 + 0.005 * i for i in range(1519)]
        if var == "Latitude" and slices == ["0:1:0", "0:1:1518"]:
            # irregular: Latitude varies along nx -> must raise
            return [44.0 + 0.001 * i for i in range(1519)]
        if var == "Longitude" and slices == ["0:1:1552", "0:1:0"]:
            return [-129.0] * 1553
        raise AssertionError(slices)
    monkeypatch.setattr(ofs_thredds, "fetch_ascii_values", fake_ascii)
    with pytest.raises(ValueError, match="non-regular"):
        fetch_grid_vectors(UP_N003, 1553, 1519)


# ---------------------------------------------------------------------------
# DAP2 binary parsing (binary vs recorded ASCII cross-check)
# ---------------------------------------------------------------------------


def _ascii_grid(name: str, path: str):
    """Parse a recorded .ascii subset into {(var): 2-D list}."""
    text = open(path).read()
    out, cur = {}, None
    for line in text.splitlines():
        m = re.match(r"^(u_eastward|v_northward|temp)\[\d+\]\[\d+\]\[(\d+)\]\[(\d+)\]$",
                     line.strip())
        if m:
            cur = (m.group(1), int(m.group(2)), int(m.group(3)))
            out[cur[0]] = [[0.0] * cur[2] for _ in range(cur[1])]
            continue
        m = re.match(r"^\[0\]\[0\]\[(\d+)\],\s*(.*)$", line.strip())
        if m and cur:
            row = int(m.group(1))
            out[cur[0]][row] = [float(x) for x in m.group(2).split(",")]
    return out


def test_parse_dods_response_matches_ascii():
    raw = _read("ofs_thredds_subset_n003.dods")
    arrays = parse_dods_response(raw)
    assert set(arrays) == {"time", "u_eastward", "v_northward", "temp"}
    assert arrays["u_eastward"].shape == (1, 1, 3, 4)
    ascii_grids = _ascii_grid("u_eastward",
                              os.path.join(FIX, "ofs_thredds_subset_n003.ascii"))
    for var in ("u_eastward", "v_northward", "temp"):
        expect = np.asarray(ascii_grids[var], dtype=np.float64)
        got = arrays[var][0, 0]
        np.testing.assert_allclose(got, expect, rtol=1e-6,
                                   err_msg=f"binary != ascii for {var}")
    # numeric time decodes to the file's valid time
    epoch = dt.datetime(2018, 1, 1, tzinfo=dt.timezone.utc)
    valid = epoch + dt.timedelta(seconds=float(arrays["time"].ravel()[0]))
    assert valid == dt.datetime(2026, 10, 1, 12, tzinfo=dt.timezone.utc)


def test_parse_dods_response_truncated():
    raw = _read("ofs_thredds_subset_n003.dods")
    marker = raw.find(b"Data:")
    with pytest.raises(UnavailableRangeError, match="truncated"):
        parse_dods_response(raw[:marker + 40])
    with pytest.raises(UnavailableRangeError, match="no DAP2 data marker"):
        parse_dods_response(b"garbage")


# ---------------------------------------------------------------------------
# index_range
# ---------------------------------------------------------------------------


def test_index_range(offline):
    lats, lons = fetch_grid_vectors(UP_N003, 1553, 1519)
    assert index_range(lats, 47.366, 47.384, "latitude", "x") == (600, 602)
    assert index_range(lons, -123.034, -123.011, "longitude", "x") == (1300, 1303)
    with pytest.raises(UnavailableRangeError, match="outside the model domain"):
        index_range(lats, 60.0, 61.0, "latitude", "x")


# ---------------------------------------------------------------------------
# End-to-end fetch (mocked HTTP)
# ---------------------------------------------------------------------------


def test_fetch_ofs_thredds_end_to_end(offline):
    field = fetch_ofs_thredds("SSCOFS", BBOX, "2026-10-01T12:00",
                              "2026-10-01T18:00", cadence_hours=6)
    assert field.source == "ofs-thredds/SSCOFS"
    assert field.times == ["2026-10-01T12:00:00+00:00",
                           "2026-10-01T18:00:00+00:00"]
    assert field.u.shape == (2, 3, 4)
    assert field.v.shape == (2, 3, 4)
    assert field.temperature.shape == (2, 3, 4)
    assert field.lats.shape == (3,) and field.lons.shape == (4,)
    assert field.lats[0] < field.lats[-1] and field.lons[0] < field.lons[-1]
    # land fill -> NaN, water values sane
    assert np.isnan(field.u).any()
    water = field.u[~np.isnan(field.u)]
    assert water.size > 0 and np.all(np.abs(water) < 5.0)
    t = field.temperature[~np.isnan(field.temperature)]
    assert np.all((t > -2.0) & (t < 40.0))
    # provenance records exact URLs
    steps = field.provenance["steps"]
    assert len(steps) == 2
    assert all(s["url"].startswith(f"{THREDDS_BASE}/dodsC/") for s in steps)
    assert "n003" in steps[0]["url"] and "f003" in steps[1]["url"]
    assert field.provenance["u_units"] == "m/s"
    assert field.forecast_hours == [-3, 3]
    assert field.model_run == "2026-10-01T15:00:00+00:00"


def test_fetch_ofs_thredds_conforms_to_current_field(offline):
    # The viz currents path reads .u/.v/.temperature/.times/.lats/.lons.
    field = fetch_ofs_thredds("SSCOFS", BBOX, "2026-10-01T12:00",
                              "2026-10-01T12:00")
    assert len(field.times) == 1
    for attr in ("u", "v", "temperature", "times", "lats", "lons",
                 "source", "provenance"):
        assert getattr(field, attr) is not None
    speed = np.hypot(np.nan_to_num(field.u), np.nan_to_num(field.v))
    assert speed.shape == field.u.shape


def test_fetch_unknown_model(offline):
    with pytest.raises(ValueError, match="unknown THREDDS OFS code"):
        fetch_ofs_thredds("NOPE", BBOX, "2026-10-01T12:00", "2026-10-01T13:00")


def test_fetch_bad_inputs(offline):
    with pytest.raises(ValueError, match="hour-aligned"):
        fetch_ofs_thredds("SSCOFS", BBOX, "2026-10-01T12:30", "2026-10-01T13:00")
    with pytest.raises(ValueError, match="bbox"):
        fetch_ofs_thredds("SSCOFS", (1.0, 2.0, 0.5, 3.0),
                          "2026-10-01T12:00", "2026-10-01T13:00")
    with pytest.raises(ValueError, match="before start"):
        fetch_ofs_thredds("SSCOFS", BBOX, "2026-10-01T13:00", "2026-10-01T12:00")
    with pytest.raises(ValueError, match="prefer"):
        fetch_ofs_thredds("SSCOFS", BBOX, "2026-10-01T12:00",
                          "2026-10-01T13:00", prefer="yesterday")


def test_fetch_bbox_outside_domain(offline):
    with pytest.raises(UnavailableRangeError, match="outside the model domain"):
        fetch_ofs_thredds("SSCOFS", (-100.0, 40.0, -99.0, 41.0),
                          "2026-10-01T12:00", "2026-10-01T12:00")


def test_fetch_missing_hour_is_honest(offline):
    with pytest.raises(UnavailableRangeError, match="no nowcast file"):
        fetch_ofs_thredds("SSCOFS", BBOX, "2026-10-01T04:00",
                          "2026-10-01T04:00")


def test_fetch_refuses_mislabeled_time(monkeypatch):
    # Serve the f003 (18:00) payload for the 12:00 request: the adapter must
    # refuse rather than mislabel.
    def fake_subset(url_path, y0, y1, x0, x1, timeout=120.0):
        raw = _read("ofs_thredds_subset_f003.dods")
        arrays = parse_dods_response(raw)
        url = dods_base_url(url_path) + ".dods?fake"
        return arrays, len(raw), url
    monkeypatch.setattr(ofs_thredds, "_http_get_bytes", _fake_get)
    monkeypatch.setattr(ofs_thredds, "fetch_subset", fake_subset)
    with pytest.raises(ValueError, match="does not match requested"):
        fetch_ofs_thredds("SSCOFS", BBOX, "2026-10-01T12:00",
                          "2026-10-01T12:00")

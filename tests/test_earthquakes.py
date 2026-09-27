"""Tests for currents.earthquakes (USGS ComCat FDSN event service).

The module is fully offline-testable: the URL builders, GeoJSON
parsers, and QuakeField model are pure, and the download seam
(``_cached_get_text``) is mocked for the fetch tests.

Live-network tests are skipped unless ``SURVEY_CURRENTS_LIVE=1``.
"""

import datetime as dt
import json as _json
import os
import urllib.parse

import numpy as np
import pytest

import currents.earthquakes as eq
from currents.earthquakes import (QuakeEvent, QuakeField, comcat_cache_dir,
                                  event_count_url, event_query_url,
                                  fetch_earthquakes)

BBOX = (-125.0, 32.0, -114.0, 42.0)  # California

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _geojson_fixture(n_events=5):
    """ComCat-style GeoJSON: n_events across two days, one without a
    magnitude, one 2-D coordinate (no depth), one quarry blast, one
    malformed feature (no geometry)."""
    feats = []
    for k in range(n_events):
        ms = int((dt.datetime(2026, 9, 20, 12, 0, tzinfo=dt.timezone.utc)
                  + dt.timedelta(hours=k * 13)).timestamp() * 1000)
        coords = [-120.0 + k * 0.5, 35.0 + k * 0.3, 10.0 + k * 5.0]
        if k == 3:
            coords = coords[:2]  # 2-D: no depth
        mag = 4.0 + k * 0.3
        if k == 4:
            mag = None
        etype = "quarry blast" if k == 2 else "earthquake"
        feats.append({
            "type": "Feature",
            "id": f"us6000ab{k:02d}",
            "properties": {
                "mag": mag,
                "magType": "ml" if mag else None,
                "place": f"{10 * k} km W of Testville, CA",
                "time": ms,
                "type": etype,
            },
            "geometry": {"type": "Point", "coordinates": coords},
        })
    feats.append({"type": "Feature", "id": "broken01", "properties": {},
                  "geometry": None})  # malformed
    return _json.dumps({
        "type": "FeatureCollection",
        "metadata": {"generated": 1790523122000, "status": 200,
                     "api": "2.7.0", "count": n_events},
        "features": feats,
    })


def _page_fixture(n, id_offset, d0, d1):
    """n synthetic events with unique ids (window + offset), spread
    over [d0, d1] — behaves like one real API page."""
    feats = []
    ndays = max(1, (d1 - d0).days + 1)
    for k in range(n):
        uid = f"{d0.isoformat()}-{id_offset + k:08d}"
        when = (dt.datetime(d0.year, d0.month, d0.day,
                            tzinfo=dt.timezone.utc)
                + dt.timedelta(days=(id_offset + k) % ndays,
                               seconds=((id_offset + k) * 977) % 86400))
        ms = int(when.timestamp() * 1000)
        feats.append({
            "type": "Feature",
            "id": f"usp{uid}",
            "properties": {
                "mag": round(2.0 + ((id_offset + k) % 50) / 10.0, 1),
                "magType": "ml",
                "place": f"Test event {uid}",
                "time": ms,
                "type": "earthquake",
            },
            "geometry": {"type": "Point",
                         "coordinates": [-120.0 + (k % 10) * 0.1,
                                         35.0 + (k % 10) * 0.1, 10.0]},
        })
    return _json.dumps({"type": "FeatureCollection",
                        "metadata": {"limit": n, "offset": id_offset},
                        "features": feats})


def _mock_comcat(monkeypatch, tmp_path, geojson=None, fixed_count=None,
                 count_per_day=3000):
    """Mock the download seam for fetch_earthquakes.

    The count endpoint returns ``fixed_count`` when given, else
    ``count_per_day * ndays`` for the queried window (so time-chunking
    recursion terminates). Query URLs return ``geojson`` when given,
    else a page fixture honoring the URL's limit/offset with unique
    ids per page — like the real API.
    """
    calls = []
    small_body = geojson

    def fake_cached_get(url, cache_dir, max_age_days, refresh=False):
        calls.append(url)
        q = dict(urllib.parse.parse_qsl(urllib.parse.urlsplit(url).query))
        d0 = dt.date.fromisoformat(q["starttime"])
        d1 = dt.date.fromisoformat(q["endtime"])
        if "/count" in url:
            if fixed_count is not None:
                total = fixed_count
            else:
                total = count_per_day * ((d1 - d0).days + 1)
            return str(total), {"url": url, "cache_hit": False,
                                "downloaded": True, "sha256": "x" * 64}
        if small_body is not None:
            return small_body, {"url": url, "cache_hit": False,
                                "downloaded": True, "sha256": "x" * 64}
        # Page fixture honoring limit/offset against the window total,
        # like the real API (last page is short). The real service's
        # offsets are 1-based: offset=1 returns the first event.
        total = (fixed_count if fixed_count is not None
                 else count_per_day * ((d1 - d0).days + 1))
        limit = int(q.get("limit", 20000))
        offset = int(q.get("offset", 1))
        first = offset - 1  # 1-based offset -> 0-based index
        remaining = max(0, total - first)
        body = _page_fixture(min(limit, remaining), first, d0, d1)
        return body, {"url": url, "cache_hit": False, "downloaded": True,
                      "sha256": "x" * 64}

    monkeypatch.setattr(eq, "_cached_get_text", fake_cached_get)
    monkeypatch.setattr(eq, "comcat_cache_dir", lambda: str(tmp_path))
    return calls


# ---------------------------------------------------------------------------
# URL builders
# ---------------------------------------------------------------------------

def test_event_query_url_params():
    url = event_query_url(BBOX, "2026-09-20", "2026-09-27",
                          min_magnitude=4.5, event_type="earthquake",
                          limit=2000, offset=4000)
    assert url.startswith("https://earthquake.usgs.gov/fdsnws/event/1/query?")
    assert "format=geojson" in url
    assert "starttime=2026-09-20" in url
    assert "endtime=2026-09-27" in url
    assert "minlatitude=32.000000" in url
    assert "maxlongitude=-114.000000" in url
    assert "minmagnitude=4.50" in url
    assert "eventtype=earthquake" in url
    assert "limit=2000" in url
    assert "offset=4000" in url


def test_event_query_url_defaults():
    url = event_query_url(BBOX, "2026-09-20", "2026-09-27")
    assert "format=geojson" in url
    assert "minmagnitude=0.00" in url
    assert "limit=" not in url and "offset=" not in url
    assert "eventtype=" not in url


def test_event_count_url():
    url = event_count_url(BBOX, "2026-09-20", "2026-09-27", min_magnitude=2.5)
    assert url.startswith("https://earthquake.usgs.gov/fdsnws/event/1/count?")
    assert "format=" not in url
    assert "minmagnitude=2.50" in url


def test_event_query_url_builders_are_pure():
    # Builders do not validate; fetch_earthquakes validates.
    url = event_query_url((0, 0, 0, 10), "2026-09-28", "2026-09-20")
    assert "starttime=2026-09-28" in url


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def test_parse_event_geojson_counts():
    events, n_malformed = eq._parse_event_geojson(_json.loads(_geojson_fixture(5)))
    assert len(events) == 5
    assert n_malformed == 1


def test_parse_event_fields():
    events, _ = eq._parse_event_geojson(_json.loads(_geojson_fixture(5)))
    e0 = events[0]
    assert e0["event_id"] == "us6000ab00"
    assert e0["lat"] == pytest.approx(35.0)
    assert e0["lon"] == pytest.approx(-120.0)
    assert e0["depth_km"] == pytest.approx(10.0)
    assert e0["magnitude"] == pytest.approx(4.0)
    assert e0["mag_type"] == "ml"
    assert "Testville" in e0["place"]
    assert e0["event_type"] == "earthquake"
    assert e0["time"].tzinfo is not None
    assert e0["time"].year == 2026


def test_parse_missing_magnitude_stays_none():
    events, _ = eq._parse_event_geojson(_json.loads(_geojson_fixture(5)))
    e4 = next(e for e in events if e["event_id"] == "us6000ab04")
    assert e4["magnitude"] is None  # never 0.0


def test_parse_2d_coords_depth_none():
    events, _ = eq._parse_event_geojson(_json.loads(_geojson_fixture(5)))
    e3 = next(e for e in events if e["event_id"] == "us6000ab03")
    assert e3["depth_km"] is None


def test_parse_quarry_blast_kept_honestly():
    events, _ = eq._parse_event_geojson(_json.loads(_geojson_fixture(5)))
    e2 = next(e for e in events if e["event_id"] == "us6000ab02")
    assert e2["event_type"] == "quarry blast"


def test_parse_skips_features_without_time():
    payload = {"features": [
        {"id": "x1", "properties": {"mag": 5.0},
         "geometry": {"coordinates": [-120.0, 35.0, 10.0]}},
    ]}
    events, n_malformed = eq._parse_event_geojson(payload)
    assert events == [] and n_malformed == 1


def test_parse_empty_body():
    events, n_malformed = eq._parse_event_geojson({})
    assert events == [] and n_malformed == 0


def test_parse_count():
    assert eq._parse_count("179\n") == 179
    assert eq._parse_count("0") == 0
    with pytest.raises(RuntimeError):
        eq._parse_count("not-a-number")


def test_dedupe_by_event_id():
    evs = [{"event_id": "a"}, {"event_id": "b"}, {"event_id": "a"}]
    out, n_dup = eq._dedupe_events(evs)
    assert [e["event_id"] for e in out] == ["a", "b"]
    assert n_dup == 1


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

def test_quake_event_round_trip():
    when = dt.datetime(2026, 9, 20, 12, 0, tzinfo=dt.timezone.utc)
    e = QuakeEvent(event_id="us1", time=when, lat=35.0, lon=-120.0,
                   depth_km=12.5, magnitude=5.2, mag_type="mww",
                   place="10 km W of Testville", event_type="earthquake")
    e2 = QuakeEvent.from_dict(e.to_dict())
    assert e2.event_id == "us1"
    assert e2.time == when
    assert e2.magnitude == pytest.approx(5.2)
    assert e2.depth_km == pytest.approx(12.5)
    assert e2.place == "10 km W of Testville"


def test_quake_event_depth_bins():
    assert QuakeEvent("a", dt.datetime.now(dt.timezone.utc),
                      0.0, 0.0, depth_km=10.0).depth_bin() == "shallow"
    assert QuakeEvent("a", dt.datetime.now(dt.timezone.utc),
                      0.0, 0.0, depth_km=150.0).depth_bin() == "intermediate"
    assert QuakeEvent("a", dt.datetime.now(dt.timezone.utc),
                      0.0, 0.0, depth_km=500.0).depth_bin() == "deep"
    assert QuakeEvent("a", dt.datetime.now(dt.timezone.utc),
                      0.0, 0.0).depth_bin() == "unknown"


def test_field_sorts_by_time():
    f = QuakeField.synthetic(seed=1)
    times = [e.time for e in f.events]
    assert times == sorted(times)


def test_field_largest():
    f = QuakeField.synthetic(n_events=20, seed=3)
    top = f.largest(3)
    mags = [e.magnitude for e in top]
    assert mags == sorted(mags, reverse=True)
    assert len(top) == 3


def test_field_select_time():
    f = QuakeField.synthetic(start="2024-01-01", end="2024-01-05",
                             n_events=30, seed=5)
    day = dt.date(2024, 1, 3)
    sel = f.select_time(day)
    assert all(e.date == day for e in sel)


def test_field_counts_by_day():
    f = QuakeField.synthetic(start="2024-01-01", end="2024-01-03",
                             n_events=30, seed=5)
    dates, counts = f.counts_by_day()
    assert len(dates) == 3
    assert int(counts.sum()) == 30


def test_field_synthetic_deterministic():
    a = QuakeField.synthetic(seed=13)
    b = QuakeField.synthetic(seed=13)
    assert [e.event_id for e in a.events] == [e.event_id for e in b.events]
    assert [e.magnitude for e in a.events] == [e.magnitude for e in b.events]
    lon_min, lat_min, lon_max, lat_max = a.bbox
    for e in a.events:
        assert lon_min <= e.lon <= lon_max
        assert lat_min <= e.lat <= lat_max


def test_field_json_round_trip(tmp_path):
    f = QuakeField.synthetic(n_events=10, seed=7)
    p = str(tmp_path / "quakes.json")
    f.to_json(p)
    g = QuakeField.from_json(p)
    assert len(g) == 10
    assert g.events[0].event_id == f.events[0].event_id
    assert g.provenance == f.provenance


def test_field_rejects_start_after_end():
    with pytest.raises(ValueError):
        QuakeField(events=[], bbox=BBOX, start="2026-09-28", end="2026-09-20")


def test_depth_bins_constant_documented():
    names = [n for n, _, _ in eq.COMCAT_DEPTH_BINS_KM]
    assert names == ["shallow", "intermediate", "deep"]


# ---------------------------------------------------------------------------
# Fetch (mocked)
# ---------------------------------------------------------------------------

def test_fetch_mocked_single_page(monkeypatch, tmp_path):
    calls = _mock_comcat(monkeypatch, tmp_path, fixed_count=5,
                         geojson=_geojson_fixture())
    f = fetch_earthquakes(BBOX, "2026-09-20", "2026-09-27",
                          min_magnitude=2.0, cache_dir=str(tmp_path))
    assert len(f) == 5
    assert f.source == "usgs"
    assert f.min_magnitude == pytest.approx(2.0)
    assert f.provenance["n_events"] == 5
    assert f.provenance["n_events_malformed"] == 1
    assert "catalog_completeness_note" in f.provenance
    assert any("/count" in u for u in calls)
    assert any("format=geojson" in u for u in calls)


def test_fetch_mocked_paginates(monkeypatch, tmp_path):
    calls = _mock_comcat(monkeypatch, tmp_path, count_per_day=643)
    # 7-day window -> 7*643 = 4501 events > page_size 2000.
    f = fetch_earthquakes(BBOX, "2026-01-01", "2026-01-07",
                          page_size=2000, cache_dir=str(tmp_path))
    page_urls = [u for u in calls if "format=geojson" in u]
    assert len(page_urls) == 3  # 2000 + 2000 + 501
    assert "offset=1" in page_urls[0]
    assert "offset=2001" in page_urls[1]
    assert "offset=4001" in page_urls[2]
    assert f.provenance["n_events_parsed"] == 4501
    assert f.provenance["n_events_duplicates"] == 0
    assert len(f) == 4501


def test_fetch_offsets_are_one_based(monkeypatch, tmp_path):
    # Regression: the ComCat FDSN service rejects offset=0 with
    # HTTP 400 ("Valid values are 1 <= offset"). The adapter must
    # never request it — verified live against service v2.7.0.
    calls = _mock_comcat(monkeypatch, tmp_path, count_per_day=643)
    fetch_earthquakes(BBOX, "2026-01-01", "2026-01-07",
                      page_size=2000, cache_dir=str(tmp_path))
    page_urls = [u for u in calls if "format=geojson" in u]
    assert page_urls, "expected at least one paged query"
    for u in page_urls:
        q = dict(urllib.parse.parse_qsl(urllib.parse.urlsplit(u).query))
        assert int(q["offset"]) >= 1, f"offset must be 1-based: {u}"


def test_fetch_mocked_empty_field(monkeypatch, tmp_path):
    _mock_comcat(monkeypatch, tmp_path, fixed_count=0,
                 geojson=_json.dumps({"type": "FeatureCollection",
                                      "metadata": {}, "features": []}))
    f = fetch_earthquakes(BBOX, "2026-09-20", "2026-09-21",
                          cache_dir=str(tmp_path))
    assert len(f) == 0
    assert "empty_reason" in f.provenance


def test_fetch_validation(monkeypatch, tmp_path):
    with pytest.raises(ValueError):
        fetch_earthquakes(BBOX, "2026-09-28", "2026-09-20")
    with pytest.raises(ValueError):
        fetch_earthquakes((0, -95, 10, 10), "2026-09-20", "2026-09-21")
    with pytest.raises(ValueError):
        fetch_earthquakes(BBOX, "1899-01-01", "1900-01-02")
    with pytest.raises(ValueError):
        fetch_earthquakes(BBOX, "2026-09-20", "2026-09-21",
                          min_magnitude=-1.0)
    with pytest.raises(ValueError):
        fetch_earthquakes(BBOX, "2026-09-20", "2026-09-21", page_size=0)
    with pytest.raises(ValueError):
        fetch_earthquakes(BBOX, "2026-09-20", "2026-09-21",
                          page_size=20001)


def test_fetch_time_chunking_over_fdsn_ceiling(monkeypatch, tmp_path):
    # 8-day window at 3000/day = 24000 > 20000 ceiling -> splits into
    # two 4-day halves of 12000 each, paged at page_size 2000.
    calls = _mock_comcat(monkeypatch, tmp_path, count_per_day=3000)
    f = fetch_earthquakes(BBOX, "2026-01-01", "2026-01-08",
                          page_size=2000, cache_dir=str(tmp_path))
    count_urls = [u for u in calls if "/count" in u]
    assert len(count_urls) == 3  # top + 2 leaves
    assert f.provenance["n_events_parsed"] == 24000
    assert len(f) == 24000


def test_fetch_chunking_single_day_still_too_big(monkeypatch, tmp_path):
    _mock_comcat(monkeypatch, tmp_path, fixed_count=21001)
    with pytest.raises(RuntimeError, match="FDSN ceiling"):
        fetch_earthquakes(BBOX, "2026-01-01", "2026-01-01",
                          cache_dir=str(tmp_path))


def test_fetch_event_type_filter_in_url(monkeypatch, tmp_path):
    calls = _mock_comcat(monkeypatch, tmp_path, fixed_count=2,
                         geojson=_geojson_fixture())
    f = fetch_earthquakes(BBOX, "2026-09-20", "2026-09-21",
                          event_type="earthquake", cache_dir=str(tmp_path))
    assert f.event_type == "earthquake"
    assert any("eventtype=earthquake" in u for u in calls)


# ---------------------------------------------------------------------------
# Cache discipline
# ---------------------------------------------------------------------------

def test_cached_get_text_writes_and_hits(tmp_path):
    data = b'{"features": []}'
    url = "https://example.com/x"
    calls = {"n": 0}

    def fake_download(u, timeout=120):
        calls["n"] += 1
        return data

    orig = eq._download_bytes
    eq._download_bytes = fake_download
    try:
        t1, p1 = eq._cached_get_text(url, str(tmp_path), 7)
        t2, p2 = eq._cached_get_text(url, str(tmp_path), 7)
    finally:
        eq._download_bytes = orig
    assert t1 == t2 == data.decode()
    assert p1["cache_hit"] is False and p1["downloaded"] is True
    assert p2["cache_hit"] is True and p2["downloaded"] is False
    assert calls["n"] == 1
    assert os.path.isfile(os.path.join(
        str(tmp_path), eq._cache_key(url) + ".txt.sha256"))


def test_cached_get_text_corrupt_sidecar_redownloads(tmp_path):
    url = "https://example.com/y"
    key = eq._cache_key(url)
    p = os.path.join(str(tmp_path), key + ".txt")
    with open(p, "w") as fh:
        fh.write("corrupt")
    with open(p + ".sha256", "w") as fh:
        fh.write("deadbeef\n")
    orig = eq._download_bytes
    eq._download_bytes = lambda u, timeout=120: b"fresh"
    try:
        text, prov = eq._cached_get_text(url, str(tmp_path), 7)
    finally:
        eq._download_bytes = orig
    assert text == "fresh"
    assert prov["downloaded"] is True


# ---------------------------------------------------------------------------
# Live (skipped without SURVEY_CURRENTS_LIVE=1)
# ---------------------------------------------------------------------------

LIVE = os.environ.get("SURVEY_CURRENTS_LIVE") == "1"


@pytest.mark.skipif(not LIVE, reason="needs network")
def test_live_small_query():
    f = fetch_earthquakes((-125.0, 32.0, -114.0, 42.0),
                          "2026-09-20", "2026-09-27", min_magnitude=4.0)
    assert f.provenance["n_events"] >= 0
    assert all(e.magnitude is None or e.magnitude >= 4.0 - 1e-9
               for e in f.events)

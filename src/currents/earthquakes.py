"""USGS earthquake catalog acquisition: the ComCat FDSN event web service.

:func:`fetch_earthquakes` queries the keyless USGS ComCat FDSN event
service (GeoJSON) for a bbox x [start, end] window, paginates large
result sets, dedupes events by ID, and returns a :class:`QuakeField`:
one :class:`QuakeEvent` per cataloged event (time, lat, lon, depth,
magnitude, magnitude type, place, event type) plus provenance.

Access truth (verified live 2026-09-27 — see docs/DATA_SOURCES.md):

* ``https://earthquake.usgs.gov/fdsnws/event/1/query`` — keyless
  anonymous HTTPS (HTTP 200, no rate-limit headers observed on
  ordinary queries).
* ``format=geojson`` + ``starttime``/``endtime`` + the four bbox
  bounds + ``minmagnitude`` selects events; ``eventtype`` restricts to
  e.g. ``earthquake`` vs ``quarry blast`` / ``explosion``.
* ``https://earthquake.usgs.gov/fdsnws/event/1/count`` returns the
  plain-text total for the same parameters — the driver for
  pagination.
* Pagination: the FDSN default limit is 20000 events per response.
  Paged responses (``limit``/``offset``) swap the metadata ``count``
  key for ``limit``/``offset`` keys (verified live). When the count
  exceeds the page size, :func:`fetch_earthquakes` pages with
  ``limit``/``offset``; when it exceeds the 20000-event FDSN ceiling,
  it recursively splits the time window until every chunk is
  page-sized — the strategy is documented and tested with mocks.
* No Earthdata login, no API key, no signup: ComCat is fully keyless.

Catalog honesty: ComCat is a catalog of OBSERVED events, not a
forecast — there is no hazard model behind it. Magnitude of
completeness varies by region and time (roughly M4.5+ globally since
~1973; M2.5+ in the contiguous US since ~2013), so small and
historical events are under-recorded — this caveat is recorded in
provenance (``catalog_completeness_note``) and surfaced by the
survey-viz renderer. Event times are UTC; magnitudes and locations
are revised by analysts, so cache entries expire after
``max_cache_age_days`` (default 7).

Cache discipline (mirrors :mod:`currents.storms`): each request URL's
payload is cached as text under ``$SURVEY_CURRENTS_CACHE/earthquakes``
(else ``~/.cache/survey-currents/earthquakes``) with atomic writes and
a ``.sha256`` sidecar; a corrupt entry is re-downloaded; entries older
than ``max_cache_age_days`` (or any entry when ``refresh=True``) are
re-fetched. Repeated renders of the same query never re-hit the
network.

Downloads use stdlib ``urllib`` only. No heavy dependencies.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json as _json
import os
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from dataclasses import field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date
from .sst_global import validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]

__all__ = [
    "COMCAT_EVENT_BASE",
    "COMCAT_COUNT_BASE",
    "COMCAT_RECORD_START",
    "COMCAT_FDSN_LIMIT",
    "COMCAT_PAGE_SIZE",
    "COMCAT_MAX_CACHE_AGE_DAYS",
    "COMCAT_DEPTH_BINS_KM",
    "QuakeEvent",
    "QuakeField",
    "comcat_cache_dir",
    "event_query_url",
    "event_count_url",
    "fetch_earthquakes",
]

#: Keyless USGS ComCat FDSN event service (verified live 2026-09-27).
COMCAT_EVENT_BASE = "https://earthquake.usgs.gov/fdsnws/event/1/query"
#: Plain-text total-count endpoint, same parameter set.
COMCAT_COUNT_BASE = "https://earthquake.usgs.gov/fdsnws/event/1/count"
#: Earliest ComCat date accepted (the catalog reaches ~1900; earlier
#: dates are rejected as absurd rather than queried).
COMCAT_RECORD_START = _dt.date(1900, 1, 1)
#: FDSN default response ceiling (events per response).
COMCAT_FDSN_LIMIT = 20000
#: Default page size for paginated queries (keeps each response and
#: each cache entry small; well under the FDSN ceiling).
COMCAT_PAGE_SIZE = 2000
#: Default cache freshness (days): magnitudes and locations are
#: revised by analysts, so catalog queries revalidate sooner than the
#: 30-day discipline used for static archives.
COMCAT_MAX_CACHE_AGE_DAYS = 7
#: Documented depth bins (km) shared with the survey-viz renderer for
#: marker coloring: shallow < 70 km, intermediate 70–300 km,
#: deep > 300 km (Wadati–Benioff convention).
COMCAT_DEPTH_BINS_KM: Tuple[Tuple[str, float, float], ...] = (
    ("shallow", 0.0, 70.0),
    ("intermediate", 70.0, 300.0),
    ("deep", 300.0, float("inf")),
)


def comcat_cache_dir() -> str:
    """User cache root for ComCat downloads.

    ``$SURVEY_CURRENTS_CACHE/earthquakes`` when set, else
    ``~/.cache/survey-currents/earthquakes``. Created on demand.
    """
    root = os.environ.get("SURVEY_CURRENTS_CACHE") or os.path.join(
        os.path.expanduser("~"), ".cache", "survey-currents")
    path = os.path.join(root, "earthquakes")
    os.makedirs(path, exist_ok=True)
    return path


def _common_params(bbox: Sequence[float], start: DateLike, end: DateLike,
                   min_magnitude: float,
                   event_type: Optional[str]) -> Dict[str, str]:
    """Shared FDSN query parameters (bbox, time, magnitude, event type)."""
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    d0, d1 = _coerce_date(start), _coerce_date(end)
    params = {
        "starttime": d0.isoformat(),
        "endtime": d1.isoformat(),
        "minlatitude": f"{lat_min:.6f}",
        "maxlatitude": f"{lat_max:.6f}",
        "minlongitude": f"{lon_min:.6f}",
        "maxlongitude": f"{lon_max:.6f}",
        "minmagnitude": f"{float(min_magnitude):.2f}",
    }
    if event_type:
        params["eventtype"] = str(event_type)
    return params


def event_query_url(bbox: Sequence[float], start: DateLike, end: DateLike,
                    min_magnitude: float = 0.0,
                    event_type: Optional[str] = None,
                    limit: Optional[int] = None,
                    offset: Optional[int] = None) -> str:
    """Keyless ComCat FDSN GeoJSON query URL (pure, offline-testable)."""
    params = _common_params(bbox, start, end, min_magnitude, event_type)
    params["format"] = "geojson"
    if limit is not None:
        params["limit"] = str(int(limit))
    if offset is not None:
        params["offset"] = str(int(offset))
    return COMCAT_EVENT_BASE + "?" + urllib.parse.urlencode(params)


def event_count_url(bbox: Sequence[float], start: DateLike, end: DateLike,
                    min_magnitude: float = 0.0,
                    event_type: Optional[str] = None) -> str:
    """Keyless ComCat plain-text total-count URL (pure, offline-testable)."""
    params = _common_params(bbox, start, end, min_magnitude, event_type)
    return COMCAT_COUNT_BASE + "?" + urllib.parse.urlencode(params)


def _tool_version() -> str:
    from . import __version__
    return __version__


# ---------------------------------------------------------------------------
# Download + cache (same discipline as currents.storms)
# ---------------------------------------------------------------------------

def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _cache_key(url: str) -> str:
    return hashlib.sha256(url.encode("utf-8")).hexdigest()[:32]


def _cached_get_text(url: str, cache_dir: str,
                     max_age_days: int, refresh: bool = False
                     ) -> Tuple[str, Dict[str, Any]]:
    """GET ``url`` with a SHA-256-verified text cache.

    Returns ``(text, provenance)`` where provenance carries ``url``,
    ``sha256``, ``cache_hit``, and ``downloaded``. Cache entries are
    written atomically with a ``.sha256`` sidecar; a corrupt entry is
    re-downloaded; entries older than ``max_age_days`` (or any entry
    when ``refresh``) are re-fetched. Mirrors the discipline of
    :func:`currents.streamgages._cached_get_text`.
    """
    os.makedirs(cache_dir, exist_ok=True)
    key = _cache_key(url)
    path = os.path.join(cache_dir, f"{key}.txt")
    sidecar = path + ".sha256"
    prov: Dict[str, Any] = {"url": url}

    def _valid_entry() -> Optional[str]:
        if not (os.path.isfile(path) and os.path.isfile(sidecar)):
            return None
        try:
            with open(sidecar, "r", encoding="utf-8") as fh:
                expected = fh.read().strip().split()[0]
            age_days = ((_dt.datetime.now(_dt.timezone.utc).timestamp()
                         - os.path.getmtime(path)) / 86400.0)
            if not refresh and age_days <= max_age_days:
                with open(path, "r", encoding="utf-8") as fh:
                    text = fh.read()
                if _sha256_bytes(text.encode("utf-8")) == expected:
                    prov.update(cache_hit=True, downloaded=False,
                                sha256=expected)
                    return text
        except (OSError, ValueError, IndexError):
            pass
        return None

    cached = _valid_entry()
    if cached is not None:
        return cached, prov

    data = _download_bytes(url, timeout=120)
    digest = _sha256_bytes(data)
    fd, tmp = tempfile.mkstemp(dir=cache_dir, suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        with open(tmp + ".sha256", "w", encoding="utf-8") as fh:
            fh.write(digest + "\n")
        os.replace(tmp, path)
        os.replace(tmp + ".sha256", sidecar)
    finally:
        for leftover in (tmp, tmp + ".sha256"):
            if os.path.exists(leftover):
                os.remove(leftover)
    prov.update(cache_hit=False, downloaded=True, sha256=digest)
    return data.decode("utf-8", errors="replace"), prov


def _download_bytes(url: str, timeout: int = 120) -> bytes:
    """Download ``url`` (mockable seam for the offline test suite)."""
    req = urllib.request.Request(url, headers={"User-Agent": "survey-currents"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"ComCat download failed for {url}: {exc}. The ComCat FDSN "
            "event service is keyless HTTPS — this is usually a "
            "connectivity issue.") from exc

# ---------------------------------------------------------------------------
# Response parsing (pure, offline-testable)
# ---------------------------------------------------------------------------

def _epoch_ms_to_utc(ms: Any) -> Optional[_dt.datetime]:
    """ComCat epoch-milliseconds -> aware UTC datetime (None when absent)."""
    if ms is None:
        return None
    try:
        return _dt.datetime.fromtimestamp(
            float(ms) / 1000.0, tz=_dt.timezone.utc)
    except (TypeError, ValueError, OverflowError, OSError):
        return None


def _parse_event_feature(feature: Dict[str, Any]
                         ) -> Optional[Dict[str, Any]]:
    """Parse one ComCat GeoJSON feature into an event dict.

    Returns None for malformed features (no id, no usable geometry or
    time) — the caller counts them honestly in provenance. Missing
    magnitudes stay None (never 0.0 — a fabricated magnitude would be
    a lie); 2-D coordinates (no depth) become depth ``None``.
    """
    if not isinstance(feature, dict):
        return None
    event_id = feature.get("id")
    props = feature.get("properties") or {}
    geom = feature.get("geometry") or {}
    coords = geom.get("coordinates") if isinstance(geom, dict) else None
    if not event_id or not isinstance(coords, (list, tuple)) or len(coords) < 2:
        return None
    try:
        lon = float(coords[0])
        lat = float(coords[1])
    except (TypeError, ValueError):
        return None
    depth_km: Optional[float] = None
    if len(coords) >= 3:
        try:
            depth_km = float(coords[2])
        except (TypeError, ValueError):
            depth_km = None
    when = _epoch_ms_to_utc(props.get("time"))
    if when is None:
        return None
    mag: Optional[float] = None
    raw_mag = props.get("mag")
    if raw_mag is not None:
        try:
            mag = float(raw_mag)
        except (TypeError, ValueError):
            mag = None
    mag_type = props.get("magType")
    return {
        "event_id": str(event_id),
        "time": when,
        "lat": lat,
        "lon": lon,
        "depth_km": depth_km,
        "magnitude": mag,
        "mag_type": str(mag_type) if mag_type is not None else None,
        "place": str(props.get("place") or ""),
        "event_type": str(props.get("type") or "unknown"),
    }


def _parse_event_geojson(payload: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], int]:
    """Parse a ComCat GeoJSON body.

    Returns ``(events, n_malformed)``: parsed event dicts (each via
    :func:`_parse_event_feature`) and the count of skipped malformed
    features — never silently dropped, recorded in provenance.
    """
    events: List[Dict[str, Any]] = []
    n_malformed = 0
    features = payload.get("features") if isinstance(payload, dict) else None
    for feature in (features or []):
        parsed = _parse_event_feature(feature)
        if parsed is None:
            n_malformed += 1
        else:
            events.append(parsed)
    return events, n_malformed


def _dedupe_events(events: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], int]:
    """Dedupe events by event_id (first occurrence wins).

    Returns ``(deduped, n_duplicates)``. Deterministic: input order is
    preserved (callers sort by time before deduping).
    """
    seen: set = set()
    out: List[Dict[str, Any]] = []
    n_dup = 0
    for ev in events:
        eid = ev["event_id"]
        if eid in seen:
            n_dup += 1
            continue
        seen.add(eid)
        out.append(ev)
    return out, n_dup


# ---------------------------------------------------------------------------
# Field model
# ---------------------------------------------------------------------------

@dataclass
class QuakeEvent:
    """One ComCat cataloged event.

    ``magnitude``/``depth_km`` are None when ComCat did not report them
    (never 0.0). ``event_type`` is the ComCat type (``"earthquake"``,
    ``"quarry blast"``, ``"explosion"``, …) — non-tectonic events are
    kept and labeled honestly, never relabeled.
    """

    event_id: str
    time: _dt.datetime
    lat: float
    lon: float
    depth_km: Optional[float] = None
    magnitude: Optional[float] = None
    mag_type: Optional[str] = None
    place: str = ""
    event_type: str = "earthquake"

    def __post_init__(self) -> None:
        self.event_id = str(self.event_id)
        if isinstance(self.time, str):
            self.time = _dt.datetime.fromisoformat(self.time)
        if self.time.tzinfo is None:
            self.time = self.time.replace(tzinfo=_dt.timezone.utc)
        self.lat = float(self.lat)
        self.lon = float(self.lon)
        if self.depth_km is not None:
            self.depth_km = float(self.depth_km)
        if self.magnitude is not None:
            self.magnitude = float(self.magnitude)
        self.place = str(self.place or "")
        self.event_type = str(self.event_type or "unknown")

    @property
    def date(self) -> _dt.date:
        """UTC calendar date of the event."""
        return self.time.date()

    def depth_bin(self) -> str:
        """Wadati–Benioff depth bin shared with the renderer."""
        if self.depth_km is None:
            return "unknown"
        for name, lo, hi in COMCAT_DEPTH_BINS_KM:
            if lo <= self.depth_km < hi:
                return name
        return "unknown"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "event_id": self.event_id,
            "time": self.time.isoformat(),
            "lat": self.lat,
            "lon": self.lon,
            "depth_km": self.depth_km,
            "magnitude": self.magnitude,
            "mag_type": self.mag_type,
            "place": self.place,
            "event_type": self.event_type,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "QuakeEvent":
        return cls(
            event_id=data["event_id"],
            time=_dt.datetime.fromisoformat(data["time"]),
            lat=data["lat"], lon=data["lon"],
            depth_km=data.get("depth_km"),
            magnitude=data.get("magnitude"),
            mag_type=data.get("mag_type"),
            place=data.get("place", ""),
            event_type=data.get("event_type", "unknown"))


@dataclass
class QuakeField:
    """A ComCat event catalog for one query.

    ``provenance`` records the exact request URLs, the count URL,
    retrieval time, event counts (queried / parsed / malformed /
    duplicates / returned), the applied ``min_magnitude`` / event-type
    filter, and the catalog-completeness caveat — following the
    :mod:`currents.storms` conventions.
    """

    events: List[QuakeEvent]
    bbox: Tuple[float, float, float, float]
    start: _dt.date
    end: _dt.date
    min_magnitude: float = 0.0
    event_type: Optional[str] = None
    source: str = "usgs"
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        self.bbox = validate_sst_bbox(self.bbox)
        self.start = _coerce_date(self.start)
        self.end = _coerce_date(self.end)
        if self.start > self.end:
            raise ValueError(
                f"QuakeField: start {self.start} is after end {self.end}")
        self.events = sorted(
            list(self.events), key=lambda e: (e.time, e.event_id))
        self.min_magnitude = float(self.min_magnitude)
        self.source = str(self.source or "usgs")

    def __len__(self) -> int:
        return len(self.events)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    # -- accessors -----------------------------------------------------

    def largest(self, n: int = 5) -> List[QuakeEvent]:
        """The ``n`` largest events by magnitude (descending).

        Events with unreported magnitudes sort last. Ties break by
        earlier time — deterministic.
        """
        return sorted(
            self.events,
            key=lambda e: (
                -1.0 if e.magnitude is None else -e.magnitude,
                e.time, e.event_id))[:max(0, int(n))]

    def select_time(self, day: DateLike) -> List[QuakeEvent]:
        """Events whose UTC date equals ``day`` (time-ordered)."""
        d = _coerce_date(day)
        return [e for e in self.events if e.date == d]

    def counts_by_day(self) -> Tuple[List[_dt.date], np.ndarray]:
        """(dates, counts): per-day event counts over [start, end]."""
        ndays = (self.end - self.start).days + 1
        dates = [self.start + _dt.timedelta(days=k) for k in range(ndays)]
        counts = np.zeros(ndays, dtype=int)
        for e in self.events:
            k = (e.date - self.start).days
            if 0 <= k < ndays:
                counts[k] += 1
        return dates, counts

    def magnitudes(self) -> np.ndarray:
        """Magnitudes as a float array (unreported -> NaN)."""
        return np.asarray(
            [e.magnitude if e.magnitude is not None else np.nan
             for e in self.events], dtype=float)

    # -- serialization ---------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "quake_events": [e.to_dict() for e in self.events],
            "bbox": list(self.bbox),
            "start": self.start.isoformat(),
            "end": self.end.isoformat(),
            "min_magnitude": self.min_magnitude,
            "event_type": self.event_type,
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "QuakeField":
        return cls(
            events=[QuakeEvent.from_dict(e)
                    for e in data.get("quake_events", [])],
            bbox=tuple(data["bbox"]),
            start=_coerce_date(data["start"]),
            end=_coerce_date(data["end"]),
            min_magnitude=float(data.get("min_magnitude", 0.0)),
            event_type=data.get("event_type"),
            source=data.get("source", "usgs"),
            provenance=dict(data.get("provenance", {})))

    def to_json(self, path: str) -> str:
        """Write the field as JSON; return ``path``."""
        with open(path, "w", encoding="utf-8") as fh:
            _json.dump(self.to_dict(), fh, indent=2)
        return path

    @classmethod
    def from_json(cls, path: str) -> "QuakeField":
        """Read a field written by :meth:`to_json`."""
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(_json.load(fh))

    # -- synthetic fixture -------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-140.0, 30.0, -110.0, 50.0),
                  start: DateLike = "2024-01-01",
                  end: DateLike = "2024-01-31",
                  n_events: int = 40, seed: int = 13,
                  source: str = "synthetic") -> "QuakeField":
        """Deterministic synthetic events (offline tests / demos).

        Magnitudes follow a Gutenberg–Richter-ish exponential
        distribution (b=1), depths cluster shallow with a deep tail, and
        times spread across the window. Two events share a place name
        so duplicate-place handling stays testable; magnitudes are
        always reported on synthetic events.
        """
        rng = np.random.default_rng(seed)
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        d0, d1 = _coerce_date(start), _coerce_date(end)
        if d1 < d0:
            raise ValueError(f"end {d1} is before start {d0}")
        ndays = (d1 - d0).days + 1
        places = ["SYNTHETIC FAULT ZONE", "MOCK TRENCH", "DEMO RIDGE",
                  "FAKE VOLCANIC ARC", "PLACEBO BASIN"]
        events: List[QuakeEvent] = []
        for k in range(n_events):
            day_off = int(rng.integers(0, ndays))
            sec_off = int(rng.integers(0, 86400))
            when = (_dt.datetime(d0.year, d0.month, d0.day,
                                 tzinfo=_dt.timezone.utc)
                    + _dt.timedelta(days=day_off, seconds=sec_off))
            # b=1: P(M > m) = 10^-(m - m_min).
            mag = 2.5 - float(np.log10(rng.random())) 
            depth = float(abs(rng.normal(25.0, 30.0)))
            if rng.random() < 0.08:
                depth = float(rng.uniform(300.0, 650.0))
            events.append(QuakeEvent(
                event_id=f"synth{k:05d}",
                time=when,
                lat=float(lat_min + (lat_max - lat_min) * rng.random()),
                lon=float(lon_min + (lon_max - lon_min) * rng.random()),
                depth_km=depth,
                magnitude=round(mag, 1),
                mag_type="ml",
                place=f"{places[k % len(places)]} — synthetic event {k}",
                event_type="earthquake"))
        return cls(events=events,
                   bbox=(lon_min, lat_min, lon_max, lat_max),
                   start=d0, end=d1, min_magnitude=0.0,
                   event_type=None, source=source,
                   provenance={"synthetic": True, "seed": seed})

# ---------------------------------------------------------------------------
# Public fetch
# ---------------------------------------------------------------------------

def _parse_count(text: str) -> int:
    """Parse the count-endpoint body (plain-text integer)."""
    try:
        total = int(text.strip().split()[0])
    except (ValueError, IndexError) as exc:
        raise RuntimeError(
            f"ComCat count endpoint returned an unparsable body: "
            f"{text[:80]!r}") from exc
    if total < 0:
        raise RuntimeError(
            f"ComCat count endpoint returned a negative total: {total}")
    return total


def _fetch_window(bbox: Sequence[float], d0: _dt.date, d1: _dt.date,
                  min_magnitude: float, event_type: Optional[str],
                  page_size: int, cdir: str, max_cache_age_days: int,
                  refresh: bool) -> Tuple[List[Dict[str, Any]], int,
                                          List[str], List[Dict[str, Any]]]:
    """Fetch one time window: page with limit/offset when needed.

    Returns ``(events, n_malformed, request_urls, page_provs)``.
    Recursively splits the window when the count exceeds the FDSN
    20000-event ceiling (each sub-window is counted and fetched
    independently, then merged) — a giant query never fails on the
    ceiling and never silently truncates.
    """
    count_url = event_count_url(bbox, d0, d1, min_magnitude, event_type)
    count_text, count_prov = _cached_get_text(
        count_url, cdir, max_cache_age_days, refresh=refresh)
    total = _parse_count(count_text)
    request_urls = [count_url]
    page_provs: List[Dict[str, Any]] = [count_prov]
    n_malformed = 0

    if total == 0:
        return [], 0, request_urls, page_provs

    if total > COMCAT_FDSN_LIMIT:
        # Recurse: split the window in half by day count (min 1-day
        # leaves) until every chunk is page-sized. Dedupe at merge so
        # boundary events are never double-counted.
        ndays = (d1 - d0).days + 1
        if ndays <= 1:
            raise RuntimeError(
                f"ComCat query for {d0.isoformat()} x bbox still "
                f"exceeds the FDSN ceiling ({total} > "
                f"{COMCAT_FDSN_LIMIT}) on a single day — raise "
                "min_magnitude or subdivide the bbox.")
        mid = d0 + _dt.timedelta(days=ndays // 2 - 1)
        left = _fetch_window(bbox, d0, mid, min_magnitude, event_type,
                             page_size, cdir, max_cache_age_days, refresh)
        right = _fetch_window(bbox, mid + _dt.timedelta(days=1), d1,
                              min_magnitude, event_type, page_size,
                              cdir, max_cache_age_days, refresh)
        merged, _ = _dedupe_events(left[0] + right[0])
        return (merged, left[1] + right[1],
                request_urls + left[2] + right[2],
                page_provs + left[3] + right[3])

    events: List[Dict[str, Any]] = []
    offset = 0
    while True:
        qurl = event_query_url(bbox, d0, d1, min_magnitude, event_type,
                               limit=min(page_size, COMCAT_FDSN_LIMIT),
                               offset=offset)
        request_urls.append(qurl)
        text, prov = _cached_get_text(
            qurl, cdir, max_cache_age_days, refresh=refresh)
        page_provs.append(prov)
        payload = _json.loads(text)
        page_events, malformed = _parse_event_geojson(payload)
        n_malformed += malformed
        events.extend(page_events)
        # Paged metadata carries limit/offset, not count (verified live
        # 2026-09-27), so stop on a short page — never on a count key.
        if len(page_events) < min(page_size, COMCAT_FDSN_LIMIT):
            break
        offset += len(page_events)
        if offset >= total:
            break
    return events, n_malformed, request_urls, page_provs


def fetch_earthquakes(bbox: Sequence[float], start: DateLike, end: DateLike,
                      min_magnitude: float = 0.0,
                      event_type: Optional[str] = None,
                      page_size: int = COMCAT_PAGE_SIZE,
                      cache_dir: Optional[str] = None,
                      max_cache_age_days: int = COMCAT_MAX_CACHE_AGE_DAYS,
                      refresh: bool = False) -> QuakeField:
    """Fetch the USGS ComCat earthquake catalog for ``bbox`` x ``[start, end]``.

    Queries the keyless ComCat FDSN event service (GeoJSON), counts
    first, then pages with ``limit``/``offset``; windows whose count
    exceeds the 20000-event FDSN ceiling are recursively split by time
    until every chunk is page-sized. Every response is cached
    (SHA-256-verified, atomic writes, ``max_cache_age_days``
    revalidation) so repeated renders never re-hit the network. Events
    are deduped by event ID and sorted by time.

    Args:
        min_magnitude: ComCat ``minmagnitude`` filter (the catalog's
            magnitude of completeness varies by region/time — see the
            provenance note).
        event_type: ComCat ``eventtype`` filter (e.g.
            ``"earthquake"``, ``"quarry blast"``); None = all types.
            Non-tectonic types are kept and labeled honestly, never
            relabeled.
        page_size: events per paged request (must be >= 1 and <= 20000).

    Event times are UTC datetimes; missing magnitudes stay None, never
    0.0; missing depths stay None. ComCat is a catalog of observed
    events, NOT a forecast — the provenance records that explicitly.
    """
    box = validate_sst_bbox(bbox)
    d0, d1 = _coerce_date(start), _coerce_date(end)
    if d1 < d0:
        raise ValueError(f"end {d1} is before start {d0}")
    if d0 < COMCAT_RECORD_START:
        raise ValueError(
            f"start {d0} predates the ComCat record floor "
            f"({COMCAT_RECORD_START})")
    if d1 > _dt.date.today() + _dt.timedelta(days=2):
        raise ValueError(f"end {d1.isoformat()} is in the future")
    min_mag = float(min_magnitude)
    if min_mag < 0.0:
        raise ValueError(f"min_magnitude must be >= 0, got {min_magnitude}")
    if event_type is not None and not str(event_type).strip():
        raise ValueError("event_type must be a non-empty string or None")
    if not (1 <= int(page_size) <= COMCAT_FDSN_LIMIT):
        raise ValueError(
            f"page_size must be in [1, {COMCAT_FDSN_LIMIT}], "
            f"got {page_size}")

    cdir = cache_dir or comcat_cache_dir()
    events, n_malformed, request_urls, page_provs = _fetch_window(
        box, d0, d1, min_mag, event_type, int(page_size),
        cdir, max_cache_age_days, refresh=refresh)
    # Sort by time (the service returns time-descending by default;
    # ascending is the deterministic field order), then dedupe.
    events.sort(key=lambda e: (e["time"], e["event_id"]))
    deduped, n_duplicates = _dedupe_events(events)

    quakes = [QuakeEvent(**e) for e in deduped]
    provenance = {
        "source": "usgs",
        "product": "USGS Earthquake Catalog (ComCat), FDSN event service",
        "service": COMCAT_EVENT_BASE,
        "count_service": COMCAT_COUNT_BASE,
        "request_urls": request_urls,
        "retrieved_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "bbox": list(box),
        "start": d0.isoformat(),
        "end": d1.isoformat(),
        "min_magnitude": min_mag,
        "event_type": event_type,
        "page_size": int(page_size),
        "n_events_queried": _parse_count_safe(request_urls, cdir),
        "n_events_parsed": len(events),
        "n_events_malformed": n_malformed,
        "n_events_duplicates": n_duplicates,
        "n_events": len(quakes),
        "catalog_completeness_note": (
            "ComCat is a catalog of OBSERVED events, not a forecast — "
            "there is no hazard model behind it. Magnitude of "
            "completeness varies by region and time: roughly M4.5+ "
            "globally since ~1973 and M2.5+ in the contiguous US since "
            "~2013; small and historical events are under-recorded. "
            "Non-tectonic event types (quarry blast, explosion, …) are "
            "included unless event_type filters them, and are labeled "
            "honestly. Magnitudes and locations are revised by "
            "analysts; cache entries revalidate after "
            f"{max_cache_age_days} days."),
        "cache": {
            "pages": [{k: p.get(k)
                       for k in ("sha256", "cache_hit", "downloaded")}
                      for p in page_provs],
        },
        "tool": f"survey-currents {_tool_version()}",
    }
    if not quakes:
        provenance["empty_reason"] = (
            "no ComCat events matched bbox "
            f"{list(box)} for {d0.isoformat()}…{d1.isoformat()} "
            f"(min_magnitude={min_mag}, event_type={event_type}). An "
            "empty catalog is a legitimate observation — quiet regions "
            "and short windows yield no events.")

    return QuakeField(events=quakes, bbox=box, start=d0, end=d1,
                      min_magnitude=min_mag, event_type=event_type,
                      source="usgs", provenance=provenance)


def _parse_count_safe(request_urls: List[str], cache_dir: str) -> Optional[int]:
    """Recover the count-endpoint total from the cached count body.

    The first request URL is always the count URL (see
    :func:`_fetch_window`); re-parse its cached payload so provenance
    records what the service claimed. None when unreadable.
    """
    if not request_urls:
        return None
    path = os.path.join(cache_dir, _cache_key(request_urls[0]) + ".txt")
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return _parse_count(fh.read())
    except (OSError, RuntimeError):
        return None

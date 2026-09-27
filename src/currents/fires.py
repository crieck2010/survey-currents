"""NASA FIRMS active-fire acquisition: detections -> FireField -> density grids.

Active-fire detections are POINTS (lat, lon, time, brightness, FRP,
confidence, instrument), not a gridded field. This module:

* :func:`fetch_firms` — downloads detections from the NASA FIRMS area
  API (verified live 2026-09-26; see docs/DATA_SOURCES.md), loops the
  requested range in <= 5-day windows per instrument, and returns a
  :class:`FireField`.
* :class:`FireField` — the canonical detection model (times/lats/lons
  plus brightness, FRP, confidence, satellite, instrument, day/night),
  with provenance, JSON round-trip, ``select_time``/``select_bbox``
  filters, a deterministic ``synthetic()`` fixture, and
  :meth:`FireField.to_density_grid`, which bins detections into daily
  fire-count (or FRP-weighted) grids on a regular lat/lon grid. The
  density grid is a plain ``times``/``lats``/``lons``/``values`` dict —
  the exact form ``survey-viz``'s renderer consumes — so fire maps
  render with zero renderer changes.

FIRMS API shape (verified live 2026-09-26):

* Endpoint: ``https://firms.modaps.eosdis.nasa.gov/api/area/csv/
  {MAP_KEY}/{PRODUCT}/{W},{S},{E},{N}/{DAY_RANGE}/{DATE}``
* ``DAY_RANGE`` is 1-5 days when a ``DATE`` window start is given
  (1-10 without it, for the most recent data only).
* A free ``MAP_KEY`` is REQUIRED — request one at
  https://firms.modaps.eosdis.nasa.gov/api/map_key/ and export it as
  ``FIRMS_MAP_KEY``. Without it a :class:`CredentialsMissing` error
  explains the setup (same pattern as the MUR/ERA5 adapters).
* Products are latency-tiered per instrument: ``*_NRT`` (near real-time,
  recent dates) and ``*_SP`` (standard processing, older dates), picked
  per date by the pure, offline-testable :func:`firms_product_for`
  rule. Detections are returned oldest-first per request.

Downloads use stdlib ``urllib`` only (via :func:`glsea._download_bytes`).
The ``MAP_KEY`` is read from the environment (or an explicit argument)
and is NEVER stored in provenance records — stored URLs carry a
``<redacted>`` placeholder instead.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import os
from dataclasses import dataclass, field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date, _download_bytes
from .sst_global import CredentialsMissing as _BaseCredentialsMissing

DateLike = Union[_dt.date, _dt.datetime, str]


# ---------------------------------------------------------------------------
# FIRMS constants (verified live against the FIRMS API 2026-09-26)
# ---------------------------------------------------------------------------

#: Base URL for the FIRMS area API.
FIRMS_API_BASE = "https://firms.modaps.eosdis.nasa.gov/api"
#: Area-API path (CSV format); the full request is
#: ``{FIRMS_API_BASE}/area/csv/{MAP_KEY}/{PRODUCT}/{W},{S},{E},{N}/{DAY_RANGE}/{DATE}``.
FIRMS_AREA_PATH = "/area/csv"
#: Environment variable carrying the free FIRMS MAP_KEY.
FIRMS_MAP_KEY_ENV = "FIRMS_MAP_KEY"
#: Where to request a free key (instant, no approval).
FIRMS_MAP_KEY_URL = "https://firms.modaps.eosdis.nasa.gov/api/map_key/"
#: Max DAY_RANGE accepted by the API when a DATE window start is given.
FIRMS_MAX_DAY_RANGE = 5

#: Per-instrument product families and (approximate) record starts.
#: NRT = near real-time (recent dates); SP = standard processing
#: (older dates). Record starts are the operational start of each
#: product family, good to month precision.
FIRMS_INSTRUMENTS: Dict[str, Dict[str, Any]] = {
    "VIIRS_SNPP": {
        "nrt": "VIIRS_SNPP_NRT",
        "sp": "VIIRS_SNPP_SP",
        "start": _dt.date(2012, 1, 1),
        "resolution_m": 375,
    },
    "VIIRS_NOAA20": {
        "nrt": "VIIRS_NOAA20_NRT",
        "sp": "VIIRS_NOAA20_SP",
        "start": _dt.date(2018, 1, 1),
        "resolution_m": 375,
    },
    "VIIRS_NOAA21": {
        "nrt": "VIIRS_NOAA21_NRT",
        "sp": "VIIRS_NOAA21_SP",
        "start": _dt.date(2023, 3, 1),
        "resolution_m": 375,
    },
    "MODIS": {
        "nrt": "MODIS_NRT",
        "sp": "MODIS_SP",
        "start": _dt.date(2000, 11, 1),
        "resolution_m": 1000,
    },
}
#: Earliest date any supported instrument covers (MODIS/Terra).
FIRMS_START = _dt.date(2000, 11, 1)

#: Dates newer than this many days before "today" use the NRT product;
#: older dates use the SP product. FIRMS transitions NRT -> SP on an
#: operational schedule (~2-3 months), so this is a documented heuristic,
#: not an exact boundary — callers can override per request with
#: ``products=`` on :func:`fetch_firms`.
FIRMS_NRT_WINDOW_DAYS = 60


class CredentialsMissing(_BaseCredentialsMissing):
    """A FIRMS MAP_KEY is required but was not found.

    Active-fire downloads need a free NASA FIRMS MAP_KEY:

    1. Request one (instant, free) at
       https://firms.modaps.eosdis.nasa.gov/api/map_key/
    2. Export it (best for scripts / CI)::

           export FIRMS_MAP_KEY="your-map-key-here"

       — or pass ``map_key=...`` directly to :func:`fetch_firms`.

    The key is never written to logs, provenance sidecars, or stored
    URLs (they are redacted there).
    """


# ---------------------------------------------------------------------------
# Pure, offline-testable helpers
# ---------------------------------------------------------------------------

def firms_map_key(explicit: Optional[str] = None) -> str:
    """Return the FIRMS MAP_KEY (explicit arg wins, else the env var).

    Raises:
        CredentialsMissing: no key found anywhere.
    """
    key = (explicit or "").strip() or os.environ.get(FIRMS_MAP_KEY_ENV, "").strip()
    if not key:
        raise CredentialsMissing(
            "No FIRMS MAP_KEY found: set the "
            f"{FIRMS_MAP_KEY_ENV} environment variable (free key at "
            f"{FIRMS_MAP_KEY_URL}) or pass map_key=... to fetch_firms()."
        )
    return key


def validate_firms_bbox(bbox: Sequence[float]) -> Tuple[float, float, float, float]:
    """Validate a (lon_min, lat_min, lon_max, lat_max) bbox.

    Antimeridian-crossing boxes (lon_min > lon_max) are allowed and are
    split into two requests by :func:`firms_area_windows`.
    """
    try:
        lon_min, lat_min, lon_max, lat_max = (float(x) for x in bbox)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"bbox must be 4 numbers, got {bbox!r}") from exc
    if not (lat_min < lat_max and -90.0 <= lat_min <= 90.0
            and -90.0 <= lat_max <= 90.0):
        raise ValueError(
            f"bbox latitudes must satisfy -90 <= lat_min < lat_max <= 90, "
            f"got {bbox!r}")
    if not (-180.0 <= lon_min <= 180.0 and -180.0 <= lon_max <= 180.0):
        raise ValueError(
            f"bbox longitudes must be within [-180, 180], got {bbox!r}")
    if lon_min == lon_max:
        raise ValueError(f"bbox has zero longitude span: {bbox!r}")
    return lon_min, lat_min, lon_max, lat_max


def firms_area_windows(bbox: Sequence[float]) -> List[Tuple[float, float, float, float]]:
    """Split ``bbox`` into FIRMS area-API windows.

    The area API takes a single W,S,E,N box, so antimeridian-crossing
    boxes (lon_min > lon_max) become two windows (lon_min..180 and
    -180..lon_max). Anything else passes through unchanged.
    """
    lon_min, lat_min, lon_max, lat_max = validate_firms_bbox(bbox)
    if lon_min < lon_max:
        return [(lon_min, lat_min, lon_max, lat_max)]
    return [(lon_min, lat_min, 180.0, lat_max),
            (-180.0, lat_min, lon_max, lat_max)]


def firms_product_for(instrument: str, day: _dt.date,
                      today: Optional[_dt.date] = None) -> str:
    """Pick the FIRMS product name for ``instrument`` on ``day``.

    Pure and offline-testable: dates within
    :data:`FIRMS_NRT_WINDOW_DAYS` of ``today`` use the ``*_NRT`` product,
    older dates the ``*_SP`` product. ``today`` defaults to the current
    UTC date; pass it explicitly in tests.
    """
    key = str(instrument).strip().upper()
    if key not in FIRMS_INSTRUMENTS:
        raise ValueError(
            f"unknown FIRMS instrument {instrument!r}; choose from "
            f"{sorted(FIRMS_INSTRUMENTS)}")
    fam = FIRMS_INSTRUMENTS[key]
    if day < fam["start"]:
        raise ValueError(
            f"{key} has no data before {fam['start']} (requested {day})")
    ref = today or _dt.datetime.now(_dt.timezone.utc).date()
    if (ref - day).days <= FIRMS_NRT_WINDOW_DAYS:
        return str(fam["nrt"])
    return str(fam["sp"])


def firms_window_chunks(start: _dt.date, end: _dt.date,
                        max_days: int = FIRMS_MAX_DAY_RANGE
                        ) -> List[Tuple[_dt.date, _dt.date]]:
    """Split ``[start, end]`` into contiguous windows of at most ``max_days``.

    The area API accepts DAY_RANGE 1-5 with a DATE window start; longer
    ranges must be looped. Returns ``(window_start, window_end)``
    inclusive pairs.
    """
    if end < start:
        raise ValueError(f"end {end} is before start {start}")
    if max_days < 1:
        raise ValueError(f"max_days must be >= 1, got {max_days}")
    chunks: List[Tuple[_dt.date, _dt.date]] = []
    cur = start
    while cur <= end:
        nxt = min(cur + _dt.timedelta(days=max_days - 1), end)
        chunks.append((cur, nxt))
        cur = nxt + _dt.timedelta(days=1)
    return chunks


def firms_area_url(map_key: str, product: str,
                   window: Tuple[float, float, float, float],
                   day_range: int, date: _dt.date) -> str:
    """Build one FIRMS area-API request URL (the key is embedded as given).

    Callers building provenance records should pass the key through
    :func:`redact_map_key` first.
    """
    w, s, e, n = window
    area = f"{w:.4f},{s:.4f},{e:.4f},{n:.4f}"
    return (f"{FIRMS_API_BASE}{FIRMS_AREA_PATH}/{map_key}/{product}/"
            f"{area}/{int(day_range)}/{date.isoformat()}")


def redact_map_key(url: str, map_key: str) -> str:
    """Replace the MAP_KEY in a FIRMS URL with ``<redacted>`` for storage."""
    return url.replace(f"/{map_key}/", "/<redacted>/")


def _looks_like_csv(text: str) -> bool:
    first = text.lstrip().splitlines()
    return bool(first) and first[0].lower().startswith("latitude")


def _parse_firms_csv(text: str) -> List[Dict[str, str]]:
    """Parse one FIRMS area-API CSV payload into row dicts.

    Handles both the VIIRS schema (``bright_ti4``/``bright_ti5``,
    single-letter confidence) and the MODIS schema (``brightness`` /
    ``bright_t31``, numeric confidence). Raises :class:`ValueError`
    with the API's message when the payload is an error string (e.g.
    ``"Invalid MAP_KEY."``).
    """
    if not _looks_like_csv(text):
        raise ValueError(
            f"FIRMS API did not return CSV: {text.strip()[:200]!r} "
            "(check the MAP_KEY and product name)")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    header = [h.strip() for h in lines[0].split(",")]
    rows: List[Dict[str, str]] = []
    for ln in lines[1:]:
        parts = [p.strip() for p in ln.split(",")]
        if len(parts) != len(header):
            continue  # tolerate ragged trailing lines
        rows.append(dict(zip(header, parts)))
    return rows


def _row_datetime(row: Dict[str, str]) -> Optional[_dt.datetime]:
    try:
        d = _dt.date.fromisoformat(row.get("acq_date", "").strip())
    except ValueError:
        return None
    t = row.get("acq_time", "").strip().zfill(4)
    try:
        hh, mm = int(t[:2]), int(t[2:4])
    except ValueError:
        hh, mm = 0, 0
    return _dt.datetime(d.year, d.month, d.day, hh, mm,
                        tzinfo=_dt.timezone.utc)


def _row_float(row: Dict[str, str], *names: str) -> float:
    for name in names:
        raw = row.get(name, "").strip()
        if raw:
            try:
                return float(raw)
            except ValueError:
                continue
    return math.nan


# ---------------------------------------------------------------------------
# FireField — canonical active-fire detection model
# ---------------------------------------------------------------------------

@dataclass
class FireField:
    """Active-fire detections from NASA FIRMS.

    Detections are points: ``times``/``lats``/``lons`` are parallel
    1-D sequences (one entry per detection), plus per-detection
    ``brightness`` (K), ``frp`` (MW), ``confidence`` (raw FIRMS value),
    ``satellite``, ``instrument``, and ``daynight`` (``"D"``/``"N"``).

    ``provenance`` records the (key-redacted) request URLs, per-request
    SHA-256 digests, the retrieval timestamp, and the tool version —
    following the :mod:`currents.sst_global` conventions.
    """

    times: List[_dt.datetime]
    lats: np.ndarray
    lons: np.ndarray
    brightness: np.ndarray
    frp: np.ndarray
    confidence: List[str]
    satellite: List[str]
    instrument: List[str]
    daynight: List[str]
    bbox: Tuple[float, float, float, float]
    instruments: Tuple[str, ...] = ("VIIRS_SNPP",)
    products: Tuple[str, ...] = ()
    source: str = "firms"
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        n = len(self.times)
        self.lats = np.asarray(self.lats, dtype=float).reshape(n)
        self.lons = np.asarray(self.lons, dtype=float).reshape(n)
        self.brightness = np.asarray(self.brightness, dtype=float).reshape(n)
        self.frp = np.asarray(self.frp, dtype=float).reshape(n)
        for name in ("confidence", "satellite", "instrument", "daynight"):
            vals = list(getattr(self, name))
            if len(vals) != n:
                raise ValueError(
                    f"FireField.{name}: length {len(vals)} != {n} detections")
            setattr(self, name, vals)
        self.bbox = validate_firms_bbox(self.bbox)
        self.instruments = tuple(str(i).strip().upper()
                                 for i in self.instruments)
        self.products = tuple(str(p) for p in self.products)

    def __len__(self) -> int:
        return len(self.times)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    @property
    def time_range(self) -> Tuple[Optional[_dt.datetime], Optional[_dt.datetime]]:
        """(earliest, latest) detection time, or (None, None) when empty."""
        if not self.times:
            return None, None
        return min(self.times), max(self.times)

    # -- filters ---------------------------------------------------------

    def _subset(self, idx: np.ndarray) -> "FireField":
        idx = np.asarray(idx)
        return FireField(
            times=[self.times[i] for i in idx],
            lats=self.lats[idx], lons=self.lons[idx],
            brightness=self.brightness[idx], frp=self.frp[idx],
            confidence=[self.confidence[i] for i in idx],
            satellite=[self.satellite[i] for i in idx],
            instrument=[self.instrument[i] for i in idx],
            daynight=[self.daynight[i] for i in idx],
            bbox=self.bbox, instruments=self.instruments,
            products=self.products, source=self.source,
            provenance=dict(self.provenance),
        )

    def select_time(self, start: DateLike, end: DateLike) -> "FireField":
        """Detections with ``start <= time <= end`` (inclusive)."""
        d0 = _coerce_date(start)
        d1 = _coerce_date(end)
        keep = np.array([d0 <= t.date() <= d1 for t in self.times])
        return self._subset(np.flatnonzero(keep))

    def select_bbox(self, bbox: Sequence[float]) -> "FireField":
        """Detections inside ``bbox`` (no antimeridian wrap here)."""
        lon_min, lat_min, lon_max, lat_max = validate_firms_bbox(bbox)
        keep = ((self.lons >= lon_min) & (self.lons <= lon_max)
                & (self.lats >= lat_min) & (self.lats <= lat_max))
        return self._subset(np.flatnonzero(keep))

    # -- density grid (the survey-viz scalar path) ------------------------

    def to_density_grid(self, resolution: float = 0.25,
                        frp_weighted: bool = False) -> Dict[str, Any]:
        """Bin detections into daily grids on a regular lat/lon mesh.

        Returns a ``times``/``lats``/``lons``/``values`` dict — the exact
        form ``survey-viz``'s renderer consumes — with one grid per UTC
        date holding fire *counts* per cell (``frp_weighted=True`` sums
        Fire Radiative Power in MW instead). Cells with no detections
        are 0. The grid spans the field's bbox; longitudes keep the
        conventional -180..180 axis.
        """
        if resolution <= 0:
            raise ValueError(f"resolution must be > 0, got {resolution}")
        lon_min, lat_min, lon_max, lat_max = self.bbox
        # Antimeridian-crossing bboxes render as a 0..360-spanning grid
        # only when the caller asked for one; the area API itself was
        # queried in two windows. Keep the grid on the caller's axis.
        span = (lon_max - lon_min) if lon_min < lon_max else (lon_max - lon_min) % 360.0
        nx = max(1, int(math.ceil(span / resolution)))
        ny = max(1, int(math.ceil((lat_max - lat_min) / resolution)))
        lons = lon_min + (np.arange(nx) + 0.5) * resolution
        lats = lat_min + (np.arange(ny) + 0.5) * resolution

        days = sorted({t.date() for t in self.times})
        times = [_dt.datetime(d.year, d.month, d.day,
                              tzinfo=_dt.timezone.utc) for d in days]
        values = np.zeros((len(days), ny, nx), dtype=float)
        day_index = {d: i for i, d in enumerate(days)}
        weights = (np.nan_to_num(self.frp, nan=0.0) if frp_weighted
                   else np.ones(len(self.times)))
        for k, t in enumerate(self.times):
            i = day_index[t.date()]
            ix = int((self.lons[k] - lon_min) // resolution)
            iy = int((self.lats[k] - lat_min) // resolution)
            ix = min(max(ix, 0), nx - 1)
            iy = min(max(iy, 0), ny - 1)
            values[i, iy, ix] += weights[k]
        return {
            "times": times,
            "lats": lats,
            "lons": lons,
            "values": values,
            "resolution": resolution,
            "weighting": "frp_MW" if frp_weighted else "count",
            "n_detections": len(self.times),
            "source": self.source,
        }

    # -- serialization ---------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "times": [t.isoformat() for t in self.times],
            "lats": [float(x) for x in self.lats],
            "lons": [float(x) for x in self.lons],
            "brightness": [float(x) for x in self.brightness],
            "frp": [float(x) for x in self.frp],
            "confidence": list(self.confidence),
            "satellite": list(self.satellite),
            "instrument": list(self.instrument),
            "daynight": list(self.daynight),
            "bbox": list(self.bbox),
            "instruments": list(self.instruments),
            "products": list(self.products),
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "FireField":
        return cls(
            times=[_dt.datetime.fromisoformat(t) for t in data["times"]],
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            brightness=np.asarray(data["brightness"], dtype=float),
            frp=np.asarray(data["frp"], dtype=float),
            confidence=list(data["confidence"]),
            satellite=list(data["satellite"]),
            instrument=list(data["instrument"]),
            daynight=list(data["daynight"]),
            bbox=tuple(data["bbox"]),
            instruments=tuple(data.get("instruments", ("VIIRS_SNPP",))),
            products=tuple(data.get("products", ())),
            source=data.get("source", "firms"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        """Write the field as JSON; return ``path``."""
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)
        return path

    @classmethod
    def from_json(cls, path: str) -> "FireField":
        """Read a field written by :meth:`to_json`."""
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture -------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-125.0, 32.0, -114.0, 42.0),
                  start: DateLike = "2024-08-01", end: DateLike = "2024-08-07",
                  n: int = 60, seed: int = 7,
                  instruments: Sequence[str] = ("VIIRS_SNPP",),
                  source: str = "synthetic") -> "FireField":
        """Deterministic synthetic detections (offline tests / demos).

        Detections cluster around two fake fire perimeters inside
        ``bbox`` so density grids show structure; one detection per day
        minimum keeps every daily grid non-empty when ``n`` is large
        enough.
        """
        rng = np.random.default_rng(seed)
        lon_min, lat_min, lon_max, lat_max = validate_firms_bbox(bbox)
        d0, d1 = _coerce_date(start), _coerce_date(end)
        if d1 < d0:
            raise ValueError(f"end {d1} is before start {d0}")
        ndays = (d1 - d0).days + 1
        # Two cluster centers, deterministic from the seed.
        centers = [
            (lon_min + 0.3 * (lon_max - lon_min),
             lat_min + 0.6 * (lat_max - lat_min)),
            (lon_min + 0.7 * (lon_max - lon_min),
             lat_min + 0.35 * (lat_max - lat_min)),
        ]
        pick = rng.integers(0, 2, size=n)
        lons = np.array([centers[p][0] for p in pick]) + rng.normal(0, 0.35, n)
        lats = np.array([centers[p][1] for p in pick]) + rng.normal(0, 0.25, n)
        lons = np.clip(lons, lon_min, lon_max)
        lats = np.clip(lats, lat_min, lat_max)
        days = [_dt.datetime(d0.year, d0.month, d0.day,
                             tzinfo=_dt.timezone.utc)
                + _dt.timedelta(days=int(i)) for i in rng.integers(0, ndays, n)]
        # Sort oldest-first like the API returns.
        order = np.argsort([t.timestamp() for t in days])
        sats = ["Suomi-NPP", "NOAA-20"]
        return cls(
            times=[days[i] for i in order],
            lats=lats[order], lons=lons[order],
            brightness=rng.normal(330.0, 25.0, n),
            frp=np.abs(rng.normal(15.0, 12.0, n)),
            confidence=[["l", "n", "h"][int(x)] for x in rng.integers(0, 3, n)],
            satellite=[sats[int(x)] for x in rng.integers(0, 2, n)],
            instrument=["VIIRS"] * n,
            daynight=[["D", "N"][int(x)] for x in rng.integers(0, 2, n)],
            bbox=(lon_min, lat_min, lon_max, lat_max),
            instruments=tuple(str(i).upper() for i in instruments),
            products=("VIIRS_SNPP_SP",),
            source=source,
            provenance={"synthetic": True, "seed": seed},
        )


# ---------------------------------------------------------------------------
# fetch_firms
# ---------------------------------------------------------------------------

def fetch_firms(bbox: Sequence[float], start: DateLike, end: DateLike,
                instruments: Sequence[str] = ("VIIRS_SNPP",),
                map_key: Optional[str] = None,
                products: Optional[Sequence[str]] = None,
                timeout: int = 300) -> FireField:
    """Fetch NASA FIRMS active-fire detections for ``bbox``/``[start, end]``.

    Args:
        bbox: (lon_min, lat_min, lon_max, lat_max); antimeridian-crossing
            boxes are split into two area-API requests.
        start/end: inclusive date range (dates, datetimes, or ISO strings).
        instruments: FIRMS instrument families, e.g.
            ``("VIIRS_SNPP", "MODIS")``. One product per instrument per
            date is requested (NRT for recent dates, SP for older ones —
            see :func:`firms_product_for`).
        map_key: explicit FIRMS MAP_KEY; defaults to the
            ``FIRMS_MAP_KEY`` environment variable.
        products: explicit product names overriding the NRT/SP tiering
            (advanced; one per instrument, order-matched).
        timeout: per-request HTTP timeout in seconds.

    Returns:
        A :class:`FireField` with key-redacted provenance.

    Raises:
        CredentialsMissing: no MAP_KEY found.
        ValueError: bad bbox / dates / instrument / API error payload.
    """
    lon_min, lat_min, lon_max, lat_max = validate_firms_bbox(bbox)
    d0, d1 = _coerce_date(start), _coerce_date(end)
    if d1 < d0:
        raise ValueError(f"end {d1} is before start {d0}")
    if d0 < FIRMS_START:
        raise ValueError(
            f"FIRMS coverage starts {FIRMS_START} (requested {d0})")
    key = firms_map_key(map_key)

    inst_list = [str(i).strip().upper() for i in instruments]
    if not inst_list:
        raise ValueError("instruments must name at least one instrument")
    for inst in inst_list:
        if inst not in FIRMS_INSTRUMENTS:
            raise ValueError(
                f"unknown FIRMS instrument {inst!r}; choose from "
                f"{sorted(FIRMS_INSTRUMENTS)}")

    windows = firms_area_windows((lon_min, lat_min, lon_max, lat_max))
    chunks = firms_window_chunks(d0, d1)

    all_rows: List[Dict[str, str]] = []
    used_products: List[str] = []
    request_records: List[Dict[str, Any]] = []
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()

    for ii, inst in enumerate(inst_list):
        for (c0, c1) in chunks:
            product = (products[ii] if products
                       else firms_product_for(inst, c0))
            if product not in used_products:
                used_products.append(product)
            day_range = (c1 - c0).days + 1
            for window in windows:
                url = firms_area_url(key, product, window, day_range, c0)
                payload = _download_bytes(url, timeout=timeout)
                digest = hashlib.sha256(payload).hexdigest()
                request_records.append({
                    "url": redact_map_key(url, key),
                    "product": product,
                    "window": list(window),
                    "day_range": day_range,
                    "date": c0.isoformat(),
                    "sha256": digest,
                    "bytes": len(payload),
                })
                text = payload.decode("utf-8", errors="replace")
                all_rows.extend(_parse_firms_csv(text))

    # Oldest-first, like the API returns; drop exact duplicates that can
    # arise at window seams.
    seen = set()
    rows: List[Dict[str, str]] = []
    for row in all_rows:
        sig = (row.get("latitude"), row.get("longitude"),
               row.get("acq_date"), row.get("acq_time"),
               row.get("satellite"))
        if sig in seen:
            continue
        seen.add(sig)
        rows.append(row)

    times: List[_dt.datetime] = []
    lats: List[float] = []
    lons: List[float] = []
    brightness: List[float] = []
    frp: List[float] = []
    confidence: List[str] = []
    satellite: List[str] = []
    instrument: List[str] = []
    daynight: List[str] = []
    for row in rows:
        t = _row_datetime(row)
        if t is None:
            continue
        try:
            la = float(row["latitude"])
            lo = float(row["longitude"])
        except (KeyError, ValueError):
            continue
        times.append(t)
        lats.append(la)
        lons.append(lo)
        brightness.append(_row_float(row, "bright_ti4", "brightness"))
        frp.append(_row_float(row, "frp"))
        confidence.append(row.get("confidence", "").strip())
        satellite.append(row.get("satellite", "").strip())
        instrument.append(row.get("instrument", "").strip())
        daynight.append(row.get("daynight", "").strip())

    order = sorted(range(len(times)), key=lambda i: times[i].timestamp())
    provenance = {
        "source": "NASA FIRMS active fires (area API)",
        "api_base": f"{FIRMS_API_BASE}{FIRMS_AREA_PATH}",
        "bbox": [lon_min, lat_min, lon_max, lat_max],
        "time_window": [d0.isoformat(), d1.isoformat()],
        "instruments": inst_list,
        "products_used": used_products,
        "requests": request_records,
        "n_requests": len(request_records),
        "n_detections": len(times),
        "retrieved_utc": retrieved_at,
        "map_key": "<redacted>",
        "tool": f"survey-currents/{_tool_version()}",
    }
    return FireField(
        times=[times[i] for i in order],
        lats=np.asarray([lats[i] for i in order]),
        lons=np.asarray([lons[i] for i in order]),
        brightness=np.asarray([brightness[i] for i in order]),
        frp=np.asarray([frp[i] for i in order]),
        confidence=[confidence[i] for i in order],
        satellite=[satellite[i] for i in order],
        instrument=[instrument[i] for i in order],
        daynight=[daynight[i] for i in order],
        bbox=(lon_min, lat_min, lon_max, lat_max),
        instruments=tuple(inst_list),
        products=tuple(used_products),
        provenance=provenance,
    )


def _tool_version() -> str:
    from . import __version__
    return __version__

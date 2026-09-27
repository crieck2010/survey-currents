"""USGS Water Services (NWIS) streamgage daily-value acquisition.

:func:`fetch_usgs` discovers stream gages inside a bbox through the
keyless USGS site service (RDB), fetches per-site daily values for
discharge (parameter 00060, ft³/s) and gage height (00065, ft) through
the keyless daily-values service (JSON), and returns a
:class:`GageField`: one :class:`GageRecord` per site with daily time
series, site metadata (site number, name, lat/lon, HUC, drainage area
when published), and honest gap accounting — sites whose record is
shorter than ``min_record_days`` are excluded (counted in
provenance), and days with no reported value are recorded as missing,
never filled.

Access truth (verified live 2026-09-27 — see docs/DATA_SOURCES.md):

* ``https://waterservices.usgs.gov/nwis/site/`` — keyless anonymous
  HTTPS (RDB); ``hasDataTypeCd=dv&parameterCd=00060`` filters to sites
  with daily discharge data; ``siteOutput=expanded`` adds drainage
  area (``drain_area_va``, sq mi), HUC, and timezone.
* ``https://waterservices.usgs.gov/nwis/dv/`` — keyless anonymous
  HTTPS (JSON, ``format=json``); ``statCd=00003`` requests daily MEAN
  values (the dv default, stated explicitly for determinism).
* The USGS missing-value sentinel is ``"-999999"`` — parsed to NaN.
  Value qualifiers (``P`` provisional, ``A`` approved, ``e``
  estimated, …) are carried on each series.
* NWIS is US-only (plus Puerto Rico / Pacific territories): a bbox
  with no USGS gages yields an empty field, honestly reported.

Cache discipline (mirrors :mod:`currents.storms`): the site inventory
(keyed by bbox + parameters) and each per-batch dv payload (keyed by
sorted site numbers + dates + parameters) are cached under
``$SURVEY_CURRENTS_CACHE/streamgages`` (else
``~/.cache/survey-currents/streamgages``) with atomic writes and
``.sha256`` sidecars; a corrupt cache entry is re-downloaded.
Entries older than ``max_cache_age_days`` (default 7 — recent daily
values are provisional and get revised) are re-fetched; ``refresh=True``
forces re-download. Repeated renders of the same query never re-hit
the network.

Units: the adapter returns NATIVE USGS units by default
(ft³/s for 00060, ft for 00065) — ``units="si"`` converts to m³/s / m
on :class:`GageField.to_si` (factors documented below). The renderer
and provenance always state which unit system a field uses.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json as _json
import os
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from dataclasses import field as _dc_field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_date
from .sst_global import validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]

__all__ = [
    "USGS_SITE_BASE",
    "USGS_DV_BASE",
    "USGS_PARAMETERS",
    "USGS_MISSING",
    "USGS_RECORD_START",
    "USGS_MAX_CACHE_AGE_DAYS",
    "FT3S_TO_M3S",
    "FT_TO_M",
    "USGS_SITE_LIMIT",
    "GageRecord",
    "GageField",
    "usgs_cache_dir",
    "site_inventory_url",
    "dv_request_url",
    "fetch_usgs",
]

#: Keyless USGS Water Services endpoints (verified live 2026-09-27).
USGS_SITE_BASE = "https://waterservices.usgs.gov/nwis/site/"
USGS_DV_BASE = "https://waterservices.usgs.gov/nwis/dv/"
#: Supported daily-value parameters: code -> (label, native unit).
USGS_PARAMETERS: Dict[str, Tuple[str, str]] = {
    "00060": ("discharge", "ft3/s"),
    "00065": ("gage height", "ft"),
}
#: USGS missing-value sentinel in dv JSON (verified live 2026-09-27).
USGS_MISSING = "-999999"
#: Earliest daily values in NWIS (19th century for a few long-record
#: gages); used only to reject absurd starts.
USGS_RECORD_START = _dt.date(1857, 1, 1)
#: Default cache freshness (days): recent daily values are provisional
#: and get revised, so dv payloads revalidate sooner than the
#: 30-day discipline used for static archives.
USGS_MAX_CACHE_AGE_DAYS = 7
#: Exact conversion factors (NIST).
FT3S_TO_M3S = 0.028316846592
FT_TO_M = 0.3048
#: Default cap on sites per fetch (deterministic: sorted by site number).
USGS_SITE_LIMIT = 200
#: Sites per dv request (the service accepts long site lists; 20 keeps
#: each request bounded and each cache entry small).
_USGS_DV_BATCH = 20


def usgs_cache_dir() -> str:
    """User cache root for USGS downloads.

    ``$SURVEY_CURRENTS_CACHE/streamgages`` when set, else
    ``~/.cache/survey-currents/streamgages``. Created on demand.
    """
    root = os.environ.get("SURVEY_CURRENTS_CACHE") or os.path.join(
        os.path.expanduser("~"), ".cache", "survey-currents")
    path = os.path.join(root, "streamgages")
    os.makedirs(path, exist_ok=True)
    return path


def site_inventory_url(bbox: Sequence[float],
                       parameters: Sequence[str] = ("00060",)) -> str:
    """Keyless site-service URL discovering dv sites in ``bbox``."""
    lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
    params = ",".join(parameters)
    return (
        f"{USGS_SITE_BASE}?format=rdb"
        f"&bBox={lon_min},{lat_min},{lon_max},{lat_max}"
        f"&hasDataTypeCd=dv&parameterCd={params}"
        f"&siteOutput=expanded"
    )


def dv_request_url(sites: Sequence[str], start: DateLike, end: DateLike,
                   parameters: Sequence[str] = ("00060",)) -> str:
    """Keyless dv-service URL for one batch of ``sites``."""
    d0, d1 = _coerce_date(start), _coerce_date(end)
    site_list = ",".join(sites)
    params = ",".join(parameters)
    return (
        f"{USGS_DV_BASE}?format=json"
        f"&sites={site_list}"
        f"&startDT={d0.isoformat()}&endDT={d1.isoformat()}"
        f"&parameterCd={params}&statCd=00003"
    )


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
    :func:`currents.storms.ensure_ibtracs_file`, adapted for text.
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
            f"USGS download failed for {url}: {exc}. The USGS Water "
            "Services endpoints are keyless HTTPS — this is usually a "
            "connectivity issue.") from exc


# ---------------------------------------------------------------------------
# Response parsing (pure, offline-testable)
# ---------------------------------------------------------------------------

def _parse_site_rdb(text: str) -> List[Dict[str, Any]]:
    """Parse a site-service RDB body into site dicts.

    Skips comment lines (``#``), the header row, and the column-width
    row; every subsequent non-empty line is a site. Missing fields
    (empty strings) become ``None`` for numerics.
    """
    sites: List[Dict[str, Any]] = []
    header: Optional[List[str]] = None
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if header is None:
            header = [c.strip() for c in cols]
            continue
        if all(set(c.strip()) <= set("0123456789s") for c in cols):
            # Column-width row (e.g. "5s 15s 50s ...") — skip.
            continue
        row = {h: (c.strip() or None) for h, c in zip(header, cols)}

        def _num(key: str) -> Optional[float]:
            v = row.get(key)
            if v in (None, ""):
                return None
            try:
                return float(v)
            except ValueError:
                return None

        sites.append({
            "agency": row.get("agency_cd"),
            "site_no": row.get("site_no"),
            "site_name": row.get("station_nm"),
            "site_type": row.get("site_tp_cd"),
            "lat": _num("dec_lat_va"),
            "lon": _num("dec_long_va"),
            "drain_area_sqmi": _num("drain_area_va"),
            "huc": row.get("huc_cd"),
            "tz": row.get("tz_cd"),
            "alt_ft": _num("alt_va"),
        })
    return [s for s in sites if s.get("site_no")]


def _parse_dv_json(payload: Dict[str, Any]) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """Parse a dv-service JSON body into {site: {param: series}}.

    Each series is ``{"dates", "values", "unit", "qualifiers"}``;
    ``values`` are floats with the USGS ``"-999999"`` sentinel mapped
    to NaN; ``qualifiers`` is the per-point qualifier list (``P``
    provisional, ``A`` approved, ``e`` estimated, …). Empty timeSeries
    (a site with no data in the window) yield no entry — the caller
    counts them honestly.
    """
    out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for ts in payload.get("value", {}).get("timeSeries", []):
        info = ts.get("sourceInfo", {})
        codes = info.get("siteCode", [])
        site_no = codes[0].get("value") if codes else None
        var = ts.get("variable", {})
        var_codes = var.get("variableCode", [])
        param = var_codes[0].get("value") if var_codes else None
        unit = (var.get("unit") or {}).get("unitCode")
        values = (ts.get("values") or [{}])[0].get("value", [])
        if not site_no or not param or not values:
            continue
        dates: List[_dt.date] = []
        vals: List[float] = []
        quals: List[List[str]] = []
        for pt in values:
            raw = str(pt.get("value", ""))
            if raw == USGS_MISSING or raw == "":
                v = float("nan")
            else:
                try:
                    v = float(raw)
                except ValueError:
                    v = float("nan")
            dtxt = str(pt.get("dateTime", ""))[:10]
            try:
                d = _dt.date.fromisoformat(dtxt)
            except ValueError:
                continue
            dates.append(d)
            vals.append(v)
            quals.append([str(q) for q in (pt.get("qualifiers") or [])])
        site = out.setdefault(site_no, {})
        site[param] = {
            "dates": dates,
            "values": np.asarray(vals, dtype=float),
            "unit": unit,
            "qualifiers": quals,
        }
    return out


# ---------------------------------------------------------------------------
# Field model
# ---------------------------------------------------------------------------

@dataclass
class GageRecord:
    """One USGS streamgage: metadata + per-parameter daily series.

    ``series`` maps parameter code (``"00060"``, ``"00065"``) to
    ``{"dates", "values", "unit", "qualifiers"}``. Days with no
    reported value are NaN in ``values`` and listed in
    :meth:`missing_days` — gaps are recorded, never filled.
    """

    site_no: str
    site_name: str
    lat: float
    lon: float
    huc: Optional[str] = None
    drain_area_sqmi: Optional[float] = None
    tz: Optional[str] = None
    series: Dict[str, Dict[str, Any]] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        self.site_no = str(self.site_no)
        self.lat = float(self.lat)
        self.lon = float(self.lon)
        self.series = dict(self.series or {})

    # -- accessors -----------------------------------------------------

    def param_series(self, parameter: str = "00060"
                     ) -> Optional[Dict[str, Any]]:
        """The daily series for ``parameter`` (None when absent)."""
        return self.series.get(str(parameter))

    def valid_mask(self, parameter: str = "00060") -> np.ndarray:
        """Boolean mask of finite values for ``parameter``."""
        s = self.param_series(parameter)
        if s is None:
            return np.zeros(0, dtype=bool)
        return np.isfinite(np.asarray(s["values"], dtype=float))

    def n_valid(self, parameter: str = "00060") -> int:
        """Count of finite daily values for ``parameter``."""
        return int(np.sum(self.valid_mask(parameter)))

    def missing_days(self, parameter: str = "00060",
                     start: Optional[DateLike] = None,
                     end: Optional[DateLike] = None) -> List[_dt.date]:
        """Dates in [start, end] with no finite value (honest gaps)."""
        s = self.param_series(parameter)
        if s is None:
            return []
        dates = s["dates"]
        vals = np.asarray(s["values"], dtype=float)
        d0 = _coerce_date(start) if start is not None else dates[0]
        d1 = _coerce_date(end) if end is not None else dates[-1]
        by_date = {d: np.isfinite(v) for d, v in zip(dates, vals)}
        missing: List[_dt.date] = []
        d = d0
        while d <= d1:
            if not by_date.get(d, False):
                missing.append(d)
            d += _dt.timedelta(days=1)
        return missing

    def latest(self, parameter: str = "00060"
               ) -> Tuple[Optional[_dt.date], Optional[float]]:
        """(date, value) of the latest finite value, or (None, None)."""
        s = self.param_series(parameter)
        if s is None:
            return None, None
        vals = np.asarray(s["values"], dtype=float)
        idx = np.flatnonzero(np.isfinite(vals))
        if idx.size == 0:
            return None, None
        k = int(idx[-1])
        return s["dates"][k], float(vals[k])

    def percentile_of_record(self, parameter: str = "00060") -> Optional[float]:
        """Percentile (0–100) of the latest value within this record.

        The documented rule behind marker coloring: the latest finite
        daily value's rank among the site's own retrieved record for
        ``parameter`` (``100 * (# below) / (n - 1)``, NaN-safe). None
        when the site has no valid value.
        """
        s = self.param_series(parameter)
        if s is None:
            return None
        vals = np.asarray(s["values"], dtype=float)
        finite = vals[np.isfinite(vals)]
        if finite.size == 0:
            return None
        if finite.size == 1:
            return 50.0
        latest_val = float(finite[-1])
        rank = float(np.sum(finite < latest_val))
        return 100.0 * rank / (finite.size - 1)

    # -- units -----------------------------------------------------------

    def to_si(self) -> "GageRecord":
        """Copy with values converted to SI (m³/s for 00060, m for 00065)."""
        factor = {"00060": FT3S_TO_M3S, "00065": FT_TO_M}
        unit = {"00060": "m3/s", "00065": "m"}
        series: Dict[str, Dict[str, Any]] = {}
        for code, s in self.series.items():
            f = factor.get(str(code), 1.0)
            series[str(code)] = {
                "dates": list(s["dates"]),
                "values": np.asarray(s["values"], dtype=float) * f,
                "unit": unit.get(str(code), s.get("unit")),
                "qualifiers": list(s.get("qualifiers", [])),
            }
        return GageRecord(
            site_no=self.site_no, site_name=self.site_name,
            lat=self.lat, lon=self.lon, huc=self.huc,
            drain_area_sqmi=self.drain_area_sqmi, tz=self.tz,
            series=series)

    # -- serialization ---------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "site_no": self.site_no,
            "site_name": self.site_name,
            "lat": self.lat,
            "lon": self.lon,
            "huc": self.huc,
            "drain_area_sqmi": self.drain_area_sqmi,
            "tz": self.tz,
            "series": {
                code: {
                    "dates": [d.isoformat() for d in s["dates"]],
                    "values": [None if v != v else float(v)
                               for v in np.asarray(s["values"], dtype=float)],
                    "unit": s.get("unit"),
                    "qualifiers": [list(q) for q in s.get("qualifiers", [])],
                }
                for code, s in self.series.items()
            },
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GageRecord":
        series: Dict[str, Dict[str, Any]] = {}
        for code, s in (data.get("series") or {}).items():
            vals = [float("nan") if v is None else float(v)
                    for v in s.get("values", [])]
            series[str(code)] = {
                "dates": [_dt.date.fromisoformat(d)
                          for d in s.get("dates", [])],
                "values": np.asarray(vals, dtype=float),
                "unit": s.get("unit"),
                "qualifiers": [list(q) for q in s.get("qualifiers", [])],
            }
        return cls(
            site_no=data["site_no"], site_name=data.get("site_name", ""),
            lat=data["lat"], lon=data["lon"], huc=data.get("huc"),
            drain_area_sqmi=data.get("drain_area_sqmi"), tz=data.get("tz"),
            series=series)


@dataclass
class GageField:
    """A set of USGS gage records for one query.

    ``provenance`` records the exact request URLs, retrieval time,
    site counts (discovered / requested / with data / no data /
    excluded by ``min_record_days``), and the unit system — following
    the :mod:`currents.storms` conventions.
    """

    records: List[GageRecord]
    bbox: Tuple[float, float, float, float]
    start: _dt.date
    end: _dt.date
    parameters: Tuple[str, ...] = ("00060",)
    units: str = "native"
    source: str = "usgs"
    provenance: Dict[str, Any] = _dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        self.bbox = validate_sst_bbox(self.bbox)
        self.start = _coerce_date(self.start)
        self.end = _coerce_date(self.end)
        if self.start > self.end:
            raise ValueError(
                f"GageField: start {self.start} is after end {self.end}")
        self.records = list(self.records)
        self.parameters = tuple(str(p) for p in self.parameters)
        self.units = str(self.units or "native")

    def __len__(self) -> int:
        return len(self.records)

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """The requested (lon_min, lat_min, lon_max, lat_max)."""
        return self.bbox

    @property
    def unit_label(self) -> str:
        """Display unit for the primary parameter."""
        native = USGS_PARAMETERS.get(self.parameters[0], ("", "ft3/s"))[1]
        if self.units == "si":
            return {"ft3/s": "m³/s", "ft": "m"}.get(native, native)
        return {"ft3/s": "ft³/s"}.get(native, native)

    def site_numbers(self) -> List[str]:
        """Site numbers in record order."""
        return [r.site_no for r in self.records]

    def select_site(self, site_no: str) -> "GageField":
        """Field with only the named site (empty when unknown)."""
        kept = [r for r in self.records if r.site_no == str(site_no)]
        return GageField(records=kept, bbox=self.bbox, start=self.start,
                         end=self.end, parameters=self.parameters,
                         units=self.units, source=self.source,
                         provenance=dict(self.provenance))

    def regional_median(self, parameter: str = "00060"
                        ) -> Tuple[List[_dt.date], np.ndarray]:
        """Daily median across sites with valid values (hydrograph rule).

        For each date in [start, end], the median of the finite values
        across all records carrying ``parameter``; dates with no valid
        value anywhere are NaN. Documented as the multi-gage
        hydrograph selection rule (see docs/STREAMFLOW.md): exactly one
        site -> that site's own series; otherwise the regional median.
        """
        d0, d1 = self.start, self.end
        ndays = (d1 - d0).days + 1
        dates = [d0 + _dt.timedelta(days=k) for k in range(ndays)]
        grid = np.full((len(self.records), ndays), np.nan)
        for i, rec in enumerate(self.records):
            s = rec.param_series(parameter)
            if s is None:
                continue
            vals = np.asarray(s["values"], dtype=float)
            for d, v in zip(s["dates"], vals):
                k = (d - d0).days
                if 0 <= k < ndays and np.isfinite(v):
                    grid[i, k] = v
        with np.errstate(all="ignore"):
            med = np.nanmedian(grid, axis=0)
        return dates, med

    def to_si(self) -> "GageField":
        """Copy with all records converted to SI units."""
        return GageField(
            records=[r.to_si() for r in self.records],
            bbox=self.bbox, start=self.start, end=self.end,
            parameters=self.parameters, units="si", source=self.source,
            provenance={**self.provenance, "units": "si"})

    # -- serialization ---------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gage_records": [r.to_dict() for r in self.records],
            "bbox": list(self.bbox),
            "start": self.start.isoformat(),
            "end": self.end.isoformat(),
            "parameters": list(self.parameters),
            "units": self.units,
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GageField":
        return cls(
            records=[GageRecord.from_dict(r)
                     for r in data.get("gage_records", [])],
            bbox=tuple(data["bbox"]),
            start=_coerce_date(data["start"]),
            end=_coerce_date(data["end"]),
            parameters=tuple(data.get("parameters", ("00060",))),
            units=data.get("units", "native"),
            source=data.get("source", "usgs"),
            provenance=dict(data.get("provenance", {})))

    def to_json(self, path: str) -> str:
        """Write the field as JSON; return ``path``."""
        with open(path, "w", encoding="utf-8") as fh:
            _json.dump(self.to_dict(), fh, indent=2)
        return path

    @classmethod
    def from_json(cls, path: str) -> "GageField":
        """Read a field written by :meth:`to_json`."""
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(_json.load(fh))

    # -- synthetic fixture -------------------------------------------------

    @classmethod
    def synthetic(cls, bbox: Sequence[float] = (-83.5, 42.0, -82.0, 43.5),
                  start: DateLike = "2024-06-01",
                  end: DateLike = "2024-08-31",
                  n_sites: int = 4, seed: int = 11,
                  parameters: Sequence[str] = ("00060",),
                  source: str = "synthetic") -> "GageField":
        """Deterministic synthetic gage records (offline tests / demos).

        Builds ``n_sites`` sites inside ``bbox`` with smooth seasonal
        discharge (base flow + freshet peak + noise) and one site with
        a gap block (NaN days) so gap accounting stays testable. The
        first site is always numbered ``"01111110"`` ("SYNTHETIC RIVER
        AT TESTVILLE") so selection tests are deterministic.
        """
        rng = np.random.default_rng(seed)
        lon_min, lat_min, lon_max, lat_max = validate_sst_bbox(bbox)
        d0, d1 = _coerce_date(start), _coerce_date(end)
        if d1 < d0:
            raise ValueError(f"end {d1} is before start {d0}")
        ndays = (d1 - d0).days + 1
        dates = [d0 + _dt.timedelta(days=k) for k in range(ndays)]
        names = ["SYNTHETIC RIVER AT TESTVILLE",
                 "MOCK CREEK NEAR DEMOTON",
                 "FAKE RIVER AT SAMPLEVILLE",
                 "PLACEBO BROOK NEAR NULLTOWN"]
        records: List[GageRecord] = []
        for s in range(n_sites):
            base = float(rng.uniform(80.0, 400.0))
            k = np.arange(ndays, dtype=float)
            flow = (base
                    + 0.6 * base * np.exp(-((k - ndays * 0.35) ** 2)
                                    / (2 * (ndays * 0.12) ** 2))
                    + 0.08 * base * rng.standard_normal(ndays))
            flow = np.clip(flow, 1.0, None)
            if s == 1:
                # A gap block: days 20..29 are missing (NaN).
                flow[20:30] = np.nan
            series = {
                "00060": {
                    "dates": list(dates),
                    "values": flow,
                    "unit": "ft3/s",
                    "qualifiers": [["A"]] * ndays,
                }
            }
            if "00065" in parameters:
                gh = 2.0 + 0.004 * flow + 0.02 * rng.standard_normal(ndays)
                series["00065"] = {
                    "dates": list(dates),
                    "values": np.where(np.isnan(flow), np.nan, gh),
                    "unit": "ft",
                    "qualifiers": [["A"]] * ndays,
                }
            records.append(GageRecord(
                site_no=f"0111111{s}",
                site_name=names[s % len(names)],
                lat=float(lat_min + 0.2 * (lat_max - lat_min)
                          + 0.6 * (lat_max - lat_min) * rng.random()),
                lon=float(lon_min + 0.2 * (lon_max - lon_min)
                          + 0.6 * (lon_max - lon_min) * rng.random()),
                huc="04090001",
                drain_area_sqmi=float(rng.uniform(50.0, 900.0)),
                tz="EST",
                series=series))
        return cls(records=records,
                   bbox=(lon_min, lat_min, lon_max, lat_max),
                   start=d0, end=d1, parameters=tuple(parameters),
                   units="native", source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# Public fetch
# ---------------------------------------------------------------------------

def fetch_usgs(bbox: Sequence[float], start: DateLike, end: DateLike,
               parameters: Sequence[str] = ("00060",),
               min_record_days: int = 30,
               site_limit: int = USGS_SITE_LIMIT,
               units: str = "native",
               cache_dir: Optional[str] = None,
               max_cache_age_days: int = USGS_MAX_CACHE_AGE_DAYS,
               refresh: bool = False) -> GageField:
    """Fetch USGS streamgage daily values for ``bbox`` x ``[start, end]``.

    Discovers sites through the keyless site service
    (``hasDataTypeCd=dv``, so every returned site has daily data of at
    least one requested parameter), then fetches daily values in
    batches through the keyless dv service. Both responses are cached
    (SHA-256-verified, atomic writes, ``max_cache_age_days``
    revalidation) so repeated renders never re-hit the network.

    Args:
        parameters: dv parameter codes — ``"00060"`` discharge
            (ft³/s), ``"00065"`` gage height (ft). Unknown codes raise
            ``ValueError``.
        min_record_days: keep only sites with at least this many
            finite daily values of the primary parameter in the window
            (excluded sites are counted in provenance, not silently
            dropped).
        site_limit: deterministic cap on sites fetched (sorted by site
            number); large bboxes should be subdivided by the caller —
            the national inventory is tens of thousands of sites.
        units: ``"native"`` (ft³/s / ft, the USGS publication units)
            or ``"si"`` (m³/s / m, converted with the exact NIST
            factors).

    Missing daily values stay NaN; :meth:`GageRecord.missing_days`
    lists them. The NWIS record is US-only: bboxes outside USGS
    coverage yield an empty field with ``provenance["empty_reason"]``.
    """
    box = validate_sst_bbox(bbox)
    d0, d1 = _coerce_date(start), _coerce_date(end)
    if d1 < d0:
        raise ValueError(f"end {d1} is before start {d0}")
    if d0 < USGS_RECORD_START:
        raise ValueError(
            f"start {d0} predates the NWIS daily record "
            f"({USGS_RECORD_START})")
    if d1 > _dt.date.today() + _dt.timedelta(days=2):
        raise ValueError(f"end {d1.isoformat()} is in the future")
    params = tuple(str(p) for p in parameters)
    for p in params:
        if p not in USGS_PARAMETERS:
            raise ValueError(
                f"unknown USGS parameter {p!r}; supported: "
                f"{sorted(USGS_PARAMETERS)}")
    if int(min_record_days) < 0:
        raise ValueError(f"min_record_days must be >= 0, got {min_record_days}")
    if int(site_limit) < 1:
        raise ValueError(f"site_limit must be >= 1, got {site_limit}")
    if units not in ("native", "si"):
        raise ValueError(f"units must be 'native' or 'si', got {units!r}")

    cdir = cache_dir or usgs_cache_dir()
    request_urls: List[str] = []

    # 1. Site discovery (cached by bbox + parameters).
    inv_url = site_inventory_url(box, params)
    request_urls.append(inv_url)
    inv_text, inv_prov = _cached_get_text(
        inv_url, cdir, max_cache_age_days, refresh=refresh)
    discovered = sorted(_parse_site_rdb(inv_text),
                        key=lambda s: s["site_no"])
    primary = params[0]

    # 2. Daily values, batched (each batch cached independently).
    site_nos = [s["site_no"] for s in discovered[:int(site_limit)]]
    by_site: Dict[str, Dict[str, Dict[str, Any]]] = {}
    batch_provs: List[Dict[str, Any]] = []
    for k in range(0, len(site_nos), _USGS_DV_BATCH):
        batch = site_nos[k:k + _USGS_DV_BATCH]
        dv_url = dv_request_url(batch, d0, d1, params)
        request_urls.append(dv_url)
        dv_text, dv_prov = _cached_get_text(
            dv_url, cdir, max_cache_age_days, refresh=refresh)
        batch_provs.append(dv_prov)
        payload = _json.loads(dv_text)
        parsed = _parse_dv_json(payload)
        for site_no, series in parsed.items():
            by_site.setdefault(site_no, {}).update(series)

    # 3. Assemble records; enforce the minimum-record rule.
    meta = {s["site_no"]: s for s in discovered}
    records: List[GageRecord] = []
    n_no_data = 0
    n_excluded = 0
    missing_days_total = 0
    for site_no in site_nos:
        s = meta.get(site_no)
        series = by_site.get(site_no, {})
        if s is None or not series:
            n_no_data += 1
            continue
        rec = GageRecord(
            site_no=site_no, site_name=s.get("site_name") or site_no,
            lat=s["lat"] if s.get("lat") is not None else 0.0,
            lon=s["lon"] if s.get("lon") is not None else 0.0,
            huc=s.get("huc"), drain_area_sqmi=s.get("drain_area_sqmi"),
            tz=s.get("tz"), series=series)
        if rec.n_valid(primary) < int(min_record_days):
            n_excluded += 1
            continue
        missing_days_total += len(rec.missing_days(primary, d0, d1))
        records.append(rec)

    provenance = {
        "source": "usgs",
        "product": "USGS Water Services (NWIS) daily values",
        "site_service_url": inv_url,
        "dv_service": USGS_DV_BASE,
        "request_urls": request_urls,
        "retrieved_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "bbox": list(box),
        "start": d0.isoformat(),
        "end": d1.isoformat(),
        "parameters": list(params),
        "parameter_meta": {p: {"label": USGS_PARAMETERS[p][0],
                               "native_unit": USGS_PARAMETERS[p][1]}
                           for p in params},
        "statistic": "00003 (daily mean, explicit)",
        "min_record_days": int(min_record_days),
        "site_limit": int(site_limit),
        "n_sites_discovered": len(discovered),
        "n_sites_requested": len(site_nos),
        "n_sites_with_data": len(records),
        "n_sites_no_data": n_no_data,
        "n_sites_excluded_short_record": n_excluded,
        "missing_days_total": missing_days_total,
        "units": units,
        "cache": {
            "inventory": {k: inv_prov.get(k)
                         for k in ("sha256", "cache_hit", "downloaded")},
            "batches": [{k: bp.get(k)
                         for k in ("sha256", "cache_hit", "downloaded")}
                        for bp in batch_provs],
        },
        "tool": f"survey-currents {_tool_version()}",
    }
    if not records:
        provenance["empty_reason"] = (
            "no USGS streamgages with daily data matched bbox "
            f"{list(box)} for {d0.isoformat()}…{d1.isoformat()} "
            f"(parameters {list(params)}). NWIS is US-only; bboxes "
            "outside USGS coverage legitimately return no sites.")

    field = GageField(records=records, bbox=box, start=d0, end=d1,
                      parameters=params, units="native", source="usgs",
                      provenance=provenance)
    if units == "si":
        field = field.to_si()
    return field

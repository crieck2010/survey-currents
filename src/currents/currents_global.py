"""Global ocean-current acquisition: NASA PODAAC OSCAR v2.0 + CMEMS.

Two complementary global surface-current sources (both verified 2026-09-26;
see docs/DATA_SOURCES.md):

* :func:`fetch_oscar` — NASA PODAAC OSCAR v2.0 (Ocean Surface Current
  Analyses Real-time): daily-averaged surface currents, 1993-present,
  0.25-degree global grid, variables ``u``/``v`` (m/s, east/north
  positive). Three latency tiers live in separate NASA CMR collections
  (Final / Interim / Near-Real-Time); the fetch picks the right tier per
  date with a pure, offline-testable rule. Granule names are
  deterministic AND verified through the keyless NASA CMR granule
  search (same pattern as the MUR adapter). Downloads go through the
  Earthdata OPeNDAP endpoint — a free Earthdata Login is required;
  without credentials a :class:`CredentialsMissing` error explains the
  setup.
* :func:`fetch_cmems_currents` — Copernicus Marine Service global ocean
  physics (``global-physics-daily`` preset, ``uo``/``vo``/``thetao``,
  1/12-degree, daily) wrapped into the standard fetch signature. Reuses
  :mod:`currents.cmems` (:func:`~currents.cmems.subset_cmems` +
  :func:`~currents.cmems.parse_cmems_netcdf`); needs the
  ``copernicusmarine`` toolbox and a free CMEMS account.

Both return :class:`~currents.models.CurrentField` (u/v in m/s,
``temperature=None`` for OSCAR, potential temperature for the CMEMS
preset — carried as-is, see the provenance note), following the
:mod:`currents.sst_global` provenance conventions (exact URLs, SHA-256
digests, retrieval time).

Verified-source corrections (documented in docs/DATA_SOURCES.md):

* OSCAR is NOT on CoastWatch ERDDAP — the old ``jplOscar_LonPM180``
  dataset id 404s; it was removed. OSCAR v2.0 lives on NASA's
  Earthdata OPeNDAP behind Earthdata Login.
* NOMADS OPeNDAP is retired (Service Change Notice 25-81).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import (_coerce_datetime_utc, _download_bytes, _require_netcdf4)
from .models import CurrentField
from .sst_global import (CMR_GRANULE_SEARCH, CredentialsMissing as
                         _BaseCredentialsMissing, earthdata_credentials,
                         validate_sst_bbox, _open_nc_bytes)

DateLike = Union[_dt.date, _dt.datetime, str]

# ---------------------------------------------------------------------------
# NASA PODAAC OSCAR v2.0 constants (verified via NASA CMR 2026-09-26)
# ---------------------------------------------------------------------------

#: Latency tiers -> NASA CMR collection concept ids. "Final" is the
#: reprocessed record (1993-01-01 -> present, ~1.5-year latency);
#: "Interim" bridges the gap (~1-month latency); "NRT" is the 2-day
#: latency stream. Each tier has its own deterministic granule naming.
OSCAR_COLLECTIONS: Dict[str, str] = {
    "final": "C2098858642-POCLOUD",
    "interim": "C2102959417-POCLOUD",
    "nrt": "C2102958977-POCLOUD",
}

#: Per-tier granule filename stems; the full title is
#: ``{stem}_YYYYMMDD.nc`` (verified live via CMR).
OSCAR_GRANULE_STEMS: Dict[str, str] = {
    "final": "oscar_currents_final",
    "interim": "oscar_currents_interim",
    "nrt": "oscar_currents_nrt",
}

#: Earthdata OPeNDAP service-URL pattern (REQUIRES Earthdata Login;
#: unauthenticated requests 302 to urs.earthdata.nasa.gov).
OSCAR_OPENDAP_PATTERN = (
    "https://opendap.earthdata.nasa.gov/collections/"
    "{collection_id}/granules/{granule_title}"
)

#: Grid geometry (verified live): 0.25-degree, lon 0..359.75 (0-360
#: convention), lat -89.75..89.75, daily averages, time units
#: "days since 1990-1-1".
OSCAR_RES = 0.25
OSCAR_LON_MIN, OSCAR_LON_MAX = 0.0, 359.75
OSCAR_LAT_MIN, OSCAR_LAT_MAX = -89.75, 89.75
OSCAR_NLON, OSCAR_NLAT = 1440, 720
OSCAR_VAR_U, OSCAR_VAR_V = "u", "v"
OSCAR_FILL = -999.0
#: NOTE the unusual dimension order: (time, longitude, latitude), not
#: the (time, latitude, longitude) most gridded products use. The parse
#: step asserts this and transposes to the canonical (nt, ny, nx).
OSCAR_DIM_ORDER = ("time", "longitude", "latitude")
#: The Final record starts here (also the floor of the collection rule).
OSCAR_START = _dt.date(1993, 1, 1)
#: Collection-pick latencies (days before "today").
OSCAR_FINAL_MAX_AGE_DAYS = 540   # older than this -> Final tier
OSCAR_INTERIM_MAX_AGE_DAYS = 45  # older than this -> Interim tier, else NRT

#: CMEMS preset + variables for the currents wrapper.
CMEMS_CURRENTS_PRESET = "global-physics-daily"
CMEMS_CURRENTS_VARIABLES = ("uo", "vo", "thetao")


class CredentialsMissing(_BaseCredentialsMissing):
    """Earthdata Login credentials are required but were not found.

    NASA's OPeNDAP endpoint (which serves OSCAR v2.0) requires a free
    Earthdata Login account — unauthenticated requests are redirected
    to the login page (HTTP 302). Provide credentials one of two ways:

    1. Environment variables (best for CI / servers)::

           export EARTHDATA_USERNAME="your_username"
           export EARTHDATA_PASSWORD="your_password"

    2. A ``~/.netrc`` entry (best for interactive use)::

           machine opendap.earthdata.nasa.gov
           login your_username
           password your_password

    Register for free at https://urs.earthdata.nasa.gov/users/new.
    """


# ---------------------------------------------------------------------------
# OSCAR collection picking (pure function — offline-testable)
# ---------------------------------------------------------------------------


def oscar_collection_for(day: DateLike,
                         today: Optional[DateLike] = None) -> str:
    """Pick the OSCAR latency tier (``"final"``/``"interim"``/``"nrt"``)
    for ``day``.

    Rule (verified against the three CMR collection records 2026-09-26):

    * ``day`` older than ``today - 540`` days -> ``"final"``
      (reprocessed record, 1993-01-01 -> present, ~1.5-year latency)
    * ``day`` older than ``today - 45`` days -> ``"interim"``
      (2020-01-01 -> present, ~1-month latency)
    * otherwise -> ``"nrt"`` (2021-01-01 -> present, ~2-day latency)

    Raises :class:`ValueError` for days before 1993-01-01 (the record
    floor) or after ``today``.
    """
    from .glsea import _coerce_date
    d = _coerce_date(day)
    t = _coerce_date(today) if today is not None else _dt.date.today()
    if d < OSCAR_START:
        raise ValueError(
            f"OSCAR v2.0 starts {OSCAR_START.isoformat()}; "
            f"day {d.isoformat()} is before the record")
    if d > t:
        raise ValueError(
            f"day {d.isoformat()} is after today ({t.isoformat()})")
    age = (t - d).days
    if age > OSCAR_FINAL_MAX_AGE_DAYS:
        return "final"
    if age > OSCAR_INTERIM_MAX_AGE_DAYS:
        return "interim"
    return "nrt"


def oscar_granule_title(collection: str, day: DateLike) -> str:
    """Deterministic OSCAR granule title, e.g.
    ``oscar_currents_final_20240115.nc`` (verified live via CMR).
    """
    from .glsea import _coerce_date
    if collection not in OSCAR_COLLECTIONS:
        raise ValueError(
            f"unknown OSCAR collection {collection!r}; "
            f"expected one of {sorted(OSCAR_COLLECTIONS)}")
    d = _coerce_date(day)
    return f"{OSCAR_GRANULE_STEMS[collection]}_{d.strftime('%Y%m%d')}.nc"


def oscar_service_url(collection: str, day: DateLike) -> str:
    """Base Earthdata OPeNDAP service URL for one OSCAR granule."""
    collection_id = OSCAR_COLLECTIONS[collection]  # KeyError -> honest
    return OSCAR_OPENDAP_PATTERN.format(
        collection_id=collection_id,
        granule_title=oscar_granule_title(collection, day))


# ---------------------------------------------------------------------------
# keyless CMR granule discovery (same pattern as the MUR adapter)
# ---------------------------------------------------------------------------


def cmr_search_oscar_granules(collection: str,
                             start: DateLike, end: DateLike,
                             page_size: int = 200) -> List[Dict[str, Any]]:
    """Discover real OSCAR granules via the public NASA CMR API.

    Keyless — discovery needs no Earthdata credentials (only the data
    download does). Each result is a dict with ``title`` (the
    authoritative granule name), ``time_start``/``time_end`` (ISO),
    ``opendap_url`` (the CMR-advertised OPeNDAP data link, if any), and
    ``links`` (the raw CMR link list).
    """
    if collection not in OSCAR_COLLECTIONS:
        raise ValueError(
            f"unknown OSCAR collection {collection!r}; "
            f"expected one of {sorted(OSCAR_COLLECTIONS)}")
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    if page_size < 1:
        raise ValueError(f"page_size must be >= 1, got {page_size}")
    params = urllib.parse.urlencode({
        "collection_concept_id": OSCAR_COLLECTIONS[collection],
        "temporal": f"{d0.date().isoformat()}T00:00:00Z,"
                    f"{d1.date().isoformat()}T23:59:59Z",
        "page_size": page_size,
        "sort_key": "start_date",
    })
    url = f"{CMR_GRANULE_SEARCH}?{params}"
    payload = _download_bytes(url)
    try:
        data = json.loads(payload.decode("utf-8"))
    except ValueError as exc:
        raise RuntimeError(
            f"CMR granule search returned non-JSON ({exc}).") from exc
    granules: List[Dict[str, Any]] = []
    for entry in data.get("feed", {}).get("entry", []):
        links = entry.get("links") or []
        opendap_url = None
        for link in links:
            href = str(link.get("href", ""))
            if ("opendap" in href.lower()
                    or "OPENDAP" in str(link.get("rel", ""))):
                opendap_url = href
                break
        granules.append({
            "title": str(entry.get("title", "")),
            "time_start": entry.get("time_start"),
            "time_end": entry.get("time_end"),
            "opendap_url": opendap_url,
            "links": links,
        })
    return granules


def oscar_match_granule(day: _dt.date,
                        discovered: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Match one sampled ``day`` to a CMR-discovered OSCAR granule.

    Matches on the granule title's ``YYYYMMDD`` stamp, falling back to
    the ``time_start`` date. Raises :class:`RuntimeError` when no
    granule matches — an honest gap, never an invented granule name.
    """
    stamp = day.strftime("%Y%m%d")
    for granule in discovered:
        title = str(granule.get("title", ""))
        if re.search(rf"^oscar_currents_(final|interim|nrt)_{stamp}\.nc$",
                     title):
            return dict(granule)
    for granule in discovered:
        time_start = str(granule.get("time_start", ""))
        if time_start[:10] == day.isoformat():
            return dict(granule)
    raise RuntimeError(
        f"CMR granule search returned no OSCAR granule for {day.isoformat()}; "
        "cannot fetch that date without inventing a granule name. "
        "Try a narrower date range or check the OSCAR collection status.")


def _oscar_service_url(granule: Dict[str, Any], collection: str,
                       day: _dt.date) -> Tuple[str, str]:
    """Base OPeNDAP URL for a CMR-discovered granule -> ``(url, source)``.

    Prefers the OPeNDAP link CMR advertises on the granule; otherwise
    falls back to the documented Earthdata OPeNDAP URL pattern for the
    collection (verified live via CMR 2026-09-26). ``source`` is
    ``"cmr-link"`` or ``"constructed"`` and is recorded in provenance.
    """
    link = granule.get("opendap_url")
    if link:
        return str(link).rstrip("/"), "cmr-link"
    return oscar_service_url(collection, day), "constructed"


# ---------------------------------------------------------------------------
# bbox -> OSCAR 0-360 longitude windows + grid index windows
# ---------------------------------------------------------------------------


def oscar_lon_windows(bbox: Sequence[float]) -> List[Tuple[float, float]]:
    """Convert a -180..180 bbox to OSCAR 0-360 longitude request windows.

    Mirrors :func:`currents.sst_global.oisst_lon_windows` (the OSCAR grid
    is the same 0-360 convention, floor 0.0, ceiling 359.75):

    * **Antimeridian-crossing bboxes** (e.g. ``(170, ., -170, .)``) wrap
      to a single 0-360 window (``(170.0, 190.0)``).
    * **Windows that cross 360° in 0-360 space** (e.g. ``(-170, ., 170, .)``
      becomes ``(190.0, 530.0)``) are split into two windows:
      ``(190.0, 359.75)`` and ``(0.0, 170.0)``.
    * **Full-globe bboxes** (span >= 359.9°) request the entire grid.
    """
    minx, _, maxx, _ = validate_sst_bbox(bbox)
    span = maxx - minx
    if span < 0:
        span += 360.0
    if span >= 359.9:
        return [(OSCAR_LON_MIN, OSCAR_LON_MAX)]
    lo, hi = minx % 360.0, maxx % 360.0
    if hi <= lo:
        hi += 360.0
    if hi <= OSCAR_LON_MAX:
        return [(max(lo, OSCAR_LON_MIN), hi)]
    return [(max(lo, OSCAR_LON_MIN), OSCAR_LON_MAX),
            (OSCAR_LON_MIN, hi - 360.0)]


def oscar_index_windows(bbox: Sequence[float]
                        ) -> List[Tuple[int, int, int, int]]:
    """Map a -180..180 bbox onto OSCAR grid index windows
    ``(lon0, lon1, lat0, lat1)`` (inclusive), one per 0-360 lon window.

    The OSCAR grid is a documented regular 0.25-degree grid
    (lon 0..359.75, lat -89.75..89.75); indices are computed
    arithmetically and clamped to the axis sizes. The parse step
    verifies the on-wire dimension order independently.
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)

    def lon_idx(lon360: float) -> int:
        return max(0, min(OSCAR_NLON - 1, int(round(lon360 / OSCAR_RES))))

    def lat_idx(lat: float) -> int:
        return max(0, min(OSCAR_NLAT - 1,
                          int(round((lat - OSCAR_LAT_MIN) / OSCAR_RES))))

    j_lo = lat_idx(max(miny, OSCAR_LAT_MIN))
    j_hi = lat_idx(min(maxy, OSCAR_LAT_MAX))
    lat_lo, lat_hi = (j_lo, j_hi) if j_lo <= j_hi else (j_hi, j_lo)
    windows = []
    for lo360, hi360 in oscar_lon_windows(bbox):
        i0, i1 = lon_idx(lo360), lon_idx(hi360)
        lon_lo, lon_hi = (i0, i1) if i0 <= i1 else (i1, i0)
        windows.append((lon_lo, lon_hi, lat_lo, lat_hi))
    return windows


def oscar_subset_urls(day: DateLike, collection: str, bbox: Sequence[float],
                      base_url: Optional[str] = None) -> List[str]:
    """Build the constrained OPeNDAP URL(s) for one OSCAR granule.

    One URL per 0-360 longitude window (antimeridian-safe); each URL
    constrains ``u`` and ``v`` to the single daily timestep and the bbox
    index window, honoring the on-wire ``(time, longitude, latitude)``
    dimension order::

        {base}?u[0:1:0][{i0}:1:{i1}][{j0}:1:{j1}],v[0:1:0][{i0}:1:{i1}][{j0}:1:{j1}]

    ``base_url`` (from :func:`_oscar_service_url` on a CMR-discovered
    granule) overrides the constructed :func:`oscar_service_url`.
    """
    from .glsea import _coerce_date
    day = _coerce_date(day)
    base = base_url or oscar_service_url(collection, day)
    urls = []
    for i0, i1, j0, j1 in oscar_index_windows(bbox):
        # NOTE: OSCAR granule titles already end in ".nc", so the DAP2
        # constraint is appended directly (unlike the MUR adapter, whose
        # titles carry no suffix).
        constraint = (
            f"{OSCAR_VAR_U}[0:1:0][{i0}:1:{i1}][{j0}:1:{j1}],"
            f"{OSCAR_VAR_V}[0:1:0][{i0}:1:{i1}][{j0}:1:{j1}]"
        )
        urls.append(f"{base}?{constraint}")
    return urls


# ---------------------------------------------------------------------------
# authenticated download + NetCDF parsing
# ---------------------------------------------------------------------------


def _oscar_opener() -> urllib.request.OpenerDirector:
    """Build an authenticated Earthdata urllib opener.

    Raises :class:`CredentialsMissing` when no credentials are found.
    """
    creds = earthdata_credentials()
    if creds is None:
        raise CredentialsMissing(
            "OSCAR v2.0 needs Earthdata Login credentials, but none were "
            "found. Set EARTHDATA_USERNAME/EARTHDATA_PASSWORD or add a "
            "~/.netrc entry for opendap.earthdata.nasa.gov (register free "
            "at https://urs.earthdata.nasa.gov/users/new).")
    user, pw = creds
    mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    mgr.add_password(None, "https://opendap.earthdata.nasa.gov", user, pw)
    auth = urllib.request.HTTPBasicAuthHandler(mgr)
    return urllib.request.build_opener(auth)


def _oscar_get(opener: urllib.request.OpenerDirector, url: str,
               context: str) -> bytes:
    req = urllib.request.Request(
        url, headers={"User-Agent": "survey-currents/0.5.0"})
    try:
        with opener.open(req, timeout=180) as resp:
            return resp.read()
    except urllib.error.HTTPError as exc:
        if exc.code == 401:
            raise CredentialsMissing(
                f"Earthdata Login rejected the credentials for {context} "
                "(HTTP 401). Check EARTHDATA_USERNAME/EARTHDATA_PASSWORD or "
                "the ~/.netrc entry for opendap.earthdata.nasa.gov.") from exc
        raise RuntimeError(
            f"OSCAR {context} failed (HTTP {exc.code}: {exc.reason}).") from exc
    except Exception as exc:
        raise RuntimeError(
            f"OSCAR {context} failed ({type(exc).__name__}: {exc}).") from exc


def _parse_oscar_bytes(payload: bytes) -> Tuple[np.ndarray, np.ndarray,
                                                np.ma.MaskedArray,
                                                np.ma.MaskedArray]:
    """Parse one OSCAR subset NetCDF payload.

    Returns ``(lats, lons, u[1, ny, nx], v[1, ny, nx])`` — the on-wire
    ``(time, longitude, latitude)`` order is asserted (an honest failure
    if PODAAC ever reorders the axes) and transposed to the canonical
    ``(time, lat, lon)``. Longitudes are normalized to -180..180 and
    sorted increasing; latitudes are flipped to increasing if needed.
    Fill value -999.0 is masked.
    """
    nc = _require_netcdf4()
    with _open_nc_bytes(payload) as ds:
        for name in (OSCAR_VAR_U, OSCAR_VAR_V):
            if name not in ds.variables:
                raise ValueError(f"OSCAR payload missing variable {name!r}")
        uvar = ds.variables[OSCAR_VAR_U]
        vvar = ds.variables[OSCAR_VAR_V]
        for var, name in ((uvar, OSCAR_VAR_U), (vvar, OSCAR_VAR_V)):
            if tuple(var.dimensions) != OSCAR_DIM_ORDER:
                raise ValueError(
                    f"OSCAR variable {name!r} has dimensions "
                    f"{tuple(var.dimensions)}; expected {OSCAR_DIM_ORDER} "
                    "(time, longitude, latitude)")
        lon360 = np.asarray(ds.variables["longitude"][:], dtype=float)
        lats = np.asarray(ds.variables["latitude"][:], dtype=float)
        # On-wire order is (time, longitude, latitude) -> transpose the
        # spatial axes to (time, lat, lon).
        u = np.ma.asarray(uvar[:], dtype=float).transpose(0, 2, 1)
        v = np.ma.asarray(vvar[:], dtype=float).transpose(0, 2, 1)
    u = np.ma.masked_values(np.ma.masked_invalid(u), OSCAR_FILL,
                            rtol=1e-5, atol=1e-3)
    v = np.ma.masked_values(np.ma.masked_invalid(v), OSCAR_FILL,
                            rtol=1e-5, atol=1e-3)
    if lats[0] > lats[-1]:
        lats = lats[::-1]
        u = u[:, ::-1, :]
        v = v[:, ::-1, :]
    lons = ((lon360 + 180.0) % 360.0) - 180.0
    order = np.argsort(lons, kind="stable")
    return lats, lons[order], u[:, :, order], v[:, :, order]


def _concat_oscar_parts(
        parts: List[Tuple[np.ndarray, np.ndarray,
                          np.ma.MaskedArray, np.ma.MaskedArray]]
) -> Tuple[np.ndarray, np.ndarray, np.ma.MaskedArray, np.ma.MaskedArray]:
    """Concatenate parsed OSCAR lon-window parts along longitude."""
    lats = parts[0][0]
    lons = np.concatenate([p[1] for p in parts])
    u = np.ma.concatenate([p[2] for p in parts], axis=2)
    v = np.ma.concatenate([p[3] for p in parts], axis=2)
    order = np.argsort(lons, kind="stable")
    lons_sorted = lons[order]
    # Antimeridian splits can share the seam meridian (180 == -180);
    # drop exact duplicates so the grid stays strictly increasing.
    _, unique_idx = np.unique(lons_sorted, return_index=True)
    keep = order[np.sort(unique_idx)]
    return lats, lons[keep], u[:, :, keep], v[:, :, keep]


def _oscar_sample_dates(d0: _dt.datetime, d1: _dt.datetime,
                        stride_days: int) -> List[_dt.date]:
    days: List[_dt.date] = []
    cursor = d0.date()
    while cursor <= d1.date():
        days.append(cursor)
        cursor += _dt.timedelta(days=stride_days)
    return days


def _validate_oscar_dates(d0: _dt.datetime, d1: _dt.datetime) -> None:
    if d0 > d1:
        raise ValueError(f"start {d0.date()} is after end {d1.date()}")
    if d0.date() < OSCAR_START:
        raise ValueError(
            f"OSCAR v2.0 starts {OSCAR_START.isoformat()}; "
            f"start {d0.date().isoformat()} is before the record")
    if d1.date() > _dt.date.today():
        raise ValueError(
            f"end {d1.date().isoformat()} is in the future")


# ---------------------------------------------------------------------------
# public fetch: OSCAR
# ---------------------------------------------------------------------------


def fetch_oscar(bbox: Sequence[float], start: DateLike, end: DateLike,
                stride_days: int = 5) -> CurrentField:
    """Fetch NASA PODAAC OSCAR v2.0 daily surface currents for ``bbox``.

    For each sampled day the latency tier is picked with
    :func:`oscar_collection_for` (Final / Interim / NRT), granules are
    discovered through the keyless NASA CMR API (never invented), and
    each granule is subset through the authenticated Earthdata OPeNDAP
    endpoint.

    Args:
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180 degrees
            (ascending; antimeridian-crossing boxes wrap internally).
        start, end: dates/datetimes/ISO strings within 1993-01-01..today;
            one granule per ``stride_days``.
        stride_days: granule sampling stride, >= 1.

    Returns:
        :class:`CurrentField` (u/v in m/s, ``temperature=None`` — OSCAR
        is a currents-only product) with provenance: exact OPeNDAP URLs,
        per-payload SHA-256, retrieval time, the three CMR collection
        ids, the pick rule, and the per-date collection used.

    Raises:
        CredentialsMissing: no Earthdata Login credentials found (or the
            server rejected them with HTTP 401).
        ValueError: invalid bbox / dates / stride.
        RuntimeError: CMR search failures, a sampled day with no
            discovered granule, or network/NetCDF failures (with context).
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_oscar_dates(d0, d1)
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")

    opener = _oscar_opener()  # raises CredentialsMissing early
    days = _oscar_sample_dates(d0, d1, stride_days)
    today = _dt.date.today()

    # Collection per date, then one keyless CMR discovery per tier.
    collections = [oscar_collection_for(day, today) for day in days]
    by_collection: Dict[str, List[int]] = {}
    for i, coll in enumerate(collections):
        by_collection.setdefault(coll, []).append(i)
    discovered: Dict[str, List[Dict[str, Any]]] = {}
    for coll, idxs in by_collection.items():
        span = [days[i] for i in idxs]
        discovered[coll] = cmr_search_oscar_granules(
            coll, min(span), max(span))

    lats_all = lons_all = None
    u_frames: List[np.ma.MaskedArray] = []
    v_frames: List[np.ma.MaskedArray] = []
    times: List[str] = []
    urls: List[str] = []
    digests: List[str] = []
    titles: List[str] = []
    url_sources: List[str] = []
    date_records: List[Dict[str, str]] = []
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()
    for day, coll in zip(days, collections):
        granule = oscar_match_granule(day, discovered[coll])
        title = str(granule["title"])
        titles.append(title)
        granule_url, url_source = _oscar_service_url(granule, coll, day)
        url_sources.append(url_source)
        date_records.append({
            "date": day.isoformat(),
            "collection": coll,
            "collection_id": OSCAR_COLLECTIONS[coll],
            "granule_title": title,
        })
        # Frame timestamp: the daily granule's own date at 00:00 UTC
        # (OSCAR granules are daily averages).
        stamp = _dt.datetime.combine(
            day, _dt.time(0, 0), tzinfo=_dt.timezone.utc).isoformat()
        day_u: List[np.ma.MaskedArray] = []
        day_v: List[np.ma.MaskedArray] = []
        day_lons: List[np.ndarray] = []
        for url in oscar_subset_urls(day, coll, (minx, miny, maxx, maxy),
                                     base_url=granule_url):
            urls.append(url)
            payload = _oscar_get(opener, url, f"granule {title}")
            digests.append(hashlib.sha256(payload).hexdigest())
            try:
                lats, lons, u, v = _parse_oscar_bytes(payload)
            except ValueError:
                raise
            except Exception as exc:
                raise RuntimeError(
                    f"OSCAR granule {title} parse failed "
                    f"({type(exc).__name__}: {exc}). netCDF4 is required "
                    "(pip install netCDF4)."
                ) from exc
            if lats_all is None:
                lats_all = lats
            day_u.append(u)
            day_v.append(v)
            day_lons.append(lons)
        _, day_lons_sorted, u_day, v_day = _concat_oscar_parts(
            [(lats_all, ln, uu, vv)
             for ln, uu, vv in zip(day_lons, day_u, day_v)])
        if lons_all is None:
            lons_all = day_lons_sorted
        u_frames.append(u_day)
        v_frames.append(v_day)
        times.append(stamp)

    u_cube = np.ma.concatenate(u_frames, axis=0)
    v_cube = np.ma.concatenate(v_frames, axis=0)
    combined = hashlib.sha256()
    for d in digests:
        combined.update(d.encode())
    tiers_used = sorted(set(collections))
    return CurrentField(
        u=np.ma.filled(u_cube, np.nan), v=np.ma.filled(v_cube, np.nan),
        temperature=None, times=times, lats=lats_all, lons=lons_all,
        source="podaac/OSCAR_L4_OC_V2.0",
        provenance={
            "dataset": "NASA PODAAC OSCAR v2.0 "
                       "(Ocean Surface Current Analyses Real-time, "
                       "daily-averaged surface currents)",
            "collections": dict(OSCAR_COLLECTIONS),
            "collections_used": tiers_used,
            "collection_pick_rule": (
                "date older than today-540d -> final "
                "(1993-01-01->present, ~1.5yr latency); older than "
                "today-45d -> interim (2020-01-01->present, ~1mo latency); "
                "else nrt (2021-01-01->present, ~2d latency)"),
            "dates": date_records,
            "cmr_search": {
                "endpoint": CMR_GRANULE_SEARCH,
                "per_collection": {
                    coll: {
                        "collection_concept_id": OSCAR_COLLECTIONS[coll],
                        "n_results": len(discovered[coll]),
                    } for coll in tiers_used
                },
            },
            "granule_titles": titles,
            "opendap_url_source": url_sources,
            "service_urls": urls,
            "granule_sha256": digests,
            "combined_sha256": combined.hexdigest(),
            "retrieved_at": retrieved_at,
            "bbox_requested": [minx, miny, maxx, maxy],
            "time_requested": [d0.date().isoformat(), d1.date().isoformat()],
            "stride_days": stride_days,
            "n_granules": len(days),
            "units": "m/s (u eastward, v northward)",
            "grid": f"{OSCAR_RES}-degree global ocean grid "
                    f"(lon {OSCAR_LON_MIN}..{OSCAR_LON_MAX}, "
                    f"lat {OSCAR_LAT_MIN}..{OSCAR_LAT_MAX})",
            "corrections": [
                "OSCAR is NOT on CoastWatch ERDDAP (the old "
                "jplOscar_LonPM180 dataset id 404s — removed); v2.0 is "
                "served from Earthdata OPeNDAP behind Earthdata Login.",
                "NOMADS OPeNDAP is retired (Service Change Notice 25-81); "
                "it is not used.",
            ],
        },
    )


# ---------------------------------------------------------------------------
# public fetch: CMEMS global physics, wrapped to the standard signature
# ---------------------------------------------------------------------------


def fetch_cmems_currents(bbox: Sequence[float], start: DateLike,
                         end: DateLike, stride_days: int = 1,
                         work_dir: Optional[str] = None) -> CurrentField:
    """Fetch CMEMS global ocean physics daily currents for ``bbox``.

    Wraps :func:`currents.cmems.subset_cmems` +
    :func:`currents.cmems.parse_cmems_netcdf` (preset
    ``global-physics-daily``: ``uo``/``vo``/``thetao``, 1/12-degree,
    daily) into the standard ``(bbox, start, end, stride)`` fetch
    signature: one NetCDF is downloaded for the whole range, then
    timesteps are stride-selected.

    Args:
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180 degrees.
        start, end: dates/datetimes/ISO strings (start <= end, end not
            in the future).
        stride_days: keep every Nth daily timestep, >= 1.
        work_dir: directory for the downloaded NetCDF (defaults to a
            fresh temp dir; the NetCDF + its provenance sidecar stay
            there so the fetch is auditable).

    Returns:
        :class:`CurrentField` with ``temperature`` = the preset's
        ``thetao`` (potential temperature in degC, carried as-is — it is
        not a foundation SST).

    Raises:
        ImportError: the ``copernicusmarine`` toolbox is not installed.
        RuntimeError: CMEMS credentials missing (free account required)
            or the subset/parse failed.
        ValueError: invalid bbox / dates / stride.
    """
    from .cmems import parse_cmems_netcdf, subset_cmems
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    if d0 > d1:
        raise ValueError(f"start {d0.date()} is after end {d1.date()}")
    if d1.date() > _dt.date.today():
        raise ValueError(
            f"end {d1.date().isoformat()} is in the future")
    if stride_days < 1:
        raise ValueError(f"stride_days must be >= 1, got {stride_days}")
    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="survey-currents-cmems_")
    os.makedirs(work_dir, exist_ok=True)

    start_iso = d0.date().isoformat() + "T00:00:00"
    end_iso = d1.date().isoformat() + "T23:59:59"
    path = subset_cmems(
        CMEMS_CURRENTS_PRESET, (minx, miny, maxx, maxy),
        start_iso, end_iso, work_dir,
        variables=list(CMEMS_CURRENTS_VARIABLES))
    field = parse_cmems_netcdf(path)
    if not field.times:
        raise RuntimeError(
            f"CMEMS subset returned no timesteps (file {path}).")

    keep = list(range(0, len(field.times), stride_days))
    temperature = (None if field.temperature is None
                   else field.temperature[keep])
    out = CurrentField(
        u=field.u[keep], v=field.v[keep], temperature=temperature,
        times=[field.times[i] for i in keep],
        lats=field.lats, lons=field.lons, crs=field.crs,
        source=f"cmems:{CMEMS_CURRENTS_PRESET}",
        provenance={
            **dict(field.provenance),
            "dataset": "Copernicus Marine Service global ocean physics "
                       "analysis+forecast, daily (1/12-degree)",
            "preset": CMEMS_CURRENTS_PRESET,
            "variables": {"u": "uo", "v": "vo",
                          "temperature": "thetao"},
            "netcdf_path": path,
            "work_dir": work_dir,
            "time_requested": [start_iso, end_iso],
            "stride_days": stride_days,
            "n_timesteps_total": len(field.times),
            "n_timesteps_kept": len(keep),
            "units": "m/s (uo eastward, vo northward); degC (thetao)",
            "temperature_note": "thetao is potential temperature (degC), "
                                "carried as-is — it is not a foundation SST",
            "retrieved_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        },
    )
    return out


# ---------------------------------------------------------------------------
# synthetic fixtures (tests, CLI demo, offline use)
# ---------------------------------------------------------------------------


def oscar_synthetic(nt: int = 4, ny: int = 6, nx: int = 8,
                    lats: Sequence[float] = (-30.0, 30.0),
                    lons: Sequence[float] = (-60.0, 60.0),
                    start: str = "2024-01-01T00:00:00+00:00",
                    step_days: int = 5,
                    seed: int = 13) -> CurrentField:
    """Deterministic synthetic OSCAR-shaped field: gyre currents, no temperature.

    Stdlib+numpy only. Used by the test suite and the offline demo.
    """
    field = CurrentField.synthetic(
        nt=nt, ny=ny, nx=nx, lats=lats, lons=lons, start=start,
        step_hours=24 * step_days, seed=seed,
        source="podaac/OSCAR_L4_OC_V2.0 (synthetic)")
    field.temperature = None  # OSCAR is a currents-only product
    field.provenance.update(
        {"synthetic": True, "seed": seed, "adapter": "fetch_oscar"})
    return field


def cmems_currents_synthetic(nt: int = 4, ny: int = 6, nx: int = 8,
                             lats: Sequence[float] = (-30.0, 30.0),
                             lons: Sequence[float] = (-60.0, 60.0),
                             start: str = "2024-01-01T00:00:00+00:00",
                             step_days: int = 1,
                             seed: int = 17) -> CurrentField:
    """Deterministic synthetic CMEMS-shaped field: gyre + temperature.

    Stdlib+numpy only. Used by the test suite and the offline demo.
    """
    field = CurrentField.synthetic(
        nt=nt, ny=ny, nx=nx, lats=lats, lons=lons, start=start,
        step_hours=24 * step_days, seed=seed,
        source=f"cmems:{CMEMS_CURRENTS_PRESET} (synthetic)")
    field.provenance.update(
        {"synthetic": True, "seed": seed, "adapter": "fetch_cmems_currents"})
    return field


# ---------------------------------------------------------------------------
# offline demo
# ---------------------------------------------------------------------------


def main_demo() -> None:
    """Print a small offline summary (no network). Mirrors sst_global.main_demo."""
    for label, make in (("oscar", oscar_synthetic),
                        ("cmems-currents", cmems_currents_synthetic)):
        f = make(nt=3)
        print(f"[{label}] synthetic {f.u.shape} m/s, "
              f"t0={f.times[0]}, mean speed={f.zonal_mean(0, 'speed'):.3f} m/s")


__all__ = [
    "CredentialsMissing",
    "DateLike",
    "OSCAR_COLLECTIONS",
    "OSCAR_GRANULE_STEMS",
    "OSCAR_OPENDAP_PATTERN",
    "OSCAR_RES",
    "OSCAR_LON_MIN",
    "OSCAR_LON_MAX",
    "OSCAR_LAT_MIN",
    "OSCAR_LAT_MAX",
    "OSCAR_DIM_ORDER",
    "OSCAR_START",
    "OSCAR_FINAL_MAX_AGE_DAYS",
    "OSCAR_INTERIM_MAX_AGE_DAYS",
    "CMEMS_CURRENTS_PRESET",
    "CMEMS_CURRENTS_VARIABLES",
    "CMR_GRANULE_SEARCH",
    "oscar_collection_for",
    "oscar_granule_title",
    "oscar_service_url",
    "cmr_search_oscar_granules",
    "oscar_match_granule",
    "oscar_lon_windows",
    "oscar_index_windows",
    "oscar_subset_urls",
    "fetch_oscar",
    "fetch_cmems_currents",
    "oscar_synthetic",
    "cmems_currents_synthetic",
    "main_demo",
]

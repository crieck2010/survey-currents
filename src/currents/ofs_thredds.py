"""NOAA Operational Forecast Systems via CO-OPS THREDDS OPeNDAP subsetting.

The direct S3 full-file OFS path is dead for reel work: hourly field files
are 62-70 MB each (gigabytes per reel). The CO-OPS THREDDS server instead
serves the same OFS ``regulargrid`` NetCDFs over OPeNDAP, so a surface-u/v
(+ temperature) subset for a small bbox and a few hourly steps is
kilobytes, not megabytes — keyless, no credentials, stdlib + numpy only.

Verified live 2026-10-01: a Puget Sound bbox (-123.2, 47.2, -122.2, 48.4,
240x201 grid points) pulls u_eastward + v_northward + temp + time for one
hourly step in a single ~0.55 MB ``.dods`` response, versus ~65 MB for the
full field file (~120x smaller). Values sanity-checked: u in
[-1.40, 1.53] m/s, v in [-1.89, 1.11] m/s, temp in [9.77, 15.23] degC,
timestamp 2026-10-01T15:00:00+00:00 matching the file's own ``Times``.

**How it works**

* The THREDDS catalog at ``/thredds/catalog/NOAA/<MODEL>/MODELS/<yyyy>/<mm>/
  <dd>/catalog.xml`` lists that day's per-hour ``regulargrid`` NetCDFs,
  e.g. ``sscofs.t15z.20261001.regulargrid.n006.nc``. Each file holds exactly
  one timestep (``time = 1``).
* Filenames encode the valid time. Verified live 2026-10-01 across SSCOFS,
  CBOFS, NGOFS2, GOMOFS, DBOFS, SFBOFS and WCOFS: forecast files
  ``f{HHH}`` are valid at ``cycle + HHH`` hours; nowcast files ``n{HHH}``
  are valid at ``cycle - span + HHH`` where ``span`` is calibrated per
  model by probing one nowcast file's own ``time`` variable (6 h for most
  models, 24 h for WCOFS). The valid time read back from each fetched file
  is cross-checked against the requested hour — a served step that
  mismatches its requested hour is refused, never mislabeled.
* The ``regulargrid`` schema is ``u_eastward`` / ``v_northward``
  (m s-1), ``temp`` (degC), on a regular lat/lon grid stored as 2-D
  ``Latitude``/``Longitude`` arrays (verified row/column-constant per
  model at fetch time; a non-regular grid raises ``ValueError``).
  ``Depth[0]`` is the surface (0.0 m, verified live 2026-10-01).
* Subsets use DAP2 constraint expressions on ``/thredds/dodsC/`` with the
  binary ``.dods`` response, parsed with a DDS-driven XDR decoder
  (stdlib ``struct`` — the server reorders coordinate variables first, so
  the parser follows the response DDS, not the request order).

**Honest limits** (verified live 2026-10-01):

* THREDDS keeps roughly the **last 31 days** of OFS output. A date with no
  day catalog, or an hour with no matching file, raises
  :class:`UnavailableRangeError` naming the exact URL — never silent,
  never padded.
* 12 of the 15 OFS models serve ``regulargrid`` files with this schema:
  SSCOFS, CBOFS, WCOFS, NGOFS2, GOMOFS, DBOFS, SFBOFS, LEOFS, LMHOFS,
  LOOFS, LSOFS, CIOFS. NYOFS and SJROFS have no ``regulargrid`` files;
  TBOFS is stale (2021). Unknown codes raise ``ValueError``.
* Each requested hour is one HTTP request (~0.2-1 MB for typical reel
  bboxes). The caller owns the frame budget.
* Land cells carry the file ``_FillValue`` and are returned as NaN;
  ``provenance["water_fraction"]`` records the share of wet points.
"""

from __future__ import annotations

import datetime as _dt
import re
import struct
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .models import CurrentField

DateLike = Union[_dt.date, _dt.datetime, str]

#: CO-OPS THREDDS root (keyless HTTPS).
THREDDS_BASE = "https://opendap.co-ops.nos.noaa.gov/thredds"

#: OFS models verified live 2026-10-01 to serve per-hour ``regulargrid``
#: NetCDFs (u_eastward / v_northward / temp on a regular lat/lon grid) on
#: CO-OPS THREDDS. NYOFS/SJROFS have no regulargrid files; TBOFS is stale.
OFS_THREDDS_MODELS = (
    "SSCOFS", "CBOFS", "WCOFS", "NGOFS2", "GOMOFS", "DBOFS",
    "SFBOFS", "LEOFS", "LMHOFS", "LOOFS", "LSOFS", "CIOFS",
)

#: Approximate THREDDS retention of OFS output, in days including today
#: (verified live 2026-10-01: 2026/09 has days 21-30, 2026/10 has day 01).
THREDDS_RETENTION_DAYS = 31

#: regulargrid variable names.
U_VAR = "u_eastward"
V_VAR = "v_northward"
TEMP_VAR = "temp"
LAT_VAR = "Latitude"
LON_VAR = "Longitude"
TIME_VAR = "time"

#: ``{code}.t{cc}z.{yyyymmdd}.regulargrid.{n|f}{hhh}.nc``
_FILENAME_RE = re.compile(
    r"^([a-z0-9]+)\.t(\d{2})z\.(\d{8})\.regulargrid\.([nf])(\d{3})\.nc$"
)

_CATALOG_NS = {"t": "http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0"}
_XLINK_HREF = "{http://www.w3.org/1999/xlink}href"

_DDS_VAR_RE = re.compile(
    r"(Float32|Float64|Int32)\s+(\w+)((?:\[\w+ = \d+\])+);"
)
_DDS_DIM_RE = re.compile(r"\[(\w+) = (\d+)\]")
_DAS_TIME_UNITS_RE = re.compile(
    r"time\s*\{(?:[^{}]|\{[^{}]*\})*?String units\s+\"([^\"]+)\"", re.S
)
_DAS_FILL_RE = re.compile(
    r"u_eastward\s*\{(?:[^{}]|\{[^{}]*\})*?(?:Float32|Float64|Int32)\s+_FillValue\s+([0-9.Ee+-]+)",
    re.S,
)
_ASCII_ROW_RE = re.compile(r"^\[\d+\],\s*(.*)$")
_TIME_UNITS_RE = re.compile(
    r"seconds since (\d{4})-(\d{2})-(\d{2})[ T](\d{2}):(\d{2}):(\d{2})"
)


class UnavailableRangeError(ValueError):
    """The requested OFS model/date/hour is not on CO-OPS THREDDS.

    Raised for dates outside the ~31-day retention window, hours with no
    matching model file, bboxes outside the model domain, and HTTP errors
    from the server. The message always names the exact URL or catalog
    that was missing — never silent, never padded.
    """


# ---------------------------------------------------------------------------
# HTTP (one indirection so tests can run fully offline)
# ---------------------------------------------------------------------------


def _http_get_bytes(url: str, timeout: float = 120.0) -> bytes:
    """GET ``url`` and return the raw response bytes.

    HTTP/URL errors are wrapped in :class:`UnavailableRangeError` naming
    the exact URL.
    """
    req = urllib.request.Request(
        url, headers={"User-Agent": "survey-currents/ofs-thredds"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except (urllib.error.HTTPError, urllib.error.URLError) as exc:
        raise UnavailableRangeError(
            f"THREDDS request failed for {url}: {exc}") from exc


# ---------------------------------------------------------------------------
# Catalog / dataset URL construction
# ---------------------------------------------------------------------------


def day_catalog_url(model: str, date: _dt.date) -> str:
    """THREDDS catalog XML URL listing one model's files for one day."""
    return (f"{THREDDS_BASE}/catalog/NOAA/{model}/MODELS/"
            f"{date:%Y/%m/%d}/catalog.xml")


def dods_base_url(url_path: str) -> str:
    """OPeNDAP base URL for a catalog ``urlPath``."""
    return f"{THREDDS_BASE}/dodsC/{url_path}"


def list_regulargrid_files(model: str, date: _dt.date,
                           timeout: float = 120.0) -> List[Dict[str, str]]:
    """Parse a day catalog; return ``regulargrid`` file entries.

    Each entry is ``{"name": ..., "url_path": ...}``. A missing day
    catalog (date outside the ~31-day retention window) raises
    :class:`UnavailableRangeError`.
    """
    url = day_catalog_url(model, date)
    try:
        raw = _http_get_bytes(url, timeout=timeout)
    except UnavailableRangeError as exc:
        raise UnavailableRangeError(
            f"no THREDDS day catalog for {model} on {date} "
            f"({url}); CO-OPS THREDDS keeps roughly the last "
            f"{THREDDS_RETENTION_DAYS} days of OFS output") from exc
    root = ET.fromstring(raw)
    out = []
    for ds in root.findall(".//t:dataset", _CATALOG_NS):
        path = ds.get("urlPath")
        name = ds.get("name") or ""
        if path and "regulargrid" in name and name.endswith(".nc"):
            out.append({"name": name, "url_path": path})
    return out


def parse_ofs_filename(name: str) -> Dict[str, Any]:
    """Split a regulargrid filename into code/cycle/kind/hour.

    >>> parse_ofs_filename("sscofs.t15z.20261001.regulargrid.n006.nc")["hour"]
    6
    """
    m = _FILENAME_RE.match(name)
    if not m:
        raise ValueError(f"not an OFS regulargrid filename: {name!r}")
    code, cycle, ymd, kind, hour = m.groups()
    cycle_dt = _dt.datetime.strptime(f"{ymd} {cycle}", "%Y%m%d %H").replace(
        tzinfo=_dt.timezone.utc)
    return {"code": code, "cycle": cycle_dt, "kind": kind,
            "hour": int(hour)}


def nominal_valid_time(parsed: Dict[str, Any],
                       nowcast_span_hours: Optional[int]) -> _dt.datetime:
    """Filename-implied valid time.

    Forecast ``f{HHH}`` files are valid at ``cycle + HHH`` hours;
    nowcast ``n{HHH}`` files at ``cycle - span + HHH`` hours, where
    ``span`` is the per-model nowcast window calibrated by probing one
    nowcast file (6 h for most models, 24 h for WCOFS).
    """
    if parsed["kind"] == "f":
        return parsed["cycle"] + _dt.timedelta(hours=parsed["hour"])
    if nowcast_span_hours is None:
        raise ValueError("no nowcast files available; cannot date n-files")
    return (parsed["cycle"]
            + _dt.timedelta(hours=parsed["hour"] - nowcast_span_hours))


# ---------------------------------------------------------------------------
# DDS / DAS / ASCII helpers
# ---------------------------------------------------------------------------


def fetch_dds(url_path: str,
              timeout: float = 120.0) -> List[Tuple[str, str, List[Tuple[str, int]]]]:
    """Parse a dataset's ``.dds`` into ``(name, dtype, [(dim, size)])``."""
    raw = _http_get_bytes(dods_base_url(url_path) + ".dds", timeout=timeout)
    text = raw.decode("utf-8", "replace")
    specs = []
    for m in _DDS_VAR_RE.finditer(text):
        dtype, name, dims = m.groups()
        specs.append((name, dtype,
                      [(d, int(n)) for d, n in _DDS_DIM_RE.findall(dims)]))
    if not specs:
        raise UnavailableRangeError(
            f"empty DDS for {dods_base_url(url_path)}")
    return specs


def fetch_das(url_path: str, timeout: float = 120.0) -> str:
    """Return a dataset's ``.das`` attribute text."""
    raw = _http_get_bytes(dods_base_url(url_path) + ".das", timeout=timeout)
    return raw.decode("utf-8", "replace")


def parse_time_units(das: str) -> _dt.datetime:
    """Epoch for the ``time`` variable from DAS units text."""
    m = _DAS_TIME_UNITS_RE.search(das)
    if not m:
        raise UnavailableRangeError("no time units found in dataset DAS")
    mu = _TIME_UNITS_RE.match(m.group(1))
    if not mu:
        raise UnavailableRangeError(
            f"unsupported time units {m.group(1)!r} (want 'seconds since ...')")
    return _dt.datetime(*map(int, mu.groups()),
                        tzinfo=_dt.timezone.utc)


def parse_fill_value(das: str) -> float:
    """``_FillValue`` of ``u_eastward`` from DAS text (fallback -99999.0)."""
    m = _DAS_FILL_RE.search(das)
    return float(m.group(1)) if m else -99999.0


def fetch_ascii_values(url_path: str, var: str,
                       slices: Sequence[str],
                       timeout: float = 120.0) -> List[float]:
    """Fetch one variable slice via the ``.ascii`` response; return floats.

    ``slices`` is one ``"start:stride:stop"`` per dimension, e.g.
    ``["0:1:1552", "0:1:0"]`` for ``Latitude[0:1:1552][0:1:0]``.
    """
    ce = (var + "".join(f"[{s}]" for s in slices)
          ).replace("[", "%5B").replace("]", "%5D")
    raw = _http_get_bytes(dods_base_url(url_path) + ".ascii?" + ce,
                          timeout=timeout)
    text = raw.decode("utf-8", "replace")
    vals: List[float] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("Dataset", "}", "---")):
            continue
        m = _ASCII_ROW_RE.match(line)
        if m:
            vals.extend(float(x) for x in m.group(1).split(",") if x.strip())
            continue
        # bare values (e.g. scalar ``time`` responses): "3.392712E8"
        try:
            vals.extend(float(x) for x in line.split(",") if x.strip())
        except ValueError:
            continue
    return vals


def fetch_valid_time(url_path: str, epoch: _dt.datetime,
                     timeout: float = 120.0) -> _dt.datetime:
    """Valid time of a one-timestep file from its numeric ``time``."""
    vals = fetch_ascii_values(url_path, TIME_VAR, ["0:1:0"], timeout=timeout)
    if not vals:
        raise UnavailableRangeError(
            f"no time value in {dods_base_url(url_path)}")
    return epoch + _dt.timedelta(seconds=vals[0])


# ---------------------------------------------------------------------------
# Grid vectors + regularity check
# ---------------------------------------------------------------------------


def fetch_grid_vectors(url_path: str, ny: int, nx: int,
                       timeout: float = 120.0) -> Tuple[np.ndarray, np.ndarray]:
    """1-D latitude/longitude vectors with a regularity guard.

    The ``regulargrid`` schema stores 2-D ``Latitude``/``Longitude``;
    the grid is regular (each constant along one axis), which is verified
    here — a non-regular grid raises ``ValueError``.
    """
    lat_col = fetch_ascii_values(url_path, LAT_VAR,
                                 [f"0:1:{ny - 1}", "0:1:0"], timeout=timeout)
    lon_row = fetch_ascii_values(url_path, LON_VAR,
                                 ["0:1:0", f"0:1:{nx - 1}"], timeout=timeout)
    lat_row0 = fetch_ascii_values(url_path, LAT_VAR,
                                  ["0:1:0", f"0:1:{nx - 1}"], timeout=timeout)
    lon_col0 = fetch_ascii_values(url_path, LON_VAR,
                                  [f"0:1:{ny - 1}", "0:1:0"], timeout=timeout)
    lats = np.asarray(lat_col, dtype=float)
    lons = np.asarray(lon_row, dtype=float)
    if lats.size != ny or lons.size != nx:
        raise UnavailableRangeError(
            f"grid vector size mismatch for {url_path}: "
            f"got ({lats.size}, {lons.size}), DDS says ({ny}, {nx})")
    if not np.allclose(lat_row0, lats[0]) or not np.allclose(lon_col0, lons[0]):
        raise ValueError(
            f"non-regular lat/lon grid in {url_path}: Latitude varies along "
            f"nx or Longitude varies along ny; regridding is not supported")
    return lats, lons


def index_range(vec: np.ndarray, lo: float, hi: float,
                what: str, url: str) -> Tuple[int, int]:
    """Inclusive index range of ``vec`` covering ``[lo, hi]``.

    Raises :class:`UnavailableRangeError` when the interval does not
    intersect the grid.
    """
    if vec[0] < vec[-1]:
        i0 = int(np.searchsorted(vec, lo, side="left"))
        i1 = int(np.searchsorted(vec, hi, side="right")) - 1
    else:
        i0 = int(np.searchsorted(-vec, -hi, side="left"))
        i1 = int(np.searchsorted(-vec, -lo, side="right")) - 1
    i0 = max(i0, 0)
    i1 = min(i1, len(vec) - 1)
    if i0 > i1:
        raise UnavailableRangeError(
            f"{what} range [{lo}, {hi}] is outside the model domain "
            f"[{min(vec[0], vec[-1])}, {max(vec[0], vec[-1])}] ({url})")
    return i0, i1


# ---------------------------------------------------------------------------
# DAP2 binary subset
# ---------------------------------------------------------------------------


def subset_constraint(y0: int, y1: int, x0: int, x1: int) -> str:
    """DAP2 constraint expression for surface u/v/temp/time (URL-escaped)."""
    box = f"[{y0}:1:{y1}][{x0}:1:{x1}]"
    parts = [
        f"{U_VAR}[0:1:0][0:1:0]{box}",
        f"{V_VAR}[0:1:0][0:1:0]{box}",
        f"{TEMP_VAR}[0:1:0][0:1:0]{box}",
        f"{TIME_VAR}[0:1:0]",
    ]
    return ",".join(parts).replace("[", "%5B").replace("]", "%5D")


_DAP2_FMT = {"Float32": "f", "Float64": "d", "Int32": "i"}


def parse_dods_response(raw: bytes) -> Dict[str, np.ndarray]:
    """Parse a ``.dods`` response into ``{variable: ndarray}``.

    The parser is driven by the response's own DDS header (the server
    reorders coordinate variables first, so request order is *not*
    trusted). Each array is framed as two big-endian uint32 element
    counts followed by XDR data (verified live 2026-10-01).
    """
    marker = raw.find(b"Data:")
    if marker < 0:
        raise UnavailableRangeError("no DAP2 data marker in response")
    header = raw[:marker].decode("utf-8", "replace")
    specs = []
    for m in _DDS_VAR_RE.finditer(header):
        dtype, name, dims = m.groups()
        nelem = 1
        shape = []
        for _d, n in _DDS_DIM_RE.findall(dims):
            nelem *= int(n)
            shape.append(int(n))
        specs.append((name, dtype, nelem, tuple(shape)))
    if not specs:
        raise UnavailableRangeError("empty DDS in .dods response")
    payload = raw[marker + len(b"Data:\n"):]
    out: Dict[str, np.ndarray] = {}
    for name, dtype, nelem, shape in specs:
        if len(payload) < 8:
            raise UnavailableRangeError(
                f"truncated .dods payload at variable {name!r}")
        n1, n2 = struct.unpack(">II", payload[:8])
        payload = payload[8:]
        if (n1, n2) != (nelem, nelem):
            raise UnavailableRangeError(
                f"DAP2 framing mismatch for {name!r}: header says "
                f"{nelem} elements, wire says ({n1}, {n2})")
        fmt = _DAP2_FMT[dtype]
        need = struct.calcsize(fmt) * nelem
        if len(payload) < need:
            raise UnavailableRangeError(
                f"truncated .dods payload for {name!r}: need {need} bytes, "
                f"have {len(payload)}")
        vals = struct.unpack(">" + fmt * nelem, payload[:need])
        payload = payload[need:]
        out[name] = np.asarray(vals, dtype=np.float64).reshape(shape)
    return out


def fetch_subset(url_path: str, y0: int, y1: int, x0: int, x1: int,
                 timeout: float = 120.0) -> Tuple[Dict[str, np.ndarray], int, str]:
    """Fetch one hourly surface subset; return (arrays, byte count, URL)."""
    url = dods_base_url(url_path) + ".dods?" + subset_constraint(y0, y1, x0, x1)
    raw = _http_get_bytes(url, timeout=timeout)
    return parse_dods_response(raw), len(raw), url


# ---------------------------------------------------------------------------
# Input coercion
# ---------------------------------------------------------------------------


def _coerce_hour_utc(value: DateLike, name: str) -> _dt.datetime:
    """Coerce to an hour-aligned UTC datetime (the model is hourly)."""
    if isinstance(value, _dt.date) and not isinstance(value, _dt.datetime):
        d = _dt.datetime.combine(value, _dt.time(0, 0),
                                 tzinfo=_dt.timezone.utc)
    elif isinstance(value, _dt.datetime):
        d = value if value.tzinfo else value.replace(tzinfo=_dt.timezone.utc)
        d = d.astimezone(_dt.timezone.utc)
    else:
        text = str(value).strip()
        if text.endswith(("Z", "z")):
            text = text[:-1] + "+00:00"
        try:
            d = _dt.datetime.fromisoformat(text)
        except ValueError as exc:
            raise ValueError(
                f"cannot parse {name}={value!r} as a datetime "
                f"(want YYYY-MM-DD or ISO)") from exc
        if d.tzinfo is None:
            d = d.replace(tzinfo=_dt.timezone.utc)
        d = d.astimezone(_dt.timezone.utc)
    if not (d.minute == 0 and d.second == 0 and d.microsecond == 0):
        raise ValueError(
            f"{name}={value!r} is not hour-aligned; OFS THREDDS files are hourly")
    return d


def _validate_bbox(bbox: Sequence[float]) -> Tuple[float, float, float, float]:
    try:
        lon0, lat0, lon1, lat1 = (float(x) for x in bbox)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "bbox must be (min_lon, min_lat, max_lon, max_lat)") from exc
    if not (lon0 < lon1 and lat0 < lat1):
        raise ValueError(f"bbox needs min < max, got {bbox!r}")
    if not (-180.0 <= lon0 <= 180.0 and -180.0 <= lon1 <= 180.0
            and -90.0 <= lat0 <= 90.0 and -90.0 <= lat1 <= 90.0):
        raise ValueError(f"bbox out of geographic range: {bbox!r}")
    return lon0, lat0, lon1, lat1


# ---------------------------------------------------------------------------
# End-to-end fetch
# ---------------------------------------------------------------------------


def _day_file_index(model: str, days: Sequence[_dt.date],
                    timeout: float) -> Dict[_dt.date, List[Dict[str, Any]]]:
    """Map each day to parsed regulargrid entries (cached per fetch call)."""
    index: Dict[_dt.date, List[Dict[str, Any]]] = {}
    for day in days:
        entries = []
        for e in list_regulargrid_files(model, day, timeout=timeout):
            try:
                parsed = parse_ofs_filename(e["name"])
            except ValueError:
                continue
            if parsed["code"] != model.lower():
                continue
            entries.append({**e, **parsed})
        index[day] = entries
    return index


def _calibrate_nowcast_span(model: str,
                            index: Dict[_dt.date, List[Dict[str, Any]]],
                            timeout: float) -> Tuple[Optional[int], _dt.datetime]:
    """Probe one nowcast file to learn the per-model nowcast span + epoch.

    Returns ``(span_hours_or_None, time_epoch)``. ``span`` satisfies
    ``valid = cycle - span + HHH`` for ``n{HHH}`` files (verified 6 h for
    SSCOFS/CBOFS/NGOFS2/GOMOFS/DBOFS, 24 h for WCOFS, live 2026-10-01).
    """
    for day in sorted(index):
        for e in index[day]:
            if e["kind"] != "n":
                continue
            das = fetch_das(e["url_path"], timeout=timeout)
            epoch = parse_time_units(das)
            valid = fetch_valid_time(e["url_path"], epoch, timeout=timeout)
            span = int(round((e["cycle"] + _dt.timedelta(hours=e["hour"])
                              - valid).total_seconds() / 3600))
            return span, epoch
    # No nowcast files at all: fall back to a forecast file for the epoch.
    for day in sorted(index):
        for e in index[day]:
            das = fetch_das(e["url_path"], timeout=timeout)
            return None, parse_time_units(das)
    raise UnavailableRangeError(
        f"no regulargrid files at all for {model} in the requested window")


def _pick_file(target: _dt.datetime,
               index: Dict[_dt.date, List[Dict[str, Any]]],
               span: Optional[int], prefer: str) -> Dict[str, Any]:
    """Choose the file whose valid time is exactly ``target``."""
    cands = []
    for day in (target.date(), target.date() - _dt.timedelta(days=1)):
        for e in index.get(day, []):
            try:
                if nominal_valid_time(e, span) == target:
                    cands.append(e)
            except ValueError:
                continue
    if not cands:
        avail = sorted({nominal_valid_time(e, span)
                        for day in index for e in index[day]
                        if e["kind"] == "f" or span is not None})
        raise UnavailableRangeError(
            f"no {prefer} file for valid time {target.isoformat()} "
            f"(nearest available hours: "
            f"{', '.join(t.isoformat() for t in avail[:6])}...)")
    # Prefer the requested kind (nowcast = analysis), then the *latest*
    # model cycle: the most recent run assimilates the most data and gives
    # the shortest forecast lead for the target hour.
    cands.sort(key=lambda e: (0 if e["kind"] == prefer[0] else 1,
                              -e["cycle"].timestamp(), e["hour"]))
    return cands[0]


def fetch_ofs_thredds(ofs_code: str,
                      bbox: Sequence[float],
                      start: DateLike,
                      end: DateLike,
                      cadence_hours: int = 6,
                      prefer: str = "nowcast",
                      timeout: float = 120.0) -> CurrentField:
    """Fetch OFS surface currents + temperature via THREDDS OPeNDAP.

    ``ofs_code`` is e.g. ``"SSCOFS"`` (see :data:`OFS_THREDDS_MODELS`);
    ``bbox`` is ``(min_lon, min_lat, max_lon, max_lat)`` and is subset
    **server-side** — only the requested window crosses the wire.
    ``start``/``end`` are hour-aligned UTC datetimes; ``end`` is
    inclusive; one timestep is fetched per ``cadence_hours``.

    Returns a :class:`CurrentField` with ``u``/``v`` in m/s,
    ``temperature`` in degC, ``source="ofs-thredds/<CODE>"``, and
    per-timestep provenance recording the exact OPeNDAP URLs.
    Missing dates/hours raise :class:`UnavailableRangeError` — the
    series is never silently padded.
    """
    code = str(ofs_code).strip().upper()
    if code not in OFS_THREDDS_MODELS:
        raise ValueError(
            f"unknown THREDDS OFS code {ofs_code!r}; known: "
            f"{', '.join(OFS_THREDDS_MODELS)}")
    lon0, lat0, lon1, lat1 = _validate_bbox(bbox)
    t0 = _coerce_hour_utc(start, "start")
    t1 = _coerce_hour_utc(end, "end")
    if t1 < t0:
        raise ValueError(f"end {t1.isoformat()} is before start {t0.isoformat()}")
    cadence = int(cadence_hours)
    if cadence < 1:
        raise ValueError(f"cadence_hours must be >= 1, got {cadence_hours!r}")
    if prefer not in ("nowcast", "forecast"):
        raise ValueError(f"prefer must be 'nowcast' or 'forecast', got {prefer!r}")

    targets = []
    t = t0
    while t <= t1:
        targets.append(t)
        t += _dt.timedelta(hours=cadence)

    days = sorted({d for tgt in targets
                   for d in (tgt.date(), tgt.date() - _dt.timedelta(days=1))})
    index = _day_file_index(code, days, timeout)
    span, epoch = _calibrate_nowcast_span(code, index, timeout)

    # Grid vectors once, from the first chosen file.
    first = _pick_file(targets[0], index, span, prefer)
    dds = fetch_dds(first["url_path"], timeout=timeout)
    dimmap = {name: dict(dd) for name, _dtype, dd in dds}
    for var in (U_VAR, V_VAR, TEMP_VAR, LAT_VAR, LON_VAR, TIME_VAR):
        if var not in dimmap:
            raise UnavailableRangeError(
                f"variable {var!r} missing in "
                f"{dods_base_url(first['url_path'])}; variables: "
                f"{sorted(dimmap)}")
    ny = dimmap[LAT_VAR].get("ny")
    nx = dimmap[LON_VAR].get("nx")
    if not ny or not nx:
        raise UnavailableRangeError(
            f"cannot find ny/nx dims in {dods_base_url(first['url_path'])}: "
            f"{dimmap[LAT_VAR]}, {dimmap[LON_VAR]}")
    lats, lons = fetch_grid_vectors(first["url_path"], ny, nx, timeout=timeout)
    das = fetch_das(first["url_path"], timeout=timeout)
    fill = parse_fill_value(das)
    y0, y1 = index_range(lats, lat0, lat1, "latitude",
                         dods_base_url(first["url_path"]))
    x0, x1 = index_range(lons, lon0, lon1, "longitude",
                         dods_base_url(first["url_path"]))
    clipped = (lons[x0] != lon0 or lons[x1] != lon1
               or lats[y0] != lat0 or lats[y1] != lat1)

    # CurrentField wants increasing 1-D coordinates.
    lat_rev = lats[0] > lats[-1]
    lon_rev = lons[0] > lons[-1]
    glats = lats[y0:y1 + 1][::-1] if lat_rev else lats[y0:y1 + 1]
    glons = lons[x0:x1 + 1][::-1] if lon_rev else lons[x0:x1 + 1]

    steps_u, steps_v, steps_t, times = [], [], [], []
    step_prov: List[Dict[str, Any]] = []
    n_bytes = 0
    model_run = ""
    forecast_hours: List[int] = []
    for tgt in targets:
        e = _pick_file(tgt, index, span, prefer)
        arrays, nbytes, url = fetch_subset(e["url_path"], y0, y1, x0, x1,
                                           timeout=timeout)
        n_bytes += nbytes
        for var in (U_VAR, V_VAR, TEMP_VAR, TIME_VAR):
            if var not in arrays:
                raise UnavailableRangeError(
                    f"variable {var!r} missing in .dods response for {url}")
        valid = epoch + _dt.timedelta(seconds=float(arrays[TIME_VAR].ravel()[0]))
        if abs((valid - tgt).total_seconds()) > 1:
            raise ValueError(
                f"served step valid time {valid.isoformat()} does not match "
                f"requested {tgt.isoformat()} ({url}); refusing to mislabel")

        def clean(a: np.ndarray) -> np.ndarray:
            # Constrained shape is (time=1, Depth=1, ny, nx); squeeze the
            # singletons down to (ny, nx), mask land, fix axis order.
            a = np.squeeze(np.asarray(a, dtype=np.float64))
            if a.ndim != 2:
                raise UnavailableRangeError(
                    f"unexpected subset shape {a.shape} for {url}")
            a = a.copy()
            a[a == fill] = np.nan
            if lat_rev:
                a = a[::-1, :]
            if lon_rev:
                a = a[:, ::-1]
            return a[None, :, :]

        steps_u.append(clean(arrays[U_VAR]))
        steps_v.append(clean(arrays[V_VAR]))
        steps_t.append(clean(arrays[TEMP_VAR]))
        times.append(valid.isoformat())
        fh = e["hour"] if e["kind"] == "f" else e["hour"] - (span or 0)
        forecast_hours.append(fh)
        if not model_run:
            model_run = e["cycle"].isoformat()
        step_prov.append({
            "url": url, "file": e["name"], "valid_time": valid.isoformat(),
            "kind": e["kind"], "hour": e["hour"], "cycle": e["cycle"].isoformat(),
            "bytes": nbytes,
        })

    field = CurrentField(
        u=np.concatenate(steps_u, axis=0),
        v=np.concatenate(steps_v, axis=0),
        temperature=np.concatenate(steps_t, axis=0),
        times=times,
        lats=np.asarray(glats, dtype=float),
        lons=np.asarray(glons, dtype=float),
        crs="EPSG:4326",
        source=f"ofs-thredds/{code}",
        model_run=model_run,
        forecast_hours=forecast_hours,
        provenance={
            "service": "NOAA CO-OPS THREDDS OPeNDAP",
            "thredds_base": THREDDS_BASE,
            "model": code,
            "requested_bbox": [lon0, lat0, lon1, lat1],
            "grid_bbox": [float(glons[0]), float(glats[0]),
                           float(glons[-1]), float(glats[-1])],
            "bbox_clipped_to_grid": bool(clipped),
            "variables": {"u": U_VAR, "v": V_VAR,
                          "temperature": TEMP_VAR, "time": TIME_VAR},
            "u_units": "m/s", "v_units": "m/s",
            "temperature_units": "degC",
            "time_units": f"seconds since {epoch.isoformat()}",
            "surface_depth_m": 0.0,
            "fill_value": fill,
            "nowcast_span_hours": span,
            "prefer": prefer,
            "cadence_hours": cadence,
            "n_timesteps": len(times),
            "n_bytes": n_bytes,
            "water_fraction": float(np.mean(~np.isnan(
                np.concatenate(steps_u, axis=0)))),
            "steps": step_prov,
        },
    )
    return field

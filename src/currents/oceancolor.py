"""Ocean color acquisition: chlorophyll-a from NOAA CoastWatch ERDDAP + NASA OBPG.

Complementary global ocean-color sources:

* :func:`fetch_oceancolor` — NOAA CoastWatch ERDDAP (keyless HTTPS),
  MODIS Aqua R2022 Science Quality chlorophyll-a, Level 3, global 4 km
  grid, variable ``chlor_a`` in mg m^-3, 2002-present (default
  ``sensor="modis-aqua"``); monthly (default), weekly, and daily
  composites. VIIRS SNPP Science Quality (2012-present) and the
  multi-sensor ESA OC-CCI v6.0 merge (1997-present) are selectable via
  ``sensor=``. The longitude axis is already -180..180, so no wrap
  conversion is needed.
* ``source="obpg"`` — NASA Ocean Biology Processing Group (OBPG)
  MODIS Aqua Level-3 mapped chlorophyll-a via Earthdata-authenticated
  HTTPS (``oceandata.sci.gsfc.nasa.gov`` direct data access; needs a
  free Earthdata Login — without it the :class:`CredentialsMissing`
  error from :mod:`sst_global` explains the setup, never an
  interactive prompt). The credentials-required fallback when
  CoastWatch is unreachable.
* ``source="cmems"`` — Copernicus Marine
  ``OCEANCOLOUR_GLO_BGC_L4_MY_009_104`` (global ocean colour, L4
  multi-year, monthly, 4 km, 1997-present, variable ``CHL`` in mg/m^3)
  via the ``copernicusmarine`` toolbox. Optional extra source (needs
  free CMEMS credentials); without them the actionable error from
  :func:`cmems.require_toolbox` explains the setup (never an
  interactive prompt).
* ``source="auto"`` — try CoastWatch, then OBPG (when Earthdata
  credentials exist), then CMEMS (when CMEMS credentials exist).

Returns :class:`OceanColorField`, which follows the :class:`SstField`
conventions in ``sst_global.py`` (and duck-types into what
``survey-viz``'s ``render_viz`` consumes: ``.times``/``.lats``/``.lons``
plus a 3-D ``.values`` grid).

Chlorophyll is cloud-sensitive by nature: cloud-covered pixels are NaN
in the source products and are KEPT as NaN here — never interpolated
or filled. Monthly composites are the default because they give the
most complete global coverage; daily/weekly frames routinely have
large cloud gaps, which are reported honestly per frame
(:meth:`OceanColorField.gap_fraction`, ``provenance["gap_fractions"]``).

Downloads use stdlib ``urllib`` only. ``netCDF4`` is a lazy import, via
:func:`sst_global._open_nc_bytes`.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import os
import re
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .cmems import REGISTER_URL as CMEMS_REGISTER_URL
from .glsea import (_coerce_date, _coerce_datetime_utc, _download_bytes,
                    _erddap_iso)
from .sst_global import _open_nc_bytes, validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]

# ---------------------------------------------------------------------------
# NOAA CoastWatch ERDDAP constants
# ---------------------------------------------------------------------------

#: Canonical CoastWatch ERDDAP host. NOTE: the legacy
#: ``coastwatch.pfeg.noaa.gov`` host (used by survey-currents <= 0.13.0
#: for OISST) was decommissioned — it returns empty replies. The
#: current canonical host is ``coastwatch.noaa.gov``.
COASTWATCH_ERDDAP_BASE = "https://coastwatch.noaa.gov/erddap"

#: Verified dataset IDs on the CoastWatch ERDDAP (confirmed from the
#: CoastWatch tutorials repo, the NOAA sanctuary caption pages, and the
#: CoastWatch WMS dataset list — see docs/DATA_SOURCES.md). Keyed by
#: (sensor, cadence). Only IDs with an independent confirmation live
#: here; anything else raises instead of guessing.
#:
#: * ``viirs-snpp`` — "Chlorophyll, NOAA S-NPP VIIRS, Science Quality,
#:   Global 4km, Level 3, 2012-present" (nesdisVHNSQchla{Daily,Weekly,Monthly}).
#: * ``modis-aqua`` — "MODIS Aqua, Chlorophyll-a, V.2022, Monthly,
#:   Global, 4km, 2002-present" (erdMH1chlamday_R2022SQ).
#: * ``multi`` — "ESA CCI Ocean Colour Dataset, v6.0, Monthly, Global,
#:   4km, 1997-Present" (pmlEsaCCI60OceanColorMonthly; multi-sensor
#:   merge: SeaWiFS/MODIS/MERIS/VIIRS) — verified 2026-09-27 in the
#:   official coastwatch-training/coastwatch-tutorials repo
#:   (Tutorial2-timeseries-compare-sensors) and an independent research
#:   repo (dtcrompton/albatross-foraging-ecology), both citing the ID
#:   against NOAA CoastWatch ERDDAP.
COASTWATCH_DATASETS: Dict[Tuple[str, str], Dict[str, Any]] = {
    ("viirs-snpp", "monthly"): {
        "id": "nesdisVHNSQchlaMonthly",
        "title": "Chlorophyll, NOAA S-NPP VIIRS, Science Quality, "
                 "Global 4km, Level 3, Monthly",
        "start": _dt.date(2012, 1, 1),
        "sensor": "VIIRS (Suomi NPP)",
        "processing": "NOAA STAR science quality, Level 3",
    },
    ("viirs-snpp", "weekly"): {
        "id": "nesdisVHNSQchlaWeekly",
        "title": "Chlorophyll, NOAA S-NPP VIIRS, Science Quality, "
                 "Global 4km, Level 3, Weekly",
        "start": _dt.date(2012, 1, 1),
        "sensor": "VIIRS (Suomi NPP)",
        "processing": "NOAA STAR science quality, Level 3",
    },
    ("viirs-snpp", "daily"): {
        "id": "nesdisVHNSQchlaDaily",
        "title": "Chlorophyll, NOAA S-NPP VIIRS, Science Quality, "
                 "Global 4km, Level 3, Daily",
        "start": _dt.date(2012, 1, 1),
        "sensor": "VIIRS (Suomi NPP)",
        "processing": "NOAA STAR science quality, Level 3",
    },
    ("modis-aqua", "monthly"): {
        "id": "erdMH1chlamday_R2022SQ",
        "title": "Chlorophyll-a, Aqua MODIS, V.2022, Monthly, "
                 "Global, 4km, Science Quality",
        "start": _dt.date(2002, 7, 1),
        "sensor": "MODIS (Aqua)",
        "processing": "NASA OBPG R2022 reprocessing, Level 3",
    },
    ("multi", "monthly"): {
        "id": "pmlEsaCCI60OceanColorMonthly",
        "title": "ESA CCI Ocean Colour v6.0, Monthly, Global, 4km "
                 "(multi-sensor merge)",
        "start": _dt.date(1997, 9, 1),
        "sensor": "multi-sensor (SeaWiFS/MODIS/MERIS/VIIRS)",
        "processing": "ESA Climate Change Initiative v6.0, Level 3",
    },
}

#: Grid geometry per the nesdisVHNSQchlaMonthly metadata (verified via
#: the CoastWatch tutorials repo). All three VIIRS products share the
#: global 4 km grid; MODIS/CCI products are the same nominal 4 km.
OCEANCOLOR_LON_MIN, OCEANCOLOR_LON_MAX = -179.9812, 179.9813
OCEANCOLOR_LAT_MIN, OCEANCOLOR_LAT_MAX = -89.75626, 89.75625
OCEANCOLOR_RES_DEG = 0.0417  # ~4 km nominal
OCEANCOLOR_VAR = "chlor_a"
OCEANCOLOR_FILL = -999.0
OCEANCOLOR_UNITS = "mg m^-3"

# ---------------------------------------------------------------------------
# CMEMS constants (authenticated fallback)
# ---------------------------------------------------------------------------

#: CMEMS global ocean-colour product (per the Nov-2025 OC QUID
#: CMEMS-OC-QUID-009-101to104-111-113-116-118): L4 multi-year,
#: monthly, 4 km, 1997-present. The exact dataset string should be
#: confirmed against the CMEMS catalogue when credentials are
#: available — the fetcher accepts ``cmems_dataset_id`` to override.
CMEMS_OCEANCOLOR_PRODUCT = "OCEANCOLOUR_GLO_BGC_L4_MY_009_104"
CMEMS_OCEANCOLOR_VARIABLE = "CHL"

# ---------------------------------------------------------------------------
# NASA OBPG constants (credentials-required fallback)
# ---------------------------------------------------------------------------

#: NASA Ocean Biology Processing Group direct data access for MODIS
#: Aqua Level-3 mapped products (Earthdata Login required).
OBPG_DIRECTACCESS_BASE = ("https://oceandata.sci.gsfc.nasa.gov"
                          "/directdataaccess/Level-3%20Mapped/Aqua-MODIS")

#: NASA CMR collection for MODIS Aqua L3 mapped chlorophyll
#: (short name MODISA_L3m_CHL, version 2022.0, provider OB_CLOUD) —
#: verified live via the keyless CMR collections API 2026-09-27
#: (time start 2002-07-04). Used for granule discovery metadata;
#: per-window files come from the OBPG direct-access filenames below.
OBPG_CMR_COLLECTION_ID = "C3380709133-OB_CLOUD"
OBPG_CMR_SHORT_NAME = "MODISA_L3m_CHL"

#: OBPG Level-3 mapped file stem (SeaDAS format), per the OBPG file
#: naming convention (verified against real CMR granule titles
#: 2026-09-27, e.g. ``AQUA_MODIS.20020704_20250228.L3m.CU.CHL.chlor_a.4km.nc``):
#: ``AQUA_MODIS.<start>_<end>.L3m.<DAY|8D|MO>.CHL.chlor_a.4km.nc``.
#: Period files are addressed by their documented date ranges —
#: monthly = first-to-last day of month, 8-day = OBPG 8-day periods
#: (day-of-year 1-8, 9-16, ...), daily = the single date.
OBPG_FILE_SENSOR = "AQUA_MODIS"
OBPG_FILE_VARIABLE = "chlor_a"
OBPG_FILE_RESOLUTION = "4km"

#: MODIS Aqua ocean-color record start (matches the CoastWatch
#: R2022 monthly dataset and the CMR collection time start).
OBPG_RECORD_START = _dt.date(2002, 7, 4)

# ---------------------------------------------------------------------------
# cache discipline (mirrors streamgages.py)
# ---------------------------------------------------------------------------

#: Cache entries older than this are re-downloaded.
OCEANCOLOR_MAX_CACHE_AGE_DAYS = 7


def oceancolor_cache_dir() -> str:
    """Cache root for ocean-color downloads.

    ``$SURVEY_CURRENTS_CACHE/oceancolor`` when set, else
    ``~/.cache/survey-currents/oceancolor``.
    """
    root = os.environ.get("SURVEY_CURRENTS_CACHE") or os.path.join(
        os.path.expanduser("~"), ".cache", "survey-currents")
    return os.path.join(root, "oceancolor")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _cache_key(parts: Sequence[str]) -> str:
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:32]


def _cached_get_bytes(url: str, cache_dir: str,
                      max_age_days: int = OCEANCOLOR_MAX_CACHE_AGE_DAYS,
                      refresh: bool = False) -> Tuple[bytes, Dict[str, Any]]:
    """GET ``url`` with an atomic SHA-256-verified cache.

    Returns ``(payload, provenance)`` with ``sha256``, ``cache_hit``,
    ``downloaded``, and ``cache_path`` keys. A corrupt cache entry is
    re-downloaded rather than trusted.
    """
    os.makedirs(cache_dir, exist_ok=True)
    key = _cache_key([url])
    path = os.path.join(cache_dir, f"{key}.nc")
    sidecar = path + ".sha256"
    prov: Dict[str, Any] = {"url": url, "cache_path": path}

    def _valid_entry() -> Optional[bytes]:
        try:
            with open(sidecar, encoding="utf-8") as fh:
                expected = fh.read().strip().split()[0]
            if (not refresh and os.path.getmtime(path) >
                    _dt.datetime.now().timestamp() - max_age_days * 86400):
                with open(path, "rb") as fh:
                    data = fh.read()
                if _sha256_bytes(data) == expected:
                    prov.update(cache_hit=True, downloaded=False,
                                sha256=expected)
                    return data
        except (OSError, ValueError, IndexError):
            pass
        return None

    cached = _valid_entry()
    if cached is not None:
        return cached, prov
    data = _download_bytes(url)
    digest = _sha256_bytes(data)
    fd, tmp = tempfile.mkstemp(dir=cache_dir, suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        with open(tmp + ".sha256", "w", encoding="utf-8") as fh:
            fh.write(digest + "  " + os.path.basename(path) + "\n")
        os.replace(tmp + ".sha256", sidecar)
        os.replace(tmp, path)
    finally:
        for leftover in (tmp, tmp + ".sha256"):
            try:
                os.unlink(leftover)
            except OSError:
                pass
    prov.update(cache_hit=False, downloaded=True, sha256=digest)
    return data, prov


# ---------------------------------------------------------------------------
# OceanColorField — canonical ocean-color field model
# ---------------------------------------------------------------------------


@dataclass
class OceanColorField:
    """Time-indexed chlorophyll-a grids (global ocean-color sources).

    ``values`` is ``(nt, ny, nx)`` chlorophyll-a concentration in
    mg m^-3 as a numpy masked array (land / cloud / missing cells
    masked — never interpolated). ``lats``/``lons`` are 1-D coordinate
    vectors, increasing, in the -180..180 convention. ``chl`` is an
    alias for ``values``. Follows the :class:`SstField` conventions in
    ``sst_global.py``.
    """

    values: np.ma.MaskedArray             # (nt, ny, nx) chlorophyll-a, mg/m^3
    times: List[str]                     # ISO-8601 timestamps, one per step
    lats: np.ndarray                     # (ny,) degrees north, increasing
    lons: np.ndarray                     # (nx,) degrees east, -180..180, increasing
    crs: str = "EPSG:4326"
    source: str = "noaa-coastwatch/erdMH1chlamday_R2022SQ"
    provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        self.values = np.ma.asarray(self.values, dtype=float)
        # Coerce coordinates BEFORE the shape check: callers may pass
        # plain lists (from_dict, duck-typed peers), which have no
        # .shape.
        self.lats = np.asarray(self.lats, dtype=float)
        self.lons = np.asarray(self.lons, dtype=float)
        if self.values.shape != (nt, self.lats.shape[0], self.lons.shape[0]):
            raise ValueError(
                f"values shape {self.values.shape} does not match "
                f"(nt, ny, nx)=({nt}, {self.lats.shape[0]}, {self.lons.shape[0]})"
            )
        self.bounds = self._bounds()

    @property
    def chl(self) -> np.ma.MaskedArray:
        """Alias for :attr:`values` (chlorophyll-a, mg m^-3)."""
        return self.values

    def _bounds(self) -> Tuple[float, float, float, float]:
        return (float(self.lons[0]), float(self.lats[0]),
                float(self.lons[-1]), float(self.lats[-1]))

    # -- derived quantities --------------------------------------------------

    def gap_fraction(self, index: int) -> float:
        """Fraction of NaN (cloud/land/missing) cells at timestep ``index``.

        Chlorophyll is cloud-sensitive: this is an honest data-quality
        record, not something to fill — see docs/OCEANCOLOR.md.
        """
        arr = np.ma.masked_invalid(self.values[index])
        total = arr.size
        if total == 0:
            return 1.0
        return float(np.ma.getmaskarray(arr).sum() / total)

    def spatial_median(self, index: int) -> float:
        """NaN-aware spatial median chlorophyll-a (mg m^-3).

        The median (not the mean) is the standard central tendency for
        lognormally distributed chlorophyll — robust to the long tail
        of bloom pixels.
        """
        arr = np.ma.masked_invalid(self.values[index]).compressed()
        if arr.size == 0:
            return float("nan")
        return float(np.median(arr))

    # -- selection -----------------------------------------------------------

    def select_time(self, index: int) -> "OceanColorField":
        """Return the single-timestep field at ``index``."""
        return OceanColorField(
            values=self.values[index:index + 1], times=[self.times[index]],
            lats=self.lats, lons=self.lons, crs=self.crs,
            source=self.source, provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "OceanColorField":
        """Spatial subset to ``(min_lon, min_lat, max_lon, max_lat)``."""
        minx, miny, maxx, maxy = validate_sst_bbox(bbox)
        ix = np.where((self.lons >= minx) & (self.lons <= maxx))[0]
        iy = np.where((self.lats >= miny) & (self.lats <= maxy))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(
                f"bbox {list(bbox)} does not overlap field bounds {self.bounds}")
        return OceanColorField(
            values=self.values[:, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1],
            times=list(self.times),
            lats=self.lats[iy[0]:iy[-1] + 1], lons=self.lons[ix[0]:ix[-1] + 1],
            crs=self.crs, source=self.source,
            provenance=dict(self.provenance),
        )

    # -- serialization -------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dict (masked cells become NaN).

        The grid travels under the ``"values"`` key — the form
        ``survey-viz``'s ``render_viz`` documents for plain dicts —
        with ``"chl"`` kept as a human-friendly alias.
        """
        grid = np.ma.filled(self.values, np.nan).tolist()
        return {
            "values": grid,
            "chl": grid,
            "times": list(self.times),
            "lats": self.lats.tolist(),
            "lons": self.lons.tolist(),
            "crs": self.crs,
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "OceanColorField":
        """Rebuild from :meth:`to_dict` (NaN -> masked). Raises on missing keys."""
        required = ("times", "lats", "lons")
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(f"OceanColorField dict missing keys: {missing}")
        grid = data.get("values", data.get("chl"))
        if grid is None:
            raise ValueError("OceanColorField dict missing 'values' (or 'chl')")
        return cls(
            values=np.ma.masked_invalid(np.asarray(grid, dtype=float)),
            times=list(data["times"]),
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            crs=data.get("crs", "EPSG:4326"),
            source=data.get("source", "noaa-coastwatch/nesdisVHNSQchlaMonthly"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "OceanColorField":
        import json
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    # -- synthetic fixture (tests, demos, offline use) -----------------------

    @classmethod
    def synthetic(
        cls,
        nt: int = 4, ny: int = 6, nx: int = 8,
        lats: Sequence[float] = (-30.0, 30.0),
        lons: Sequence[float] = (-60.0, 60.0),
        start: str = "2020-01-01T00:00:00+00:00",
        step_days: int = 30,
        seed: int = 12,
        source: str = "synthetic",
        gap_fraction: float = 0.25,
    ) -> "OceanColorField":
        """Deterministic synthetic chlorophyll-a field (lognormal).

        Stdlib+numpy only. Chlorophyll is lognormally distributed in
        the real ocean, so the fixture draws from a lognormal
        (median ~0.3 mg/m^3 offshore, higher toward the coasts) with a
        deterministic cloud mask covering ``gap_fraction`` of cells —
        the honest NaN treatment the real fetch applies.
        """
        rng = np.random.default_rng(seed)
        lat = np.linspace(lats[0], lats[1], ny)
        lon = np.linspace(lons[0], lons[1], nx)
        yy, xx = np.meshgrid(np.linspace(-1, 1, ny),
                             np.linspace(-1, 1, nx), indexing="ij")
        # Higher chlorophyll near coasts/equator, oligotrophic gyres
        # offshore: log10 median from ~-0.5 (0.3 mg/m^3) offshore to
        # ~0.7 (5 mg/m^3) near the domain edges.
        log_med = -0.5 + 1.2 * np.clip(np.abs(xx) + 0.3 * (1 - yy ** 2),
                                       0.0, 1.0)
        values = np.ma.empty((nt, ny, nx))
        t0 = _dt.datetime.fromisoformat(start)
        times = []
        for k in range(nt):
            field_k = 10.0 ** (log_med + rng.normal(0.0, 0.25, size=(ny, nx)))
            mask = rng.random((ny, nx)) < gap_fraction
            values[k] = np.ma.array(field_k, mask=mask)
            times.append((t0 + _dt.timedelta(days=k * step_days)).isoformat())
        return cls(values=values, times=times, lats=lat, lons=lon,
                   source=source,
                   provenance={"synthetic": True, "seed": seed,
                               "gap_fraction": gap_fraction})


# ---------------------------------------------------------------------------
# CoastWatch ERDDAP fetch
# ---------------------------------------------------------------------------


def _dataset_entry(sensor: str, cadence: str) -> Dict[str, Any]:
    """Look up the CoastWatch dataset for ``(sensor, cadence)``.

    Raises :class:`ValueError` listing the valid combinations instead
    of guessing a dataset ID.
    """
    key = (str(sensor).strip().lower(), str(cadence).strip().lower())
    try:
        return COASTWATCH_DATASETS[key]
    except KeyError as exc:
        known = ", ".join(f"{s}/{c}" for s, c in sorted(COASTWATCH_DATASETS))
        raise ValueError(
            f"unknown ocean-color (sensor, cadence) {key!r}; "
            f"known combinations: {known}") from exc


def _validate_oceancolor_dates(d0: _dt.datetime, d1: _dt.datetime,
                               earliest: _dt.date, label: str) -> None:
    if d0 > d1:
        raise ValueError(f"start {d0.date()} is after end {d1.date()}")
    if d0.date() < earliest:
        raise ValueError(
            f"{label} starts {earliest.isoformat()}; "
            f"start {d0.date().isoformat()} is before the record")
    if d1.date() > _dt.date.today() + _dt.timedelta(days=2):
        raise ValueError(f"end {d1.date().isoformat()} is in the future")


def oceancolor_urls(bbox: Sequence[float], start: DateLike, end: DateLike,
                    sensor: str = "modis-aqua", cadence: str = "monthly",
                    stride: int = 1,
                    erddap_base: str = COASTWATCH_ERDDAP_BASE) -> List[str]:
    """Build the CoastWatch griddap NetCDF URL(s) (exact bytes land in provenance).

    The longitude axis is -180..180 on the wire (per the dataset
    metadata), so the requested bbox is used as-is. Each query indexes
    the singleton ``altitude`` axis explicitly with ``[(0.0)]`` — the
    same requirement the OISST adapter documents for ``zlev``.
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    entry = _dataset_entry(sensor, cadence)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_oceancolor_dates(d0, d1, entry["start"], entry["id"])
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}")
    t0, t1 = _erddap_iso(d0), _erddap_iso(d1)
    query = (
        f"{OCEANCOLOR_VAR}[({t0}):{stride}:({t1})][(0.0)]"
        f"[({miny}):({maxy})][({minx}):({maxx})]"
    )
    return [f"{erddap_base.rstrip('/')}/griddap/{entry['id']}.nc?{query}"]


def _parse_coastwatch_bytes(payload: bytes) -> Tuple[np.ndarray, np.ndarray,
                                                     List[str],
                                                     np.ma.MaskedArray]:
    """Parse one CoastWatch chlorophyll griddap NetCDF payload.

    Returns (lats, lons, times, chl). Lons are normalized to -180..180
    and sorted increasing; the chl array has shape (nt, ny, nx) in
    mg m^-3 with fill/missing values masked (never filled).
    """
    nc_mod = None
    try:
        from .sst_global import _require_netcdf4
        nc_mod = _require_netcdf4()
    except ImportError:
        pass
    with _open_nc_bytes(payload) as ds:
        if OCEANCOLOR_VAR not in ds.variables:
            raise ValueError(
                f"ocean-color payload missing variable {OCEANCOLOR_VAR!r}")
        cvar = ds.variables[OCEANCOLOR_VAR]
        dims = cvar.dimensions
        # (time, altitude, latitude, longitude); tolerate a missing
        # altitude axis on future regriddings.
        data = np.ma.asarray(cvar[:], dtype=float)
        if data.ndim == 4:
            data = data[:, 0, :, :]
        elif data.ndim != 3:
            raise ValueError(
                f"unexpected {OCEANCOLOR_VAR} dimensions {dims}; "
                "expected (time, altitude, latitude, longitude)")
        lats = np.asarray(ds.variables["latitude"][:], dtype=float)
        lons_raw = np.asarray(ds.variables["longitude"][:], dtype=float)
        time_var = ds.variables["time"]
        try:
            dts = nc_mod.num2date(time_var[:], time_var.units) \
                if nc_mod is not None else time_var[:]
        except Exception as exc:
            raise ValueError(
                f"could not decode ocean-color time axis: {exc}") from exc
    times = [_coerce_datetime_utc(d).isoformat() for d in dts]
    lons = ((lons_raw + 180.0) % 360.0) - 180.0
    order = np.argsort(lons)
    lons = lons[order]
    chl = np.ma.masked_values(data[:, :, order], OCEANCOLOR_FILL,
                              rtol=1e-5, atol=1e-3)
    chl = np.ma.masked_invalid(chl)
    return lats, lons, times, chl


def fetch_oceancolor_coastwatch(
        bbox: Sequence[float], start: DateLike, end: DateLike,
        sensor: str = "modis-aqua", cadence: str = "monthly",
        stride: int = 1,
        erddap_base: str = COASTWATCH_ERDDAP_BASE,
        work_dir: Optional[str] = None,
        refresh: bool = False) -> OceanColorField:
    """Fetch chlorophyll-a from the keyless CoastWatch ERDDAP.

    Args:
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180
            degrees (ascending). The server subsets by bbox, so the
            download stays proportional to the requested region —
            never the full 4320x8640 global grid.
        start, end: dates/datetimes/ISO strings.
        sensor: ``"modis-aqua"`` (default), ``"viirs-snpp"``,
            ``"multi"`` (ESA OC-CCI merge).
        cadence: ``"monthly"`` (default — the most cloud-complete),
            ``"weekly"``, ``"daily"``.
        stride: time-axis index stride (>= 1).
        erddap_base: override the ERDDAP host (keyless mirrors).
        work_dir: cache directory (default: the shared cache).
        refresh: ignore the cache and re-download.

    Returns:
        :class:`OceanColorField` in mg m^-3 with provenance (dataset,
        exact URLs, SHA-256, byte counts, retrieval time, sensor,
        cadence, units, per-frame gap fractions, cache state).

    Raises:
        ValueError: invalid bbox / dates / sensor / cadence / stride.
        RuntimeError: network or NetCDF failures (with context) — e.g.
            the CoastWatch ERDDAP outage observed 2026-09-27.
    """
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    entry = _dataset_entry(sensor, cadence)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_oceancolor_dates(d0, d1, entry["start"], entry["id"])
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}")
    cache_dir = work_dir or oceancolor_cache_dir()

    urls = oceancolor_urls((minx, miny, maxx, maxy), d0, d1,
                           sensor=sensor, cadence=cadence, stride=stride,
                           erddap_base=erddap_base)
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()
    lats = lons = None
    times: List[str] = []
    frames: List[np.ma.MaskedArray] = []
    url_provs: List[Dict[str, Any]] = []
    n_payload_bytes = 0
    payload_for_hash = hashlib.sha256()
    for url in urls:
        try:
            payload, cprov = _cached_get_bytes(url, cache_dir,
                                              refresh=refresh)
        except Exception as exc:
            raise RuntimeError(
                f"CoastWatch ocean-color request failed "
                f"({type(exc).__name__}: {exc}) — URL: {url}. "
                "The CoastWatch ERDDAP was unreachable on 2026-09-27 "
                "(HTTP 502/503 on griddap); retry later or use "
                "source=\"cmems\" with CMEMS credentials."
            ) from exc
        try:
            pl_lats, pl_lons, pl_times, chl = _parse_coastwatch_bytes(payload)
        except ValueError:
            raise
        except Exception as exc:
            raise RuntimeError(
                f"CoastWatch ocean-color payload parse failed "
                f"({type(exc).__name__}: {exc}). netCDF4 is required "
                "(pip install netCDF4)."
            ) from exc
        lats, lons = pl_lats, pl_lons
        times.extend(pl_times)
        frames.append(chl)
        n_payload_bytes += len(payload)
        payload_for_hash.update(payload)
        url_provs.append(cprov)
    cube = np.ma.concatenate(frames, axis=0)
    gap_fractions = [float(np.ma.masked_invalid(cube[k]).mask.mean())
                     if cube[k].size else 1.0
                     for k in range(cube.shape[0])]
    return OceanColorField(
        values=cube, times=times, lats=lats, lons=lons,
        source=f"noaa-coastwatch/{entry['id']}",
        provenance={
            "dataset": entry["title"],
            "dataset_id": entry["id"],
            "erddap_base": erddap_base.rstrip("/"),
            "variable": OCEANCOLOR_VAR,
            "units": OCEANCOLOR_UNITS,
            "sensor": entry["sensor"],
            "processing_level": entry["processing"],
            "cadence": cadence,
            "urls": urls,
            "sha256": payload_for_hash.hexdigest(),
            "n_bytes": n_payload_bytes,
            "retrieved_at": retrieved_at,
            "bbox_requested": [minx, miny, maxx, maxy],
            "time_requested": [_erddap_iso(d0), _erddap_iso(d1)],
            "stride": stride,
            "grid": f"~{OCEANCOLOR_RES_DEG}-degree global "
                    f"({OCEANCOLOR_LON_MIN}..{OCEANCOLOR_LON_MAX} lon, "
                    f"{OCEANCOLOR_LAT_MIN}..{OCEANCOLOR_LAT_MAX} lat)",
            "gap_fractions": gap_fractions,
            "cloud_gaps": "NaN cells are cloud/land/missing — "
                          "never interpolated or filled",
            "cache": url_provs,
        },
    )


# ---------------------------------------------------------------------------
# CMEMS authenticated fallback
# ---------------------------------------------------------------------------


def fetch_oceancolor_cmems(
        bbox: Sequence[float], start: DateLike, end: DateLike,
        work_dir: Optional[str] = None,
        cmems_dataset_id: str = CMEMS_OCEANCOLOR_PRODUCT,
        variable: str = CMEMS_OCEANCOLOR_VARIABLE,
        refresh: bool = False) -> OceanColorField:
    """Fetch chlorophyll-a from CMEMS global ocean colour (authenticated).

    Uses :func:`cmems.require_toolbox` (fail-fast credential check —
    never an interactive prompt), then ``copernicusmarine.subset`` on
    the global ocean-colour product. Monthly, 4 km, 1997-present.

    The exact ``cmems_dataset_id`` string should be confirmed against
    the CMEMS catalogue when credentials are available (taken here
    from the Nov-2025 OC QUID); pass the catalogue's dataset id via
    ``cmems_dataset_id`` if it differs.
    """
    from .cmems import require_toolbox
    cm = require_toolbox()
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_oceancolor_dates(d0, d1, _dt.date(1997, 9, 1),
                               cmems_dataset_id)
    out_dir = work_dir or oceancolor_cache_dir()
    os.makedirs(out_dir, exist_ok=True)
    out_name = (f"oceancolor_cmems_{d0.date().isoformat()}_"
                f"{d1.date().isoformat()}_{minx}_{miny}_{maxx}_{maxy}.nc"
                ).replace(" ", "_")
    out_path = os.path.join(out_dir, out_name)
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()
    if refresh or not os.path.exists(out_path):
        cm.subset(
            dataset_id=cmems_dataset_id,
            variables=[variable],
            minimum_longitude=minx, maximum_longitude=maxx,
            minimum_latitude=miny, maximum_latitude=maxy,
            start_datetime=d0.isoformat(), end_datetime=d1.isoformat(),
            output_filename=os.path.basename(out_path),
            output_directory=out_dir,
        )
    try:
        from .sst_global import _require_netcdf4
        nc = _require_netcdf4()
    except ImportError as exc:
        raise ImportError(
            "parsing the CMEMS ocean-color NetCDF needs netCDF4 "
            "(pip install 'survey-currents[noaa]')") from exc
    with nc.Dataset(out_path, mode="r") as ds:
        if variable not in ds.variables:
            raise ValueError(
                f"CMEMS ocean-color file missing variable {variable!r}")
        var = ds.variables[variable]
        lat_name = "lat" if "lat" in ds.variables else "latitude"
        lon_name = "lon" if "lon" in ds.variables else "longitude"
        lats = np.asarray(ds.variables[lat_name][:], dtype=float)
        lons_raw = np.asarray(ds.variables[lon_name][:], dtype=float)
        data = np.ma.asarray(var[:], dtype=float)
        if data.ndim == 4:
            data = data[:, 0, :, :]
        time_var = ds.variables["time"]
        try:
            dts = nc.num2date(time_var[:], time_var.units)
        except Exception as exc:
            raise ValueError(
                f"could not decode CMEMS time axis: {exc}") from exc
    times = [_coerce_datetime_utc(d).isoformat() for d in dts]
    lons = ((lons_raw + 180.0) % 360.0) - 180.0
    order = np.argsort(lons)
    data = np.ma.masked_invalid(data[:, :, order])
    with open(out_path, "rb") as fh:
        digest = _sha256_bytes(fh.read())
    gap_fractions = [float(np.ma.getmaskarray(data[k]).mean())
                     if data[k].size else 1.0
                     for k in range(data.shape[0])]
    return OceanColorField(
        values=data, times=times, lats=lats, lons=lons[order],
        source=f"cmems:{cmems_dataset_id}",
        provenance={
            "dataset": "Copernicus Marine Global Ocean Colour "
                       "(OCEANCOLOUR_GLO_BGC_L4_MY_009_104)",
            "cmems_dataset_id": cmems_dataset_id,
            "variable": variable,
            "units": "mg m^-3",
            "cadence": "monthly",
            "local_path": out_path,
            "sha256": digest,
            "retrieved_at": retrieved_at,
            "bbox_requested": [minx, miny, maxx, maxy],
            "time_requested": [d0.date().isoformat(), d1.date().isoformat()],
            "grid": "4 km global",
            "gap_fractions": gap_fractions,
            "cloud_gaps": "NaN cells are cloud/land/missing — "
                          "never interpolated or filled",
        },
    )


# ---------------------------------------------------------------------------
# NASA OBPG fallback (credentials-required)
# ---------------------------------------------------------------------------


def _obpg_opener() -> urllib.request.OpenerDirector:
    """Build an Earthdata-authenticated opener; raise CredentialsMissing."""
    from .sst_global import CredentialsMissing, earthdata_credentials
    creds = earthdata_credentials()
    if creds is None:
        raise CredentialsMissing(
            "Ocean color via NASA OBPG needs Earthdata Login credentials, "
            "but none were found. Set EARTHDATA_USERNAME/EARTHDATA_PASSWORD "
            "or add a ~/.netrc entry for urs.earthdata.nasa.gov (register "
            "free at https://urs.earthdata.nasa.gov/users/new). "
            "The keyless NOAA CoastWatch source (source=\"coastwatch\") "
            "needs no credentials.")
    user, pw = creds
    # OBPG direct data access accepts Earthdata Login via HTTP Basic
    # Auth (the documented curl -u / wget --user approach).
    mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    mgr.add_password(None, "oceandata.sci.gsfc.nasa.gov", user, pw)
    return urllib.request.build_opener(
        urllib.request.HTTPBasicAuthHandler(mgr))


def _obpg_8day_periods(d0: _dt.date, d1: _dt.date):
    """OBPG 8-day bin periods overlapping [d0, d1] as (start, end) dates.

    OBPG 8-day periods are fixed day-of-year bins: 1-8, 9-16, ...,
    353-360, 361-365/366.
    """
    periods = []
    year = d0.year
    while _dt.date(year, 1, 1) <= d1:
        leap = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
        for start_doy in range(1, leap + 1, 8):
            end_doy = min(start_doy + 7, leap)
            p0 = _dt.date(year, 1, 1) + _dt.timedelta(days=start_doy - 1)
            p1 = _dt.date(year, 1, 1) + _dt.timedelta(days=end_doy - 1)
            if p1 >= d0 and p0 <= d1:
                periods.append((p0, p1))
        year += 1
    return periods


def _obpg_months(d0: _dt.date, d1: _dt.date):
    """Calendar months overlapping [d0, d1] as (first, last) dates."""
    months = []
    y, m = d0.year, d0.month
    while (y, m) <= (d1.year, d1.month):
        first = _dt.date(y, m, 1)
        last = (_dt.date(y + 1, 1, 1) if m == 12 else _dt.date(y, m + 1, 1)) \
            - _dt.timedelta(days=1)
        months.append((first, last))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return months


def obpg_file_names(start: DateLike, end: DateLike,
                    cadence: str = "monthly") -> List[str]:
    """OBPG Level-3 mapped filenames covering [start, end].

    Follows the documented OBPG file naming convention
    (``AQUA_MODIS.<start>_<end>.L3m.<DAY|8D|MO>.CHL.chlor_a.4km.nc`` —
    the stem was verified against real CMR granule titles 2026-09-27).
    Pure function: no network, no credentials — unit-testable.
    """
    cad = str(cadence).strip().lower()
    if cad not in ("monthly", "weekly", "daily"):
        raise ValueError(
            f"unknown OBPG cadence {cadence!r}; expected "
            "'monthly', 'weekly', or 'daily'")
    d0 = _coerce_datetime_utc(start).date()
    d1 = _coerce_datetime_utc(end).date()
    names = []
    if cad == "monthly":
        for p0, p1 in _obpg_months(d0, d1):
            names.append(
                f"{OBPG_FILE_SENSOR}.{p0.strftime('%Y%m%d')}_"
                f"{p1.strftime('%Y%m%d')}.L3m.MO.CHL."
                f"{OBPG_FILE_VARIABLE}.{OBPG_FILE_RESOLUTION}.nc")
    elif cad == "weekly":
        for p0, p1 in _obpg_8day_periods(d0, d1):
            names.append(
                f"{OBPG_FILE_SENSOR}.{p0.strftime('%Y%m%d')}_"
                f"{p1.strftime('%Y%m%d')}.L3m.8D.CHL."
                f"{OBPG_FILE_VARIABLE}.{OBPG_FILE_RESOLUTION}.nc")
    else:
        day = d0
        while day <= d1:
            names.append(
                f"{OBPG_FILE_SENSOR}.{day.strftime('%Y%m%d')}.L3m.DAY.CHL."
                f"{OBPG_FILE_VARIABLE}.{OBPG_FILE_RESOLUTION}.nc")
            day += _dt.timedelta(days=1)
    return names


def _parse_obpg_bytes(payload: bytes, url: str):
    """Parse a SeaDAS L3 mapped NetCDF into (lats, lons, times, chl).

    Defensive against the documented OBPG L3 mapped structure:
    ``chlor_a`` may be 2-D ``(lat, lon)`` or 3-D ``(time, lat, lon)``;
    ``scale_factor``/``add_offset`` are applied when present;
    ``_FillValue`` (or NaN) cells become masked. Times come from the
    file's ``time`` variable via ``num2date``, falling back to ``None``
    (the caller stamps the period midpoint) when absent.
    """
    from .sst_global import _require_netcdf4
    nc = _require_netcdf4()
    lats = lons = None
    times: List[str] = []
    chl = None
    with _open_nc_bytes(payload) as ds:
        if OBPG_FILE_VARIABLE not in ds.variables:
            raise ValueError(
                f"OBPG file {url} missing variable "
                f"{OBPG_FILE_VARIABLE!r}")
        var = ds.variables[OBPG_FILE_VARIABLE]
        lat_name = next((n for n in ("lat", "latitude") if n in ds.variables),
                        None)
        lon_name = next((n for n in ("lon", "longitude") if n in ds.variables),
                        None)
        if lat_name is None or lon_name is None:
            raise ValueError(f"OBPG file {url} missing lat/lon variables")
        lats = np.asarray(ds.variables[lat_name][:], dtype=float)
        lons = np.asarray(ds.variables[lon_name][:], dtype=float)
        data = np.ma.asarray(var[:], dtype=float)
        # SeaDAS packing: apply scale/offset when the file uses them.
        try:
            scale = float(var.getncattr("scale_factor"))
            offset = float(var.getncattr("add_offset"))
            data = data * scale + offset
        except Exception:
            pass
        fill = var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() \
            else None
        if fill is not None:
            try:
                data = np.ma.masked_values(data, float(fill),
                                           rtol=1e-5, atol=1e-3)
            except Exception:
                pass
        data = np.ma.masked_invalid(data)
        if data.ndim == 2:
            data = data[np.newaxis, :, :]
        elif data.ndim == 3:
            pass
        else:
            raise ValueError(
                f"OBPG chlor_a has unexpected shape {data.shape} in {url}")
        chl = data
        if "time" in ds.variables:
            tvar = ds.variables["time"]
            try:
                units = tvar.getncattr("units")
                dts = nc.num2date(tvar[:], units)
                times = [_coerce_datetime_utc(d).isoformat() for d in dts]
            except Exception:
                times = []
    if lats is None or lons is None or chl is None:
        raise ValueError(f"could not parse OBPG file {url}")
    return lats, lons, times, chl


def fetch_oceancolor_obpg(bbox: Sequence[float], start: DateLike,
                          end: DateLike, cadence: str = "monthly",
                          work_dir: Optional[str] = None,
                          refresh: bool = False) -> OceanColorField:
    """Fetch MODIS Aqua chlorophyll-a from NASA OBPG (Earthdata-authenticated).

    Downloads the documented OBPG Level-3 mapped period files
    (``AQUA_MODIS.<start>_<end>.L3m.<DAY|8D|MO>.CHL.chlor_a.4km.nc``)
    from ``oceandata.sci.gsfc.nasa.gov`` direct data access, subsets
    each to ``bbox``, and stacks the frames. Cloud/land/missing cells
    stay NaN — never interpolated.

    The exact download URLs land in provenance. This path was
    implemented against the documented OBPG file convention and the
    verified CMR collection, but could not be live-verified in this
    environment (no Earthdata credentials available and the
    CoastWatch outage meant the fallback chain was never exercised
    end-to-end) — see docs/DATA_SOURCES.md.

    Raises:
        CredentialsMissing: no Earthdata Login credentials found (or
            the server rejected them with HTTP 401).
        ValueError: invalid bbox / dates / cadence.
        RuntimeError: network/NetCDF failures (with context).
    """
    from .sst_global import CredentialsMissing
    opener = _obpg_opener()  # raises CredentialsMissing early
    minx, miny, maxx, maxy = validate_sst_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_oceancolor_dates(d0, d1, OBPG_RECORD_START, "NASA OBPG")
    names = obpg_file_names(d0, d1, cadence)
    out_dir = work_dir or oceancolor_cache_dir()
    os.makedirs(out_dir, exist_ok=True)
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()

    frames: List[np.ma.MaskedArray] = []
    times: List[str] = []
    urls: List[str] = []
    digests: List[str] = []
    lats = lons = None
    for name in names:
        url = f"{OBPG_DIRECTACCESS_BASE}/{name}"
        urls.append(url)
        cache_path = os.path.join(out_dir, f"obpg_{name}")
        if refresh or not os.path.exists(cache_path):
            try:
                with opener.open(url, timeout=120) as resp:
                    payload = resp.read()
            except urllib.error.HTTPError as exc:
                if exc.code == 401:
                    raise CredentialsMissing(
                        "NASA OBPG rejected the Earthdata credentials "
                        f"(HTTP 401) for {url}. Check EARTHDATA_USERNAME/"
                        "EARTHDATA_PASSWORD or the ~/.netrc entry.") from exc
                raise RuntimeError(
                    f"OBPG download failed for {url} "
                    f"({type(exc).__name__}: {exc})") from exc
            except Exception as exc:
                raise RuntimeError(
                    f"OBPG download failed for {url} "
                    f"({type(exc).__name__}: {exc})") from exc
            with open(cache_path, "wb") as fh:
                fh.write(payload)
        else:
            with open(cache_path, "rb") as fh:
                payload = fh.read()
        with open(cache_path, "rb") as fh:
            digests.append(_sha256_bytes(fh.read()))
        f_lats, f_lons, f_times, chl = _parse_obpg_bytes(payload, url)
        # Latitude may run north-to-south in SeaDAS files; normalize.
        if f_lats[0] > f_lats[-1]:
            f_lats = f_lats[::-1]
            chl = chl[:, ::-1, :]
        order = np.argsort(f_lons)
        f_lons = ((f_lons[order] + 180.0) % 360.0) - 180.0
        chl = chl[:, :, order]
        iy = np.where((f_lats >= miny) & (f_lats <= maxy))[0]
        ix = np.where((f_lons >= minx) & (f_lons <= maxx))[0]
        if iy.size == 0 or ix.size == 0:
            raise ValueError(
                f"bbox {list(bbox)} does not overlap OBPG grid of {url}")
        if lats is None:
            lats, lons = f_lats[iy], f_lons[ix]
        sub = chl[:, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
        # Period midpoint from the filename as the fallback stamp when
        # a file carries no decodable time variable (defensive — OBPG
        # L3 files normally have one).
        m = re.search(r"(\d{8})(?:_(\d{8}))?", name)
        if m:
            p0 = _dt.datetime.strptime(m.group(1), "%Y%m%d")
            p1 = _dt.datetime.strptime(m.group(2) or m.group(1), "%Y%m%d")
            fallback_stamp = (p0 + (p1 - p0) / 2).replace(
                tzinfo=_dt.timezone.utc).isoformat()
        else:
            fallback_stamp = d0.isoformat()
        for k in range(sub.shape[0]):
            stamp = f_times[k] if k < len(f_times) else fallback_stamp
            frames.append(np.ma.masked_invalid(sub[k]))
            times.append(stamp)
    if not frames:
        raise RuntimeError("OBPG fetch produced no frames")
    values = np.ma.stack(frames)
    gap_fractions = [float(np.ma.getmaskarray(values[k]).mean())
                     if values[k].size else 1.0
                     for k in range(values.shape[0])]
    return OceanColorField(
        values=values, times=times, lats=lats, lons=lons,
        source="nasa-obpg/AQUA_MODIS_L3m_CHL_chlor_a_4km",
        provenance={
            "dataset": "NASA OBPG MODIS Aqua Level-3 mapped "
                       "chlorophyll-a (direct data access)",
            "cmr_collection": f"{OBPG_CMR_SHORT_NAME} "
                              f"({OBPG_CMR_COLLECTION_ID})",
            "urls": urls,
            "sha256": digests,
            "retrieved_at": retrieved_at,
            "cadence": str(cadence).strip().lower(),
            "units": "mg m^-3",
            "bbox_requested": [minx, miny, maxx, maxy],
            "gap_fractions": gap_fractions,
            "cloud_gaps": "NaN cells are cloud/land/missing — "
                          "never interpolated or filled",
            "live_verified": False,
            "live_verified_note": "OBPG network path not live-verified "
                                  "in this environment (no Earthdata "
                                  "credentials available).",
        },
    )


# ---------------------------------------------------------------------------
# unified entry point
# ---------------------------------------------------------------------------


def fetch_oceancolor(bbox: Sequence[float], start: DateLike, end: DateLike,
                     product: str = "chlorophyll-a",
                     cadence: str = "monthly",
                     sensor: str = "modis-aqua",
                     source: str = "coastwatch",
                     work_dir: Optional[str] = None,
                     refresh: bool = False,
                     stride: int = 1,
                     erddap_base: str = COASTWATCH_ERDDAP_BASE,
                     **kwargs: Any) -> OceanColorField:
    """Fetch ocean color (chlorophyll-a) for ``bbox`` over [start, end].

    Args:
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180
            degrees (ascending).
        start, end: dates/datetimes/ISO strings.
        product: ``"chlorophyll-a"`` (the only product in v0.14.0;
            kd490/PAR are future work — anything else raises).
        cadence: ``"monthly"`` (default — the most cloud-complete),
            ``"weekly"``, ``"daily"``.
        sensor: ``"modis-aqua"`` (default — MODIS Aqua R2022,
            2002-present), ``"viirs-snpp"`` (VIIRS Suomi NPP science
            quality, 2012-
            present), ``"multi"`` (ESA OC-CCI v6.0 merge, 1997-present).
        source: ``"coastwatch"`` (default — keyless NOAA CoastWatch
            ERDDAP), ``"obpg"`` (NASA OBPG direct data access — needs
            a free Earthdata Login; the credentials-required
            fallback), ``"cmems"`` (authenticated Copernicus Marine —
            optional extra source), ``"auto"`` (try CoastWatch, then
            OBPG when Earthdata credentials exist, then CMEMS when
            CMEMS credentials exist).
        work_dir: cache directory override.
        refresh: ignore the cache and re-download.
        stride: time-axis index stride (>= 1).
        erddap_base: override the CoastWatch ERDDAP host (keyless
            mirrors; documented in docs/DATA_SOURCES.md).
        **kwargs: forwarded to the CMEMS path (``cmems_dataset_id``,
            ``variable``).

    Returns:
        :class:`OceanColorField` in mg m^-3 (masked array — cloud/land
        gaps stay NaN, never filled) with full provenance.

    Raises:
        ValueError: invalid bbox / dates / product / sensor / cadence.
        CredentialsMissing: ``source="obpg"`` without Earthdata Login
            credentials.
        RuntimeError: the CoastWatch ERDDAP is unreachable (with the
            2026-09-27 outage noted) or no working source exists.
    """
    if str(product).strip().lower() != "chlorophyll-a":
        raise ValueError(
            f"unknown ocean-color product {product!r}; v0.14.0 supports "
            "'chlorophyll-a' only (kd490/PAR are future work)")
    src = str(source).strip().lower()
    if src == "coastwatch":
        return fetch_oceancolor_coastwatch(
            bbox, start, end, sensor=sensor, cadence=cadence,
            stride=stride, erddap_base=erddap_base,
            work_dir=work_dir, refresh=refresh)
    if src == "obpg":
        return fetch_oceancolor_obpg(
            bbox, start, end, cadence=cadence,
            work_dir=work_dir, refresh=refresh)
    if src == "cmems":
        return fetch_oceancolor_cmems(bbox, start, end, work_dir=work_dir,
                                      refresh=refresh, **kwargs)
    if src == "auto":
        failures = []
        try:
            return fetch_oceancolor_coastwatch(
                bbox, start, end, sensor=sensor, cadence=cadence,
                stride=stride, erddap_base=erddap_base,
                work_dir=work_dir, refresh=refresh)
        except Exception as exc:
            failures.append(f"CoastWatch: {type(exc).__name__}: {exc}")
        try:
            from .sst_global import CredentialsMissing, earthdata_credentials
            if earthdata_credentials() is None:
                raise CredentialsMissing("no Earthdata credentials")
            return fetch_oceancolor_obpg(
                bbox, start, end, cadence=cadence,
                work_dir=work_dir, refresh=refresh)
        except Exception as exc:
            failures.append(f"OBPG: {type(exc).__name__}: {exc}")
        try:
            from .cmems import require_toolbox
            require_toolbox()
            return fetch_oceancolor_cmems(bbox, start, end,
                                          work_dir=work_dir,
                                          refresh=refresh, **kwargs)
        except Exception as exc:
            failures.append(f"CMEMS: {type(exc).__name__}: {exc}")
        raise RuntimeError(
            "ocean-color fetch failed on every source "
            f"({'; '.join(failures)}). The keyless CoastWatch ERDDAP "
            "was unreachable; for the authenticated paths register free "
            "at https://urs.earthdata.nasa.gov/users/new (OBPG) or "
            f"{CMEMS_REGISTER_URL} (CMEMS), or retry CoastWatch later.")
    raise ValueError(
        f"unknown ocean-color source {source!r}; expected "
        "'coastwatch', 'obpg', 'cmems', or 'auto'")


# ---------------------------------------------------------------------------
# offline demo
# ---------------------------------------------------------------------------


def main_demo() -> None:
    """Print a small offline summary (no network). Mirrors sst_global.main_demo."""
    f = OceanColorField.synthetic(nt=3)
    print(f"[oceancolor] synthetic {f.values.shape} mg/m^3, "
          f"t0={f.times[0]}, median={f.spatial_median(0):.3f} mg/m^3, "
          f"gap={f.gap_fraction(0):.0%}")


__all__ = [
    "DateLike",
    "COASTWATCH_ERDDAP_BASE",
    "COASTWATCH_DATASETS",
    "OCEANCOLOR_LON_MIN",
    "OCEANCOLOR_LON_MAX",
    "OCEANCOLOR_LAT_MIN",
    "OCEANCOLOR_LAT_MAX",
    "OCEANCOLOR_VAR",
    "OCEANCOLOR_UNITS",
    "OCEANCOLOR_MAX_CACHE_AGE_DAYS",
    "CMEMS_OCEANCOLOR_PRODUCT",
    "CMEMS_OCEANCOLOR_VARIABLE",
    "OBPG_DIRECTACCESS_BASE",
    "OBPG_CMR_COLLECTION_ID",
    "OBPG_CMR_SHORT_NAME",
    "OBPG_RECORD_START",
    "obpg_file_names",
    "fetch_oceancolor_obpg",
    "OceanColorField",
    "oceancolor_cache_dir",
    "oceancolor_urls",
    "fetch_oceancolor",
    "fetch_oceancolor_coastwatch",
    "fetch_oceancolor_cmems",
    "main_demo",
]

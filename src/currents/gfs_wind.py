"""NOAA GFS 10-m winds + 2-m air temperature via the NOMADS GRIB filter.

:func:`fetch_gfs_wind` pulls 10-m wind and 2-m temperature fields from
the 0.25-degree Global Forecast System (``00``/``06``/``12``/``18`` UTC
cycles) via NCEP's NOMADS server — completely **keyless** (plain HTTPS,
no account, no token), which is what makes it usable from unattended
automation where credentialed sources (OSCAR, CMEMS, ERA5/CDS, FIRMS)
cannot run.

The request goes through the GRIB filter CGI
(``filter_gfs_0p25.pl``)::

    https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl
        ?dir=/gfs.YYYYMMDD/HH/atmos
        &file=gfs.t{HH}z.pgrb2.0p25.f{HHH}
        &subregion=on
        &leftlon=..&rightlon=..&toplat=..&bottomlat=..
        &lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on
        &lev_2_m_above_ground=on&var_TMP=on

``subregion=on`` is what makes the filter actually honor the bounding
box (verified 2026-10-01: without it the filter silently returns the
full 1440x721 global grid). The adapter still **crops locally with
numpy after the read** and verifies the returned grid covers the
requested bbox — defense in depth against the filter's quirks.

.. note::
    NOMADS retired the DODS/OpenDAP interface in 2025 (NWS notice
    scn25-81). The GRIB filter CGI above is the working keyless path;
    do not "upgrade" this module to an OpenDAP URL.

Fetched fields (one snapshot per ``(day, forecast_hour)`` step):

* ``"u10"`` / ``"v10"`` — 10-m wind components, m/s (``10u``/``10v`` on
  the wire)
* ``"t2m"`` — 2-m air temperature, **°C here** (Kelvin on the wire)

``forecast_hours=(0,)`` is the f000 **analysis** snapshot per sampled
day. The GFS 0.25° is hourly through f120, so e.g.
``forecast_hours=(0, 1, 6)`` assembles a sub-daily time series whose
timestamps carry the forecast hour (``field.times`` renders ``HH:MM`` in
viz). The valid time of each step is taken from the GRIB message's own
``validityDate``/``validityTime`` and cross-checked against
``cycle + forecast_hour`` — a served step that does not match its
requested hour is refused, never mislabeled.

:class:`GfsWindField` matches the :class:`Era5Field` shape
(``grids``/``times``/``lats``/``lons`` plus ``values``/``overlay_grids``)
so it flows through the ``survey-viz`` ``wind`` variable path
unchanged, and additionally exposes ``air_temperature`` (3-D, **°F** —
the warming.watch strand-color convention ``survey-viz``'s
``dark_strands`` preset reads) and ``temperature_unit``. A
:class:`CurrentField` was *not* reused: viz's wind path keys off
``grids["u10"]``/``grids["v10"]`` and ``spec.variable == "wind"``, which
the ocean-current shape does not satisfy (documented choice).

**Honest limits** (verified live 2026-10-01):

* NOMADS keeps roughly the **last 10 days** of 0.25° GFS
  (``GFS_RETENTION_DAYS``). Dates outside ``[today-9, today]`` raise
  :class:`UnavailableRangeError` — never silent, never padded.
* A cycle that has not posted yet (e.g. today's ``18`` at 09:00 UTC)
  returns HTTP 404 from NOMADS, which surfaces as
  :class:`UnavailableRangeError` naming the exact URL. The same
  applies per forecast hour: a missing ``(day, hour)`` is never
  skipped or padded.
* f000 is the **analysis**, f001..f120 are hourly **forecasts**;
  ``stride_days >= 1`` only (days, not hours, are the coarse axis).
* **Download cost is honest and per-step**: each ``(day, hour,
  window)`` is one ~2.4 MB GRIB for a North-America box. A full
  ``forecast_hours=tuple(range(121))`` day is ~290 MB — the caller
  (reel pipeline, CLI) owns the frame budget, not this adapter.

``cfgrib`` is a lazy import — the module imports cleanly without it,
and :func:`fetch_gfs_wind` raises an actionable ``ImportError``
(``pip install survey-currents[gfs]``) when it is missing. Note
``cfgrib.open_file`` **fails** on these payloads (mixed 10-m / 2-m
levels break its dataset builder); the parser iterates
``cfgrib.messages.FileStream`` and assembles the numpy arrays manually.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import os
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .glsea import _coerce_datetime_utc
from .sst_global import validate_sst_bbox

DateLike = Union[_dt.date, _dt.datetime, str]

# ---------------------------------------------------------------------------
# GFS / NOMADS constants
# ---------------------------------------------------------------------------

#: NOMADS GRIB-filter CGI for the 0.25° GFS (keyless HTTPS).
NOMADS_FILTER_BASE = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl"
#: GFS analysis cycles served under /gfs.YYYYMMDD/HH/atmos.
GFS_CYCLES = ("00", "06", "12", "18")
#: Forecast hour requested: f000 = the analysis. The default — a daily
#: reel of analyses is what the wind-strand use case needs. Pass other
#: hours via ``fetch_gfs_wind(..., forecast_hours=(0, 1, 6))``.
GFS_FORECAST_HOUR = "000"
#: Highest forecast hour served hourly by the GFS 0.25° (f000..f120).
#: Verified live 2026-10-01: f001 and f006 both returned HTTP 200 through
#: the same GRIB filter path.
GFS_FORECAST_HOUR_MAX = 120
#: Native grid spacing, degrees.
GFS_RES = 0.25
#: NOMADS retention of the 0.25° GFS, in days including today.
#: Verified live 2026-10-01: the f000 00z subset returned HTTP 200 for
#: 2026-09-22..2026-10-01 (10 days) and HTTP 404 for 2026-09-21.
GFS_RETENTION_DAYS = 10

#: Wire shortNames we request -> canonical grid keys, on-wire units,
#: and the converted units stored on the field.
GFS_GRIDS: Dict[str, Dict[str, str]] = {
    "10u": {"grid": "u10", "wire_units": "m s**-1", "units": "m/s"},
    "10v": {"grid": "v10", "wire_units": "m s**-1", "units": "m/s"},
    "2t": {"grid": "t2m", "wire_units": "K", "units": "\u00b0C"},
}

#: Kelvin -> Celsius offset applied to 2t on ingest.
_KELVIN_OFFSET = 273.15


class UnavailableRangeError(ValueError):
    """The requested GFS date/cycle is not on NOMADS (honest, not silent).

    NOMADS keeps roughly the last :data:`GFS_RETENTION_DAYS` days of the
    0.25° GFS, and a cycle that has not posted yet returns HTTP 404.
    Either way the caller gets this error naming the exact request URL —
    the adapter never pads, interpolates, or silently substitutes.
    """


# ---------------------------------------------------------------------------
# GfsWindField — canonical multi-variable GFS wind field model
# ---------------------------------------------------------------------------


@dataclass
class GfsWindField:
    """Time-indexed GFS 10-m wind + 2-m temperature grids.

    ``grids`` maps canonical grid names (``"u10"``, ``"v10"`` in m/s,
    ``"t2m"`` in °C) to ``(nt, ny, nx)`` numpy masked arrays. Follows the
    :class:`Era5Field` conventions (``times``/``lats``/``lons``,
    ``values``, ``overlay_grids``) so ``survey-viz``'s ``wind`` variable
    path consumes it unchanged, and additionally exposes
    ``air_temperature`` (3-D, °F — the strand-color scalar
    ``survey-viz``'s ``dark_strands`` preset reads on wind) and
    ``temperature_unit``.
    """

    grids: Dict[str, np.ma.MaskedArray]  # grid name -> (nt, ny, nx)
    times: List[str]                     # ISO-8601 timestamps, one per step
    lats: np.ndarray                     # (ny,) degrees north, increasing
    lons: np.ndarray                     # (nx,) degrees east, -180..180, increasing
    cycle: str = "00"                    # GFS analysis cycle ("00"/"06"/"12"/"18")
    units: Dict[str, str] = field(default_factory=dict)
    crs: str = "EPSG:4326"
    source: str = "noaa-nomads/gfs-0p25"
    provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        ny = self.lats.shape[0]
        nx = self.lons.shape[0]
        for name, grid in self.grids.items():
            arr = np.ma.asarray(grid, dtype=float)
            if arr.shape != (nt, ny, nx):
                raise ValueError(
                    f"grid {name!r} shape {arr.shape} does not match "
                    f"(nt, ny, nx)=({nt}, {ny}, {nx})")
            self.grids[name] = arr
        self.lats = np.asarray(self.lats, dtype=float)
        self.lons = np.asarray(self.lons, dtype=float)
        missing = [g for g in ("u10", "v10", "t2m") if g not in self.grids]
        if missing:
            raise ValueError(
                f"GfsWindField needs grids {missing} "
                "(u10/v10 10-m wind + t2m 2-m temperature)")
        if self.cycle not in GFS_CYCLES:
            raise ValueError(
                f"cycle must be one of {GFS_CYCLES}, got {self.cycle!r}")
        if not self.units:
            self.units = {spec["grid"]: spec["units"]
                          for spec in GFS_GRIDS.values()}
        self.bounds = (float(self.lons[0]), float(self.lats[0]),
                       float(self.lons[-1]), float(self.lats[-1]))

    # -- derived quantities --------------------------------------------------

    @property
    def wind_speed(self) -> np.ma.MaskedArray:
        """Wind speed magnitude ``sqrt(u10^2 + v10^2)`` in m/s, (nt, ny, nx)."""
        return np.ma.sqrt(self.grids["u10"] ** 2 + self.grids["v10"] ** 2)

    @property
    def values(self) -> np.ma.MaskedArray:
        """The rendered base grid: wind speed in m/s, (nt, ny, nx).

        This is what generic scalar-grid consumers (``survey-viz``'s
        ``render_viz`` documented dict form, the reel-studio cache
        fingerprint) pick up.
        """
        return self.wind_speed

    @property
    def air_temperature(self) -> np.ma.MaskedArray:
        """2-m air temperature in **°F**, (nt, ny, nx).

        The warming.watch convention ``survey-viz``'s ``dark_strands``
        preset colors wind strands by (not by wind speed). The field
        stores t2m in °C; this property converts to the reference unit.
        """
        return self.grids["t2m"] * 9.0 / 5.0 + 32.0

    @property
    def temperature_unit(self) -> str:
        """Unit of :attr:`air_temperature`: always ``"°F"``."""
        return "\u00b0F"

    @property
    def overlay_grids(self) -> Dict[str, np.ma.MaskedArray]:
        """Non-wind grids as render-ready 3-D grids: ``{"t2m": ...}`` in °C."""
        return {"t2m": self.grids["t2m"]}

    def overlay_grid(self, name: str, index: int) -> np.ma.MaskedArray:
        """2-D overlay grid for ``name`` at timestep ``index``."""
        grids = self.overlay_grids
        if name not in grids:
            raise KeyError(
                f"no overlay {name!r} on this field; available: "
                f"{sorted(grids)}")
        return grids[name][index]

    def spatial_mean(self, variable: str, index: int) -> float:
        """Masked/NaN-aware spatial mean of ``variable`` at timestep ``index``.

        ``variable`` is ``"wind"`` (m/s), ``"t2m"`` (°C), or
        ``"air_temperature"`` (°F).
        """
        if variable == "wind":
            grid = self.wind_speed
        elif variable == "t2m":
            grid = self.grids["t2m"]
        elif variable == "air_temperature":
            grid = self.air_temperature
        else:
            raise KeyError(
                f"unknown variable {variable!r}; "
                "choose 'wind', 't2m', or 'air_temperature'")
        return float(np.ma.masked_invalid(grid[index]).mean())

    # -- selection -----------------------------------------------------------

    def select_time(self, index: int) -> "GfsWindField":
        """Return the single-timestep field at ``index``."""
        return GfsWindField(
            grids={k: v[index:index + 1] for k, v in self.grids.items()},
            times=[self.times[index]], lats=self.lats, lons=self.lons,
            cycle=self.cycle, units=dict(self.units), crs=self.crs,
            source=self.source, provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "GfsWindField":
        """Spatial subset to ``(min_lon, min_lat, max_lon, max_lat)``."""
        minx, miny, maxx, maxy = (float(x) for x in bbox)
        ix = np.where((self.lons >= minx) & (self.lons <= maxx))[0]
        iy = np.where((self.lats >= miny) & (self.lats <= maxy))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(
                f"bbox {list(bbox)} does not overlap field bounds {self.bounds}")
        sl = (slice(None), slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1))
        return GfsWindField(
            grids={k: v[sl] for k, v in self.grids.items()},
            times=list(self.times),
            lats=self.lats[iy[0]:iy[-1] + 1], lons=self.lons[ix[0]:ix[-1] + 1],
            cycle=self.cycle, units=dict(self.units), crs=self.crs,
            source=self.source, provenance=dict(self.provenance),
        )

    # -- serialization -------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dict (masked cells become NaN).

        Carries ``grids`` (u10/v10 m/s, t2m °C), ``air_temperature``
        (°F — what ``survey-viz``'s ``dark_strands`` reads on wind),
        ``temperature_unit``, and ``values`` (wind speed, m/s).
        """
        return {
            "grids": {k: np.ma.filled(v, np.nan).tolist()
                      for k, v in self.grids.items()},
            "air_temperature": np.ma.filled(self.air_temperature,
                                           np.nan).tolist(),
            "temperature_unit": self.temperature_unit,
            "values": np.ma.filled(self.values, np.nan).tolist(),
            "times": list(self.times),
            "lats": self.lats.tolist(),
            "lons": self.lons.tolist(),
            "cycle": self.cycle,
            "units": dict(self.units),
            "crs": self.crs,
            "source": self.source,
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GfsWindField":
        """Rebuild from :meth:`to_dict` (NaN -> masked). Raises on missing keys."""
        required = ("grids", "times", "lats", "lons")
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(f"GfsWindField dict missing keys: {missing}")
        return cls(
            grids={k: np.ma.masked_invalid(np.asarray(v, dtype=float))
                   for k, v in data["grids"].items()},
            times=list(data["times"]),
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            cycle=data.get("cycle", "00"),
            units=dict(data.get("units", {})),
            crs=data.get("crs", "EPSG:4326"),
            source=data.get("source", "noaa-nomads/gfs-0p25"),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        import json
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "GfsWindField":
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
        start: str = "2024-01-01T00:00:00+00:00",
        step_days: int = 1,
        cycle: str = "00",
        seed: int = 11,
        source: str = "synthetic",
    ) -> "GfsWindField":
        """Deterministic synthetic GFS wind field: westerlies + a warm pool.

        Stdlib+numpy only. Used by the test suite, the offline demo, and
        the reel-studio render path test.
        """
        rng = np.random.default_rng(seed)
        lat = np.linspace(lats[0], lats[1], ny)
        lon = np.linspace(lons[0], lons[1], nx)
        lon2d, lat2d = np.meshgrid(lon, lat)
        t0 = _dt.datetime.fromisoformat(start)
        times = [(t0 + _dt.timedelta(days=k * step_days)).isoformat()
                 for k in range(nt)]
        # Zonal westerlies strengthening poleward + a cyclonic swirl.
        u = 8.0 + 0.25 * np.abs(lat2d) + rng.normal(0, 1.5, (ny, nx))
        v = 4.0 * np.exp(-((lon2d / 40.0) ** 2 + (lat2d / 25.0) ** 2))
        v = v * np.sign(lon2d + 1e-9) + rng.normal(0, 1.0, (ny, nx))
        base_t = 28.0 - 0.6 * np.abs(lat2d)  # °C
        grids = {
            "u10": np.ma.array(
                np.repeat(u[None, :, :], nt, axis=0)
                + rng.normal(0, 0.5, (nt, ny, nx)),
                mask=np.zeros((nt, ny, nx), dtype=bool)),
            "v10": np.ma.array(
                np.repeat(v[None, :, :], nt, axis=0)
                + rng.normal(0, 0.5, (nt, ny, nx)),
                mask=np.zeros((nt, ny, nx), dtype=bool)),
            "t2m": np.ma.array(
                np.repeat(base_t[None, :, :], nt, axis=0)
                + rng.normal(0, 0.5, (nt, ny, nx)),
                mask=np.zeros((nt, ny, nx), dtype=bool)),
        }
        return cls(grids=grids, times=times, lats=lat, lons=lon,
                   cycle=cycle, source=source,
                   provenance={"synthetic": True, "seed": seed})


# ---------------------------------------------------------------------------
# bbox / cycle / date validation
# ---------------------------------------------------------------------------


def validate_gfs_bbox(bbox: Sequence[float]) -> Tuple[float, float, float, float]:
    """Validate ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180 degrees.

    Antimeridian-crossing boxes (``max_lon < min_lon``) are valid and are
    fetched as two NOMADS subregion requests internally (the GFS grid is
    global).
    """
    return validate_sst_bbox(bbox)


def gfs_lon_windows(bbox: Sequence[float]
                    ) -> List[Tuple[float, float, float, float]]:
    """Split a validated bbox into NOMADS subregion windows.

    * Antimeridian-crossing bboxes (``max_lon < min_lon``) become two
      windows: ``(min_lon, 180]`` and ``[-180, max_lon)``, concatenated
      along longitude after the read.
    * Otherwise a single window.
    """
    minx, miny, maxx, maxy = validate_gfs_bbox(bbox)
    if maxx < minx:  # antimeridian crossing
        return [(minx, miny, 180.0, maxy), (-180.0, miny, maxx, maxy)]
    return [(minx, miny, maxx, maxy)]


def _validate_forecast_hours(
        forecast_hours: Sequence[int]) -> Tuple[int, ...]:
    """Validate ``forecast_hours`` and normalize to sorted, deduped ints.

    GFS 0.25° is hourly f000..f120 (:data:`GFS_FORECAST_HOUR_MAX`).
    Raises ``ValueError`` for non-ints (bools included), empties, or
    hours outside 0..120. Returns the hours sorted ascending and
    deduplicated, so every fetch assembles its per-day time series in
    chronological order regardless of the order the caller passed.
    """
    if isinstance(forecast_hours, bool):
        raise ValueError(
            f"forecast_hours must be a sequence of ints 0.."
            f"{GFS_FORECAST_HOUR_MAX}, got {forecast_hours!r}")
    try:
        hours = tuple(forecast_hours)
    except TypeError:
        raise ValueError(
            f"forecast_hours must be a sequence of ints 0.."
            f"{GFS_FORECAST_HOUR_MAX}, got {forecast_hours!r}") from None
    if not hours:
        raise ValueError(
            f"forecast_hours must be a non-empty sequence of ints 0.."
            f"{GFS_FORECAST_HOUR_MAX}")
    for h in hours:
        if isinstance(h, bool) or not isinstance(h, int):
            raise ValueError(
                f"forecast_hours must be ints 0..{GFS_FORECAST_HOUR_MAX}, "
                f"got {h!r}")
        if not 0 <= h <= GFS_FORECAST_HOUR_MAX:
            raise ValueError(
                f"forecast hour {h} is outside 0..{GFS_FORECAST_HOUR_MAX} "
                "(GFS 0.25° is hourly f000..f120)")
    return tuple(sorted(set(hours)))


def _validate_gfs_cycle(cycle: str) -> str:
    c = str(cycle).strip()
    if c not in GFS_CYCLES:
        raise ValueError(
            f"cycle must be one of {GFS_CYCLES}, got {cycle!r}")
    return c


def gfs_retention_window(today: Optional[_dt.date] = None
                         ) -> Tuple[_dt.date, _dt.date]:
    """Oldest and newest fetchable GFS dates (inclusive).

    NOMADS keeps roughly the last :data:`GFS_RETENTION_DAYS` days of the
    0.25° GFS (verified live 2026-10-01).
    """
    today = today or _dt.date.today()
    oldest = today - _dt.timedelta(days=GFS_RETENTION_DAYS - 1)
    return oldest, today


def _validate_gfs_dates(d0: _dt.datetime, d1: _dt.datetime) -> None:
    if d0 > d1:
        raise ValueError(f"start {d0.date()} is after end {d1.date()}")
    oldest, newest = gfs_retention_window()
    if d0.date() < oldest or d1.date() > newest:
        raise UnavailableRangeError(
            f"GFS 10-m winds are only on NOMADS for the last "
            f"{GFS_RETENTION_DAYS} days ({oldest.isoformat()}.."
            f"{newest.isoformat()}); requested "
            f"{d0.date().isoformat()}..{d1.date().isoformat()}. NOMADS "
            "rolls old cycles off — this is a retention limit, not a "
            "transient failure.")


def gfs_sample_plan(d0: _dt.datetime, d1: _dt.datetime,
                    stride_days: int) -> List[_dt.date]:
    """Sample ``[d0, d1]`` every ``stride_days`` days (dates, ascending)."""
    if not isinstance(stride_days, int) or stride_days < 1:
        raise ValueError(f"stride_days must be an int >= 1, got {stride_days!r}")
    plan: List[_dt.date] = []
    cursor = d0.date()
    while cursor <= d1.date():
        plan.append(cursor)
        cursor += _dt.timedelta(days=stride_days)
    return plan


# ---------------------------------------------------------------------------
# NOMADS GRIB-filter request construction
# ---------------------------------------------------------------------------


def gfs_filter_url(day: _dt.date, cycle: str,
                   window: Sequence[float],
                   forecast_hour: int = 0) -> str:
    """Build the NOMADS GRIB-filter URL for one GFS forecast step.

    ``window`` is ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180
    (ascending — antimeridian windows are split by
    :func:`gfs_lon_windows` before this is called). ``subregion=on`` is
    what makes the filter honor the box; the adapter still verifies and
    crops locally after the read. ``forecast_hour`` selects the file
    ``gfs.t{HH}z.pgrb2.0p25.f{HHH}`` (f000 = analysis, f001..f120 =
    hourly forecasts); it defaults to ``0`` so existing callers keep
    getting the analysis URL.
    """
    minx, miny, maxx, maxy = (float(x) for x in window)
    c = _validate_gfs_cycle(cycle)
    hours = _validate_forecast_hours((forecast_hour,))
    fff = f"{hours[0]:03d}"
    datestr = day.strftime("%Y%m%d")
    params = (
        f"dir=/gfs.{datestr}/{c}/atmos"
        f"&file=gfs.t{c}z.pgrb2.0p25.f{fff}"
        "&subregion=on"
        f"&leftlon={minx}&rightlon={maxx}&toplat={maxy}&bottomlat={miny}"
        "&lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on"
        "&lev_2_m_above_ground=on&var_TMP=on"
    )
    return f"{NOMADS_FILTER_BASE}?{params}"


# ---------------------------------------------------------------------------
# download + cache (SHA-256-verified, atomic writes)
# ---------------------------------------------------------------------------


def gfs_cache_dir(work_dir: Optional[str] = None) -> str:
    """Cache root for GFS GRIB payloads.

    ``work_dir`` when given, else ``$SURVEY_CURRENTS_CACHE/gfs-wind``,
    else ``~/.cache/survey-currents/gfs-wind``. Created on demand.
    """
    if work_dir:
        root = os.path.abspath(os.path.expanduser(work_dir))
    else:
        root = os.environ.get("SURVEY_CURRENTS_CACHE") or os.path.join(
            os.path.expanduser("~"), ".cache", "survey-currents")
        root = os.path.join(root, "gfs-wind")
    os.makedirs(root, exist_ok=True)
    return root


def _cache_key(url: str) -> str:
    return hashlib.sha256(url.encode("utf-8")).hexdigest()


def _cached_get_bytes(url: str, cache_dir: str) -> Tuple[bytes, bool]:
    """GET ``url`` with a SHA-256-verified byte cache.

    Returns ``(payload, from_cache)``. Cache entries are keyed by the
    SHA-256 of the URL (one entry per day/cycle/window — GFS analyses are
    immutable once posted, so entries never expire). Writes are atomic
    (temp file + rename) with a JSON sidecar recording the URL, payload
    SHA-256, byte count, and retrieval time.

    NOMADS HTTP errors surface as :class:`UnavailableRangeError` naming
    the URL: a 404 means the date/cycle has rolled off NOMADS or the
    cycle has not posted yet — never a silent skip.
    """
    key = _cache_key(url)
    data_path = os.path.join(cache_dir, key + ".grib2")
    meta_path = os.path.join(cache_dir, key + ".json")
    if os.path.exists(data_path) and os.path.exists(meta_path):
        with open(data_path, "rb") as fh:
            return fh.read(), True
    req = urllib.request.Request(
        url, headers={"User-Agent": "survey-currents/gfs-wind"})
    try:
        with urllib.request.urlopen(req, timeout=300) as resp:
            payload = resp.read()
    except urllib.error.HTTPError as exc:
        raise UnavailableRangeError(
            f"NOMADS returned HTTP {exc.code} for {url}. "
            "The requested GFS date/cycle has likely rolled off NOMADS "
            f"(retention ~{GFS_RETENTION_DAYS} days) or the cycle has not "
            "posted yet."
        ) from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"GFS download failed ({type(exc).__name__}: {exc}); "
            "check the network connection and that "
            "nomads.ncep.noaa.gov is reachable."
        ) from exc
    digest = hashlib.sha256(payload).hexdigest()
    import json
    tmp_fd, tmp_path = tempfile.mkstemp(dir=cache_dir, suffix=".grib2")
    try:
        with os.fdopen(tmp_fd, "wb") as fh:
            fh.write(payload)
        os.replace(tmp_path, data_path)
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
    meta = {
        "url": url,
        "sha256": digest,
        "n_bytes": len(payload),
        "retrieved_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
    }
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)
    return payload, False


# ---------------------------------------------------------------------------
# GRIB2 parsing (cfgrib, lazy import; message iteration, not open_file)
# ---------------------------------------------------------------------------


def _require_cfgrib():
    """Lazy import of cfgrib with an actionable error."""
    try:
        from cfgrib.messages import FileStream  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "fetch_gfs_wind needs the cfgrib package "
            "(pip install 'survey-currents[gfs]'); the survey-currents "
            "engine itself has no hard dependency on it."
        ) from exc
    return FileStream


def _parse_grib_payload(payload: bytes, url: str, day: _dt.date,
                        cycle: str, forecast_hour: int = 0
                        ) -> Dict[str, Dict[str, Any]]:
    """Parse one NOMADS GRIB-filter payload into per-variable grids.

    Iterates ``cfgrib.messages.FileStream`` (``cfgrib.open_file`` fails
    on these payloads — the mixed 10-m / 2-m levels break its dataset
    builder) and returns, per canonical grid key (``"u10"``, ``"v10"``,
    ``"t2m"``), ``{"values": (Ny, Nx) float array, "lats": (Ny,),
    "lons": (Nx,) (degrees east, -180..180, increasing), "data_date": str,
    "data_time": str, "valid_time": str}``.

    ``forecast_hour`` is the requested file step (f000 = analysis,
    f001..f120 = hourly forecasts). Each message's ``stepRange`` must
    equal it, and the message's ``validityDate``/``validityTime`` must
    equal ``cycle + forecast_hour`` — the valid time is what lands on
    the field's ``times`` axis, so a served step that does not match
    its requested hour is refused rather than mislabeled.

    Raises :class:`RuntimeError` when a requested variable is missing,
    when a message's data date/cycle or forecast step does not match
    the request (the filter must never silently serve the wrong day or
    the wrong hour).
    """
    FileStream = _require_cfgrib()
    with tempfile.NamedTemporaryFile(suffix=".grib2", delete=False) as tmp:
        tmp.write(payload)
        tmp_path = tmp.name
    try:
        # NOTE: FileStream.items() cannot be list()-ed (cfgrib's
        # FileStreamItems lacks the mapping protocol list() probes);
        # iterate it instead.
        messages = [(off, m) for off, m in FileStream(tmp_path).items()]
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
    if not messages:
        raise RuntimeError(
            f"NOMADS payload for {url} contained no GRIB messages.")
    want = {short: spec["grid"] for short, spec in GFS_GRIDS.items()}
    cycle_dt = _dt.datetime(day.year, day.month, day.day,
                            int(cycle), tzinfo=_dt.timezone.utc)
    expected_valid = cycle_dt + _dt.timedelta(hours=forecast_hour)
    found: Dict[str, Dict[str, Any]] = {}
    for _offset, m in messages:
        short = m.get("shortName")
        if short not in want:
            continue
        key = want[short]
        if key in found:
            raise RuntimeError(
                f"duplicate {short!r} message in NOMADS payload for {url}; "
                "refusing to guess which one is the analysis")
        data_date = str(m.get("dataDate"))
        data_time = str(m.get("dataTime")).zfill(4)[:2]
        if data_date != day.strftime("%Y%m%d") or data_time != cycle:
            raise RuntimeError(
                f"NOMADS served the wrong day/cycle for {url}: message "
                f"dataDate={data_date} dataTime={data_time}, requested "
                f"{day.strftime('%Y%m%d')}/{cycle}. Refusing to mislabel.")
        step = str(m.get("stepRange"))
        if step != str(forecast_hour):
            raise RuntimeError(
                f"NOMADS served the wrong forecast step for {url}: message "
                f"stepRange={step!r}, requested f{forecast_hour:03d}. "
                "Refusing to mislabel.")
        # The valid time is authoritative from the payload itself; it
        # must agree with cycle + forecast_hour (it crosses midnight
        # for e.g. 18z f006 — that is expected, not an error).
        validity = _dt.datetime.strptime(
            f"{int(m.get('validityDate')):08d}"
            f"{int(m.get('validityTime')):04d}",
            "%Y%m%d%H%M").replace(tzinfo=_dt.timezone.utc)
        if validity != expected_valid:
            raise RuntimeError(
                f"NOMADS payload for {url} carries valid time "
                f"{validity.isoformat()} but cycle {cycle}z + "
                f"f{forecast_hour:03d} is {expected_valid.isoformat()}. "
                "Refusing to mislabel.")
        ny, nx = int(m["Ny"]), int(m["Nx"])
        vals = np.asarray(m["values"], dtype=float).reshape(ny, nx)
        lat0 = float(m["latitudeOfFirstGridPointInDegrees"])
        lat1 = float(m["latitudeOfLastGridPointInDegrees"])
        lon0 = float(m["longitudeOfFirstGridPointInDegrees"])
        lon1 = float(m["longitudeOfLastGridPointInDegrees"])
        lats = np.linspace(lat0, lat1, ny)
        lons = np.linspace(lon0, lon1, nx)
        # Normalize to -180..180, increasing.
        lons = ((lons + 180.0) % 360.0) - 180.0
        order = np.argsort(lons, kind="stable")
        lons = lons[order]
        vals = vals[:, order]
        if lats[0] > lats[-1]:
            lats = lats[::-1]
            vals = vals[::-1, :]
        found[key] = {"values": vals, "lats": lats, "lons": lons,
                      "data_date": data_date, "data_time": data_time,
                      "valid_time": validity.isoformat()}
    missing = [k for k in ("u10", "v10", "t2m") if k not in found]
    if missing:
        raise RuntimeError(
            f"NOMADS payload for {url} is missing variables {missing}; "
            "the filter request may have been truncated.")
    return found


def _verify_coverage(parsed: Dict[str, Dict[str, Any]],
                    window: Sequence[float], url: str) -> None:
    """Verify the returned grid covers the requested window (honest check).

    The filter has silently returned the wrong region in the past (missing
    ``subregion=on``); if the parsed grid does not cover the requested
    window within one grid step, fail loudly instead of mislabeling data.
    """
    minx, miny, maxx, maxy = (float(x) for x in window)
    tol = GFS_RES
    for key, part in parsed.items():
        lats, lons = part["lats"], part["lons"]
        if not (lats[0] <= miny + tol and lats[-1] >= maxy - tol
                and lons[0] <= minx + tol and lons[-1] >= maxx - tol):
            raise RuntimeError(
                f"NOMADS served a grid for {url} that does not cover the "
                f"requested window {list(window)}: grid spans lon "
                f"{lons[0]:.2f}..{lons[-1]:.2f}, lat {lats[0]:.2f}.."
                f"{lats[-1]:.2f}. Refusing to crop/mislabel.")


def _crop_to_window(parsed: Dict[str, Dict[str, Any]],
                    window: Sequence[float]) -> Tuple[np.ndarray, np.ndarray,
                                                     Dict[str, np.ndarray]]:
    """Crop parsed grids to the exact requested window (local, numpy).

    Returns ``(lats, lons, {grid: (ny, nx)})``. The filter usually honors
    the subregion request; this crop is defense in depth (and snaps
    off-grid request edges to the 0.25° grid).
    """
    minx, miny, maxx, maxy = (float(x) for x in window)
    ref = parsed["u10"]
    iy = np.where((ref["lats"] >= miny - GFS_RES / 2)
                  & (ref["lats"] <= maxy + GFS_RES / 2))[0]
    ix = np.where((ref["lons"] >= minx - GFS_RES / 2)
                  & (ref["lons"] <= maxx + GFS_RES / 2))[0]
    if iy.size == 0 or ix.size == 0:  # pragma: no cover - guarded by _verify_coverage
        raise RuntimeError(
            f"requested window {list(window)} has no grid points in the "
            "returned GFS grid")
    lats = ref["lats"][iy[0]:iy[-1] + 1]
    lons = ref["lons"][ix[0]:ix[-1] + 1]
    out = {}
    for key, part in parsed.items():
        if part["lats"].shape != ref["lats"].shape or \
                part["lons"].shape != ref["lons"].shape:
            raise RuntimeError(
                f"GFS variables have mismatched grids in one payload "
                f"({key}: {part['lats'].shape}/{part['lons'].shape} vs "
                f"u10: {ref['lats'].shape}/{ref['lons'].shape})")
        out[key] = part["values"][iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
    return lats, lons, out


# ---------------------------------------------------------------------------
# public fetch
# ---------------------------------------------------------------------------


def fetch_gfs_wind(bbox: Sequence[float],
                   start: DateLike, end: DateLike,
                   stride_days: int = 1,
                   cycle: str = "00",
                   forecast_hours: Sequence[int] = (0,),
                   work_dir: Optional[str] = None) -> GfsWindField:
    """Fetch NOAA GFS 10-m winds + 2-m air temperature (keyless, NOMADS).

    One snapshot per ``(sampled day, forecast hour)`` step, from the
    0.25° GFS via the NOMADS GRIB filter — no account, no API key.

    Args:
        bbox: ``(min_lon, min_lat, max_lon, max_lat)`` in -180..180
            degrees. Antimeridian-crossing boxes (``max_lon < min_lon``)
            are fetched as two subregion requests and concatenated.
        start, end: dates/datetimes/ISO strings within the last
            :data:`GFS_RETENTION_DAYS` days (NOMADS retention).
        stride_days: day sampling stride, int >= 1 (default 1 = every day).
        cycle: GFS analysis cycle — ``"00"``, ``"06"``, ``"12"`` or
            ``"18"`` (default ``"00"``).
        forecast_hours: tuple/list of forecast hours to fetch per
            sampled day, ints 0..120 (default ``(0,)`` = the f000
            analysis only — backwards compatible). ``(0, 1, 6)`` gives
            a sub-daily series; hours are normalized to sorted order so
            the time series is chronological. Each hour costs one
            ~2.4 MB download per window for a North-America box, so
            ``range(121)`` is ~290 MB/day — the caller owns the frame
            budget.
        work_dir: cache directory for the downloaded GRIB payloads
            (default ``$SURVEY_CURRENTS_CACHE/gfs-wind``).

    Returns:
        :class:`GfsWindField` with ``u10``/``v10`` in m/s and ``t2m``
        in °C, one timestep per ``(day, forecast_hour)`` in
        chronological order; ``times`` carry the forecast hour (valid
        time from the GRIB message itself). Provenance records a
        per-timestep list: exact request URLs, per-file SHA-256, byte
        counts, retrieval time.

    Raises:
        ImportError: ``cfgrib`` is not installed.
        UnavailableRangeError: a date is outside the NOMADS retention
            window, or NOMADS returned HTTP 404 for a (day, hour)
            request — never a silent skip or silent padding.
        ValueError: invalid bbox / stride / cycle / date order /
            forecast hour.
        RuntimeError: download or GRIB parse failures.
    """
    minx, miny, maxx, maxy = validate_gfs_bbox(bbox)
    d0 = _coerce_datetime_utc(start)
    d1 = _coerce_datetime_utc(end)
    _validate_gfs_dates(d0, d1)
    c = _validate_gfs_cycle(cycle)
    hours = _validate_forecast_hours(forecast_hours)
    if not isinstance(stride_days, int) or isinstance(stride_days, bool) \
            or stride_days < 1:
        raise ValueError(f"stride_days must be an int >= 1, got {stride_days!r}")

    plan = gfs_sample_plan(d0, d1, stride_days)
    windows = gfs_lon_windows((minx, miny, maxx, maxy))
    cache = gfs_cache_dir(work_dir)

    step_grids: List[Dict[str, np.ndarray]] = []
    step_times: List[str] = []
    step_lats: Optional[np.ndarray] = None
    step_lons: Optional[np.ndarray] = None
    requests: List[Dict[str, Any]] = []
    n_payload_bytes = 0
    combined = hashlib.sha256()
    n_cached = 0
    retrieved_at = _dt.datetime.now(_dt.timezone.utc).isoformat()

    for day in plan:
        for hour in hours:
            win_parts: List[Tuple[np.ndarray, np.ndarray,
                                  Dict[str, np.ndarray]]] = []
            valid_time: Optional[str] = None
            for window in windows:
                url = gfs_filter_url(day, c, window, forecast_hour=hour)
                payload, from_cache = _cached_get_bytes(url, cache)
                if from_cache:
                    n_cached += 1
                n_payload_bytes += len(payload)
                combined.update(hashlib.sha256(payload).digest())
                requests.append({
                    "url": url,
                    "date": day.isoformat(),
                    "cycle": c,
                    "forecast_hour": hour,
                    "window": list(window),
                    "sha256": hashlib.sha256(payload).hexdigest(),
                    "n_bytes": len(payload),
                    "from_cache": from_cache,
                })
                parsed = _parse_grib_payload(payload, url, day, c,
                                             forecast_hour=hour)
                _verify_coverage(parsed, window, url)
                if valid_time is None:
                    valid_time = parsed["u10"]["valid_time"]
                win_parts.append(_crop_to_window(parsed, window))
            # Concatenate antimeridian windows along longitude (each part
            # is already -180..180 increasing; sort + dedupe the seam
            # meridian).
            lats = win_parts[0][0]
            lons = np.concatenate([p[1] for p in win_parts])
            order = np.argsort(lons, kind="stable")
            lons_sorted = lons[order]
            _, unique_idx = np.unique(lons_sorted, return_index=True)
            keep = order[np.sort(unique_idx)]
            merged = {key: np.concatenate([p[2][key] for p in win_parts],
                                          axis=1)[:, keep]
                      for key in ("u10", "v10", "t2m")}
            step_lats, step_lons = lats, lons[keep]
            step_grids.append(merged)
            assert valid_time is not None  # parser records it on every part
            step_times.append(valid_time)

    nt = len(step_times)
    lats = np.asarray(step_lats, dtype=float)
    lons = np.asarray(step_lons, dtype=float)
    grids = {
        "u10": np.ma.array(np.stack([g["u10"] for g in step_grids]),
                           mask=False),
        "v10": np.ma.array(np.stack([g["v10"] for g in step_grids]),
                           mask=False),
        # 2t is Kelvin on the wire -> °C here.
        "t2m": np.ma.array(np.stack([g["t2m"] for g in step_grids])
                           - _KELVIN_OFFSET, mask=False),
    }
    times = step_times

    return GfsWindField(
        grids=grids, times=times, lats=lats, lons=lons, cycle=c,
        source="noaa-nomads/gfs-0p25",
        provenance={
            "dataset": "NOAA GFS 0.25-degree analyses + hourly forecast "
                       "steps (f000..f120) via the NCEP NOMADS GRIB filter "
                       "(keyless HTTPS)",
            "nomads_filter": NOMADS_FILTER_BASE,
            "requests": requests,
            "n_requests": len(requests),
            "n_cached": n_cached,
            "sha256": combined.hexdigest(),
            "n_bytes": n_payload_bytes,
            "retrieved_at": retrieved_at,
            "bbox_requested": [minx, miny, maxx, maxy],
            "lon_windows": [list(w) for w in windows],
            "time_requested": [d0.date().isoformat(), d1.date().isoformat()],
            "cycle": c,
            "forecast_hours": list(hours),
            "forecast_hours_note": (
                "f000 = analysis; f001..f120 = hourly forecast steps. "
                "Each step is one NOMADS GRIB-filter request (~2.4 MB "
                "for a North-America box at 0.25°); a full 120-step day "
                "is ~290 MB — the caller owns the frame budget."),
            "stride_days": stride_days,
            "retention_days": GFS_RETENTION_DAYS,
            "retention_note": (
                f"NOMADS keeps roughly the last {GFS_RETENTION_DAYS} days "
                "of the 0.25° GFS (verified 2026-10-01); older dates "
                "raise UnavailableRangeError."),
            "subset_note": (
                "Requested with subregion=on; the returned grid was "
                "verified to cover the bbox and cropped locally with "
                "numpy (the filter has silently ignored subregions "
                "without subregion=on)."),
            "units": {spec["grid"]: spec["units"]
                      for spec in GFS_GRIDS.values()},
            "unit_notes": {
                "t2m": "Kelvin on the wire, converted to degC",
                "air_temperature": "2-m air temperature in degF "
                                   "(warming.watch strand-color convention)",
            },
            "grid": f"{GFS_RES}-degree global (per-step analysis/forecast "
                     "snapshots)",
            "access": "keyless",
        },
    )


# ---------------------------------------------------------------------------
# offline demo
# ---------------------------------------------------------------------------


def main_demo() -> None:
    """Print a small offline summary (no network). Mirrors era5.main_demo."""
    f = GfsWindField.synthetic(nt=3)
    print(f"[gfs-wind] synthetic grids "
          f"{ {k: v.shape for k, v in f.grids.items()} }, "
          f"t0={f.times[0]}, mean wind={f.spatial_mean('wind', 0):.2f} m/s, "
          f"mean 2m air temp={f.spatial_mean('air_temperature', 0):.1f} F")


__all__ = [
    "DateLike",
    "GfsWindField",
    "GFS_CYCLES",
    "GFS_FORECAST_HOUR",
    "GFS_FORECAST_HOUR_MAX",
    "GFS_GRIDS",
    "GFS_RES",
    "GFS_RETENTION_DAYS",
    "NOMADS_FILTER_BASE",
    "UnavailableRangeError",
    "fetch_gfs_wind",
    "gfs_cache_dir",
    "gfs_filter_url",
    "gfs_lon_windows",
    "gfs_retention_window",
    "gfs_sample_plan",
    "validate_gfs_bbox",
    "main_demo",
]

"""Canonical data model: the OFS registry and the CurrentField.

Everything in survey-currents converges on :class:`CurrentField` — a
regular-grid, time-ordered stack of surface-current vectors (u/v, m/s)
and water temperature (degC). It is deliberately boring: plain numpy
arrays plus metadata, so survey-flow (particle advection), survey-animate
(frames), and the survey-suite engines can consume it without knowing
which model or service produced it.
"""

from __future__ import annotations

import datetime as _dt
import json
import math
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# NOAA Operational Forecast System registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OfsModel:
    """One NOAA Operational Forecast System entry."""

    code: str            # e.g. "LMHOFS" — used in S3 path construction
    name: str            # human-readable name
    region: str          # geographic coverage
    resolution: str      # nominal grid resolution
    horizon_hours: int   # forecast horizon in hours
    cycles: Tuple[str, ...] = ("00", "06", "12", "18")
    notes: str = ""


#: Registry of NOAA OFS models served from the public S3 buckets.
#: Resolution / horizon values follow the NOAA OFS open-data docs
#: (verified 2026-09-25); see docs/DATA_SOURCES.md.
OFS_REGISTRY: Dict[str, OfsModel] = {
    m.code: m
    for m in (
        OfsModel("GLOFS", "Great Lakes Operational Forecast System",
                 "All five Great Lakes", "5 km", 60,
                 notes="Basin-wide GLOFS grid; use LMHOFS/LEOFS for finer detail."),
        OfsModel("LMHOFS", "Lake Michigan-Huron Operational Forecast System",
                 "Lakes Michigan and Huron", "50 m - 2.5 km", 120,
                 notes="Finest public grid for Lake Michigan currents."),
        OfsModel("LEOFS", "Lake Erie Operational Forecast System",
                 "Lake Erie", "400 m - 4 km", 120),
        OfsModel("CBOFS", "Chesapeake Bay Operational Forecast System",
                 "Chesapeake Bay", "50 m - 3 km", 48),
        OfsModel("CIOFS", "Cook Inlet Operational Forecast System",
                 "Cook Inlet, Alaska", "10 m - 3.5 km", 48),
        OfsModel("CREOFS", "Columbia River Estuary Operational Forecast System",
                 "Columbia River Estuary", "100 m - 4 km", 48),
        OfsModel("DBOFS", "Delaware Bay Operational Forecast System",
                 "Delaware Bay", "100 m - 3 km", 48),
        OfsModel("GoMOFS", "Gulf of Maine Operational Forecast System",
                 "Gulf of Maine", "700 m", 72),
        OfsModel("NGOFS2", "Northern Gulf of Mexico Operational Forecast System",
                 "Northern Gulf of Mexico", "45 m - 300 m", 48),
        OfsModel("SFBOFS", "San Francisco Bay Operational Forecast System",
                 "San Francisco Bay", "100 m - 4 km", 48),
        OfsModel("TBOFS", "Tampa Bay Operational Forecast System",
                 "Tampa Bay", "100 m - 1.2 km", 48),
        OfsModel("WCOFS", "West Coast Operational Forecast System",
                 "US West Coast", "4 km", 72),
        OfsModel("SSCOFS", "Salish Sea and Columbia River Operational Forecast System",
                 "Salish Sea and Columbia River", "100 m - 3 km", 72),
    )
}


def get_ofs_model(code: str) -> OfsModel:
    """Return the registry entry for ``code`` (case-insensitive)."""
    try:
        return OFS_REGISTRY[code.strip().upper()]
    except KeyError as exc:
        known = ", ".join(sorted(OFS_REGISTRY))
        raise KeyError(f"unknown OFS code {code!r}; known codes: {known}") from exc


def list_ofs_models() -> List[OfsModel]:
    """All registered OFS models, sorted by code."""
    return [OFS_REGISTRY[k] for k in sorted(OFS_REGISTRY)]


# ---------------------------------------------------------------------------
# CurrentField
# ---------------------------------------------------------------------------


@dataclass
class CurrentField:
    """Surface currents + water temperature on a regular grid over time.

    Arrays are ``(nt, ny, nx)`` with ``nt == len(times)``. ``lats``/``lons``
    are 1-D coordinate vectors (increasing). Units: u/v in m/s (east/north
    positive), temperature in degrees Celsius. Missing values are NaN.
    """

    u: np.ndarray                       # (nt, ny, nx) eastward velocity, m/s
    v: np.ndarray                       # (nt, ny, nx) northward velocity, m/s
    temperature: Optional[np.ndarray]   # (nt, ny, nx) water temperature, degC
    times: List[str]                    # ISO-8601 timestamps, one per step
    lats: np.ndarray                    # (ny,) degrees north, increasing
    lons: np.ndarray                    # (nx,) degrees east, increasing
    crs: str = "EPSG:4326"
    source: str = ""                    # e.g. "noaa-ofs/LMHOFS"
    model_run: str = ""                 # ISO timestamp of the model cycle
    forecast_hours: List[int] = field(default_factory=list)
    provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        nt = len(self.times)
        for name, arr in (("u", self.u), ("v", self.v)):
            arr = np.asarray(arr, dtype=float)
            if arr.shape != (nt, self.lats.shape[0], self.lons.shape[0]):
                raise ValueError(
                    f"{name} shape {arr.shape} does not match "
                    f"(nt, ny, nx)=({nt}, {self.lats.shape[0]}, {self.lons.shape[0]})"
                )
            object.__setattr__(self, name, arr)
        if self.temperature is not None:
            t = np.asarray(self.temperature, dtype=float)
            if t.shape != self.u.shape:
                raise ValueError(
                    f"temperature shape {t.shape} does not match u/v shape {self.u.shape}"
                )
            self.temperature = t
        self.lats = np.asarray(self.lats, dtype=float)
        self.lons = np.asarray(self.lons, dtype=float)
        if self.forecast_hours and len(self.forecast_hours) != nt:
            raise ValueError("forecast_hours length must match times length")
        self.bounds = self._bounds()

    def _bounds(self) -> Tuple[float, float, float, float]:
        return (float(self.lons[0]), float(self.lats[0]),
                float(self.lons[-1]), float(self.lats[-1]))

    # -- derived quantities ------------------------------------------------

    def speed(self) -> np.ndarray:
        """Current speed magnitude in m/s, shape (nt, ny, nx)."""
        with np.errstate(invalid="ignore"):
            return np.sqrt(self.u ** 2 + self.v ** 2)

    def direction_deg(self) -> np.ndarray:
        """Current direction in degrees clockwise from north, (nt, ny, nx)."""
        with np.errstate(invalid="ignore"):
            return (90.0 - np.degrees(np.arctan2(self.v, self.u))) % 360.0

    # -- selection ----------------------------------------------------------

    def select_time(self, index: int) -> "CurrentField":
        """Return the single-timestep field at ``index``."""
        temp = None if self.temperature is None else self.temperature[index:index + 1]
        hours = [self.forecast_hours[index]] if self.forecast_hours else []
        return CurrentField(
            u=self.u[index:index + 1], v=self.v[index:index + 1],
            temperature=temp, times=[self.times[index]],
            lats=self.lats, lons=self.lons, crs=self.crs,
            source=self.source, model_run=self.model_run,
            forecast_hours=hours, provenance=dict(self.provenance),
        )

    def select_bbox(self, bbox: Sequence[float]) -> "CurrentField":
        """Spatial subset to ``(min_lon, min_lat, max_lon, max_lat)``."""
        minx, miny, maxx, maxy = (float(x) for x in bbox)
        ix = np.where((self.lons >= minx) & (self.lons <= maxx))[0]
        iy = np.where((self.lats >= miny) & (self.lats <= maxy))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(f"bbox {list(bbox)} does not overlap field bounds {self.bounds}")
        temp = None if self.temperature is None else self.temperature[:, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
        return CurrentField(
            u=self.u[:, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1],
            v=self.v[:, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1],
            temperature=temp, times=list(self.times),
            lats=self.lats[iy[0]:iy[-1] + 1], lons=self.lons[ix[0]:ix[-1] + 1],
            crs=self.crs, source=self.source, model_run=self.model_run,
            forecast_hours=list(self.forecast_hours),
            provenance=dict(self.provenance),
        )

    def zonal_mean(self, index: int, stat: str = "speed") -> float:
        """NaN-aware spatial mean of ``stat`` at timestep ``index``.

        ``stat`` is one of ``"speed"``, ``"u"``, ``"v"``,
        ``"temperature"``, ``"direction"``.
        """
        arr = {
            "speed": self.speed(), "u": self.u, "v": self.v,
            "direction": self.direction_deg(),
        }.get(stat)
        if arr is None and stat == "temperature":
            if self.temperature is None:
                raise ValueError("field has no temperature grid")
            arr = self.temperature
        if arr is None:
            raise ValueError(f"unknown stat {stat!r}; want speed/u/v/temperature/direction")
        with np.errstate(invalid="ignore"):
            mean = float(np.nanmean(arr[index]))
        return mean

    # -- serialization -------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dict (arrays become nested lists)."""
        return {
            "u": self.u.tolist(),
            "v": self.v.tolist(),
            "temperature": None if self.temperature is None else self.temperature.tolist(),
            "times": list(self.times),
            "lats": self.lats.tolist(),
            "lons": self.lons.tolist(),
            "crs": self.crs,
            "source": self.source,
            "model_run": self.model_run,
            "forecast_hours": list(self.forecast_hours),
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CurrentField":
        """Rebuild from :meth:`to_dict`. Raises on missing keys."""
        required = ("u", "v", "times", "lats", "lons")
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(f"CurrentField dict missing keys: {missing}")
        return cls(
            u=np.asarray(data["u"], dtype=float),
            v=np.asarray(data["v"], dtype=float),
            temperature=None if data.get("temperature") is None
            else np.asarray(data["temperature"], dtype=float),
            times=list(data["times"]),
            lats=np.asarray(data["lats"], dtype=float),
            lons=np.asarray(data["lons"], dtype=float),
            crs=data.get("crs", "EPSG:4326"),
            source=data.get("source", ""),
            model_run=data.get("model_run", ""),
            forecast_hours=list(data.get("forecast_hours", [])),
            provenance=dict(data.get("provenance", {})),
        )

    def to_json(self, path: str) -> str:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh)
        return path

    @classmethod
    def from_json(cls, path: str) -> "CurrentField":
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    def to_netcdf(self, path: str) -> str:
        """Write CF-ish NetCDF (needs the ``noaa`` extra: xarray/netCDF4)."""
        try:
            import xarray as xr
        except ImportError as exc:
            raise ImportError(
                "to_netcdf needs xarray/netCDF4 (pip install 'survey-currents[noaa]')"
            ) from exc
        coords = {
            "time": ("time", self.times),
            "lat": ("lat", self.lats, {"units": "degrees_north"}),
            "lon": ("lon", self.lons, {"units": "degrees_east"}),
        }
        data_vars = {
            "u": (("time", "lat", "lon"), self.u,
                  {"units": "m s-1", "long_name": "eastward surface current"}),
            "v": (("time", "lat", "lon"), self.v,
                  {"units": "m s-1", "long_name": "northward surface current"}),
        }
        if self.temperature is not None:
            data_vars["temperature"] = (
                ("time", "lat", "lon"), self.temperature,
                {"units": "degC", "long_name": "surface water temperature"})
        ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs={
            "source": self.source, "model_run": self.model_run,
            "crs": self.crs, "created": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        })
        ds.to_netcdf(path)
        return path

    # -- synthetic fixture (tests, demos, offline use) ------------------------

    @classmethod
    def synthetic(
        cls,
        nt: int = 4, ny: int = 6, nx: int = 8,
        lats: Sequence[float] = (41.5, 46.5),
        lons: Sequence[float] = (-92.5, -84.5),
        start: str = "2026-09-16T00:00:00+00:00",
        step_hours: int = 1,
        seed: int = 7,
        source: str = "synthetic",
    ) -> "CurrentField":
        """Deterministic synthetic field: rotating gyre + warm core.

        Stdlib+numpy only. Used by the test suite and the offline demo.
        """
        rng = np.random.default_rng(seed)
        lat = np.linspace(lats[0], lats[1], ny)
        lon = np.linspace(lons[0], lons[1], nx)
        yy, xx = np.meshgrid(np.linspace(-1, 1, ny), np.linspace(-1, 1, nx),
                             indexing="ij")
        base_u = -yy * 0.35   # clockwise gyre, m/s
        base_v = xx * 0.35
        base_t = 18.0 - 4.0 * (xx ** 2 + yy ** 2)  # warm core, degC
        u = np.empty((nt, ny, nx)); v = np.empty((nt, ny, nx)); t = np.empty((nt, ny, nx))
        t0 = _dt.datetime.fromisoformat(start)
        times = []
        for k in range(nt):
            jitter = rng.normal(0.0, 0.02, size=(ny, nx))
            u[k] = base_u + jitter
            v[k] = base_v + jitter * 0.5
            t[k] = base_t + rng.normal(0.0, 0.1, size=(ny, nx)) + 0.2 * k
            times.append((t0 + _dt.timedelta(hours=k * step_hours)).isoformat())
        return cls(u=u, v=v, temperature=t, times=times, lats=lat, lons=lon,
                   source=source, model_run=start,
                   forecast_hours=[k * step_hours for k in range(nt)])

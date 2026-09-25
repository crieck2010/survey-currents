"""Interoperability with the survey-suite / earthwatch engines.

* :class:`CurrentsPassProvider` implements survey-monitor's
  ``PassProvider`` interface (``list_passes`` / ``metrics``) without
  importing survey-monitor — the pipeline duck-types providers, and the
  import stays lazy/optional exactly like ``ImageryPassProvider`` does.
  Each forecast hour becomes one "pass"; metrics are NaN-aware spatial
  means over the site AOI (``speed_mean``, ``u_mean``, ``v_mean``,
  ``temp_mean`` when the field carries temperature).

* :func:`align_to_thermal_zone` clips a :class:`CurrentField` to a
  survey-thermal zone dict and returns per-timestep mean water
  temperature, so model SST can be compared against Landsat LST passes
  from ``thermal.acquire.acquire_lst_passes``.

Neither helper imports the sibling package; both accept plain dicts /
fields, which keeps this engine dependency-light and testable offline.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from .models import CurrentField


def _pass_info_cls():
    """survey-monitor's PassInfo when installed, else a local twin."""
    try:
        from monitor.scheduler import PassInfo  # type: ignore
        return PassInfo
    except ImportError:
        @dataclass
        class _LocalPassInfo:
            pass_id: str
            date: str
            cloud_cover: Optional[float] = None
            source: str = ""
        return _LocalPassInfo


def _bbox_of_geometry(geometry: Dict[str, Any]) -> Sequence[float]:
    """Accept a GeoJSON geometry dict or ``{"bbox": [...]}`` (monitor style)."""
    if "bbox" in geometry:
        return [float(x) for x in geometry["bbox"]]
    coords = geometry.get("coordinates")
    if coords is None:
        raise ValueError("geometry has neither 'bbox' nor 'coordinates'")
    flat = []

    def walk(c):
        if isinstance(c[0], (int, float)):
            flat.append((float(c[0]), float(c[1])))
        else:
            for sub in c:
                walk(sub)

    walk(coords)
    xs = [p[0] for p in flat]
    ys = [p[1] for p in flat]
    return [min(xs), min(ys), max(xs), max(ys)]


class CurrentsPassProvider:
    """survey-monitor pass provider over a CurrentField.

    Parameters
    ----------
    field:
        The current field (already fetched, e.g. via
        :func:`currents.noaa_ofs.fetch_currents`).
    site_id:
        Monitor site id, used in pass ids.
    geometry:
        GeoJSON geometry dict or ``{"bbox": [minx, miny, maxx, maxy]}``
        in the monitor's ``SiteConfig.geometry`` shape.
    source:
        Provenance label recorded on each pass.
    """

    def __init__(self, field: CurrentField, site_id: str,
                 geometry: Dict[str, Any], source: str = "") -> None:
        self.field = field
        self.site_id = site_id
        self.bbox = _bbox_of_geometry(geometry)
        self.source = source or field.source
        clipped = field.select_bbox(self.bbox)
        self._clipped = clipped

    # -- PassProvider interface (duck-typed) -------------------------------

    def list_passes(self, since: Optional[str] = None) -> List[Any]:
        """One pass per forecast hour, oldest first.

        Returns survey-monitor ``PassInfo`` objects when survey-monitor is
        installed, otherwise a local equivalent with the same fields.
        """
        PassInfo = _pass_info_cls()
        out = []
        for k, iso in enumerate(self._clipped.times):
            date = iso[:10]
            if since is not None and date <= since:
                continue
            hour = (self._clipped.forecast_hours[k]
                    if self._clipped.forecast_hours else k)
            out.append(PassInfo(
                pass_id=f"currents-{self.site_id}-f{hour:03d}-{date}",
                date=date, cloud_cover=0.0, source=self.source))
        return out

    def metrics(self, pass_id: str) -> Dict[str, float]:
        """NaN-aware spatial means for one pass."""
        try:
            hour = int(pass_id.split("-f")[1].split("-")[0])
        except (IndexError, ValueError) as exc:
            raise KeyError(f"unknown currents pass_id {pass_id!r}") from exc
        hours = self._clipped.forecast_hours or list(range(len(self._clipped.times)))
        if hour not in hours:
            raise KeyError(f"unknown currents pass_id {pass_id!r}")
        k = hours.index(hour)
        out = {
            "speed_mean": self._clipped.zonal_mean(k, "speed"),
            "u_mean": self._clipped.zonal_mean(k, "u"),
            "v_mean": self._clipped.zonal_mean(k, "v"),
        }
        if self._clipped.temperature is not None:
            out["temp_mean"] = self._clipped.zonal_mean(k, "temperature")
        return out


def align_to_thermal_zone(field: CurrentField,
                          zone: Dict[str, Any],
                          start: Optional[str] = None,
                          end: Optional[str] = None) -> Dict[str, Any]:
    """Clip ``field`` to a survey-thermal zone and summarize water temp.

    ``zone`` is a thermal ``Zone`` as a dict (``id`` + ``geometry``);
    ``start``/``end`` are ISO dates bounding the model timesteps.
    Returns ``{"zone_id", "bbox", "timesteps": [...], "temp_mean_series":
    [...]}`` — directly comparable with the per-pass LST CSV rows from
    ``thermal.acquire.acquire_lst_passes``.
    """
    bbox = _bbox_of_geometry(zone.get("geometry", {}))
    clipped = field.select_bbox(bbox)
    steps: List[Dict[str, Any]] = []
    series: List[Optional[float]] = []
    for k, iso in enumerate(clipped.times):
        if start and iso[:10] < start:
            continue
        if end and iso[:10] > end:
            continue
        mean = (clipped.zonal_mean(k, "temperature")
                if clipped.temperature is not None else None)
        steps.append({"time": iso,
                      "forecast_hour": (clipped.forecast_hours[k]
                                       if clipped.forecast_hours else k),
                      "temp_mean_c": mean})
        series.append(mean)
    return {
        "zone_id": zone.get("id", ""),
        "bbox": list(bbox),
        "source": field.source,
        "timesteps": steps,
        "temp_mean_series": series,
    }

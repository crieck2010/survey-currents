"""GeoTIFF / COG export: one file per timestep, four bands.

Band order is the stable interchange contract survey-flow will consume:

1. ``speed`` — current speed magnitude, m/s
2. ``u`` — eastward velocity, m/s
3. ``v`` — northward velocity, m/s
4. ``temperature`` — water temperature, degC (NaN band when the field
   has no temperature grid)

Filenames: ``{prefix}_{YYYYMMDDTHHMM}.tif``. Needs rasterio
(``pip install 'survey-currents[raster]'``).
"""

from __future__ import annotations

import datetime as _dt
import os
from typing import List

import numpy as np

from .models import CurrentField

BAND_NAMES = ("speed", "u", "v", "temperature")


def _require_rasterio():
    try:
        import rasterio
        from rasterio.transform import from_origin
    except ImportError as exc:
        raise ImportError(
            "COG export needs rasterio (pip install 'survey-currents[raster]')"
        ) from exc
    return rasterio, from_origin


def _timestep_stamp(iso: str) -> str:
    try:
        return _dt.datetime.fromisoformat(iso).strftime("%Y%m%dT%H%M")
    except ValueError:
        return "".join(c for c in iso if c.isalnum())[:13]


def export_cogs(field: CurrentField, out_dir: str,
                prefix: str = "currents") -> List[str]:
    """Write one 4-band GeoTIFF per timestep; return the file paths."""
    rasterio, from_origin = _require_rasterio()
    os.makedirs(out_dir, exist_ok=True)
    ny, nx = field.u.shape[1], field.u.shape[2]
    res_x = (field.lons[-1] - field.lons[0]) / max(nx - 1, 1)
    res_y = (field.lats[-1] - field.lats[0]) / max(ny - 1, 1)
    transform = from_origin(float(field.lons[0]) - res_x / 2,
                            float(field.lats[-1]) + res_y / 2,
                            res_x, res_y)
    speed = field.speed()
    paths = []
    for k, iso in enumerate(field.times):
        stack = np.stack([
            speed[k],
            field.u[k],
            field.v[k],
            field.temperature[k] if field.temperature is not None
            else np.full((ny, nx), np.nan),
        ]).astype("float32")
        # rasterio origin is upper-left: flip latitude axis.
        stack = stack[:, ::-1, :]
        path = os.path.join(out_dir, f"{prefix}_{_timestep_stamp(iso)}.tif")
        profile = {
            "driver": "GTiff", "height": ny, "width": nx, "count": 4,
            "dtype": "float32", "crs": field.crs, "transform": transform,
            "nodata": float("nan"), "compress": "deflate", "tiled": True,
        }
        try:
            with rasterio.open(path, "w", **{**profile, "driver": "COG"}) as dst:
                dst.write(stack)
                dst.descriptions = BAND_NAMES
        except Exception:
            with rasterio.open(path, "w", **profile) as dst:
                dst.write(stack)
                dst.descriptions = BAND_NAMES
        paths.append(path)
    return paths

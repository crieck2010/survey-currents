"""survey-currents: surface-current and water-temperature acquisition engine.

Acquire hourly current-vector (u/v) and water-temperature fields from
operational hydrodynamic forecast models — NOAA Operational Forecast
Systems (anonymous AWS S3) and the Copernicus Marine Service — plus
NOAA GLSEA satellite-derived daily sea-surface temperature and
lake-average temperature series, into one canonical
:class:`~currents.models.CurrentField` model (and the SST-focused
:class:`~currents.glsea.GlseaField`) with provenance, COG export, and
survey-suite interoperability.

Engine core is stdlib+numpy. Heavy dependencies (xarray, boto3,
netCDF4, copernicusmarine, rasterio) are lazy optional imports; every
public function that needs one raises an actionable ImportError naming
the ``pip install`` extra.
"""

from __future__ import annotations

__version__ = "0.2.0"

from .glsea import (LAKE_COLUMNS, GlseaField, LakeSeries,
                    fetch_glsea_lake_averages, fetch_glsea_sst,
                    glsea_averages_url, glsea_sst_url, validate_glsea_bbox)
from .models import OFS_REGISTRY, CurrentField, OfsModel, get_ofs_model
from .provenance import PROVENANCE_VERSION, read_provenance, verify_provenance, write_provenance

__all__ = [
    "CurrentField",
    "OfsModel",
    "OFS_REGISTRY",
    "get_ofs_model",
    "GlseaField",
    "LakeSeries",
    "LAKE_COLUMNS",
    "fetch_glsea_sst",
    "fetch_glsea_lake_averages",
    "glsea_sst_url",
    "glsea_averages_url",
    "validate_glsea_bbox",
    "PROVENANCE_VERSION",
    "read_provenance",
    "verify_provenance",
    "write_provenance",
    "__version__",
]

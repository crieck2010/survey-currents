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

__version__ = "0.3.0"

from .glsea import (LAKE_COLUMNS, GlseaField, LakeSeries,
                    fetch_glsea_lake_averages, fetch_glsea_sst,
                    glsea_averages_url, glsea_sst_url, validate_glsea_bbox)
from .models import OFS_REGISTRY, CurrentField, OfsModel, get_ofs_model
from .provenance import PROVENANCE_VERSION, read_provenance, verify_provenance, write_provenance
from .sst_global import (CMR_GRANULE_SEARCH, CredentialsMissing, SstField,
                         cmr_search_mur_granules, earthdata_credentials,
                         fetch_mur, fetch_oisst, mur_granule_title,
                         mur_match_granules, mur_opendap_url, mur_subset_urls,
                         oisst_lon_windows, oisst_sst_urls,
                         validate_sst_bbox)

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
    "SstField",
    "CredentialsMissing",
    "CMR_GRANULE_SEARCH",
    "cmr_search_mur_granules",
    "mur_match_granules",
    "fetch_oisst",
    "fetch_mur",
    "oisst_sst_urls",
    "oisst_lon_windows",
    "mur_granule_title",
    "mur_opendap_url",
    "mur_subset_urls",
    "earthdata_credentials",
    "validate_sst_bbox",
    "PROVENANCE_VERSION",
    "read_provenance",
    "verify_provenance",
    "write_provenance",
    "__version__",
]

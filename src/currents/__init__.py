"""survey-currents: surface-current and water-temperature acquisition engine.

Acquire hourly current-vector (u/v) and water-temperature fields from
operational hydrodynamic forecast models — NOAA Operational Forecast
Systems (anonymous AWS S3) and the Copernicus Marine Service — into one
canonical :class:`~currents.models.CurrentField` model with provenance,
COG export, and survey-suite interoperability.

Engine core is stdlib+numpy. Heavy dependencies (xarray, boto3,
copernicusmarine, rasterio) are lazy optional imports; every public
function that needs one raises an actionable ImportError naming the
``pip install`` extra.
"""

from __future__ import annotations

__version__ = "0.1.0"

from .models import OFS_REGISTRY, CurrentField, OfsModel, get_ofs_model
from .provenance import PROVENANCE_VERSION, read_provenance, verify_provenance, write_provenance

__all__ = [
    "CurrentField",
    "OfsModel",
    "OFS_REGISTRY",
    "get_ofs_model",
    "PROVENANCE_VERSION",
    "read_provenance",
    "verify_provenance",
    "write_provenance",
    "__version__",
]

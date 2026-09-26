# survey-currents

Surface-current and water-temperature acquisition engine for surveying and remote sensing — the first module of the earthwatch-suite **flow-field animation program** (the pipeline that produces mapped.earth-style animated current/temperature reels).

`survey-currents` pulls hourly current-vector (u/v) and water-temperature fields from **operational hydrodynamic forecast models** — not satellites — into one canonical `CurrentField` model with provenance, COG export, and survey-suite interoperability. Since v0.2.0 it also pulls **satellite-derived daily sea-surface temperature** from NOAA's GLSEA analysis (Great Lakes only) into a companion `GlseaField` model, plus lake-wide average temperature series. Since v0.3.0 it pulls **global satellite SST** from two more sources into a companion `SstField` model: NOAA OISST v2.1 (keyless, 0.25°, 1981–present) and NASA JPL MUR v4.1 (free Earthdata Login, ~1 km, 2002–present).

## Data sources

| Source | Access | Coverage | Resolution | Horizon |
|---|---|---|---|---|
| NOAA OFS via anonymous AWS S3 (`noaa-ofs-pds`, `noaa-nos-ofs-pds`) | No signup, unsigned requests | US coasts + Great Lakes | 50 m – 5 km | 48–120 h |
| Copernicus Marine Service (`copernicusmarine` toolbox) | Free account | Global ocean | 1/12° (~9 km) | NRT + forecast |
| NOAA GLSEA via ERDDAP griddap (`GLSEA_ACSPO_GCS`) | No signup | **Great Lakes only** | ~1.5 km | daily analysis, 2006–present |
| NOAA GLSEA via ERDDAP tabledap (`glsea_avgtemps_3`) | No signup | Great Lakes (per-lake daily averages) | lake-wide | daily, 2006–present |
| NOAA OISST v2.1 via CoastWatch ERDDAP (`ncdcOisst21Agg`) | No signup | **Global ocean** | 0.25° (~28 km) | daily analysis, 1981–present |
| NASA JPL MUR v4.1 via Earthdata OPeNDAP | Free Earthdata Login | **Global ocean** | ~0.01° (~1 km) | daily analysis, 2002–present |

Registered NOAA models include **GLOFS** (Great Lakes, 5 km, 60 h), **LMHOFS** (Lake Michigan/Huron, 50 m–2.5 km, 120 h — the source of the Lake Michigan reel), LEOFS, CBOFS, DBOFS, GoMOFS, WCOFS, NGOFS2, SFBOFS, TBOFS, CIOFS, CREOFS, SSCOFS. Full table in `docs/DATA_SOURCES.md`.

## Install

```bash
pip install survey-currents                  # engine core (stdlib + numpy)
pip install "survey-currents[noaa]"          # NetCDF parsing (xarray/netCDF4)
pip install "survey-currents[s3]"            # faster S3 via boto3 (optional)
pip install "survey-currents[cmems]"         # Copernicus Marine toolbox
pip install "survey-currents[raster]"        # GeoTIFF/COG export (rasterio)
pip install "survey-currents[full]"          # everything
```

## Quickstart

```python
from currents.noaa_ofs import fetch_currents

# One week of Lake Michigan surface currents, hour by hour —
# the same LMHOFS feed behind the "Lake Michigan never sits still" reel.
field = fetch_currents(
    "LMHOFS", date="2026-09-16", cycle="00",
    hours=range(0, 24),                       # forecast hours
    bbox=(-92.5, 41.5, -84.5, 46.5),          # min_lon, min_lat, max_lon, max_lat
    work_dir="lmhofs_week",
)
print(field.u.shape)        # (24, ny, nx) eastward velocity, m/s
print(field.zonal_mean(0, "speed"))   # mean current speed at hour 0

from currents.convert import export_cogs
paths = export_cogs(field, "cogs/")   # one 4-band GeoTIFF per timestep:
                                      # bands = speed, u, v, temperature
```

No network? The engine ships a deterministic synthetic field:

```python
from currents.models import CurrentField
field = CurrentField.synthetic(nt=24)   # rotating gyre + warm core, stdlib+numpy
```

CMEMS (needs a free account — see `docs/DATA_SOURCES.md`):

```python
from currents.cmems import subset_cmems, parse_cmems_netcdf
nc = subset_cmems("global-physics-daily",
                  bbox=(-10, 35, 40, 60),
                  start="2026-09-16T00:00:00", end="2026-09-18T00:00:00",
                  work_dir="cmems_data")
field = parse_cmems_netcdf(nc)
```

GLSEA satellite SST (Great Lakes only — no signup, needs `netCDF4`):

```python
from currents.glsea import fetch_glsea_sst, fetch_glsea_lake_averages

# Daily SST grids over western Lake Superior, every 30 days in 2025.
# NOTE: the GLSEA grid is clipped to the lakes region — its longitude
# floor is exactly -92.4199507342304; requests west of that raise
# ValueError before any download.
sst = fetch_glsea_sst(
    bbox=(-92.0, 46.5, -87.0, 48.0),   # min_lon, min_lat, max_lon, max_lat
    start="2025-01-01", end="2025-12-31", stride_days=30,
)
print(sst.sst.shape)          # (nt, ny, nx) masked array, degC
print(sst.spatial_mean(0))    # masked/NaN-aware mean SST at step 0

# Lake-wide daily average temperature series (stdlib only, no netCDF4).
series = fetch_glsea_lake_averages("superior", "2025-01-01", "2025-12-31")
print(series.n, f"{series.mean():.2f} degC mean")
```

## CLI

```bash
survey-currents list-models
survey-currents fetch-noaa --ofs LMHOFS --date 2026-09-16 --cycle 00 \
    --hours 0-23 --bbox -92.5,41.5,-84.5,46.5 --out lmhofs_data
survey-currents fetch-cmems --preset global-physics-daily \
    --bbox -10,35,40,60 --start 2026-09-16T00:00:00 --end 2026-09-18T00:00:00
survey-currents info lmhofs_data/nos.lmhofs.fields.f000.20260916.t00z.nc
survey-currents export-cogs field.json --out cogs/
survey-currents synthetic --nt 24 --out demo_field
survey-currents fetch-glsea-sst --bbox -92.0,46.5,-87.0,48.0 \
    --start 2025-01-01 --end 2025-12-31 --stride-days 30 --out superior_sst
survey-currents fetch-glsea-averages --lake superior \
    --start 2025-01-01 --end 2025-12-31 --out superior_avg
survey-currents glsea-synthetic --nt 4 --out glsea_demo
survey-currents fetch-oisst --bbox -80,20,-60,40 \
    --start 2020-01-01 --end 2020-12-31 --stride-days 30 --out atlantic_sst
survey-currents fetch-mur --bbox -80,20,-60,40 \
    --start 2024-01-01 --end 2024-01-31 --stride-days 7 --out atlantic_mur
survey-currents sst-synthetic --nt 4 --out sst_demo
```

## The canonical model

Everything converges on `CurrentField` (`src/currents/models.py`):

- `u`, `v`: (nt, ny, nx) velocity grids, m/s — east/north positive
- `temperature`: (nt, ny, nx) water temperature, °C (optional)
- `times`: ISO-8601 timestamps; `lats`/`lons`: increasing coordinate vectors
- `speed()`, `direction_deg()`, `select_time()`, `select_bbox()`, `zonal_mean()`
- `to_dict()`/`from_dict()` JSON round-trip, `to_netcdf()`, `to_json()`/`from_json()`
- Every download writes a `.provenance.json` sidecar (source URL, bbox, time window, SHA-256)

COG band order is the stable interchange contract for **survey-flow** (next module): `speed, u, v, temperature`.

Satellite SST converges on `GlseaField` (`src/currents/glsea.py`):

- `sst`: (nt, ny, nx) numpy masked array, °C — land/missing cells masked
- `times`: ISO-8601 timestamps (daily, 12:00 UTC); `lats`/`lons`: increasing vectors
- `spatial_mean()`, `select_time()`, `select_bbox()`, JSON round-trip, `GlseaField.synthetic()`
- `LakeSeries`: lake-wide daily average temps (dates + temps lists, °C), parsed with stdlib `csv` — no netCDF4 needed
- Provenance dict carries the exact ERDDAP URL, SHA-256 of the downloaded bytes, and retrieval timestamp

## Interoperability

- **survey-monitor**: `CurrentsPassProvider` in `currents/interop.py` implements the `PassProvider` interface (`list_passes`/`metrics`) — each forecast hour becomes a monitored pass with `speed_mean`/`u_mean`/`v_mean`/`temp_mean` metrics.
- **survey-thermal**: `align_to_thermal_zone(field, zone, start, end)` clips a field to a thermal zone and returns a per-timestep mean water-temperature series comparable with Landsat LST passes.
- **survey-sites**: proposed `currents:` site-config block is documented in `docs/INTEROP.md` (no changes made to that repo).
- **survey-flow** (next): consumes `CurrentField` + the 4-band COGs.

## Layout

```
src/currents/
    models.py      # CurrentField + NOAA OFS registry
    noaa_ofs.py    # anonymous S3 listing/download, NetCDF parsing
    glsea.py       # NOAA GLSEA satellite SST (griddap) + lake averages (tabledap)
    cmems.py       # copernicusmarine subset wrapper + parser
    convert.py     # per-timestep 4-band GeoTIFF/COG export
    provenance.py  # SHA-256 provenance sidecars
    interop.py     # survey-monitor provider, thermal zone alignment
    cli.py         # thin CLI adapter
docs/
    ARCHITECTURE.md  DATA_SOURCES.md  INTEROP.md
examples/
    glsea_demo.py  # offline GLSEA demo (GlseaField.synthetic only)
```

## Limitations

- **GLSEA is Great Lakes only.** The `GLSEA_ACSPO_GCS` longitude axis is clipped to the lakes region with a minimum of exactly **-92.4199507342304** — bboxes west of that floor (or outside the lat/lon ranges in `docs/DATA_SOURCES.md`) raise a clear `ValueError` before any download. Coverage is daily 2006–present; timesteps are stamped 12:00 UTC.
- **Regular grids only (v0.1.0).** The NetCDF parser handles regular lat/lon grids (e.g. GLOFS). Unstructured FVCOM triangular meshes (native LMHOFS/LEOFS output) raise a clear error — regridding is planned. Many NOAA S3 holdings are regridded; check `info` output.
- **Model data, not observations.** OFS/CMEMS fields are hydrodynamic model output (assimilated, but still model). Treat as guidance; validate against in-situ or satellite SST where it matters.
- **Filename conventions vary per OFS.** The engine lists the date prefix and filters keys by OFS-code/date/cycle/hour tokens rather than assuming exact names.
- Landsat-overpass-style diurnal caveats from survey-thermal apply in reverse: model hours are forecast hours from the cycle time, not local solar time — align windows deliberately when comparing with satellite passes.

## License

MIT — see `LICENSE`.

# survey-currents

Surface-current and water-temperature acquisition engine for surveying and remote sensing — the first module of the earthwatch-suite **flow-field animation program** (the pipeline that produces mapped.earth-style animated current/temperature reels).

`survey-currents` pulls hourly current-vector (u/v) and water-temperature fields from **operational hydrodynamic forecast models** — not satellites — into one canonical `CurrentField` model with provenance, COG export, and survey-suite interoperability. Since v0.2.0 it also pulls **satellite-derived daily sea-surface temperature** from NOAA's GLSEA analysis (Great Lakes only) into a companion `GlseaField` model, plus lake-wide average temperature series. Since v0.3.0 it pulls **global satellite SST** from two more sources into a companion `SstField` model: NOAA OISST v2.1 (keyless, 0.25°, 1981–present) and NASA JPL MUR v4.1 (free Earthdata Login, ~1 km, 2002–present). Since v0.4.0 it pulls **global atmospheric reanalysis** (Copernicus ERA5: wind, pressure, air temperature, precipitation) into a companion `Era5Field` model. Since v0.5.0 it pulls **global surface currents** — NASA PODAAC OSCAR v2.0 (daily, 1993–present, free Earthdata Login) and the CMEMS `global-physics-daily` preset (1/12°, free CMEMS account) — into the canonical `CurrentField` model. Since v0.7.0 it pulls **passive-microwave sea-ice concentration** from the NOAA/NSIDC Sea Ice Index (G02135 v4.0, keyless, 1978–present) into a companion `IceField` model. Since v0.8.0 it pulls **satellite-observed precipitation** from NASA GPM IMERG V07 (half-hourly, 0.1°, 2000–present, free Earthdata Login, Early/Late/Final latency tiers) into a companion `RainField` model.

## Data sources

| Source | Access | Coverage | Resolution | Horizon |
|---|---|---|---|---|
| NOAA OFS via anonymous AWS S3 (`noaa-ofs-pds`, `noaa-nos-ofs-pds`) | No signup, unsigned requests | US coasts + Great Lakes | 50 m – 5 km | 48–120 h |
| Copernicus Marine Service (`copernicusmarine` toolbox) | Free account | Global ocean | 1/12° (~9 km) | NRT + forecast |
| NOAA GLSEA via ERDDAP griddap (`GLSEA_ACSPO_GCS`) | No signup | **Great Lakes only** | ~1.5 km | daily analysis, 2006–present |
| NOAA GLSEA via ERDDAP tabledap (`glsea_avgtemps_3`) | No signup | Great Lakes (per-lake daily averages) | lake-wide | daily, 2006–present |
| NOAA OISST v2.1 via CoastWatch ERDDAP (`ncdcOisst21Agg`) | No signup | **Global ocean** | 0.25° (~28 km) | daily analysis, 1981–present |
| NASA JPL MUR v4.1 via Earthdata OPeNDAP | Free Earthdata Login | **Global ocean** | ~0.01° (~1 km) | daily analysis, 2002–present |
| Copernicus ERA5 via CDS API (`reanalysis-era5-single-levels`) | Free CDS account | **Global atmosphere** | 0.25° (~28 km) | hourly reanalysis, 1940–present |
| NASA PODAAC OSCAR v2.0 via Earthdata OPeNDAP | Free Earthdata Login | **Global ocean currents** | 0.25° (~28 km) | daily averages, 1993–present |
| CMEMS global ocean physics (`global-physics-daily` preset) | Free CMEMS account | **Global ocean currents** | 1/12° (~9 km) | daily analysis+forecast |
| NASA GPM IMERG V07 via GES DISC HTTPS | Free Earthdata Login | **Global precipitation** | 0.1° (~10 km) | half-hourly, 2000–present (Early/Late/Final runs) |
| NASA Black Marble VNP46A2 V002 via LAADS DAAC | Free Earthdata Login | **Global night lights** | 15″ native, 0.05° default output | daily DNB radiance, 2012–present |
| GEBCO 2024 via BODC/CEDA open download (per-tile ranged extraction) | No signup | **Global topography/bathymetry** | 15″ native; block-averaged outputs | static 2024 compilation |
| Natural Earth via anonymous S3 | No signup | **Global coastline/country vectors** | 110m/50m/10m cartographic scales | static |
| NOAA IBTrACS v04r01 via NCEI HTTPS (download-once, cached) | No signup | **Global tropical-cyclone best tracks** | 3-hourly fixes | 1980–present (1842– with `--full-archive`) |

Registered NOAA models include **GLOFS** (Great Lakes, 5 km, 60 h), **LMHOFS** (Lake Michigan/Huron, 50 m–2.5 km, 120 h — the source of the Lake Michigan reel), LEOFS, CBOFS, DBOFS, GoMOFS, WCOFS, NGOFS2, SFBOFS, TBOFS, CIOFS, CREOFS, SSCOFS. Full table in `docs/DATA_SOURCES.md`.

## Install

```bash
pip install survey-currents                  # engine core (stdlib + numpy)
pip install "survey-currents[noaa]"          # NetCDF parsing (xarray/netCDF4)
pip install "survey-currents[s3]"            # faster S3 via boto3 (optional)
pip install "survey-currents[cmems]"         # Copernicus Marine toolbox
pip install "survey-currents[raster]"        # GeoTIFF/COG export (rasterio)
pip install "survey-currents[imerg]"         # IMERG HDF5 parsing (h5py)
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

Copernicus ERA5 reanalysis (global atmosphere — free CDS account, needs `cdsapi` + `netCDF4`):

```python
from currents.currents_global import fetch_oscar, fetch_cmems_currents

# Daily surface currents over the Gulf Stream, one granule every 5 days —
# free Earthdata Login required (see docs/DATA_SOURCES.md).
field = fetch_oscar(bbox=(-81.0, 25.0, -55.0, 43.0),
                    start="2024-01-01", end="2024-12-31", stride_days=5)
print(field.u.shape)          # (74, ny, nx) eastward current, m/s
print(field.provenance["collections_used"])  # ['final'] — the pick rule

# Higher-resolution currents via the CMEMS global-physics-daily preset
# (free CMEMS account + copernicusmarine toolbox).
field = fetch_cmems_currents(bbox=(-81.0, 25.0, -55.0, 43.0),
                             start="2024-01-01", end="2024-01-31")
print(field.u.shape, field.temperature.shape)  # (31, ny, nx) each
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
survey-currents fetch-era5 --variables wind,msl --bbox=-98,18,-80,31 \
    --start 2024-01-01 --end 2024-01-31 --stride-hours 24 --out gom_era5
survey-currents era5-synthetic --variables wind,msl --nt 4 --out era5_demo
survey-currents glsea-synthetic --nt 4 --out glsea_demo
survey-currents fetch-oisst --bbox -80,20,-60,40 \
    --start 2020-01-01 --end 2020-12-31 --stride-days 30 --out atlantic_sst
survey-currents fetch-mur --bbox -80,20,-60,40 \
    --start 2024-01-01 --end 2024-01-31 --stride-days 7 --out atlantic_mur
survey-currents sst-synthetic --nt 4 --out sst_demo
survey-currents fetch-oscar --bbox -81,25,-55,43 \
    --start 2024-01-01 --end 2024-12-31 --stride-days 5 --out gulfstream_oscar
survey-currents fetch-cmems-currents --bbox -81,25,-55,43 \
    --start 2024-01-01 --end 2024-01-31 --out gulfstream_cmems
survey-currents currents-synthetic --nt 4 --out currents_demo
# NASA FIRMS active fires need a free MAP_KEY: export FIRMS_MAP_KEY=...
survey-currents fetch-firms --bbox -125,32,-114,42 \
    --start 2024-08-01 --end 2024-08-07 --instruments VIIRS_SNPP --out ca_fires
survey-currents fires-synthetic --n 60 --out fires_demo
# NASA GPM IMERG precipitation needs a free Earthdata Login:
# export EARTHDATA_USERNAME=... EARTHDATA_PASSWORD=...
survey-currents fetch-imerg --bbox=-125,25,-66,49 \
    --start 2024-01-01 --end 2024-01-07 --run late --out conus_rain
survey-currents rain-synthetic --out rain_demo
# GEBCO 2024 topography/bathymetry + Natural Earth vectors (no account):
survey-currents fetch-gebco --bbox=-125,25,-66,49 --resolution 0.25 \
    --out conus_topo
survey-currents fetch-naturalearth --bbox=-125,25,-66,49 \
    --scale 110m --layers coastline,countries --out conus_ne
survey-currents basemaps-synthetic --out topo_demo
# NOAA IBTrACS tropical-cyclone best tracks (no account):
survey-currents fetch-ibtracs --bbox=-100,10,-60,40 \
    --start 2024-08-01 --end 2024-11-30 --storm-name milton --out milton
survey-currents storms-synthetic --out storms_demo
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

Atmospheric reanalysis converges on `Era5Field` (`src/currents/era5.py`):

- `grids`: dict of (nt, ny, nx) numpy masked arrays (`u10`/`v10` in m/s, `msl` in hPa, `t2m` in °C, `tp` in mm)
- `times`: ISO-8601 timestamps (hourly); `lats`/`lons`: increasing vectors (−180..180)
- `values` property: the rendered base grid (wind speed for `wind`, else the single grid); `overlay_grids`: non-base variables for contour overlays (e.g. isobars)
- `spatial_mean()`, `select_time()`, `select_bbox()`, JSON round-trip, `Era5Field.synthetic()`
- `fetch_era5(variables, bbox, start, end, stride_hours=6)` via `cdsapi` (lazy import, free CDS account); requests chunked by calendar month; provenance carries the CDS request dicts, SHA-256, and retrieval timestamp

Active-fire detections converge on `FireField` (`src/currents/fires.py`):

- Point detections: `times`/`lats`/`lons` plus `brightness` (K), `frp` (MW), `confidence`, `satellite`, `instrument`, `daynight`
- `fetch_firms(bbox, start, end, instruments=("VIIRS_SNPP",))` via the NASA FIRMS area API (free MAP_KEY, `FIRMS_MAP_KEY` env; verified live 2026-09-26); requests looped in ≤5-day windows with per-date NRT/SP product tiering (`firms_product_for`, pure and offline-testable); antimeridian boxes split
- `to_density_grid(resolution=0.25, frp_weighted=False)` bins detections into daily fire-count (or FRP-weighted MW) grids — the plain `times`/`lats`/`lons`/`values` dict `survey-viz` renders with zero changes
- `select_time()`, `select_bbox()`, JSON round-trip, `FireField.synthetic()`; provenance carries the key-**redacted** request URLs, per-request SHA-256, and retrieval timestamp (the MAP_KEY is never stored)

Daily sea-ice concentration converges on `IceField` (`src/currents/sea_ice.py`):

- Gridded percent (0–100): `times`/`lats`/`lons`/`values` on a regular lat/lon grid; pole-hole/coast/land/missing cells are NaN
- `fetch_nsidc_sic(bbox, start, end, hemisphere="auto", stride_days=1, resolution=0.25)` — downloads the G02135 v4.0 daily concentration GeoTIFFs over **keyless anonymous HTTPS** (no account; verified live 2026-09-26), decodes the ×10-scaled uint16 values, and nearest-neighbor reprojects the native 25 km NSIDC polar-stereographic grids (EPSG:3411/3412, Hughes 1980) onto the target grid; record 1978-11-01–present, missing days skipped with a provenance note
- `select_time()`, `select_bbox()`, JSON round-trip, `IceField.synthetic()`; provenance carries the exact file URLs, per-file SHA-256, reprojection method, and retrieval timestamp
- CLI: `fetch-nsidc`, `ice-synthetic`

Daily satellite precipitation converges on `RainField` (`src/currents/imerg.py`):

- Gridded rates/totals: `times`/`lats`/`lons`/`values` on a regular lat/lon grid; mm/hr for `accumulate="native"` (48 half-hourly steps/day), mm/day for `accumulate="daily"` (default)
- `fetch_imerg(bbox, start, end, accumulate="daily", run="late", stride_days=1)` — downloads the GPM IMERG V07 half-hourly HDF5 granules over **authenticated HTTPS** (free Earthdata Login; verified live 2026-09-27 — OPeNDAP catalog/`.dds`/`.das` metadata are keyless, data download 302s to URS), subsets to the bbox locally, and accumulates; runs: `early` (~4 h), `late` (~14 h, default), `final` (~3.5 months, gauge-adjusted); record 2000-06-01–present; 404 granules skipped with a `skipped_slots` provenance note and per-day `{"expected": 48, "retrieved": n}` coverage
- `select_time()`, `select_bbox()`, JSON round-trip, `RainField.synthetic()`; provenance carries the exact download URLs, per-file SHA-256, run, accumulation mode, and retrieval timestamp
- CLI: `fetch-imerg`, `rain-synthetic`

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
    sst_global.py  # NOAA OISST v2.1 + NASA JPL MUR v4.1 global SST
    era5.py        # Copernicus ERA5 hourly reanalysis (wind/msl/t2m/tp) via cdsapi
    currents_global.py  # NASA PODAAC OSCAR v2.0 + CMEMS global-physics-daily currents
    fires.py       # NASA FIRMS active-fire detections (area API) -> FireField + density grids
    sea_ice.py     # NSIDC G02135 v4.0 daily sea-ice concentration (keyless HTTPS) -> IceField
    imerg.py       # NASA GPM IMERG V07 half-hourly precipitation (Earthdata HTTPS) -> RainField
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

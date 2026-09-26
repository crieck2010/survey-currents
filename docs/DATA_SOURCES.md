# Data sources

All endpoints below were verified on 2026-09-25. Model products are
**forecast-model output** (assimilated, but still model) — treat as
guidance and validate against in-situ or satellite SST where it matters.

## 1. NOAA Operational Forecast Systems (no account)

Public S3 buckets, unsigned requests (`--no-sign-request` equivalent —
this engine signs nothing):

| Bucket | Contents | Prefix template |
|---|---|---|
| `noaa-ofs-pds` | NOMADS production, last 30 days | `OFS.YYYYMMDD` |
| `noaa-nos-ofs-pds` | CO-OPS operational archive | `OFS/netcdf/YYYYMM/` |

### Registered models

| Code | System | Region | Resolution | Horizon | Cycles |
|---|---|---|---|---|---|
| GLOFS | Great Lakes OFS | All five Great Lakes | 5 km | 60 h | 00/06/12/18 |
| LMHOFS | Lake Michigan-Huron OFS | Lakes Michigan + Huron | 50 m – 2.5 km | 120 h | 00/06/12/18 |
| LEOFS | Lake Erie OFS | Lake Erie | 400 m – 4 km | 120 h | 00/06/12/18 |
| CBOFS | Chesapeake Bay OFS | Chesapeake Bay | 50 m – 3 km | 48 h | 00/06/12/18 |
| CIOFS | Cook Inlet OFS | Cook Inlet, AK | 10 m – 3.5 km | 48 h | 00/06/12/18 |
| CREOFS | Columbia River Estuary OFS | Columbia River Estuary | 100 m – 4 km | 48 h | 00/06/12/18 |
| DBOFS | Delaware Bay OFS | Delaware Bay | 100 m – 3 km | 48 h | 00/06/12/18 |
| GoMOFS | Gulf of Maine OFS | Gulf of Maine | 700 m | 72 h | 00/06/12/18 |
| NGOFS2 | Northern Gulf of Mexico OFS | N. Gulf of Mexico | 45 m – 300 m | 48 h | 00/06/12/18 |
| SFBOFS | San Francisco Bay OFS | San Francisco Bay | 100 m – 4 km | 48 h | 00/06/12/18 |
| TBOFS | Tampa Bay OFS | Tampa Bay | 100 m – 1.2 km | 48 h | 00/06/12/18 |
| WCOFS | West Coast OFS | US West Coast | 4 km | 72 h | 00/06/12/18 |
| SSCOFS | Salish Sea + Columbia River OFS | Salish Sea / Columbia R. | 100 m – 3 km | 72 h | 00/06/12/18 |

Filename conventions differ per model, so the engine lists the date
prefix and filters keys by OFS-code / date / cycle (`tCCz`) /
forecast-hour (`fHHH`) tokens rather than assuming exact names.

### Variables

The parser resolves variable names tolerantly. Typical NOS OFS
NetCDFs carry eastward/northward velocity (`u`/`v`), water temperature
(`temp`), and water level (`zeta`); only the surface level is kept.
Known spellings are listed in `noaa_ofs.U_CANDIDATES` /
`V_CANDIDATES` / `TEMP_CANDIDATES`; explicit names can be passed to
`parse_ofs_netcdf()` and are recorded in provenance.

### The Lake Michigan reel recipe

The mapped.earth "Lake Michigan never sits still" reel (16–23 Sep 2026,
hourly currents colored by water temperature) is reproducible as:

```bash
survey-currents fetch-noaa --ofs LMHOFS --date 2026-09-16 --cycle 00 \
    --hours 0-167 --bbox -92.5,41.5,-84.5,46.5 --out lmhofs_week
survey-currents export-cogs lmhofs_week.json --out cogs/
# -> survey-flow renders frames -> survey-animate encodes the reel
```

Note: native LMHOFS NetCDFs are FVCOM unstructured grids; v0.1.0 parses
regular-grid holdings and raises a clear error otherwise (see README
limitations). GLOFS (regular 5 km grid) works out of the box for
basin-scale reels.

## 2. Copernicus Marine Service (free account)

- Register: <https://marine.copernicus.eu/register> (2-minute signup,
  instant access, no quotas on the toolbox path).
- Authenticate once: `copernicusmarine login`, or set
  `COPERNICUSMARINE_SERVICE_USERNAME` /
  `COPERNICUSMARINE_SERVICE_PASSWORD`.
- The engine wraps the toolbox `subset()` (server-side bbox / time /
  variable / depth slicing) — it never walks the raw archive.

### Presets

| Preset | Dataset | Variables | Cadence |
|---|---|---|---|
| `global-physics-daily` | `cmems_mod_glo_phy_anfc_0.083deg_P1D-m` | `uo`, `vo`, `thetao` | daily |
| `global-physics-hourly` | `cmems_mod_glo_phy_anfc_0.083deg_PT1H-m` | `uo`, `vo`, `thetao` | hourly |

Resolution 1/12° (~9 km); surface level only. Any other dataset id the
toolbox recognises can be used via a custom `CmemsPreset`.

## 3. NOAA GLSEA — Great Lakes satellite SST (no account)

Sea Surface Temperature from the Great Lakes Surface Environmental
Analysis (ACSPO GLSEA), served by NOAA GLERL's ERDDAP. This is
**satellite-derived analysis** (ACSPO L3S-LEO SST from NPP, NOAA-20,
MetOp A/B/C), not hydrodynamic model output — the complement to the
OFS/CMEMS model fields above. Verified live 2026-09-26.

### 3a. Gridded daily SST — `GLSEA_ACSPO_GCS` (griddap)

- Endpoint: `https://apps.glerl.noaa.gov/erddap/griddap/GLSEA_ACSPO_GCS.nc`
- Variable `sst` (float, °C, `_FillValue` −99999.0 → masked); dims `(time, latitude, longitude)`
- Daily timesteps stamped **12:00 UTC**; coverage **2006–present**
- Grid ~0.014° (~1.5 km); the engine samples the time axis with a
  day-stride (default 30) and full-resolution lat/lon

Request shape (built by `currents.glsea.glsea_sst_url`):

```
sst[(<start>ISO):<stride>:(<end>ISO)][(<lat_min>):1:(<lat_max>)][(<lon_min>):1:(<lon_max>)]
```

with ISO like `2016-01-01T12:00:00Z`.

**Grid quirk — the −92.42 longitude floor.** The longitude axis is
clipped to the lakes region: `actual_range = -92.4199507342304,
-75.8816402880531` (latitude `38.8749871947297, 50.6059751976539`).
Any bbox west of the floor raises a clear `ValueError` *before* any
download (`validate_glsea_bbox`); latitude/longitude constraints must
be ascending.

### 3b. Lake-average daily temperature — `glsea_avgtemps_3` (tabledap)

- Endpoint: `https://apps.glerl.noaa.gov/erddap/tabledap/glsea_avgtemps_3.csv?Year,Day,<Col>&Year>=<y0>`
  (`>=` is percent-encoded as `%3E` in the request — raw `>` gets
  dropped by some proxies)
- Columns: `Year, Day` (day-of-year), plus `Sup, Mich, Huron, Erie,
  Ont` — one lake-average SST (°C) per day
- Parsed with stdlib `csv` only (no netCDF4 needed); the ERDDAP units
  row and rows with missing temperatures are skipped, and rows are
  filtered to the requested `[start, end]` window
- Lake name mapping (`currents.glsea.LAKE_COLUMNS`):
  `superior→Sup`, `michigan→Mich`, `huron→Huron`, `erie→Erie`,
  `ontario→Ont`; unknown names raise `ValueError`

## Choosing a source

- **US Great Lakes / coasts, no signup:** NOAA OFS (LMHOFS for Lake
  Michigan detail, GLOFS basin-wide).
- **Great Lakes satellite SST, no signup:** NOAA GLSEA — observed
  analysis rather than model output; pairs with OFS fields for
  model-vs-satellite temperature comparison.
- **Any other coastline:** CMEMS global physics.
- **Blending with satellites:** use `align_to_thermal_zone()` to compare
  model water temperature against survey-thermal Landsat LST passes over
  the same zone and window.

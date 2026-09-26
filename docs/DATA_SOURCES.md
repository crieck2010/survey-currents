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

## 4. NOAA OISST v2.1 — global satellite SST (no account)

Daily global sea-surface temperature analysis (AVHRR-only final
product), served by NOAA CoastWatch ERDDAP. The keyless global-SST
companion to GLSEA. Verified live 2026-09-26.

- Endpoint: `https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg.nc`
- Variable `sst` (float, °C, `_FillValue` −9.99 → masked); dims
  `(time, zlev, latitude, longitude)` — the singleton `zlev` axis
  (0.0 m) **must** be indexed explicitly as `[(0.0)]`; ERDDAP 404s
  without it
- Daily timesteps stamped **12:00 UTC**; coverage **1981-09-01–present**
- Grid 0.25°: latitude −89.875..89.875, longitude **0.125..359.875
  (0–360 convention)**

Request shape (built by `currents.sst_global.oisst_sst_urls`):

```
sst[(<start>ISO):<stride>:(<end>ISO)][(0.0)][(<lat_min>):(<lat_max>)][(<lon360_min>):(<lon360_max>)]
```

with ISO like `2020-01-01T12:00:00Z`.

**Longitude handling.** The engine accepts conventional −180..180
bboxes and converts to 0–360 internally (`oisst_lon_windows`);
downstream fields are normalized back to −180..180, sorted
increasing. Three cases:

- Normal boxes → one request, e.g. `(−80, ., −60, .)` → `(280.0, 300.0)`.
- **Antimeridian-crossing** boxes, e.g. `(170, ., −170, .)` → one
  wrapped window `(170.0, 190.0)`.
- Boxes whose 0–360 window **crosses 360°** (e.g. `(−170, ., 170, .)`
  → `(190.0, 530.0)`) cannot be expressed in one ERDDAP query: the
  engine issues **two requests** and concatenates along longitude —
  never silently truncated. Full-globe boxes (span ≥ 359.9°) request
  the whole grid.

**Time chunking.** Long windows are split into ≤5-year requests
(`OISST_MAX_YEARS_PER_REQUEST`) because ERDDAP drops very long time
ranges; results are concatenated along time.

## 5. NASA JPL MUR v4.1 — global ultra-high-resolution SST (free Earthdata account)

The GHRSST Level 4 MUR global foundation SST analysis (~0.01°,
~1 km — the highest-resolution global SST in this engine), via the
NASA Earthdata OPeNDAP service. Verified live via NASA CMR
2026-09-26 (collection `C1996881146-POCLOUD`).

**Granule discovery (not construction).** Before any data request,
`fetch_mur` calls the **public NASA CMR granule search API**
(`cmr.earthdata.nasa.gov/search/granules.json`, keyless) for
`C1996881146-POCLOUD` over the requested window
(`cmr_search_mur_granules`), then matches each sampled day to a real
granule (`mur_match_granules`, by title stamp then `time_start`). The
OPeNDAP service URL for each granule is the OPeNDAP link CMR
advertises on it when present, else the documented Earthdata URL
pattern for the collection (`_mur_service_url` — recorded per
granule as `"cmr-link"` / `"constructed"` in provenance). A sampled
day with no discovered granule raises `RuntimeError` — an honest gap,
never an invented granule name. `mur_granule_title` remains as a
tested naming/reference helper (the CMR-verified naming pattern),
not as the discovery mechanism.

- Service: `https://opendap.earthdata.nasa.gov/collections/C1996881146-POCLOUD/granules/<granule-title>`
- Granule titles embed the analysis time (always 09:00 UTC), e.g.
  `20260925090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1`
- Variable `analysed_sst` (**Kelvin** on the wire, usually packed as
  scaled shorts; converted to °C by the engine); dims `(time, lat,
  lon)` — one analysis per granule
- Coverage **2002-06-01–present**; one granule per sampled day
  (`stride_days`, default 30)

Request shape (built by `currents.sst_global.mur_subset_urls`):

```
<granule>.nc?analysed_sst[0:1:0][<j0>:1:<j1>][<i0>:1:<i1>]
```

Grid index windows come from the granule's own `.das`/`.dds`
metadata (parsed by `_parse_mur_grid`) — no hardcoded grid, so the
engine survives upstream grid revisions. MUR's longitude axis is
−180..180; antimeridian-crossing bboxes become two subset requests
concatenated along longitude.

**Authentication.** OPeNDAP requires a free Earthdata Login account.
`earthdata_credentials()` checks `EARTHDATA_USERNAME` /
`EARTHDATA_PASSWORD` first, then a `~/.netrc` entry for
`opendap.earthdata.nasa.gov`. With no credentials (or an HTTP 401),
`fetch_mur` raises `CredentialsMissing` with setup instructions —
the keyless OISST source keeps working.

## 6. Copernicus ERA5 — global atmospheric reanalysis (free CDS account)

ERA5 hourly data on single levels (`reanalysis-era5-single-levels` on the
Copernicus Climate Data Store): hourly 1940–present, 0.25° global grid.
Fetched through the `cdsapi` package (lazy import — the engine imports
cleanly without it). Needs a **free** CDS account: register at
https://cds.climate.copernicus.eu/, accept the ERA5 Terms of Use, and
create `~/.cdsapirc` (`url:` + `key:`) or set `CDSAPI_URL` /
`CDSAPI_KEY`. Without working credentials `fetch_era5` raises
`CredentialsMissing` with these exact steps.

### Variables

| Short key | CDS variable name(s) | On-wire grid(s) | Converted units |
|---|---|---|---|
| `wind` | `10m_u_component_of_wind`, `10m_v_component_of_wind` | `u10`, `v10` | m/s (rendered as wind speed) |
| `msl` | `mean_sea_level_pressure` | `msl` | hPa (Pa ÷ 100) |
| `t2m` | `2m_temperature` | `t2m` | °C (K − 273.15) |
| `tp` | `total_precipitation` | `tp` | mm per hourly step (m × 1000) |

`fetch_era5(variables, bbox, start, end, stride_hours=6)` — `variables`
is one key or a list; `"wind"` fetches the u/v pair. `stride_hours=6`
samples 00/06/12/18 UTC; `24` gives daily 12:00 UTC; `168` weekly, etc.

Request shape (built by `currents.era5.era5_request`, one per calendar
month per longitude window — CDS request-size hygiene):

```python
{
    "product_type": "reanalysis",
    "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
    "year": "2024", "month": "01",
    "day": ["06", "07"], "time": ["00:00", "06:00", "12:00", "18:00"],
    "area": [31.0, -98.0, 29.0, -97.0],  # N, W, S, E
    "grid": "0.25/0.25",
    "data_format": "netcdf", "download_format": "unarchived",
}
```

Notes:

- bboxes use the conventional −180..180 convention; antimeridian-crossing
  boxes split into two CDS requests and are concatenated (with the shared
  180°/−180° seam column deduplicated) — see `era5_area_windows`.
- The CDS latitude axis arrives descending (N→S) and longitudes 0–360;
  both are normalized to increasing / −180..180 on ingest
  (`_parse_era5_bytes`). Time comes from `valid_time` (new CDS) or `time`.
- ERA5 is a reanalysis with a ~5-day release lag — requests past today
  are refused before any download.
- Provenance records the CDS dataset id, the exact request dicts,
  per-payload SHA-256 (combined), byte counts, and retrieval time.

## Choosing a source

- **US Great Lakes / coasts, no signup:** NOAA OFS (LMHOFS for Lake
  Michigan detail, GLOFS basin-wide).
- **Great Lakes satellite SST, no signup:** NOAA GLSEA — observed
  analysis rather than model output; pairs with OFS fields for
  model-vs-satellite temperature comparison.
- **Global satellite SST, no signup:** NOAA OISST v2.1 — daily 0.25°
  analysis, 1981–present; the global answer to GLSEA.
- **Global ultra-high-resolution SST, free Earthdata account:**
  NASA JPL MUR v4.1 — daily ~1 km analysis, 2002–present; coastal
  detail OISST cannot resolve.
- **Global atmosphere (wind, pressure, air temperature, precipitation),
  free CDS account:** Copernicus ERA5 — hourly 0.25° reanalysis,
  1940–present; the answer for "winds", "storm", "heat", "rain"
  visualizations (see `currents.era5`).
- **Any other coastline:** CMEMS global physics.
- **Blending with satellites:** use `align_to_thermal_zone()` to compare
  model water temperature against survey-thermal Landsat LST passes over
  the same zone and window.

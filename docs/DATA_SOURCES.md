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

## 7. NASA PODAAC OSCAR v2.0 — global surface currents (free Earthdata account)

OSCAR (Ocean Surface Current Analyses Real-time) v2.0: daily-averaged
surface currents, **1993–present**, 0.25° global grid, variables `u`/`v`
(m/s, east/north positive). This is the answer for "currents"
visualizations outside the NOAA OFS footprints.

Access goes through the **Earthdata OPeNDAP** endpoint and **requires a
free Earthdata Login** — unauthenticated requests are redirected to the
login page (HTTP 302). Set `EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` or
add a `~/.netrc` entry for `opendap.earthdata.nasa.gov` (register free
at https://urs.earthdata.nasa.gov/users/new). Without credentials
`fetch_oscar` raises `CredentialsMissing` with these exact steps (a 401
from the server maps to the same error).

Three latency tiers live in separate NASA CMR collections; the fetch
picks the right tier per date with a pure, offline-testable rule
(`oscar_collection_for`):

| Tier | CMR collection | Coverage | Latency | Picked when |
|---|---|---|---|---|
| `final` | `C2098858642-POCLOUD` | 1993-01-01 → present | ~1.5 yr | date older than today − 540 d |
| `interim` | `C2102959417-POCLOUD` | 2020-01-01 → present | ~1 mo | date older than today − 45 d |
| `nrt` | `C2102958977-POCLOUD` | 2021-01-01 → present | ~2 d | otherwise |

Granule names are deterministic *and* verified through the keyless NASA
CMR granule search (same pattern as the MUR adapter) — the fetch never
invents a granule name:

```
oscar_currents_final_YYYYMMDD.nc        # e.g. oscar_currents_final_20240115.nc
oscar_currents_interim_YYYYMMDD.nc
oscar_currents_nrt_YYYYMMDD.nc
```

Request shape (one constrained OPeNDAP URL per date per longitude
window, built by `currents.currents_global.oscar_subset_urls`):

```
https://opendap.earthdata.nasa.gov/collections/C2098858642-POCLOUD/granules/oscar_currents_final_20240115.nc
    ?u[0:1:0][1116:1:1220][459:1:531],v[0:1:0][1116:1:1220][459:1:531]
```

Notes:

- On-wire dimension order is the unusual **(time, longitude, latitude)** —
  `u[time][lon][lat]` — not (time, lat, lon). The parse step asserts this
  and transposes to the canonical (nt, ny, nx); a reordered product
  fails loudly instead of silently transposing.
- The grid is the 0–360 convention (lon 0..359.75, lat −89.75..89.75);
  −180..180 bboxes are converted internally (same approach as OISST),
  antimeridian-crossing boxes wrap into a single 0–360 window, and
  360°-crossing windows split into two requests and are concatenated
  with the seam deduplicated.
- Fill value −999.0 is masked. Frame timestamps are the daily granule
  dates at 00:00 UTC (OSCAR granules are daily averages).
- `fetch_oscar(bbox, start, end, stride_days=5)` returns a
  `CurrentField` with `temperature=None` (OSCAR is currents-only).
  Provenance records the exact OPeNDAP URLs, per-payload SHA-256
  (combined), retrieval time, the three CMR collection ids, the pick
  rule, and the per-date collection + granule title used.
- The CMEMS `global-physics-daily` preset (`uo`/`vo`/`thetao`, 1/12°,
  daily) is also wrapped to the standard fetch signature as
  `fetch_cmems_currents(bbox, start, end, stride_days=1)` (see
  `src/currents/cmems.py` for the toolbox path): one NetCDF is
  downloaded for the whole range, then timesteps are stride-selected.
  `thetao` is potential temperature (°C), carried as-is — recorded in
  provenance; it is not a foundation SST.
- **Verified-source corrections (2026-09-26):** OSCAR is **not** on
  CoastWatch ERDDAP — the old `jplOscar_LonPM180` dataset id 404s
  (removed); v2.0 is served from Earthdata OPeNDAP behind Earthdata
  Login. NOMADS OPeNDAP is retired (Service Change Notice 25-81) and is
  not used. Both corrections are recorded in OSCAR provenance.

## 8. NASA FIRMS — global active-fire detections (free MAP_KEY)

Fire Information for Resource Management System: thermal-anomaly
detections from MODIS (Terra/Aqua, 1 km, Nov 2000–present) and VIIRS
(Suomi-NPP / NOAA-20 / NOAA-21, 375 m, 2012/2018/2023–present).

| Item | Value |
|---|---|
| API | `https://firms.modaps.eosdis.nasa.gov/api/area/csv/{MAP_KEY}/{PRODUCT}/{W},{S},{E},{N}/{DAY_RANGE}/{DATE}` |
| Auth | **Free MAP_KEY required** — request at https://firms.modaps.eosdis.nasa.gov/api/map_key/; export as `FIRMS_MAP_KEY` (verified live 2026-09-26: a bad key returns `Invalid MAP_KEY.`) |
| Window | `DAY_RANGE` 1–5 days with a `DATE` start; longer ranges are looped by the engine |
| Products | Per-instrument latency tiers picked per date by the pure, offline-testable `firms_product_for` (dates within 60 days of today → `*_NRT`; older → `*_SP`): `VIIRS_SNPP_NRT/_SP`, `VIIRS_NOAA20_NRT/_SP`, `VIIRS_NOAA21_NRT/_SP`, `MODIS_NRT/_SP` |
| Record | VIIRS columns `latitude,longitude,bright_ti4,bright_ti5,frp,scan,track,acq_date,acq_time,satellite,instrument,confidence,version,daynight`; MODIS columns `latitude,longitude,brightness,…,bright_t31,frp,daynight` (numeric confidence) — both parsed by `currents.fires` |
| Model | `FireField` (point detections) → `to_density_grid(resolution=0.25)` bins detections into daily fire-count (or FRP-weighted MW) grids — the plain `times`/`lats`/`lons`/`values` dict `survey-viz` renders |

Operational notes:

- The MAP_KEY is **never stored** in provenance, logs, or error messages —
  stored request URLs carry a `<redacted>` placeholder, and the
  provenance record keeps `map_key: "<redacted>"`.
- Antimeridian-crossing bboxes are split into two area-API requests.
- Active-fire *detections* only: FIRMS says "something hot is here now",
  not how much area burned. Burned-area / burn-severity mapping is the
  `survey-burn` module's domain (satellite imagery) — see the
  `survey-viz` parser's honest refusal for "burn scar" phrases.

## 9. NSIDC G02135 v4.0 — global sea-ice concentration (no account)

NOAA/NSIDC Sea Ice Index: daily passive-microwave sea-ice concentration,
November 1978–present, both hemispheres.

| Item | Value |
|---|---|
| Archive | `https://noaadata.apps.nsidc.org/NOAA/G02135/{north,south}/daily/geotiff/{YYYY}/{MM_Mon}/{N,S}_YYYYMMDD_concentration_v4.0.tif` |
| Auth | **None — keyless anonymous HTTPS** (verified live 2026-09-26: directory listings and real files return HTTP 200 with no credentials) |
| Grid | NSIDC polar stereographic, Hughes 1980 ellipsoid, 25 km — North EPSG:3411 (304×448), South EPSG:3412 (316×332) |
| Encoding | Unsigned 16-bit; concentration scaled ×10 (divide by 10 → percent 0–100); 2510 = Arctic pole hole, 2530 = coast, 2540 = land, 2550 = missing → NaN |
| Record | 1978-11-01–present; SMMR era (1978–1987) is every other day; no data 1987-12-03 → 1988-01-13 |
| Model | `IceField` — regular lat/lon `times`/`lats`/`lons`/`values` (percent), flags masked as NaN, per-file provenance |

Operational notes:

- The engine reprojects the native grids onto the caller's lat/lon grid
  with **nearest-neighbor** sampling (documented choice — it preserves
  the 15%-threshold ice edge exactly; bilinear would smear it).
  The formulas are verified against a live file: the pole maps to (0, 0)
  m, i.e. pixel (154, 234) on the north grid, where the pole-hole flag
  sits.
- Missing days (SMMR off-days, the 1987/88 outage) are skipped with a
  `provenance["skipped_days"]` note; if every day is missing the fetch
  raises instead of returning an empty field.
- GeoTIFF reading is stdlib + numpy only (the files are uncompressed
  single-band 16-bit) — the keyless NSIDC path needs **no new
  dependencies**. Georeferencing comes from the file's own
  ModelPixelScaleTag/ModelTiepointTag, falling back to the documented
  grid constants.
- Values 1–15% are kept as-is but are statistically irrelevant per the
  G02135 user guide (passive-microwave uncertainty below 15%); the
  guide's 15% cutoff defines the ice *extent* edge.
- `survey-viz` routes the `sea-ice` variable to NSIDC **only** in polar
  regions. Land ice (glaciers, ice sheets, icebergs) is a different
  physical product — see the `survey-viz` parser's honest refusal for
  those phrases.

## 10. NASA GPM IMERG V07 — global half-hourly precipitation (free Earthdata Login)

Integrated Multi-satellitE Retrievals for GPM: half-hourly global
precipitation rates, 2000-06-01–present, three latency runs.

| Item | Value |
|---|---|
| Archive | `https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/{GPM_3IMERGHH.07, GPM_3IMERGHHE.07, GPM_3IMERGHHL.07}/{YYYY}/{DOY}/` (the `/data/` path serves the same granules as `/opendap/` without the DAP overhead) |
| File naming | `3B-HHR[-E\|-L].MS.MRG.3IMERG.{YYYYMMDD}-S{HHMMSS}-E{HHMMSS}.{HHMM}.V07B.HDF5` — Early/Late carry the `-E`/`-L` infix; the 4th dot-field is the half-hour start (`0000`, `0030`, …, `2330`) |
| Auth | **Free Earthdata Login** (`EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` or `.netrc`, same as the MUR path). Verified live 2026-09-27: the OPeNDAP catalog and `.dds`/`.das` metadata are **keyless** (HTTP 200, no credentials), but data download redirects (302) to `urs.earthdata.nasa.gov` for the OAuth handshake and returns 401 without a login |
| Runs | `early` (~4 h latency, collection `GPM_3IMERGHHE`), `late` (~14 h, `GPM_3IMERGHHL`, **the default**), `final` (~3.5 months, `GPM_3IMERGHH`, gauge-adjusted) |
| Grid | 0.1° global, 3600 × 1800 (`Grid/lon`, `Grid/lat`); on-wire HDF5 dataspace is `(time=1, lon, lat)` per the `.das` — the adapter resolves the axis order from the actual vector sizes and raises on layout drift |
| Encoding | `precipitation` in mm/hr, `_FillValue` ≈ −9999.9 → NaN |
| Record | 2000-06-01–present for all runs (TRMM June 2000 – May 2014, GPM June 2014 – present). The live CMR archive shows V07 reprocessing reaching 1998-01-01 for Early/Final, but the adapter enforces the documented floor: pre-2000 TRMM-era estimates come from a much smaller constellation and are lower quality, especially at high latitudes |
| Model | `RainField` — `times`/`lats`/`lons`/`values` (mm/day for `accumulate="daily"`, mm/hr for `accumulate="native"`), per-file SHA-256 provenance with exact download URLs |

Operational notes:

- Each fetch downloads **complete global HDF5 granules (48 per day)**
  and subsets to the bbox locally — simple and robust, but heavy for
  long windows (each granule is tens of MB). A constrained OPeNDAP
  range request would be cheaper per byte, but server behavior for
  authenticated `.dods` queries was unreliable during verification
  (connection resets), so full-file HTTPS is the documented path.
  Prefer `run="late"` (better quality than Early for a modest latency
  cost) and keep windows short, or use `stride_days` for sampling.
- `accumulate="daily"` (default) sums the 48 half-hourly mm/hr rates
  into mm/day totals (each slot contributes rate × 0.5 h; a NaN slot
  contributes 0, an all-NaN day stays NaN). Days with missing granules
  (404s are skipped, e.g. the most recent day for the Final run)
  accumulate only the retrieved slots — check
  `provenance["day_coverage"][day]` (`{"expected": 48, "retrieved": n}`)
  and `provenance["skipped_slots"]` before trusting a daily total.
- IMERG is the **observed** precipitation source: it complements the
  ERA5 reanalysis `tp` (1940–present, model physics). `survey-viz`
  routes recent/event precipitation phrases to IMERG and long-record
  or trend phrases to ERA5.
- HDF5 parsing needs `h5py`: `pip install "survey-currents[imerg]"`
  (also in the `full` extra). The engine core stays stdlib+numpy.

## 11. NASA Black Marble VNP46A2 V002 — global daily night lights (free Earthdata Login)

Moonlight-adjusted nighttime lights: daily DNB radiance, 2012-01-19–present.

| Item | Value |
|---|---|
| Archive | LAADS DAAC: `https://data.laadsdaac.earthdatacloud.nasa.gov/prod-lads/VNP46A2/...` (Earthdata Cloud) and `https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/VNP46A2/...` (on-prem) |
| File naming | `VNP46A2.A{YYYY}{DOY}.h{HH}v{VV}.002.{production-stamp}.h5` — the production stamp is not predictable, so the adapter **discovers** exact download URLs per (day, tile) through NASA CMR granule search (`https://cmr.earthdata.nasa.gov/search/granules.json`, collection `C3365931269-LAADS`) |
| Auth | **Free Earthdata Login** (`EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` or `.netrc`, same as the MUR/IMERG path). Verified live 2026-09-26: CMR granule *discovery* is **keyless** (HTTP 200, no credentials), but both LAADS download endpoints redirect unauthenticated requests to `urs.earthdata.nasa.gov` (OAuth) — `CredentialsMissing` is raised otherwise |
| Tiles | V002 keeps the `hHHvVV` naming but the tiles are **10°×10° lat/lon tiles** (the V2 "15 arc-second linear lat/lon grid"), not sinusoidal: `h = floor((lon + 180) / 10)` (0..35), `v = floor((90 - lat_top) / 10)` (0..17). Verified against live CMR footprints (`h07v10` = lon −110..−100, lat −20..−10) |
| Grid | Native 15 arc-second (2400×2400 per tile); the adapter mosaics the tiles covering the bbox and NaN-aware block-averages to the requested `resolution` (default 0.05°) |
| Encoding | `Gap_Filled_DNB_BRDF_Corrected_NTL` (fallback `DNB_BRDF_Corrected_NTL`) in nW/cm²/sr; `_FillValue` → NaN; `scale_factor`/`add_offset` attributes honored when present |
| Record | 2012-01-19–present (daily; enforced with an honest error for earlier dates) |
| Model | `LightsField` — `times`/`lats`/`lons`/`values` (nW/cm²/sr, NaN for unlit/missing), per-file SHA-256 provenance with exact download URLs, tile list, and skipped tiles/days |

Operational notes:

- The VNP46A2 monthly (VNP46A3) and annual (VNP46A4) composites exist but
  are **not wired**: daily is the primary product because the reel
  factory needs a time axis. `product="monthly"`/`"annual"` raise a
  `ValueError` naming the future work instead of silently misbehaving.
- Each fetch downloads **one HDF5 tile per (day, tile)** (~tens of MB
  each, ~384 tiles/day globally) and mosaics locally — keep windows and
  bboxes tight, or use `stride_days` for sampling. `survey-viz` parses
  night-lights descriptions to daily cadence (one frame per day).
- HDF5 parsing needs `h5py`: `pip install "survey-currents[blackmarble]"`
  (also in the `full` extra). The engine core stays stdlib+numpy.
- A single-epoch Black Marble map is **not a change-detection product**:
  `survey-viz` parses "power outage"/"blackout" to a `power-outage`
  variable with no adapter and an honest refusal, rather than
  misrouting to the daily lights product.

## 12. GEBCO 2024 + Natural Earth — global topography/bathymetry and cartographic vectors (no account)

Basemap layers: GEBCO 2024 global terrain (15 arc-second grid, elevations
in metres, positive up) and Natural Earth coastline/country vectors for
context underlays. Neither needs an account.

| Item | Value |
|---|---|
| GEBCO archive | `https://www.bodc.ac.uk/data/open_download/gebco/gebco_2024/geotiff/` → redirects to a CEDA-hosted `gebco_2024_geotiff.zip` (~4.26 GB). Verified live 2026-09-27: anonymous HTTPS, byte-range requests honored (206). No login, no API key |
| GEBCO layout | Zip holds eight 90°×90° stripped, uncompressed int16 GeoTIFFs (`gebco_2024_n90.0_s0.0_w-180.0_e-90.0.tif`, …), each ~495–577 MB compressed / 933,255,948 bytes uncompressed (21600×21600 px, 15 arc-second; 21600 one-row strips) |
| GEBCO access | **Per-tile ranged extraction**: the adapter parses the zip central directory with a few KB of range requests, then downloads only the tile entries intersecting the bbox (no 4.26 GB download). Extracted tiles are cached locally with SHA-256 sidecars, verified on every hit, re-fetched on corruption |
| GEBCO grid | 15 arc-second lattice (tag layout verified against a real tile entry 2026-09-27: stripped uncompressed int16, positive Y pixel scale per the GeoTIFF spec, ModelPixelScaleTag, tiepoint, SampleFormat checked; window reads verified on the live tile); output grids are snapped to it and mosaicked with exact integer index math, so adjacent tiles butt-join without seams. `resolution="15s"` (native) or a degree value that is an integer multiple of 15 arc-seconds **and** divides the 90° tiles (e.g. `0.25`, `0.5`, `1.0`); coarser grids are block-averaged (mean). Native resolution has a 8192 px safety cap per axis — request a coarser grid for continental extents |
| GEBCO record | Static compilation (2024 release; GEBCO 2025 exists but has no confirmed keyless open-download path). The field carries one structural timestamp so renderers keep their 3D-field contract; `provenance["static_compilation"]` is `True` |
| Model | `TopoField` — `times` (single structural step) / `lats` / `lons` / `values` (metres, positive up), `select_bbox`, `to_json`/`from_json`, deterministic `synthetic()` |
| Natural Earth | `https://naturalearth.s3.amazonaws.com/{110m,50m,10m}_{physical,cultural}/ne_*.zip` — anonymous S3, verified live 2026-09-27. Scales `110m`/`50m`/`10m`; layers `coastline` (polylines) and `countries` (polygons). Zips cached with the same SHA-256 discipline; shapefiles parsed with a minimal stdlib reader (no fiona/geopandas) and clipped to the bbox |
| CLI | `fetch-gebco --bbox … --resolution 0.25`, `fetch-naturalearth --bbox … --scale 110m --layers coastline,countries`, `basemaps-synthetic` |

Operational notes:

- The reader is file-backed: only the IFD and the compressed internal
  tiles intersecting the window are read, so a full 90° tile never has to
  sit in RAM.
- `fetch_gebco` returns deterministic output for any bbox; cells are
  NaN only where the source grid has no data (not expected — GEBCO
  2024 has full global coverage, and every land/ocean point on Earth
  is inside some tile).
- Natural Earth vectors are cartographic context, not a data variable:
  `survey-viz` keeps them as an underlay and refuses "country borders"
  alone with an honest no-product response.

## 14. CSR GRACE / GRACE-FO RL06.3 — terrestrial water storage anomalies (no account)

The canonical satellite-gravimetry record of total water storage over
land (CSR mascon RL06.3). No account, keyless anonymous HTTPS from
the University of Texas Center for Space Research. The adapter
downloads two NetCDFs once and subsets locally.

| Item | Value |
|---|---|
| Directory | `https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/` — verified live 2026-09-27: anonymous HTTPS, `Accept-Ranges: bytes`, `Last-Modified: 2026-08-24` |
| Solution file | `CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc` (112,611,569 bytes); NetCDF-4/HDF5. `time`=258 months (2002-04 – 2026-06), `lat`/`lon`=720×1440 at 0.25°; main variable `lwe_thickness` in cm. The grid is already 0.25° — no mascon mapping needed — but the native mascon resolving power is coarser than 0.25°, so do not claim true 0.25° resolving power |
| Land-mask file | `CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc` (4,174,189 bytes, `Last-Modified: 2024-09-18`); variable `LO_val`, 0/1 on the same 720×1440 grid. The solution grid does **not** mask oceans itself (a Pacific sample held a valid number), so this file is mandatory: ocean cells are NaN |
| Anomaly baseline | `time_mean_removed = "2004.000 to 2009.999"` — every value is an anomaly against the 2004–2009 time-mean; units cm liquid-water-equivalent thickness |
| Access pattern | `ensure_grace_files()` — download-once into `$SURVEY_CURRENTS_CACHE/grace`, atomic writes, SHA-256 sidecars verified on every hit, corruption-triggered redownload, stale-file `Last-Modified` revalidation (CSR re-releases the file as new months arrive), `refresh=True` forces re-download |
| Query | `fetch_grace(bbox, start, end)` — monthly frames in `[start, end]`; longitudes normalized 0..360 → −180..180 (antimeridian and prime-meridian-crossing bboxes handled as one or two windows on the 0..360 axis) |
| Missing months | The file carries only observed solution months. Missing months are diffed from the file's actual time axis against the requested calendar — do **not** trust the file's `months_missing` attribute alone: it under-reports (2011-11 and 2015-05 are missing from the time axis but absent from the attribute). All 35 observed gaps (2017-07 … 2018-05 inter-mission gap plus 24 isolated months) become all-NaN frames listed in `WaterField.gap_months`; never interpolated, never silently dropped |
| Provenance | Product/version (CSR GRACE and GRACE-FO MASCON RL06.3M), solution+mask URLs, SHA-256 digests, anomaly baseline, requested bbox/time range, gap months, cache-hit flag, `retrieved_at`, tool version |
| Model | `WaterField` — `times` (month-start dates), `lats`, `lons` (−180..180 ascending), `values` (ntime, nlat, nlon) float; `gap_months`, `n_gap_months`, `is_gap()`, `spatial_mean()`, `select_time`, `select_bbox`, `to_dict`/`from_dict`, `to_json`/`from_json`, deterministic `synthetic()` (drying trend + seasonal cycle, ~35% ocean-masked cells, configurable gap months) |
| CLI | `fetch-grace --bbox … --start … --end …`, `grace-synthetic` |
| Optional dependency | netCDF4 is lazy/optional (the CSR archive is NetCDF-4/HDF5); parsing is pure over duck-typed datasets so tests and `synthetic()` work without it |

Operational notes:

- Record bounds: requests before 2002-04-01 raise; future end dates
  raise. CSR publishes new months roughly monthly — the 30-day
  `Last-Modified` revalidation keeps the cache current without
  re-downloading the ~107 MB file on every call.
- Interpretation: positive anomalies = wetter than the 2004–2009
  mean (blue), negative = drier (brown); GRACE measures *total* water
  storage (groundwater + soil moisture + snow + surface water), not
  groundwater alone — say so in captions. The `survey-viz`
  `water-storage` variable renders this field with a zero-centered
  diverging colormap and marks gap months as "NO GRACE OBSERVATION".

## 13. NOAA IBTrACS v4 — global tropical-cyclone best tracks (no account)

The canonical global best-track archive of tropical and subtropical
storms (IBTrACS v04r01, `since1980` subset 1980–present plus the full
1842–present `ALL` file). No account, keyless HTTPS from NCEI. The
adapter downloads the archive NetCDF once and subsets locally.

| Item | Value |
|---|---|
| Directory | `https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/netcdf/` — verified live 2026-09-27: anonymous HTTPS, `Accept-Ranges: bytes` |
| Files | `IBTrACS.since1980.v04r01.nc` (~10.8 MB, default) and `IBTrACS.ALL.v04r01.nc` (~23.4 MB, `--full-archive`); NetCDF-4/HDF5. Storm dimension 4991, per-storm fix slots 360; fixes nominally 3-hourly but cadence varies by agency/era |
| Fix variables | `usa_wind` → `wmo_wind` (kt, 1-min sustained), `usa_pres` → `wmo_pres` (hPa), `lat`, `lon`, `iso_time`, `sid`, `name`, `season`, `basin`, `usa_sshs`; integer fill −9999, `usa_sshs` fill −15 |
| Access pattern | `ensure_ibtracs_file()` — download-once, atomic write, SHA-256 sidecar verified on every hit, corruption-triggered redownload, stale-file `Last-Modified` revalidation (the archive is republished as storms are added), `refresh=True` forces re-download |
| Query | `fetch_ibtracs(bbox, start, end, min_wind=None, storm_name=None, full_archive=False)` — storms kept when ≥1 fix falls in bbox; fixes clipped to `[start, end]`; `min_wind` filters by lifetime max sustained wind (kt); `storm_name` selects one named storm (case-insensitive exact) |
| Dateline | Longitudes normalized to [−180, 180); tracks keep continuous values; renderers must break polylines on |Δlon| > 180° |
| Provenance | IBTrACS version, file URL, SHA-256, subset parameters, wind/pressure agency priority, `retrieved_at`, cache-hit flag |
| Model | `StormField` (`StormTrack` list) — `rank_by_intensity()`, `select_time`, `select_bbox`, `to_dict`/`from_dict`, `to_json`/`from_json`, deterministic `synthetic()` (includes a `TESTALPHA` track) |
| CLI | `fetch-ibtracs --bbox … --start … --end … [--storm-name katrina] [--min-wind 100] [--full-archive]`, `storms-synthetic` |

Operational notes:

- Record bounds: the default file starts 1980-01-01 (a request before
  that date raises and suggests `--full-archive`); the full file starts
  1842-01-01.
- Agency wind averaging conventions differ globally; Saffir-Simpson
  categories are formally based on 1-minute winds, so category mapping
  is most defensible against the USA-agency values — documented in
  `docs/STORMS.md` on the `survey-viz` side.

## 15. USGS Water Services (NWIS) — US streamgage daily values (no account)

The canonical US streamflow source: the USGS Water Services site
service discovers daily-data gages in a bbox (RDB), and the
daily-values service fetches per-site daily MEAN series for discharge
(00060, ft³/s) and gage height (00065, ft). Keyless anonymous HTTPS,
US-only coverage (Puerto Rico / Pacific territories included).

| Item | Value |
|---|---|
| Site service | `https://waterservices.usgs.gov/nwis/site/` — verified live 2026-09-27: keyless HTTPS, `?format=rdb&bBox=minLon,minLat,maxLon,maxLat&hasDataTypeCd=dv&parameterCd=00060&siteOutput=expanded`; `siteOutput=expanded` adds drainage area (`drain_area_va`, sq mi), HUC, timezone, altitude |
| Daily-values service | `https://waterservices.usgs.gov/nwis/dv/` — verified live 2026-09-27: keyless HTTPS, `?format=json&sites=…&startDT=…&endDT=…&parameterCd=00060,00065&statCd=00003`; daily MEAN stated explicitly; qualifiers `P` (provisional) / `A` (approved) / `e` (estimated) carried per value |
| Missing values | USGS sentinel `"-999999"` → NaN; days with no reported value are NaN in the series and listed by `GageRecord.missing_days()` — never filled, never interpolated |
| Access pattern | Per-bbox inventory + per-batch (20 sites) dv payloads cached under `$SURVEY_CURRENTS_CACHE/streamgages` with atomic writes, SHA-256 sidecars verified on every hit, corruption-triggered redownload, 7-day freshness revalidation (recent daily values are provisional and get revised), `refresh=True` forces re-download — repeated renders never re-hit the network |
| Query | `fetch_usgs(bbox, start, end, parameters=("00060",), min_record_days=30, site_limit=200, units="native")` — sites sorted by site number (deterministic); sites with fewer than `min_record_days` finite daily values of the primary parameter are excluded but counted (`n_sites_excluded_short_record`); sites returning no series are counted (`n_sites_no_data`, e.g. inventory sites with no record in the window) |
| Empty results | NWIS is US-only: bboxes outside USGS coverage (or with no qualifying gages) yield an empty `GageField` with `provenance["empty_reason"]` — never fabricated data |
| Units | Native USGS units by default (ft³/s, ft); `units="si"` (or `GageField.to_si()`) converts with exact NIST factors `1 ft³/s = 0.028316846592 m³/s`, `1 ft = 0.3048 m`; the unit system is recorded in provenance and on every `unit_label` |
| Provenance | Exact request URLs (site + each dv batch), `retrieved_at`, bbox, dates, parameters, statistic (`00003`), per-class site counts, cache status per entry, tool version |
| Model | `GageField` (`GageRecord` list) — `percentile_of_record()` (latest value's rank in its own record, the survey-viz marker-color rule), `regional_median()` (daily median across sites, NaN-aware, the regional-hydrograph rule), `select_site()`, `to_si()`, `to_dict`/`from_dict`, `to_json`/`from_json`, deterministic `synthetic()` (4 sites incl. an engineered 10-day gap block) |
| CLI | `fetch-usgs --bbox … --start … --end … [--parameters 00060,00065] [--min-record-days 30] [--site-limit 200] [--units si]`, `usgs-synthetic` |
| Interop | Consumed by the `usgs` source in survey-viz (streamflow planning); `GageField.to_dict()` feeds the dedicated gage-point renderer — point/time-series data, never converted through a scalar grid |

Operational notes:

- Site inventory can include stations with no values in a requested
  window (verified 2026-09-27: station 04160050 returns zero series for
  tested 2020/2025 windows) — those are counted, not errors.
- Recent daily values are provisional (`P`): cache freshness defaults
  to 7 days, not the 30-day discipline of static archives.
- The 20-site dv batches keep each cache entry small; large bboxes
  should be subdivided by the caller (the national inventory is tens of
  thousands of sites; default `site_limit=200`).

## 16. Ocean color (chlorophyll-a) — NOAA CoastWatch ERDDAP + NASA OBPG (mixed access)

The ocean-color source family: chlorophyll-a concentration in mg m^-3
from Level-3 satellite products, cloud-masked (NaN cells are kept as
NaN — never interpolated or filled). The longitude axis is -180..180
on the wire, so no wrap conversion is needed.

| Item | Value |
|---|---|
| Default source | NOAA CoastWatch ERDDAP (keyless HTTPS): `erdMH1chlamday_R2022SQ` — "Chlorophyll-a, Aqua MODIS, V.2022, Monthly, Global, 4km, Science Quality, 2002-present" (`sensor="modis-aqua"`, `cadence="monthly"`) |
| VIIRS SNPP | `nesdisVHNSQchlaMonthly/Weekly/Daily` — "Chlorophyll, NOAA S-NPP VIIRS, Science Quality, Global 4km, Level 3, 2012-present" (`sensor="viirs-snpp"`, `cadence="monthly"`; weekly and daily IDs also verified) |
| Multi-sensor merge | `pmlEsaCCI60OceanColorMonthly` — "ESA CCI Ocean Colour v6.0, Monthly, Global, 4km, 1997-present" (`sensor="multi"`) — the longest climate-grade record; ID verified 2026-09-27 against the official coastwatch-training/coastwatch-tutorials repo and an independent research repo, both citing it against the NOAA CoastWatch ERDDAP |
| NASA OBPG fallback (required) | `source="obpg"` — MODIS Aqua Level-3 mapped chlorophyll-a from the OBPG direct data access (`https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-3%20Mapped/Aqua-MODIS`, documented filenames `AQUA_MODIS.<start>_<end>.L3m.<DAY|8D|MO>.CHL.chlor_a.4km.nc`; monthly = calendar month, 8-day = OBPG day-of-year bins 1-8/9-16/…, daily = single date). Needs a free Earthdata Login — without it the adapter raises the credentials-gated `CredentialsMissing` (never an interactive prompt). CMR collection `C3380709133-OB_CLOUD` (`MODISA_L3m_CHL`, v2022.0, 2002-07-04–present) verified live via the keyless CMR collections API 2026-09-27. `source="auto"` tries CoastWatch, then OBPG when Earthdata credentials exist, then CMEMS when CMEMS credentials exist |
| Griddap URL | `https://coastwatch.noaa.gov/erddap/griddap/{dataset}.nc?chlor_a[({start}):{stride}:({end})][(0.0)][({miny}):({maxy})][({minx}):({maxx})]` — the singleton `altitude` axis is pinned at `(0.0)` exactly like the OISST `zlev` axis; bbox is server-side subset |
| Grid geometry | ~0.0417° nominal (~4 km); lon -179.9812..179.9813, lat -89.75626..89.75625 (≈ 4320×8640 cells); variable `chlor_a`, units `mg m^-3`, `_FillValue -999.0` → masked |
| Download-size math | A bbox subset pulls only the requested cells at 4 bytes/cell/frame: a 20°×20° monthly bbox is ~480×480 cells ≈ 0.9 MB/frame (≈ 11 MB for a year of months). The adapter never downloads the full 4320×8640 global frame (~150 MB float32/frame) unless the bbox spans the globe |
| Cloud gaps | Chlorophyll is cloud-sensitive: daily/weekly frames routinely have large NaN regions; monthly is the default cadence because it is the most cloud-complete. Per-frame gap fractions are recorded (`provenance["gap_fractions"]`, `OceanColorField.gap_fraction()`); NaN cells are rendered as "no observation", never filled |
| Access pattern | Payloads cached under `$SURVEY_CURRENTS_CACHE/oceancolor` with atomic writes, SHA-256 sidecars verified on every hit, corruption-triggered redownload, 7-day freshness (`refresh=True` forces re-download); cache key covers the full request URL (dataset + bbox + dates + cadence + stride) |
| Query | `fetch_oceancolor(bbox, start, end, product="chlorophyll-a", cadence="monthly", sensor="modis-aqua", source="coastwatch")`; `oceancolor_urls(...)` builds the exact ERDDAP URLs offline; `obpg_file_names(...)` builds the exact OBPG filenames offline |
| Model | `OceanColorField` (`.values`/`chl` `(nt, ny, nx)` masked array in mg/m³, `.times`/`.lats`/`.lons`, `spatial_median()` (NaN-aware median — the standard central tendency for lognormal chlorophyll), `select_time()`, `select_bbox()`, `to_dict`/`from_dict`, `to_json`/`from_json`, deterministic lognormal `synthetic()` with an engineered cloud mask) |
| CLI | `fetch-oceancolor --bbox … --start … --end … [--product chlorophyll-a] [--cadence monthly|weekly|daily] [--sensor modis-aqua|viirs-snpp|multi] [--source coastwatch|obpg|cmems|auto] [--stride N] [--refresh]`, `oceancolor-synthetic` |
| Interop | Consumed by the `oceancolor` source in survey-viz (`ocean-color` variable, log-scaled renderer); `OceanColorField.values` duck-types into the scalar-grid path |

Access status (recorded 2026-09-27):

- The CoastWatch ERDDAP was unreachable from the build environment
  during verification: `https://coastwatch.noaa.gov/erddap/` returned
  HTTP 503 and griddap probes returned 502/503 (the legacy
  `coastwatch.pfeg.noaa.gov` host is decommissioned — empty replies).
  Dataset IDs, grid geometry, and the griddap query grammar above are
  corroborated from the CoastWatch tutorials repo, the CoastWatch WMS
  dataset list, and NOAA sanctuary caption pages; the identical
  `.das` → griddap-subset request grammar was exercised successfully
  against a live ERDDAP server (NCEI) in the same session.
- A successful end-to-end CoastWatch data pull still needs a re-run
  once the service recovers — the adapter is written so it does
  exactly that.
- The OBPG fallback could not be live-verified in the build
  environment: no Earthdata credentials were available, and the
  CoastWatch outage meant the fallback chain was never exercised
  end-to-end. The filename convention was verified against real CMR
  granule titles (e.g. `AQUA_MODIS.20020704_20250228.L3m.CU.CHL.chlor_a.4km.nc`)
  and the CMR collection record; the parser is defensive against the
  documented SeaDAS L3 mapped structure (2-D/3-D `chlor_a`,
  `scale_factor`/`add_offset`, `_FillValue`) and is unit-tested against
  synthetic SeaDAS-shaped files. Provenance records `live_verified:
  false` until a credentialed run succeeds.
- The CMEMS path was not live-tested (no CMEMS credentials in the
  build environment); its dataset ID is taken from the Nov-2025
  CMEMS ocean-colour QUID
  (CMEMS-OC-QUID-009-101to104-111-113-116-118) and can be overridden
  via `cmems_dataset_id`; it is an optional extra source, not the
  required fallback.

## 17. USGS Earthquake Catalog (ComCat) — global seismic events (no account)

The USGS ComCat FDSN event service: a catalog of **observed**
seismic events (earthquakes plus non-tectonic types like quarry
blasts), served as GeoJSON over keyless HTTPS. It is **not** an
earthquake forecast or hazard model — say so explicitly wherever it
is surfaced.

| Item | Value |
|---|---|
| Endpoints | Query: `https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&starttime=…&endtime=…&minlatitude=…&maxlatitude=…&minlongitude=…&maxlongitude=…&minmagnitude=…[&eventtype=…][&limit=…&offset=…]`; count: `https://earthquake.usgs.gov/fdsnws/event/1/count` (same params minus `format`; returns a plain-text integer) |
| Access | Keyless HTTPS; no token, no signup. Live-verified 2026-09-27: global M4+ for 2026-09-20…2026-09-27 returned 179 events; M0+ returned 1,906; a California M4+ window returned a legitimate empty catalog; paged queries (`limit=50&offset=1`) returned 50 features. No numeric rate limit observed (no rate/quota/retry headers seen); be a polite client anyway |
| Pagination | Paged-response metadata carries `limit`/`offset`, not `count` (verified live) — the adapter stops on a short page, never on a count key. Queries whose count exceeds the 20000-event FDSN ceiling are recursively time-split until every chunk is page-sized; a single day that still exceeds the ceiling raises a clear error |
| Event fields | `id`, epoch-ms `time` (UTC), `mag` (None when unreported — never 0.0), `magType`, `place`, `type` (event type), geometry coordinates lon/lat/depth-km (2-D coordinates → depth None). Events deduped by `id`, sorted by time |
| Magnitude completeness | Varies by region and time: roughly M4.5+ globally since ~1973 and M2.5+ in the contiguous US since ~2013; small and historical events are under-recorded. Record floor 1900-01-01 |
| Access pattern | Payloads cached under `$SURVEY_CURRENTS_CACHE/comcat` (atomic writes, SHA-256 sidecars, corruption-triggered redownload, 7-day freshness, `refresh=True` forces re-download); cache key covers the full request URL |
| Query | `fetch_earthquakes(bbox, start, end, min_magnitude=0.0, event_type=None, page_size=2000)`; `event_query_url(...)` / `event_count_url(...)` build the exact URLs offline |
| Model | `QuakeField` (`.events`, `largest(n)`, `select_time()`, `counts_by_day()` for daily cumulative frames, `magnitudes()`, `to_dict`/`from_dict`, `to_json`/`from_json`, deterministic `synthetic()` with a Gutenberg–Richter-ish magnitude taper); `QuakeEvent.depth_bin()` = shallow <70 km / intermediate 70–300 km / deep >300 km |
| CLI | `fetch-earthquakes --bbox … --start … --end … [--min-magnitude 0.0] [--event-type earthquake] [--page-size 2000] [--out …]`, `earthquakes-synthetic` |
| Interop | Consumed by the `comcat` source in survey-viz (`earthquakes` variable, magnitude-scaled/depth-colored event renderer with cumulative daily frames) |

## 18. NOAA GFS 0.25° — 10-m winds + 2-m air temperature (no account)

The Global Forecast System f000 **analysis** snapshots via NCEP's
NOMADS GRIB filter (`filter_gfs_0p25.pl`) — **keyless HTTPS, no
account, no token**. This is the program's keyless wind source for
unattended automation (OSCAR, CMEMS, ERA5/CDS, and FIRMS all need
credentials).

| Item | Value |
|---|---|
| Endpoint | `https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=/gfs.YYYYMMDD/HH/atmos&file=gfs.t{HH}z.pgrb2.0p25.f000&subregion=on&leftlon=..&rightlon=..&toplat=..&bottomlat=..&lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on&lev_2_m_above_ground=on&var_TMP=on` |
| Access | Keyless HTTPS. Live-verified 2026-10-01: HTTP 200 for 2026-09-22…2026-10-01 (10 days), HTTP 404 for 2026-09-21. NOMADS retired the DODS/OpenDAP interface in 2025 (NWS SCN 25-81) — the GRIB filter CGI is the working keyless path |
| Cycles | `00`/`06`/`12`/`18` UTC analyses; `f000` is the analysis, not a forecast |
| Variables | `10u`/`10v` → `u10`/`v10` 10-m wind components (m/s); `2t` → `t2m` 2-m air temperature (Kelvin on the wire, °C on the field). `GfsWindField.air_temperature` exposes °F (the warming.watch strand-color convention survey-viz's `dark_strands` reads) |
| Subsetting | `subregion=on` makes the filter honor the bbox (verified 2026-10-01: without it the filter silently returns the full 1440×721 global grid). The adapter still verifies the returned grid covers the requested window and crops locally with numpy; antimeridian-crossing bboxes split into two subregion requests and concatenate along longitude |
| Parsing | `cfgrib.messages.FileStream` iteration — `cfgrib.open_file` fails on these payloads (mixed 10-m/2-m levels break its dataset builder). Lazy import: `pip install 'survey-currents[gfs]'` |
| Access pattern | Payloads cached under `$SURVEY_CURRENTS_CACHE/gfs-wind` (atomic writes, SHA-256 sidecars; analyses are immutable once posted, so entries never expire); cache key covers the full request URL |
| Retention | ~10 days (`GFS_RETENTION_DAYS`, verified 2026-10-01). Dates outside the window and not-yet-posted cycles raise `UnavailableRangeError` — never silent padding. A message whose data date/cycle mismatches the request is refused rather than mislabeled |
| Per-day size | KB-scale for regional bboxes (a North-America box is ~100–200 KB/day — the subregion request, not the 2.4 MB full-globe fallback) |
| Query | `fetch_gfs_wind(bbox, start, end, stride_days=1, cycle="00", work_dir=...)`; `gfs_filter_url(day, cycle, window)` builds the exact URL offline |
| Model | `GfsWindField` (Era5Field-shaped: `grids`/`times`/`lats`/`lons`, `values` = wind speed, `overlay_grids` = `{"t2m": …}`, plus `air_temperature`/`temperature_unit` for the viz wind path; `select_time()`, `select_bbox()`, JSON round-trip, deterministic `synthetic()`) |
| CLI | `fetch-gfs-wind --bbox … --start … --end … [--stride-days 1] [--cycle 00] [--work-dir …] [--out …]`, `gfs-wind-synthetic` |
| Interop | Consumed by the `gfs-wind` source in reel-studio (`wind` variable, keyless); survey-viz's `dark_strands`/`dark_flow` presets read `grids["u10"]`/`grids["v10"]` and color wind strands by `air_temperature` (°F) |

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
- **Global precipitation — observed satellite, free Earthdata account:**
  NASA GPM IMERG V07 — half-hourly 0.1° precipitation, 2000–present,
  Early/Late/Final latency tiers; the answer for "recent rain",
  "storm", "hurricane rainfall" visualizations. Prefer Late (default)
  over Early for quality; use the Final run for research-grade,
  gauge-adjusted totals (see `currents.imerg`). For climatology,
  multi-decade trends, or pre-2000 analysis, ERA5 `tp` remains the
  answer (1940–present reanalysis).
- **Global atmosphere (wind, pressure, air temperature, long-record
  precipitation), free CDS account:** Copernicus ERA5 — hourly 0.25°
  reanalysis, 1940–present; the answer for "winds", "heat", and
  multi-decade precipitation trends (see `currents.era5`). For
  *observed* recent precipitation, IMERG above is the better source.
- **Global surface currents, free Earthdata account:** NASA PODAAC
  OSCAR v2.0 — daily 0.25° currents, 1993–present, with Final/Interim/NRT
  latency tiers picked per date; the answer for "currents" visualizations
  outside the NOAA OFS footprints (see `currents.currents_global`).
  **Corrections:** not on CoastWatch ERDDAP (old `jplOscar_LonPM180`
  id 404s — removed); NOMADS OPeNDAP is retired (SCN 25-81).
- **Global high-resolution currents, free CMEMS account:** CMEMS
  `global-physics-daily` preset wrapped as `fetch_cmems_currents`
  (1/12°, daily; `thetao` carried as-is).
- **Any other coastline:** CMEMS global physics.
- **Tropical-cyclone tracks, identity, and season rankings, no
  signup:** NOAA IBTrACS v04r01 — 1980–present best tracks
  (1842–present with `--full-archive`); the answer for named-storm
  tracks, "strongest hurricanes of the season" rankings, and storm
  identity. For the *ambient* storm-environment fields (wind and
  pressure patterns around a storm), ERA5 remains the answer
  (see `currents.storms`).
- **Global 10-m winds + 2-m air temperature, no signup (recent ~10
  days):** NOAA GFS 0.25° f000 analyses via the keyless NOMADS GRIB
  filter — the answer for unattended wind-strand automation, where the
  credentialed sources above cannot run. For long-record or
  multi-decade wind climatology, ERA5 remains the answer (1940–present
  reanalysis).
- **Terrestrial water storage / groundwater anomalies, no signup:**
  CSR GRACE/GRACE-FO RL06.3 — monthly TWS anomalies in cm LWE,
  2002–present, land-only; the answer for "water storage",
  "groundwater", "aquifer", "drought", and "total water storage"
  descriptions. For *rainfall/precipitation* use GPM IMERG or ERA5
  (`tp`) — never GRACE. **For river discharge or streamflow the answer
  is USGS NWIS streamgages** (see section 15): daily observed discharge
  (00060, ft³/s) and gage height (00065, ft) from the per-site gage
  network — use GRACE for basin-scale storage context alongside it, not
  as a substitute. For sea level there is no adapter — it is a
  different observable from GRACE TWS.
- **Blending with satellites:** use `align_to_thermal_zone()` to compare
  model water temperature against survey-thermal Landsat LST passes over
  the same zone and window.

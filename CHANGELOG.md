# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.10.0] - 2026-09-27

### Added
- GEBCO 2024 + Natural Earth basemap adapter (`src/currents/basemaps.py`):
  `fetch_gebco(bbox, resolution="15s")` -> `TopoField` and
  `fetch_naturalearth(bbox, scale="110m", layers=("coastline", "countries"))`
  -> dict of GeoJSON FeatureCollections. **Verified live 2026-09-27:**
  the keyless open-download route is
  `https://www.bodc.ac.uk/data/open_download/gebco/gebco_2024/geotiff/`
  (redirects to a CEDA-hosted `gebco_2024_geotiff.zip`, ~4.26 GB,
  anonymous HTTPS, byte-range requests honored). The zip holds eight
  90°×90° stripped uncompressed int16 GeoTIFFs (~495–577 MB
  compressed / 933,255,948 bytes uncompressed each, 21600×21600 px at
  15 arc-second; layout verified against a real tile 2026-09-27:
  21600 one-row strips, positive Y pixel scale per the GeoTIFF spec).
  The adapter parses the zip central directory with a
  few KB of range requests and downloads **only the tile entries
  intersecting the bbox** — no 4.26 GB download. Tiles are cached
  locally with SHA-256 sidecars (verified on every hit, re-fetched on
  corruption). A minimal stdlib GeoTIFF reader (file-backed:
  only the IFD and intersecting strips/tiles are read; handles the
  official stripped-uncompressed layout as well as tiled deflate)
  subsets natively and block-averages to coarser grids; output grids
  snap to the 15 arc-second lattice and mosaic with exact integer
  index math (antimeridian wraps handled in unwrapped longitude
  space), so tiles butt-join without seams. `resolution="15s"` or a
  degree value that is a multiple of 15 arc-seconds and divides the
  90° tiles (`0.25`, `0.5`, `1.0`, …); native resolution has an 8192
  px safety cap per axis. `TopoField` carries one structural timestamp
  (`2024-01-01T00:00:00Z`, the grid's release year) so existing
  renderers keep their 3D-field contract, with
  `provenance["static_compilation"] = True`, the grid version, the
  requested subset bbox, and the snapped output grid.
  Natural Earth zips come from anonymous S3
  (`https://naturalearth.s3.amazonaws.com`, verified live 2026-09-27),
  cached with the same SHA-256 discipline; shapefiles are parsed with
  a minimal stdlib reader (no fiona/geopandas) and clipped to the
  bbox. 38 offline tests (`tests/test_basemaps.py`).
- CLI: `fetch-gebco`, `fetch-naturalearth`, `basemaps-synthetic`
  (deterministic synthetic topography, offline).

## [0.9.0] - 2026-09-27

### Added
- NASA Black Marble night-lights adapter (`src/currents/blackmarble.py`):
  `fetch_blackmarble(bbox, start, end, product="daily", stride_days=1,
  resolution=0.05)` -> `LightsField`. **Verified live 2026-09-26:**
  VNP46A2 V002 (daily gap-filled lunar BRDF-adjusted nighttime lights,
  CMR collection `C3365931269-LAADS`, record starts 2012-01-19). V002
  keeps the `hHHvVV` tile naming but the tiles are **10°×10° lat/lon**
  tiles (the V2 "15 arc-second linear lat/lon grid"), not sinusoidal:
  `h = floor((lon + 180) / 10)`, `v = floor((90 - lat_top) / 10)` —
  verified against live CMR footprints (`h07v10` = lon −110..−100, lat
  −20..−10). Filenames carry an unpredictable production stamp, so the
  adapter **discovers** exact download URLs per (day, tile) through the
  **keyless** NASA CMR granule search, then downloads from LAADS with
  Earthdata Login (`EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` or
  `.netrc`; `CredentialsMissing` otherwise — both LAADS endpoints were
  verified to redirect unauthenticated requests to URS OAuth).
  - Science dataset `Gap_Filled_DNB_BRDF_Corrected_NTL` (fallback
    `DNB_BRDF_Corrected_NTL`) in nW/cm²/sr; `_FillValue` → NaN;
    `scale_factor`/`add_offset` honored; tiles mosaicked and NaN-aware
    block-averaged to the requested `resolution` (default 0.05°).
  - `LightsField` follows the `IceField`/`RainField` conventions
    (`times`/`lats`/`lons`/`values`, `select_time`/`select_bbox`,
    JSON round-trip, deterministic `synthetic()`, per-file SHA-256
    provenance with exact URLs, tile list, skipped tiles/days).
  - CLI: `fetch-blackmarble`, `lights-synthetic`; extra
    `survey-currents[blackmarble]` (h5py; also in `full`).
  - 39 new tests, fully offline (mocked CMR/discovery/downloads,
    in-memory HDF5 fixtures). 417 passed, 22 pre-existing failures
    (missing netCDF4/xarray in this environment — identical on the
    pristine v0.8.0 tree), 3 skipped.
- `docs/DATA_SOURCES.md` §11 documents the product, the verified
  access truth, the tile math, and why monthly/annual composites stay
  unwired (daily is the reel factory's time axis).

## [0.8.0] - 2026-09-27

### Added
- NASA GPM IMERG precipitation adapter (`src/currents/imerg.py`):
  `fetch_imerg(bbox, start, end, accumulate="daily", run="late",
  stride_days=1)` -> `RainField`. **Verified live 2026-09-27:**
  the V07 half-hourly granules live under
  `https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/{GPM_3IMERGHH.07,
  GPM_3IMERGHHE.07, GPM_3IMERGHHL.07}/{YYYY}/{DOY}/` with filenames like
  `3B-HHR[-E|-L].MS.MRG.3IMERG.{YYYYMMDD}-S{HHMMSS}-E{HHMMSS}.{HHMM}.V07B.HDF5`
  (Early/Late carry the `-E`/`-L` infix); the on-wire HDF5 dataspace is
  `(time=1, lon=3600, lat=1800)` per the `.dds`/`.das`, units mm/hr,
  `_FillValue` ≈ -9999.9, 0.1° global grid; lat/lon axes are resolved
  from the actual `Grid/lat`/`Grid/lon` vector sizes (honest `ValueError`
  on layout drift, never a silent transpose). OPeNDAP catalog and
  `.dds`/`.das` metadata are **keyless**, but data download redirects
  (302 → 401) to a free Earthdata Login — the adapter authenticates via
  `EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` (or `.netrc`) exactly like
  the MUR path and raises `CredentialsMissing` otherwise.
  - Three latency runs: `early` (~4 h, `IMERGHHE`), `late` (~14 h,
    `IMERGHHL`, the default), `final` (~3.5 months, `IMERGHH`,
    gauge-adjusted); the documented V07 record starts 2000-06-01 for
    all three runs (TRMM June 2000 - May 2014, GPM June 2014 - present;
    enforced per run with honest errors). The live CMR archive shows
    V07 reprocessing reaching 1998-01-01 for Early/Final, but the
    adapter keeps the documented floor.
  - `accumulate="daily"` (default) sums the 48 half-hourly mm/hr rates
    into mm/day totals (rate × 0.5 h per slot; NaN treated as 0 for that
    slot, all-NaN days stay NaN); `accumulate="native"` keeps the 48
    half-hourly mm/hr rate steps. Incomplete days accumulate the slots
    that were retrieved and record
    `provenance["day_coverage"][day] = {"expected": 48, "retrieved": n}`;
    404 granules are skipped with a `skipped_slots` note; an empty
    fetch raises instead of returning an empty field.
  - `RainField`: gridded `times`/`lats`/`lons`/`values` (mm/day or
    mm/hr), `select_time`/`select_bbox`, JSON round-trip, deterministic
    `synthetic()`, per-file SHA-256 provenance with exact download
    URLs, retrieval timestamp, run, and accumulation mode; antimeridian
    bboxes wrap into two index windows and come back sorted -180..180.
  - CLI: `fetch-imerg`, `rain-synthetic`.
  - New `imerg` extra: `pip install "survey-currents[imerg]"` (h5py);
    also included in `full`. The engine core stays stdlib+numpy.

## [0.7.0] - 2026-09-26

### Added
- NSIDC Sea Ice Index adapter (`src/currents/sea_ice.py`):
  `fetch_nsidc_sic(bbox, start, end, hemisphere="auto", stride_days=1,
  resolution=0.25)` -> `IceField`. **Verified live 2026-09-26:**
  `https://noaadata.apps.nsidc.org/NOAA/G02135/{north,south}/daily/geotiff/{YYYY}/{MM_Mon}/{N,S}_YYYYMMDD_concentration_v4.0.tif`
  serves daily concentration GeoTIFFs over **keyless anonymous HTTPS**
  (no Earthdata login, no key — directory listings and real files
  return HTTP 200 anonymously). **Correction from the plan:** G02135
  v4 uses the NSIDC polar stereographic grid (EPSG:3411 North,
  EPSG:3412 South, Hughes 1980 ellipsoid, 25 km), not EASE-Grid.
  - stdlib+numpy only: a minimal reader for the uncompressed
    single-band 16-bit GeoTIFFs (`read_concentration_geotiff`, tags
    with fallback to the documented grid constants); concentration is
    ×10-scaled (÷10 → percent 0-100); flags 2510/2530/2540/2550
    (pole hole/coast/land/missing) → NaN.
  - Nearest-neighbor reprojection of the native grids onto the
    caller's regular lat/lon grid (documented choice — preserves the
    15%-threshold ice edge; the formulas are verified against a live
    file: the pole maps to pixel (154, 234) where the pole-hole flag
    sits). Southern target lats are masked off the north grid and vice
    versa (the |lat|-symmetric formulas would otherwise alias rings
    across hemispheres).
  - `IceField`: gridded percent model, `select_time`/`select_bbox`,
    JSON round-trip, deterministic `synthetic()` fixture, per-file
    SHA-256 provenance; missing days (SMMR every-other-day era,
    1987-12-03 → 1988-01-13 outage) are skipped with a provenance
    note — record starts 1978-11-01.
  - CLI: `fetch-nsidc`, `ice-synthetic`.

## [0.6.0] - 2026-09-26

### Added
- NASA FIRMS active-fire adapter (`src/currents/fires.py`):
  `fetch_firms(bbox, start, end, instruments=("VIIRS_SNPP",))` ->
  `FireField`. **Verified live 2026-09-26:** the area API shape is
  `https://firms.modaps.eosdis.nasa.gov/api/area/csv/{MAP_KEY}/{PRODUCT}/{W},{S},{E},{N}/{DAY_RANGE}/{DATE}`
  (DAY_RANGE 1–5 with a DATE start), and a free MAP_KEY is **required**
  (request at https://firms.modaps.eosdis.nasa.gov/api/map_key/; a bad
  key returns `Invalid MAP_KEY.`). Without a key, `CredentialsMissing`
  explains the setup (same pattern as the MUR/ERA5 adapters).
  - Per-date latency-tier picking via the pure, offline-testable
    `firms_product_for` (dates within 60 days of today -> `*_NRT`;
    older -> `*_SP`): `VIIRS_SNPP_NRT/_SP` (2012–present, 375 m),
    `VIIRS_NOAA20_NRT/_SP` (2018–present), `VIIRS_NOAA21_NRT/_SP`
    (2023–present), `MODIS_NRT/_SP` (Nov 2000–present). Long ranges
    loop in <=5-day windows (`firms_window_chunks`); antimeridian
    boxes split into two requests; seam duplicates dropped.
  - Both CSV schemas parsed (VIIRS `bright_ti4` + letter confidence;
    MODIS `brightness`/`bright_t31` + numeric confidence), normalized
    to brightness (K), FRP (MW), confidence, satellite, instrument,
    day/night. API error payloads raise `ValueError` with the API's
    message.
  - `FireField`: point-detection model (parallel `times`/`lats`/`lons`
    + `brightness`/`frp`/`confidence`/`satellite`/`instrument`/
    `daynight`), `select_time`/`select_bbox`, JSON round-trip,
    deterministic `synthetic()` fixture. `to_density_grid(resolution,
    frp_weighted)` bins detections into daily fire-count (or
    FRP-weighted MW) grids — the plain `times`/`lats`/`lons`/`values`
    dict `survey-viz` renders with zero changes.
  - Provenance: the MAP_KEY is **never stored** — request URLs are
    recorded with a `<redacted>` placeholder, plus per-request SHA-256,
    retrieval timestamp, products used, and request count.
  - CLI: `fetch-firms`, `fires-synthetic`.
- 40 new fully-offline tests (`tests/test_fires.py`: sample VIIRS/MODIS
  CSV fixtures, faked HTTP, explicit `today` for the tiering rule).

### Docs
- `docs/DATA_SOURCES.md` §8 (FIRMS endpoint, key setup, product
  families, operational notes), README (FireField model, CLI examples,
  layout), CHANGELOG.

## [0.5.0] - 2026-09-26

### Added
- OSCAR v2.0 global surface-current adapter (`src/currents/currents_global.py`):
  NASA PODAAC OSCAR v2.0 (Ocean Surface Current Analyses Real-time) —
  daily-averaged surface currents, 1993–present, 0.25° global grid,
  variables `u`/`v` (m/s). **Verified-source corrections (2026-09-26,
  also recorded in provenance):** OSCAR is NOT on CoastWatch ERDDAP
  (the old `jplOscar_LonPM180` dataset id 404s — removed) and NOMADS
  OPeNDAP is retired (Service Change Notice 25-81); v2.0 is served
  from Earthdata OPeNDAP behind a free Earthdata Login (unauthenticated
  requests 302 to the login page).
  - `fetch_oscar(bbox, start, end, stride_days=5)` -> `CurrentField`
    (`temperature=None` — OSCAR is currents-only). Per-date latency-tier
    picking via the pure, offline-testable `oscar_collection_for`
    (date < today−540d → Final `C2098858642-POCLOUD`; < today−45d →
    Interim `C2102959417-POCLOUD`; else NRT `C2102958977-POCLOUD`;
    floor 1993-01-01). Granule names are deterministic
    (`oscar_currents_{final,interim,nrt}_YYYYMMDD.nc`) AND verified via
    the keyless NASA CMR granule search (`cmr_search_oscar_granules` /
    `oscar_match_granule` — a missing granule is an honest
    `RuntimeError`, never an invented name).
  - On-wire dimension order is the unusual `(time, longitude, latitude)`
    — asserted at parse time and transposed to (nt, ny, nx). Grid is
    0–360 (lon 0..359.75, lat −89.75..89.75); −180..180 bboxes convert
    internally (same approach as OISST), antimeridian boxes wrap, and
    360°-crossing windows split into two OPeNDAP requests concatenated
    with the seam deduplicated (`oscar_lon_windows`,
    `oscar_index_windows`, `oscar_subset_urls`).
  - `CredentialsMissing` (OSCAR-specific text) reuses the shared
    Earthdata env/netrc credential lookup; a 401 from the server maps
    to the same error. Fill value −999.0 masked; frame timestamps are
    the daily granule dates at 00:00 UTC.
  - Provenance: exact OPeNDAP URLs, per-payload SHA-256 (combined),
    retrieval time, the three CMR collection ids, the pick rule, and
    the per-date collection + granule title used.
- `fetch_cmems_currents(bbox, start, end, stride_days=1, work_dir=None)`:
  the CMEMS `global-physics-daily` preset (`uo`/`vo`/`thetao`, 1/12°,
  daily) wrapped to the standard fetch signature via
  `currents.cmems.subset_cmems` + `parse_cmems_netcdf` — one NetCDF
  downloaded for the whole range, timesteps stride-selected.
  `thetao` (potential temperature, °C) is carried as-is and documented
  in provenance (it is not a foundation SST).
- Synthetic fixtures `oscar_synthetic()` / `cmems_currents_synthetic()`
  (deterministic, stdlib+numpy) and CLI commands `fetch-oscar`,
  `fetch-cmems-currents`, `currents-synthetic`.
- 45 new offline tests (fake NetCDF payloads with the real
  (time, longitude, latitude) dim order; mocked CMR search, download,
  credentials, and toolbox — no live network in the suite).
- Live authenticated OSCAR fetch not yet verified — needs his Earthdata
  Login credentials; live CMEMS subset not verified either (needs the
  copernicusmarine toolbox + CMEMS account).

## [0.4.0] - 2026-09-26

### Added
- ERA5 adapter (`src/currents/era5.py`): Copernicus ERA5 hourly single-level
  reanalysis (1940–present, 0.25° global) via the `cdsapi` package (lazy
  import — the engine imports cleanly without it; free CDS account required,
  `CredentialsMissing` carries the exact setup steps).
  - `fetch_era5(variables, bbox, start, end, stride_hours=6)` with short
    variable keys: `wind` (10-m u/v pair, rendered as wind speed), `msl`
    (mean sea-level pressure, Pa→hPa), `t2m` (2-m air temperature, K→°C),
    `tp` (total precipitation, m→mm per hourly step). Requests chunked by
    calendar month (CDS limits); antimeridian-crossing bboxes split into two
    requests and concatenated with the seam meridian deduplicated.
  - New `Era5Field` model: `grids` dict of `(nt, ny, nx)` masked arrays,
    `times`/`lats`/`lons` (−180..180, increasing), `values` property (the
    render-ready base grid), `overlay_grids` (non-base variables for contour
    overlays, e.g. isobars), `select_time`/`select_bbox`, JSON round-trip,
    `Era5Field.synthetic()`. Duck-types into `survey-viz`'s `render_viz`.
  - Provenance: CDS dataset id, exact request dicts, per-payload SHA-256
    (combined), byte counts, retrieval time, unit-conversion notes.
  - CLI: `fetch-era5` (`--variables`, `--bbox`, `--start`, `--end`,
    `--stride-hours`), `era5-synthetic`.
  - 45 new offline tests (fake `cdsapi` module via `sys.modules`; no live
    CDS calls in the suite). Live CDS fetch not yet verified — needs the
    user's CDS credentials.
- Docs: `docs/DATA_SOURCES.md` gained §6 (ERA5 request shape, variables,
  credentials, antimeridian handling); README adapter table + `Era5Field`
  section + Python/CLI examples.

## [0.3.0] - 2026-09-26

### Added
- Global-SST adapter (`src/currents/sst_global.py`): two new satellite
  sea-surface temperature sources returning the new `SstField` model
  (`sst` as `(nt, ny, nx)` °C masked array, `times`/`lats`/`lons`,
  provenance with exact URLs + SHA-256 + retrieval time — mirrors the
  `GlseaField` conventions and duck-types into `survey-viz`'s
  `render_viz`).
  - `fetch_oisst(bbox, start, end, stride_days=30)` — NOAA OISST v2.1
    (`ncdcOisst21Agg` griddap on CoastWatch ERDDAP): daily SST,
    1981–present, 0.25° global grid. Keyless. Handles the 0–360
    longitude axis internally (bboxes stay conventional −180..180;
    antimeridian-crossing boxes wrap; boxes whose 0–360 window crosses
    360° are fetched as two requests and concatenated — never silently
    truncated). Long windows are chunked into ≤5-year requests
    (ERDDAP drops very long time ranges).
  - `fetch_mur(bbox, start, end, stride_days=30)` — NASA JPL MUR v4.1
    (GHRSST L4, CMR collection `C1996881146-POCLOUD`) via the
    Earthdata OPeNDAP endpoint: daily SST, 2002–present, ~0.01° (~1 km)
    global grid. **Granules are discovered, not constructed:**
    `cmr_search_mur_granules` queries the public NASA CMR granule API
    (keyless) and `mur_match_granules` pairs each sampled day with a
    real granule; service URLs prefer the CMR-advertised OPeNDAP link
    (`_mur_service_url`, recorded per granule). A day with no
    discovered granule raises `RuntimeError` — never an invented name.
    Grid index windows come from each granule's own
    `.das`/`.dds` metadata (no hardcoded grid); `analysed_sst` is
    converted Kelvin→°C. Needs a free Earthdata Login account
    (`EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` or `~/.netrc`); missing
    or rejected credentials raise `CredentialsMissing` with setup
    instructions, while OISST stays fully usable without credentials.
- CLI: `fetch-oisst`, `fetch-mur`, `sst-synthetic` (offline
  deterministic fixture, `--source` label).
- 48 new fully-offline tests (mocked HTTP, fabricated NetCDF
  payloads): URL construction, 0–360 window logic, antimeridian
  handling, bbox/date/stride validation, credentials handling,
  OISST depth-dimension indexing, parsing/masking, Kelvin→Celsius,
  serialization round-trips.

## [0.2.0] - 2026-09-26

### Added
- NOAA GLSEA adapter (`src/currents/glsea.py`): satellite-derived daily
  sea-surface temperature from `GLSEA_ACSPO_GCS` (ERDDAP griddap,
  2006–present, Great Lakes only) and lake-wide daily average
  temperatures from `glsea_avgtemps_3` (ERDDAP tabledap).
- `GlseaField`: time-indexed SST grids (°C, numpy masked array),
  `lats`/`lons`/`times`, source metadata, provenance (`url`, `sha256`
  of downloaded bytes, `retrieved_at`); `select_time()`,
  `select_bbox()`, masked/NaN-aware `spatial_mean()`, JSON round-trip,
  deterministic `synthetic()` fixture — mirrors the `CurrentField`
  conventions.
- `LakeSeries`: lake-wide average temperature series (`dates` + `temps`
  lists, °C), parsed with stdlib `csv`; unknown lake names raise
  `ValueError` (mapping: superior→Sup, michigan→Mich, huron→Huron,
  erie→Erie, ontario→Ont).
- Grid-extent validation up front: the GLSEA longitude axis is clipped
  to the lakes region (floor −92.4199507342304); out-of-extent bboxes
  and non-ascending constraints raise clear `ValueError` before any
  download.
- Downloads use stdlib `urllib` with retry on transient drops;
  `netCDF4` is a lazy import with an actionable `pip install netCDF4`
  error raised only when parsing is actually attempted.
- CLI: `fetch-glsea-sst`, `fetch-glsea-averages`, `glsea-synthetic`
  (offline); `examples/glsea_demo.py` offline demo.
- Docs: README, `docs/ARCHITECTURE.md`, `docs/DATA_SOURCES.md`,
  `docs/INTEROP.md`; 53 new tests (fully offline except two
  network-guarded live tests that skip without network/netCDF4).

## [0.1.0] - 2026-09-25

### Added
- Canonical `CurrentField` model: (nt, ny, nx) u/v velocity grids (m/s) + water
  temperature (°C), timestamps, CRS/bounds, source metadata, provenance.
  `speed()`, `direction_deg()`, `select_time()`, `select_bbox()`, `zonal_mean()`,
  JSON dict round-trip, NetCDF export, deterministic `synthetic()` fixture.
- NOAA OFS registry: 13 models (GLOFS, LMHOFS, LEOFS, CBOFS, CIOFS, CREOFS,
  DBOFS, GoMOFS, NGOFS2, SFBOFS, TBOFS, WCOFS, SSCOFS) with resolution/horizon/cycles.
- NOAA acquisition via anonymous AWS S3 (`noaa-ofs-pds`, `noaa-nos-ofs-pds`):
  stdlib-only unsigned listing (S3 ListObjectsV2 XML) and download with
  optional boto3 fast path; token-based key filtering per OFS/date/cycle/
  forecast hour; tolerant NetCDF variable-name resolution; surface-level
  selection; regular-grid guard with a clear error for unstructured meshes.
- CMEMS via the `copernicusmarine` toolbox: curated presets
  (`global-physics-daily`, `global-physics-hourly`; uo/vo/thetao), server-side
  subset, surface parsing. Missing credentials raise an actionable error —
  never an interactive prompt. Fully optional: the NOAA path needs no account.
- Per-timestep 4-band GeoTIFF/COG export (`speed, u, v, temperature`) — the
  stable interchange contract for survey-flow.
- SHA-256 provenance sidecars (`.provenance.json`) for every download, with
  read/verify helpers.
- Interop: `CurrentsPassProvider` implementing survey-monitor's `PassProvider`
  interface (forecast hours as passes; speed/u/v/temp mean metrics), and
  `align_to_thermal_zone()` for comparison with survey-thermal LST passes.
  Proposed survey-sites `currents:` config block documented in `docs/INTEROP.md`.
- Thin CLI: `list-models`, `fetch-noaa`, `fetch-cmems`, `info`,
  `export-cogs`, `synthetic`.
- Docs: README, `docs/ARCHITECTURE.md`, `docs/DATA_SOURCES.md`,
  `docs/INTEROP.md`; 69 tests passing (2 rasterio-dependent tests skip when
  rasterio is not installed; no network in tests).

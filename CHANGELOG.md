# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

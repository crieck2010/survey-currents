# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

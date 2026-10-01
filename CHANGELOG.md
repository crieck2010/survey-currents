# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.18.0] - 2026-10-01

### Added
- NOAA OFS surface currents via **keyless CO-OPS THREDDS OPeNDAP subsetting**
  (`src/currents/ofs_thredds.py`):
  `fetch_ofs_thredds(ofs_code, bbox, start, end, cadence_hours=6,
  prefer="nowcast", timeout=120.0)` -> `CurrentField` (`u`/`v` in m/s,
  `temperature` in °C, `source="ofs-thredds/<CODE>"`). The direct S3
  full-file OFS path is dead for reel work (hourly field files are
  62–70 MB each); the same models' `regulargrid` NetCDFs on the CO-OPS
  THREDDS server subset server-side through DAP2 constraint expressions
  on `/thredds/dodsC/`, parsed with a DDS-driven XDR decoder (stdlib
  `struct` + numpy — the server reorders coordinate variables first, so
  the parser follows the response DDS, never the request order). Grid
  vectors come from tiny `.ascii` probes with a regularity guard (a
  non-regular grid raises `ValueError`); land cells (`_FillValue`) are
  NaN; `Depth[0]` = 0.0 m surface (verified live 2026-10-01).
  **12 of 15 OFS models verified live 2026-10-01** with the
  `u_eastward`/`v_northward`/`temp` schema: SSCOFS, CBOFS, WCOFS, NGOFS2,
  GOMOFS, DBOFS, SFBOFS, LEOFS, LMHOFS, LOOFS, LSOFS, CIOFS (NYOFS/SJROFS
  have no `regulargrid` files; TBOFS is stale). Filename→valid-time rules
  verified live across 7 models (forecast `f{HHH}` = cycle+HHH; nowcast
  `n{HHH}` = cycle−span+HHH with the per-model span calibrated by probing
  one nowcast file's own `time`: 6 h for most models, 24 h for WCOFS);
  every served step's valid time is cross-checked against the requested
  hour and refused on mismatch — never mislabeled. **Honest cost:**
  measured live 2026-10-01, a Puget Sound bbox (240×201 grid) is
  **579,229 bytes per hourly step** vs ~65 MB for the full field file
  (~120× smaller) — the caller owns the frame budget. **Honest limits:**
  THREDDS keeps roughly the last **31 days** of OFS output; a date with
  no day catalog or an hour with no matching file raises
  `UnavailableRangeError` naming the exact URL — never silent, never
  padded. Per-step provenance records the exact OPeNDAP URLs, valid
  times, kind/hour/cycle, byte counts, and `water_fraction`.
- CLI: `fetch-ofs-thredds --ofs SSCOFS --bbox … --start … --end …
  [--cadence-hours 6] [--prefer nowcast] [--out …]` (writes field JSON).
- Tests: 29 offline tests (`tests/test_ofs_thredds.py`) against recorded
  real fixtures (day catalogs, DAS/DDS, grid-vector and 3×4 subset
  `.dods`/`.ascii` pairs, time probes — all recorded live 2026-10-01;
  the binary parser is cross-validated against the ASCII form of the
  same subset). No live network in the suite.
- `docs/DATA_SOURCES.md` section 19 and README document the THREDDS
  endpoint, the 12-model coverage table, filename conventions, per-step
  size, and the caller-owned frame budget.

## [0.17.0] - 2026-10-01

### Added
- Sub-daily (hourly) GFS winds (`src/currents/gfs_wind.py`):
  `fetch_gfs_wind(..., forecast_hours=(0,))` — per sampled day the
  adapter now fetches `gfs.t{CC}z.pgrb2.0p25.f{HHH}` for each hour in
  the tuple (GFS 0.25° is hourly f000–f120; f001/f006 verified live
  2026-10-01 through the same keyless NOMADS GRIB-filter path).
  Default `(0,)` is the f000 analysis and behaves exactly as in
  v0.16.0 (backwards compatible). Hours are validated as ints in
  0–120 (`GFS_FORECAST_HOUR_MAX`; `ValueError` otherwise) and
  normalized to sorted/deduped order, so each day assembles one
  chronological time series. `field.times` carries the forecast hour:
  the valid time is read from each GRIB message's own
  `validityDate`/`validityTime` and cross-checked against
  `cycle + forecast_hour` (a served step that mismatches its requested
  hour is refused, never mislabeled). Provenance `requests` is now a
  per-timestep list (`url`, `date`, `cycle`, `forecast_hour`,
  `window`, per-file SHA-256, byte counts, cache hits). **Honest
  cost:** each `(day, hour, window)` is one download — ~110 KB per
  hourly file for a North-America box when the subregion is honored
  (measured live 2026-10-01: 109,278 bytes, 101×261 grid), ~2.4 MB at
  the full-globe fallback the filter silently returns when it ignores
  the subregion; a full 121-step day is ~13 MB, ~290 MB worst-case —
  the caller (reel pipeline, CLI) owns the frame budget. A `(day,
  hour)` the server 404s raises `UnavailableRangeError` naming the
  exact URL — never skipped, never padded.
- CLI: `fetch-gfs-wind --forecast-hours 0,1,6` (comma-separated,
  default `0`); new recorded fixture
  `tests/fixtures/gfs_wind_20261001_f001_5x5.grib2` (real 2026-10-01
  00z f001 NOMADS response).

### Changed
- `gfs_filter_url(day, cycle, window)` gains `forecast_hour=0`;
  `gfs_wind._parse_grib_payload` gains `forecast_hour=0` and verifies
  `stepRange` plus valid time.
- `docs/DATA_SOURCES.md` section 18 and README document the hourly
  window (f000–f120, ~10-day NOMADS retention), the per-hour file
  size, and the caller-owned frame budget.

## [0.16.0] - 2026-10-01

### Added
- NOAA GFS 10-m winds adapter (`src/currents/gfs_wind.py`):
  `fetch_gfs_wind(bbox, start, end, stride_days=1, cycle="00",
  work_dir=...)` -> `GfsWindField` via the **keyless NCEP NOMADS GRIB
  filter** (`filter_gfs_0p25.pl`, plain HTTPS — no account, no token).
  One f000 **analysis** snapshot per sampled day (cycles 00/06/12/18
  UTC): `u10`/`v10` 10-m wind components (m/s) and `t2m` 2-m air
  temperature (°C here, Kelvin on the wire). Requests carry
  `subregion=on` — verified live 2026-10-01 that without it the filter
  silently returns the full 1440x721 global grid (~2.4 MB) instead of
  the requested box; the adapter verifies the returned grid covers the
  bbox and crops locally with numpy regardless. Antimeridian-crossing
  bboxes split into two subregion requests and concatenate along
  longitude (seam meridian deduped). Parsing iterates
  `cfgrib.messages.FileStream` (`cfgrib.open_file` fails on the mixed
  10-m/2-m levels — documented in the module); lazy import with an
  actionable error (`pip install 'survey-currents[gfs]'`, new
  `gfs` extra, also in `full`). `GfsWindField` matches the `Era5Field`
  shape (`grids`/`times`/`lats`/`lons`, `values` = wind speed,
  `overlay_grids` = `{"t2m": ...}`) so survey-viz's `wind` variable
  path consumes it unchanged, and adds `air_temperature` (°F — the
  warming.watch strand-color convention survey-viz's `dark_strands`
  reads) and `temperature_unit`. Provenance carries the exact request
  URLs, per-file SHA-256, byte counts, cache hits, and retrieval time;
  payloads cached under `$SURVEY_CURRENTS_CACHE/gfs-wind` (atomic
  writes, SHA-256 sidecars; analyses are immutable once posted, so
  entries never expire). **Honesty contract:** NOMADS keeps roughly
  the last 10 days of the 0.25° GFS (`GFS_RETENTION_DAYS`, verified
  live 2026-10-01: HTTP 200 for 2026-09-22…2026-10-01, HTTP 404 for
  2026-09-21) — dates outside the window and not-yet-posted cycles
  raise `UnavailableRangeError`, never silent padding; a message whose
  data date/cycle mismatches the request is refused rather than
  mislabeled; f000 is the analysis, not a forecast. Live-verified
  2026-10-01 (North-America box: KB-scale per-day downloads).
  **Design note:** a `CurrentField` was not reused — viz's wind path
  keys off `grids["u10"]`/`grids["v10"]` with `spec.variable == "wind"`,
  which the ocean-current shape does not satisfy (documented in
  `docs/DATA_SOURCES.md`).
- CLI: `fetch-gfs-wind` (real fetch), `gfs-wind-synthetic` (offline).
- Tests: 35 new tests (`tests/test_gfs_wind.py` + CLI additions), fully
  offline — a 614-byte recorded NOMADS fixture
  (`tests/fixtures/gfs_wind_20261001_5x5.grib2`: real 2026-10-01 00z
  f000 2t/10u/10v over -100…-99, 40…41) covers parsing; the HTTP layer
  is mocked for fetch/URL/cache tests; cfgrib-dependent tests skip
  cleanly without it.

## [0.15.2] - 2026-09-27

### Fixed
- `currents.__version__` is now read from the installed distribution
  metadata (written from `pyproject.toml` at install time) instead of a
  hand-edited string that had gone stale at `"0.14.0"` through the
  0.15.x releases. No behavior change.

## [0.15.1] - 2026-09-27

### Fixed
- USGS earthquakes adapter: start `limit`/`offset` paging at `offset=1`.
  The ComCat FDSN event service (v2.7.0) rejects `offset=0` with
  HTTP 400 ("Bad offset value \"0\". Valid values are 1 <= offset"),
  which broke every earthquake fetch, including
  `"earthquakes in Japan over the past 10 years"` in reel-studio.
  Verified live: the Japan 7-day M4+ query now returns events.
  Added `test_fetch_offsets_are_one_based` regression test.

## [0.15.0] - 2026-09-27

### Added
- USGS earthquakes adapter (`src/currents/earthquakes.py`):
  `fetch_earthquakes(bbox, start, end, min_magnitude=0.0,
  event_type=None, page_size=2000)` -> `QuakeField` via the **keyless
  USGS ComCat FDSN event service** (GeoJSON). Count-first, then
  `limit`/`offset` pages (paged metadata carries limit/offset, not
  count — verified live 2026-09-27); windows whose count exceeds the
  20000-event FDSN ceiling are recursively time-split; responses
  cached under `$SURVEY_CURRENTS_CACHE/comcat` (atomic writes, SHA-256
  sidecars, 7-day revalidation). Per-event records (`event_id`, UTC
  time, lat/lon, depth, magnitude, magnitude type, place, event type);
  missing magnitudes stay None, never 0.0; non-tectonic event types
  kept and labeled honestly. Depth bins for rendering (shallow <70 km,
  intermediate 70–300 km, deep >300 km); `largest()`, `select_time()`,
  `counts_by_day()`, JSON round-trip, deterministic
  `QuakeField.synthetic()`. Provenance carries exact request URLs,
  count/parsed/malformed/duplicate tallies, retrieval timestamp,
  `empty_reason`, and a catalog-completeness note. **Honesty contract:**
  ComCat is a catalog of observed events — never a forecast or hazard
  model; magnitude completeness varies by region and time; record
  floor 1900-01-01. Live-verified 2026-09-27 (global M4+ week: 179
  events; M0+: 1,906; California M4+ window: legitimate empty catalog).
- CLI: `fetch-earthquakes` (real fetch), `earthquakes-synthetic`
  (offline deterministic field).
- 32 fully-offline tests (`tests/test_earthquakes.py` — URL
  builders, parsing incl. malformed features/missing magnitude/depth,
  pagination, time-chunking over the FDSN ceiling, dedupe, cache
  integrity, model methods, mocked fetch; live test behind
  `SURVEY_CURRENTS_LIVE=1`).

### Changed
- README: ComCat section + honesty-contract limitation; docs layout
  listing; `docs/DATA_SOURCES.md` gains §17 (ComCat).

## [0.14.0] - 2026-09-27

### Added
- Ocean color adapter (`src/currents/oceancolor.py`):
  `fetch_oceancolor(bbox, start, end, product="chlorophyll-a",
  cadence="monthly", sensor="modis-aqua", source="coastwatch")` ->
  `OceanColorField` (chlorophyll-a in mg/m³, `(nt, ny, nx)` masked array
  — cloud/land gaps stay NaN, never interpolated). Sensors: `modis-aqua`
  (`erdMH1chlamday_R2022SQ`, 2002–present, default), `viirs-snpp`
  (`nesdisVHNSQchlaMonthly/Weekly/Daily`, 2012–present),
  `multi` (ESA OC-CCI v6.0 `pmlEsaCCI60OceanColorMonthly`,
  1997–present; ID verified 2026-09-27 against the official
  coastwatch-training/coastwatch-tutorials repo). `source="obpg"` fetches
  MODIS Aqua L3 mapped chlorophyll-a from the NASA OBPG direct data
  access (`oceandata.sci.gsfc.nasa.gov`, documented
  `AQUA_MODIS.<dates>.L3m.<DAY|8D|MO>.CHL.chlor_a.4km.nc` filenames;
  CMR collection `C3380709133-OB_CLOUD`/`MODISA_L3m_CHL` verified live
  2026-09-27) — needs a free Earthdata Login and raises the
  credentials-gated `CredentialsMissing` without it (never an
  interactive prompt); this is the required fallback when CoastWatch is
  unreachable. `source="cmems"` uses the authenticated
  `OCEANCOLOUR_GLO_BGC_L4_MY_009_104` product via the copernicusmarine
  toolbox (free CMEMS credentials; optional extra source);
  `source="auto"` tries CoastWatch, then OBPG when Earthdata
  credentials exist, then CMEMS when credentials exist, else raises
  with an actionable message naming all three sources. Monthly is
  the default cadence (most cloud-complete). Provenance carries the
  exact URLs, SHA-256 of the payloads, byte counts, retrieval
  time, sensor, processing level, cadence, units, grid geometry, and
  per-frame NaN gap fractions (OBPG provenance also records the CMR
  collection and an honest `live_verified: false` until its network
  path is exercised with real credentials). Downloads are cached
  (`$SURVEY_CURRENTS_CACHE/oceancolor`, 7-day freshness) with atomic
  writes + SHA-256 sidecars and corruption recovery. CLI:
  `fetch-oceancolor`, `oceancolor-synthetic`. 48 fully offline tests
  (monkeypatched transport; netCDF4-gated parse fixtures). **Access
  truth 2026-09-27:** the CoastWatch ERDDAP was unreachable from the
  build environment during verification (HTTP 502/503 on griddap;
  `coastwatch.noaa.gov/erddap` itself intermittently 503). Dataset IDs,
  grid geometry, and query grammar are corroborated against the
  CoastWatch tutorials repo and sanctuary caption pages; the request
  path was exercised against a live ERDDAP server (NCEI) with the same
  `.das`/griddap grammar. A successful end-to-end CoastWatch data pull
  still needs a re-run once the service recovers — the adapter is
  written so it does exactly that.
- `docs/DATA_SOURCES.md`: ocean-color source table (MODIS Aqua R2022
  monthly default, VIIRS SNPP monthly/weekly/daily, ESA OC-CCI v6.0
  monthly, NASA OBPG fallback, CMEMS 009_104 optional), resolution/
  download-size math, and the 2026-09-27 CoastWatch outage record.

## [0.13.0] - 2026-09-27

### Added
- USGS Water Services (NWIS) streamgage daily-value adapter
  (`src/currents/streamgages.py`): `fetch_usgs(bbox, start, end,
  parameters=("00060",), min_record_days=30, site_limit=200,
  units="native")` -> `GageField` (list of `GageRecord`: site number,
  station name, lat/lon, HUC, drainage area when published, and
  per-parameter daily series for 00060 discharge in ft³/s and 00065
  gage height in ft — the native USGS publication units; the USGS
  missing-value sentinel `"-999999"` becomes NaN and is listed by
  `GageRecord.missing_days()`, never filled). **Verified live
  2026-09-27:** keyless anonymous HTTPS — site service
  `https://waterservices.usgs.gov/nwis/site/` (RDB;
  `hasDataTypeCd=dv&parameterCd=…&siteOutput=expanded` adds drainage
  area, HUC, timezone) and dv service
  `https://waterservices.usgs.gov/nwis/dv/` (JSON, `statCd=00003`
  daily MEAN stated explicitly, qualifiers P/A/e carried). Sites are
  discovered in the bbox (deterministic order by site number), daily
  values fetched in batches of 20, and sites below `min_record_days`
  finite daily values are excluded but counted in provenance
  (`n_sites_no_data`, `n_sites_excluded_short_record`); sites with no
  data in the window yield an honest empty field with `empty_reason`
  (NWIS is US-only). Cache discipline mirrors the storms/GRACE
  adapters: per-bbox inventory and per-batch dv payloads are cached
  under `$SURVEY_CURRENTS_CACHE/streamgages` with atomic writes,
  SHA-256 sidecars verified on every hit, corruption-triggered
  redownload, `max_cache_age_days=7` revalidation (recent daily values
  are provisional and get revised), `refresh=True` forces re-download.
  `units="si"` converts to m³/s / m on exact NIST factors
  (`1 ft³/s = 0.028316846592 m³/s`). `GageRecord.percentile_of_record()`
  (rank of the latest value within its own record) and
  `GageField.regional_median()` (daily median across sites, NaN-aware)
  document the survey-viz rendering rules; `GageField.to_si()`,
  `select_site()`, `to_dict`/`from_dict`, `to_json`/`from_json`, and a
  deterministic `synthetic()` fixture (4 sites incl. an engineered
  10-day gap block) round out the model. CLI: `fetch-usgs` and
  `usgs-synthetic`. 36 new offline tests (2 live opt-in tests guarded
  by `SURVEY_CURRENTS_LIVE=1`).

## [0.12.0] - 2026-09-27

### Added
- CSR GRACE / GRACE-FO RL06.3 terrestrial water storage adapter
  (`src/currents/grace.py`): `fetch_grace(bbox, start, end)` ->
  `WaterField` (monthly terrestrial water storage anomalies, cm
  liquid-water-equivalent thickness, land-only, 0.25° grid,
  longitudes normalized to -180..180). **Verified live
  2026-09-27:** keyless CSR HTTPS
  `https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/`;
  `CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc` (~107 MB,
  258 months, 2002-04 – 2026-06, already gridded at 0.25° — no mascon
  mapping needed; native mascon resolving power is coarser than the
  0.25° output grid) and the separate
  `CSR_GRACE_GRACE-FO_RL06_Mascons_v02_LandMask.nc` (~4 MB,
  `LO_val` 0/1 on the same grid). Units are cm; the 2004–2009
  time-mean has been removed (`time_mean_removed` attribute) so every
  value is an anomaly against that baseline. `ensure_grace_files()`:
  download-once into `$SURVEY_CURRENTS_CACHE/grace`, atomic writes,
  SHA-256 sidecars verified on every hit, corruption-triggered
  redownload, `Last-Modified`-based stale-file revalidation (CSR
  re-releases the file as new months arrive), `refresh=True` forces
  re-download. Months with no solution in the requested window become
  all-NaN gap frames and are listed in `WaterField.gap_months` —
  never interpolated, never silently dropped; gap detection diffs the
  file's actual time axis against the calendar (the file's
  `months_missing` attribute under-reports: 2011-11 and 2015-05 are
  missing from the time axis but absent from the attribute — the 2017-07
  … 2018-05 inter-mission gap plus 24 isolated months are all
  handled). `WaterField` has `spatial_mean()`, `is_gap()`,
  `select_time`, `select_bbox`, dict/JSON round-trips, and a
  deterministic `synthetic()` fixture (drying trend + seasonal cycle,
  ~35% ocean-masked cells, configurable gap months). netCDF4 stays
  lazy/optional (pure parsing works over duck-typed datasets). CLI:
  `fetch-grace` and `grace-synthetic`. 42 new offline tests (1 live
  opt-in test guarded by `SURVEY_CURRENTS_LIVE=1`).

## [0.11.0] - 2026-09-27

### Added
- NOAA IBTrACS v4 tropical-cyclone best-track adapter
  (`src/currents/storms.py`): `fetch_ibtracs(bbox, start, end,
  min_wind=None, storm_name=None, full_archive=False)` -> `StormField`
  (list of `StormTrack`: time, lat/lon, max sustained wind in kt,
  min central pressure in hPa, storm name, SID, basin). **Verified
  live 2026-09-27:** keyless NCEI HTTPS
  `https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/netcdf/`;
  `IBTrACS.since1980.v04r01.nc` (~10.8 MB, 1980–present, default) and
  `IBTrACS.ALL.v04r01.nc` (~23.4 MB, 1842–present) are NetCDF-4/HDF5
  (storm dim 4991, 360 fix slots; time units days since 1858-11-17;
  longitudes mixed 0..360 / -180..180). Per-fix wind priority
  `usa_wind` then `wmo_wind`, pressure `usa_pres` then `wmo_pres`
  (integer fill -9999). `ensure_ibtracs_file()`: download-once,
  atomic writes, SHA-256 sidecars verified on every hit,
  corruption-triggered redownload, `Last-Modified`-based stale-file
  revalidation (the archive is republished as storms are added),
  `refresh=True` forces re-download. Longitudes normalized to
  [-180, 180); storms kept when any fix falls in bbox, fixes clipped
  to `[start, end]`, `min_wind` filters by lifetime max sustained
  wind, `storm_name` selects one named storm (case-insensitive).
  `StormField` has `rank_by_intensity()`, `select_time`,
  `select_bbox`, JSON round-trips, and a deterministic `synthetic()`.
  Documented Saffir-Simpson mapping (TD <34, TS 34–63, C1 64–82,
  C2 83–95, C3 96–112, C4 113–136, C5 ≥137 kt) with per-category
  render colors. CLI: `fetch-ibtracs`, `storms-synthetic`.
  Docs: `docs/DATA_SOURCES.md` §13.
- 42 tests in `tests/test_storms.py` (41 offline, 1 live Katrina-2005
  verification skipped unless `SURVEY_CURRENTS_LIVE=1`).

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

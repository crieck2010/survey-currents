# Interoperability

How `survey-currents` plugs into the survey-suite / earthwatch-suite
engines. All bridges are duck-typed — this package never imports a
sibling, so install graphs stay acyclic and every bridge is testable
offline.

## survey-monitor — `CurrentsPassProvider`

`currents.interop.CurrentsPassProvider` implements the `PassProvider`
interface from `survey-monitor` (`src/monitor/acquire.py`):
`list_passes(since)` / `metrics(pass_id)`.

- Each **forecast hour becomes one pass**: pass ids look like
  `currents-<site_id>-f006-2026-09-16`.
- `metrics()` returns NaN-aware spatial means over the site AOI:
  `speed_mean`, `u_mean`, `v_mean` (m/s) and `temp_mean` (°C, when the
  field carries temperature).
- When `survey-monitor` is installed, `list_passes()` returns its real
  `PassInfo`; otherwise a local twin with identical fields.
- Geometry accepts the monitor's `SiteConfig.geometry` shape: a GeoJSON
  geometry dict or `{"bbox": [min_lon, min_lat, max_lon, max_lat]}`.

This lets scheduled monitoring watch e.g. "mean surface current speed
in the shipping channel" with the same alerting machinery as the
satellite indices — model hours simply appear as passes.

## survey-thermal — `align_to_thermal_zone`

`currents.interop.align_to_thermal_zone(field, zone, start, end)`:

- `zone`: a thermal `Zone` as a dict (`id` + `geometry`).
- Clips the field to the zone bbox, filters timesteps to
  `[start, end]`, and returns per-timestep mean water temperature —
  directly comparable with the per-pass LST CSV rows produced by
  `thermal.acquire.acquire_lst_passes`.
- Caveat: model forecast hours are hours-since-cycle, not local solar
  time like Landsat's ~10:00 overpass — align windows deliberately.

## survey-sites — proposed `currents:` config block

No changes were made to the survey-sites repo. Proposed schema
extension for a future `site-config` version (for the survey-sites
maintainer):

```yaml
currents:
  source: noaa-ofs            # or "cmems"
  ofs_code: LMHOFS            # NOAA path
  # preset: global-physics-daily   # CMEMS path
  cycle: "00"
  hours: [0, 6, 12, 18]       # forecast hours per run (monitor passes)
  bbox: [-92.5, 41.5, -84.5, 46.5]
  metrics: [speed_mean, temp_mean]
  cmems_username_env: CMEMS_USER   # optional; only for the CMEMS path
```

A monitor run with this block would fetch the listed forecast hours
through `CurrentsPassProvider` and store the resulting COGs as run
artifacts for survey-flow.

## survey-flow (next module) — the interchange contract

survey-flow consumes, in order of preference:

1. `CurrentField` objects (in-process).
2. `CurrentField` JSON (`to_json`/`from_json`) — small fields, debugging.
3. Per-timestep 4-band COGs from `convert.export_cogs` — **band order
   `(speed, u, v, temperature)` is the stable contract**; units m/s,
   m/s, m/s, °C; NaN for missing data/temperature.
4. NetCDF from `to_netcdf()` — variables `u`, `v`, `temperature`
   (optional), coords `time`, `lat`, `lon`.

These formats are versioned: breaking changes require a minor version
bump and a CHANGELOG entry.

## survey-qgis / survey-alerts (later)

- A future survey-qgis Processing algorithm can style the 4-band COGs
  (speed colormap + temperature overlay) reusing this package's band
  contract.
- survey-alerts can attach finished reels; the per-hour metrics from
  `CurrentsPassProvider` already flow through the monitor's
  `AlertEvent` machinery unchanged.

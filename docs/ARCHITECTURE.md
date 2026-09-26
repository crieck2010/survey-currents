# Architecture

`survey-currents` is a pure-Python **acquisition engine**: it turns
operational hydrodynamic forecast products into one canonical in-memory
model (`CurrentField`) plus on-disk artifacts (NetCDF, COGs, provenance
sidecars). Since v0.2.0 it also acquires NOAA GLSEA satellite
sea-surface temperature into the companion `GlseaField` model (plus
`LakeSeries` lake-average series). It does no rendering and no
scheduling — those belong to survey-flow, survey-animate, and
survey-monitor respectively.

## Module map

```
src/currents/
    models.py      Canonical model + OFS registry. Only numpy + stdlib.
                   CurrentField is the single type the rest of the
                   program passes around.
    noaa_ofs.py    Anonymous S3 access + OFS NetCDF parsing.
                   Network: s3_list / s3_download (urllib stdlib path,
                   boto3 fast path). Parsing: parse_ofs_netcdf (xarray,
                   lazy). Orchestrator: fetch_currents -> stack_steps.
    glsea.py       NOAA GLSEA satellite SST via ERDDAP (stdlib urllib;
                   netCDF4 lazy). fetch_glsea_sst -> GlseaField (masked
                   SST grids + url/sha256/retrieved_at provenance);
                   fetch_glsea_lake_averages -> LakeSeries (stdlib csv).
                   Bboxes are validated against the grid's clipped
                   lakes-region extent (notably the -92.41995 longitude
                   floor) before any download.
    cmems.py       copernicusmarine wrapper. subset_cmems downloads one
                   NetCDF; parse_cmems_netcdf -> CurrentField.
                   require_toolbox() fails fast with setup instructions
                   when credentials are absent.
    convert.py     CurrentField -> per-timestep 4-band GeoTIFF/COG
                   (rasterio, lazy). Band order is the survey-flow
                   contract: speed, u, v, temperature.
    provenance.py  .provenance.json sidecars: source, bbox, time window,
                   SHA-256, tool version. write/read/verify.
    interop.py     survey-suite bridges. No sibling imports: duck-typed
                   CurrentsPassProvider (monitor PassProvider interface)
                   and align_to_thermal_zone (thermal zone dicts).
    cli.py         Thin argparse adapter. Zero UI-framework imports.
```

## Data flow

```
NOAA OFS (S3, anonymous) ──┐
                           ├─> NetCDF per forecast hour ──> parse ──┐
CMEMS (toolbox, account) ──┘                                        │
                                                              stack_steps
                                                                    │
NOAA GLSEA (ERDDAP, no account) ──> SST subset NetCDF ──> parse ──>  │
                                                             CurrentField / GlseaField
                                                              (nt,ny,nx)
                                                                    ├──> to_netcdf / to_json
                                                                    ├──> export_cogs  -> survey-flow
                                                                    ├──> CurrentsPassProvider -> survey-monitor
                                                                    └──> align_to_thermal_zone -> survey-thermal
```

Every download writes a provenance sidecar next to the file (S3/CMEMS
path); the GLSEA path carries url/sha256/retrieved_at in the field's
`provenance` dict. The `CurrentField.provenance` dict carries the
parse-time record (variable names chosen, grid type) forward.

## Key design decisions

- **Core is stdlib+numpy.** `models.py`, `provenance.py`, `interop.py`
  import nothing heavier. xarray, boto3, copernicusmarine, and rasterio
  are lazy optional imports with actionable error messages naming the
  `pip install` extra. The engine is fully usable offline via
  `CurrentField.synthetic()`.
- **Token filtering, not filename assumptions.** NOAA OFS filename
  conventions differ per model, so `filter_keys` matches on
  OFS-code/date/cycle/hour tokens inside the documented date prefix
  instead of constructing an exact filename.
- **Validate before downloading.** The GLSEA grid is clipped to the
  lakes region, so `validate_glsea_bbox` rejects out-of-extent boxes
  (notably the −92.41995 longitude floor) before any bytes move.
  Failing fast on the client keeps ERDDAP load down and errors
  actionable.
- **Regular grids in v0.1.0.** The parser requires u/v to be indexed by
  lat/lon dims and raises a clear `ValueError` for unstructured
  (FVCOM) meshes. Regridding is deferred to a later release — the
  `CurrentField` contract does not change.
- **No sibling imports.** Interop is duck-typed: `CurrentsPassProvider`
  returns survey-monitor's real `PassInfo` when that package is
  installed and a local twin otherwise. This keeps install graphs
  acyclic and every module testable in isolation.
- **Stable interchange.** COG band order `(speed, u, v, temperature)`
  and the `CurrentField` JSON schema are versioned contracts —
  survey-flow and survey-animate will be built against them, so
  breaking changes require a minor version bump and a CHANGELOG entry.

## Scaling notes

- Fields are (nt, ny, nx) float64 in memory. A 7-day hourly LMHOFS
  subset at full resolution is large — clip with `bbox` early
  (`select_bbox` / the `fetch_currents` bbox argument) and prefer
  GLOFS for basin-scale work.
- `s3_list` paginates; `s3_download` streams in 1 MiB chunks.
- For production pipelines, persist `CurrentField` via `to_netcdf()`
  (float32 on write is a caller-side option) rather than the JSON
  form, which is intended for small fields and debugging.

#!/usr/bin/env python3
"""Offline GLSEA demo: synthetic SST field -> JSON -> summary.

Uses GlseaField.synthetic() only — no network, no netCDF4
(stdlib + numpy only). Mirrors examples/run_demo.sh for the
CurrentField path.
"""

from currents.glsea import (GlseaField, LakeSeries, fetch_glsea_lake_averages,
                            glsea_averages_url, glsea_sst_url,
                            validate_glsea_bbox)

# 1. Build a deterministic synthetic SST field (warm core + cooling trend).
field = GlseaField.synthetic(nt=6, ny=8, nx=10, seed=7)
field.to_json("examples/glsea_demo_field.json")
print(f"timesteps : {len(field.times)}")
print(f"bounds    : {field.bounds}")
print(f"shape     : {field.sst.shape} (masked cells: {field.sst.mask.sum()})")
print(f"mean SST  : {field.spatial_mean(0):.2f} degC at step 0")

# 2. Clip to a site AOI and re-read the JSON round-trip.
site = field.select_bbox((-90.0, 44.0, -87.0, 46.0))
back = GlseaField.from_json("examples/glsea_demo_field.json")
print(f"site clip : bounds {site.bounds}, step-0 mean {site.spatial_mean(0):.2f} degC")
print(f"roundtrip : {len(back.times)} steps, source={back.source!r}")

# 3. Show the URL the live fetch would request (no download here).
bbox = (-92.0, 46.5, -87.0, 48.0)
validate_glsea_bbox(bbox)  # raises before download if out of the lakes grid
print("griddap   :", glsea_sst_url(bbox, "2025-01-01", "2025-03-01"))
print("tabledap  :", glsea_averages_url("superior", 2025))

# 4. Parse a LakeSeries from captured real CSV lines (offline fixture).
csv_text = """Year,Day,Sup
,,
2025,001,3.21
2025,002,3.10
2025,003,2.98
"""
from currents.glsea import parse_glsea_averages_csv
import datetime as _dt
series = parse_glsea_averages_csv(
    csv_text, "superior", _dt.date(2025, 1, 1), _dt.date(2025, 1, 31))
print(f"lake avg  : {series.lake} n={series.n} "
      f"mean={series.mean():.2f} degC over {series.dates[0]}..{series.dates[-1]}")
series.to_json("examples/glsea_demo_series.json")
assert LakeSeries.from_json("examples/glsea_demo_series.json").temps == series.temps
print("demo OK — no network used")

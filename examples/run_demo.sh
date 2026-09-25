#!/usr/bin/env bash
# Offline demo: synthetic current field -> JSON -> summary.
# No network, no heavy dependencies (stdlib + numpy only).
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHONPATH=src python3 -m currents.cli synthetic --nt 24 --out examples/demo_field
PYTHONPATH=src python3 - <<'EOF'
from currents.models import CurrentField
from currents.interop import CurrentsPassProvider, align_to_thermal_zone

field = CurrentField.from_json("examples/demo_field.json")
print(f"timesteps : {len(field.times)}")
print(f"bounds    : {field.bounds}")
print(f"mean speed: {field.zonal_mean(0, 'speed'):.3f} m/s")

provider = CurrentsPassProvider(field, "demo-site",
                                {"bbox": [-92.5, 41.5, -84.5, 46.5]})
passes = provider.list_passes()
print(f"passes    : {len(passes)} (first: {passes[0].pass_id})")
print("metrics   :", provider.metrics(passes[0].pass_id))

zone = {"id": "demo-zone",
        "geometry": {"bbox": [-92.0, 42.0, -86.0, 45.0]}}
aligned = align_to_thermal_zone(field, zone)
print(f"zone temps: {len(aligned['timesteps'])} timesteps, "
      f"first mean {aligned['temp_mean_series'][0]:.2f} C")
EOF

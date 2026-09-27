"""Thin CLI adapter over the currents engine (no UI-framework imports)."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .models import list_ofs_models


def _parse_bbox(s: str):
    parts = [float(x) for x in s.split(",")]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("bbox needs min_lon,min_lat,max_lon,max_lat")
    return parts


def _parse_hours(s: str):
    hours = set()
    for part in s.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-", 1)
            hours.update(range(int(a), int(b) + 1))
        elif part:
            hours.add(int(part))
    return sorted(hours)


def cmd_list_models(_a: argparse.Namespace) -> int:
    for m in list_ofs_models():
        print(f"{m.code:8s} {m.region:42s} {m.resolution:16s} {m.horizon_hours}h  "
              f"cycles {','.join(m.cycles)}")
    return 0


def cmd_fetch_noaa(a: argparse.Namespace) -> int:
    from .noaa_ofs import fetch_currents
    field = fetch_currents(a.ofs, a.date, cycle=a.cycle, hours=_parse_hours(a.hours),
                           work_dir=a.out, bbox=_parse_bbox(a.bbox) if a.bbox else None)
    print(f"fetched {len(field.times)} timesteps -> {a.out}")
    print(f"source={field.source} model_run={field.model_run} bounds={field.bounds}")
    out_json = f"{a.out.rstrip('/')}.json"
    field.to_json(out_json)
    print(f"wrote {out_json}")
    return 0


def cmd_fetch_cmems(a: argparse.Namespace) -> int:
    from .cmems import parse_cmems_netcdf, subset_cmems
    path = subset_cmems(a.preset, _parse_bbox(a.bbox), a.start, a.end, a.out)
    field = parse_cmems_netcdf(path)
    print(f"downloaded {path}: {len(field.times)} timesteps, bounds={field.bounds}")
    return 0


def cmd_info(a: argparse.Namespace) -> int:
    from .models import CurrentField
    if a.input.endswith(".json"):
        field = CurrentField.from_json(a.input)
    else:
        from .noaa_ofs import parse_ofs_netcdf
        field = parse_ofs_netcdf(a.input)
    print(json.dumps({
        "source": field.source, "model_run": field.model_run,
        "times": field.times, "forecast_hours": field.forecast_hours,
        "shape": list(field.u.shape), "bounds": list(field.bounds),
        "crs": field.crs,
        "speed_mean_m_s": float(__import__("numpy").nanmean(field.speed())),
    }, indent=2))
    return 0


def cmd_export_cogs(a: argparse.Namespace) -> int:
    from .models import CurrentField
    from .convert import export_cogs
    field = CurrentField.from_json(a.input)
    paths = export_cogs(field, a.out, prefix=a.prefix)
    for p in paths:
        print(p)
    return 0


def cmd_synthetic(a: argparse.Namespace) -> int:
    from .models import CurrentField
    field = CurrentField.synthetic(nt=a.nt, seed=a.seed)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic field ({len(field.times)} steps) -> {path}")
    return 0


def cmd_fetch_glsea_sst(a: argparse.Namespace) -> int:
    from .glsea import fetch_glsea_sst
    field = fetch_glsea_sst(_parse_bbox(a.bbox), a.start, a.end,
                            stride_days=a.stride_days)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched {len(field.times)} daily SST timesteps -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"sha256={field.provenance['sha256']}")
    return 0


def cmd_fetch_glsea_averages(a: argparse.Namespace) -> int:
    from .glsea import fetch_glsea_lake_averages
    series = fetch_glsea_lake_averages(a.lake, a.start, a.end)
    path = f"{a.out.rstrip('/')}.json"
    series.to_json(path)
    print(f"fetched {series.n} daily lake-average temps ({series.lake}) -> {path}")
    print(f"mean={series.mean():.2f} degC over {series.dates[0]}..{series.dates[-1]}")
    return 0


def cmd_glsea_synthetic(a: argparse.Namespace) -> int:
    from .glsea import GlseaField
    field = GlseaField.synthetic(nt=a.nt, seed=a.seed)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic GLSEA SST field ({len(field.times)} steps) -> {path}")
    return 0


def cmd_fetch_oisst(a: argparse.Namespace) -> int:
    from .sst_global import fetch_oisst
    field = fetch_oisst(_parse_bbox(a.bbox), a.start, a.end,
                        stride_days=a.stride_days)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched {len(field.times)} OISST v2.1 SST timesteps -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"sha256={field.provenance['sha256']}")
    return 0


def cmd_fetch_mur(a: argparse.Namespace) -> int:
    from .sst_global import fetch_mur
    field = fetch_mur(_parse_bbox(a.bbox), a.start, a.end,
                      stride_days=a.stride_days)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched {len(field.times)} MUR v4.1 SST timesteps -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"combined_sha256={field.provenance['combined_sha256']}")
    return 0


def cmd_fetch_era5(a: argparse.Namespace) -> int:
    from .era5 import fetch_era5
    variables = [v.strip() for v in a.variables.split(",") if v.strip()]
    field = fetch_era5(variables, _parse_bbox(a.bbox), a.start, a.end,
                       stride_hours=a.stride_hours)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched {len(field.times)} ERA5 timesteps "
          f"(variables={','.join(field.variables)}) -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"sha256={field.provenance['sha256']}")
    return 0


def cmd_era5_synthetic(a: argparse.Namespace) -> int:
    from .era5 import Era5Field
    variables = [v.strip() for v in a.variables.split(",") if v.strip()]
    field = Era5Field.synthetic(variables=variables, nt=a.nt, seed=a.seed,
                                source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic ERA5 field ({len(field.times)} steps, "
          f"variables={','.join(field.variables)}) -> {path}")
    return 0


def cmd_fetch_oscar(a: argparse.Namespace) -> int:
    from .currents_global import fetch_oscar
    field = fetch_oscar(_parse_bbox(a.bbox), a.start, a.end,
                        stride_days=a.stride_days)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched {len(field.times)} OSCAR v2.0 current timesteps -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"collections_used={field.provenance.get('collections_used')}")
    print(f"combined_sha256={field.provenance.get('combined_sha256')}")
    return 0


def cmd_fetch_cmems_currents(a: argparse.Namespace) -> int:
    from .currents_global import fetch_cmems_currents
    field = fetch_cmems_currents(_parse_bbox(a.bbox), a.start, a.end,
                                 stride_days=a.stride_days,
                                 work_dir=a.work_dir)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched {len(field.times)} CMEMS current timesteps -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"netcdf={field.provenance.get('netcdf_path', '(in-memory)')}")
    return 0


def cmd_currents_synthetic(a: argparse.Namespace) -> int:
    from .currents_global import cmems_currents_synthetic, oscar_synthetic
    for label, make in (("oscar", oscar_synthetic),
                        ("cmems-currents", cmems_currents_synthetic)):
        field = make(nt=a.nt, seed=a.seed)
        path = f"{a.out.rstrip('/')}_{label}.json"
        field.to_json(path)
        print(f"wrote synthetic {label} field "
              f"({len(field.times)} steps) -> {path}")
    return 0


def cmd_fetch_firms(a: argparse.Namespace) -> int:
    from .fires import fetch_firms
    instruments = [i.strip() for i in a.instruments.split(",") if i.strip()]
    field = fetch_firms(_parse_bbox(a.bbox), a.start, a.end,
                        instruments=instruments)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched {len(field)} FIRMS active-fire detections -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"products_used={field.provenance.get('products_used')}")
    print(f"n_requests={field.provenance.get('n_requests')}")
    return 0


def cmd_fires_synthetic(a: argparse.Namespace) -> int:
    from .fires import FireField
    field = FireField.synthetic(bbox=_parse_bbox(a.bbox), start=a.start,
                                end=a.end, n=a.n, seed=a.seed,
                                source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic FIRMS field ({len(field)} detections) -> {path}")
    return 0


def cmd_fetch_nsidc(a: argparse.Namespace) -> int:
    from .sea_ice import fetch_nsidc_sic
    field = fetch_nsidc_sic(_parse_bbox(a.bbox), a.start, a.end,
                            hemisphere=a.hemisphere,
                            stride_days=a.stride_days,
                            resolution=a.resolution)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched NSIDC sea-ice field ({len(field)} days, "
          f"hemisphere={field.hemisphere}) -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"n_files={field.provenance.get('n_files')} "
          f"skipped={len(field.provenance.get('skipped_days', []))}")
    return 0


def cmd_ice_synthetic(a: argparse.Namespace) -> int:
    from .sea_ice import IceField
    field = IceField.synthetic(bbox=_parse_bbox(a.bbox), start=a.start,
                               end=a.end, resolution=a.resolution,
                               seed=a.seed, source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic NSIDC sea-ice field ({len(field)} days) -> {path}")
    return 0


def cmd_fetch_imerg(a: argparse.Namespace) -> int:
    from .imerg import fetch_imerg
    field = fetch_imerg(_parse_bbox(a.bbox), a.start, a.end,
                        accumulate=a.accumulate, run=a.run,
                        stride_days=a.stride_days)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched GPM IMERG rain field ({len(field)} steps, run={field.run}, "
          f"accumulate={field.accumulate}) -> {path}")
    print(f"source={field.source} units={field.units} bounds={field.bounds}")
    print(f"n_files={field.provenance.get('n_files')} "
          f"skipped={len(field.provenance.get('skipped_slots', []))}")
    return 0


def cmd_rain_synthetic(a: argparse.Namespace) -> int:
    from .imerg import RainField
    field = RainField.synthetic(bbox=_parse_bbox(a.bbox), start=a.start,
                                end=a.end, resolution=a.resolution,
                                accumulate=a.accumulate,
                                seed=a.seed, source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic IMERG rain field ({len(field)} steps, "
          f"accumulate={field.accumulate}) -> {path}")
    return 0


def cmd_fetch_blackmarble(a: argparse.Namespace) -> int:
    from .blackmarble import fetch_blackmarble
    field = fetch_blackmarble(_parse_bbox(a.bbox), a.start, a.end,
                              product=a.product,
                              stride_days=a.stride_days,
                              resolution=a.resolution)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched Black Marble lights field ({len(field)} days, "
          f"product={a.product}) -> {path}")
    print(f"source={field.source} units={field.units} bounds={field.bounds}")
    print(f"n_files={field.provenance.get('n_files')} "
          f"skipped={len(field.provenance.get('skipped', []))}")
    return 0


def cmd_lights_synthetic(a: argparse.Namespace) -> int:
    from .blackmarble import LightsField
    field = LightsField.synthetic(bbox=_parse_bbox(a.bbox), start=a.start,
                                  end=a.end, resolution=a.resolution,
                                  seed=a.seed, source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic Black Marble lights field ({len(field)} days) "
          f"-> {path}")
    return 0


def cmd_fetch_grace(a: argparse.Namespace) -> int:
    from .grace import fetch_grace
    field = fetch_grace(_parse_bbox(a.bbox), a.start, a.end)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched GRACE water-storage field ({len(field)} months, "
          f"{field.n_gap_months} gap months) -> {path}")
    print(f"source={field.source} units={field.units} bounds={field.bounds}")
    print(f"baseline: {field.anomaly_baseline}")
    print(f"sha256={field.provenance.get('solution_sha256')}")
    return 0


def cmd_grace_synthetic(a: argparse.Namespace) -> int:
    from .grace import WaterField
    field = WaterField.synthetic(bbox=_parse_bbox(a.bbox), start=a.start,
                                 end=a.end, resolution=a.resolution,
                                 seed=a.seed, source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic GRACE water-storage field ({len(field)} months, "
          f"{field.n_gap_months} gap months) -> {path}")
    return 0


def cmd_fetch_ibtracs(a: argparse.Namespace) -> int:
    from .storms import fetch_ibtracs
    field = fetch_ibtracs(_parse_bbox(a.bbox), a.start, a.end,
                          min_wind=a.min_wind, storm_name=a.storm_name,
                          full_archive=a.full_archive)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"fetched IBTrACS storm field ({len(field)} storms, "
          f"{field.n_obs} fixes) -> {path}")
    print(f"source={field.source} bounds={field.bounds}")
    print(f"sha256={field.provenance.get('sha256')}")
    return 0


def cmd_storms_synthetic(a: argparse.Namespace) -> int:
    from .storms import StormField
    field = StormField.synthetic(bbox=_parse_bbox(a.bbox), start=a.start,
                                 end=a.end, n_storms=a.n_storms,
                                 seed=a.seed, source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic IBTrACS storm field ({len(field)} storms, "
          f"{field.n_obs} fixes) -> {path}")
    return 0


def cmd_sst_synthetic(a: argparse.Namespace) -> int:
    from .sst_global import SstField
    field = SstField.synthetic(nt=a.nt, seed=a.seed, source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic global SST field ({len(field.times)} steps) -> {path}")
    return 0


def cmd_fetch_gebco(a: argparse.Namespace) -> int:
    from .basemaps import fetch_gebco
    field = fetch_gebco(_parse_bbox(a.bbox), resolution=a.resolution)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    ny, nx = field.elevation.shape[0], field.elevation.shape[1]
    print(f"fetched GEBCO topography field ({ny}x{nx} cells, "
          f"resolution={field.resolution}°) -> {path}")
    print(f"source={field.source} units={field.units} bounds={field.bounds}")
    print(f"tiles={len(field.provenance.get('tiles', []))}")
    return 0


def cmd_fetch_naturalearth(a: argparse.Namespace) -> int:
    from .basemaps import fetch_naturalearth
    layers = tuple(s.strip() for s in a.layers.split(",") if s.strip())
    collections = fetch_naturalearth(_parse_bbox(a.bbox), scale=a.scale,
                                     layers=layers)
    import json as _json
    for layer, fc in collections.items():
        path = f"{a.out.rstrip('/')}_{layer}.geojson"
        with open(path, "w", encoding="utf-8") as fh:
            _json.dump(fc, fh)
        print(f"wrote Natural Earth {layer} "
              f"({len(fc['features'])} features) -> {path}")
    return 0


def cmd_basemaps_synthetic(a: argparse.Namespace) -> int:
    from .basemaps import TopoField
    field = TopoField.synthetic(bbox=_parse_bbox(a.bbox),
                                resolution=a.resolution, seed=a.seed,
                                source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    ny, nx = field.elevation.shape[0], field.elevation.shape[1]
    print(f"wrote synthetic topography field ({ny}x{nx} cells) -> {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="survey-currents",
                                description="Surface-current / water-temperature acquisition engine")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("list-models", help="list registered NOAA OFS models")
    s.set_defaults(func=cmd_list_models)

    s = sub.add_parser("fetch-noaa", help="fetch NOAA OFS fields (needs network)")
    s.add_argument("--ofs", required=True, help="OFS code, e.g. LMHOFS")
    s.add_argument("--date", required=True, help="model run date YYYY-MM-DD")
    s.add_argument("--cycle", default="00", help="model cycle, e.g. 00/06/12/18")
    s.add_argument("--hours", default="0", help="forecast hours, e.g. 0-23 or 0,6,12")
    s.add_argument("--bbox", default=None, help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--out", default="currents_data", help="work directory")
    s.set_defaults(func=cmd_fetch_noaa)

    s = sub.add_parser("fetch-cmems", help="fetch CMEMS subset (needs account + network)")
    s.add_argument("--preset", default="global-physics-daily")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--start", required=True, help="ISO start datetime")
    s.add_argument("--end", required=True, help="ISO end datetime")
    s.add_argument("--out", default="currents_data", help="work directory")
    s.set_defaults(func=cmd_fetch_cmems)

    s = sub.add_parser("info", help="summarise a NetCDF or field JSON")
    s.add_argument("input", help=".nc file or field .json")
    s.set_defaults(func=cmd_info)

    s = sub.add_parser("export-cogs", help="field JSON -> per-timestep 4-band GeoTIFFs")
    s.add_argument("input", help="field .json")
    s.add_argument("--out", required=True, help="output directory")
    s.add_argument("--prefix", default="currents")
    s.set_defaults(func=cmd_export_cogs)

    s = sub.add_parser("synthetic", help="write a deterministic synthetic field (offline)")
    s.add_argument("--nt", type=int, default=4)
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--out", default="synthetic_field", help="output path prefix")
    s.set_defaults(func=cmd_synthetic)

    s = sub.add_parser("fetch-glsea-sst",
                       help="fetch NOAA GLSEA daily SST grids (needs network)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(must lie inside the GLSEA lakes-region grid)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--stride-days", type=int, default=30,
                   help="time stride in days (default 30)")
    s.add_argument("--out", default="glsea_sst", help="output path prefix")
    s.set_defaults(func=cmd_fetch_glsea_sst)

    s = sub.add_parser("fetch-glsea-averages",
                       help="fetch NOAA GLSEA lake-average temps (needs network)")
    s.add_argument("--lake", required=True,
                   help="superior|michigan|huron|erie|ontario")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--out", default="glsea_averages", help="output path prefix")
    s.set_defaults(func=cmd_fetch_glsea_averages)

    s = sub.add_parser("glsea-synthetic",
                       help="write a deterministic synthetic GLSEA SST field (offline)")
    s.add_argument("--nt", type=int, default=4)
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--out", default="glsea_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_glsea_synthetic)

    s = sub.add_parser("fetch-oisst",
                       help="fetch NOAA OISST v2.1 daily global SST grids (needs network)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--stride-days", type=int, default=30,
                   help="time stride in days (default 30)")
    s.add_argument("--out", default="oisst_sst", help="output path prefix")
    s.set_defaults(func=cmd_fetch_oisst)

    s = sub.add_parser("fetch-mur",
                       help="fetch NASA JPL MUR v4.1 daily global SST grids "
                            "(needs network + Earthdata Login)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--stride-days", type=int, default=30,
                   help="granule sampling stride in days (default 30)")
    s.add_argument("--out", default="mur_sst", help="output path prefix")
    s.set_defaults(func=cmd_fetch_mur)

    s = sub.add_parser("fetch-nsidc",
                       help="fetch NSIDC G02135 v4.0 daily sea-ice concentration "
                            "GeoTIFFs (needs network; no account)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--hemisphere", default="auto",
                   help="north/south/auto (default auto from bbox latitudes)")
    s.add_argument("--stride-days", type=int, default=1,
                   help="keep every Nth day (default 1)")
    s.add_argument("--resolution", type=float, default=0.25,
                   help="target lat/lon grid resolution in degrees (default 0.25)")
    s.add_argument("--out", default="nsidc_ice", help="output path prefix")
    s.set_defaults(func=cmd_fetch_nsidc)

    s = sub.add_parser("ice-synthetic",
                       help="write a deterministic synthetic sea-ice field (offline)")
    s.add_argument("--bbox", default="-180,66,180,90",
                   help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--start", default="2024-02-01", help="start date YYYY-MM-DD")
    s.add_argument("--end", default="2024-02-05", help="end date YYYY-MM-DD")
    s.add_argument("--resolution", type=float, default=1.0,
                   help="target lat/lon grid resolution in degrees (default 1.0)")
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="ice_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_ice_synthetic)

    s = sub.add_parser("fetch-imerg",
                       help="fetch NASA GPM IMERG half-hourly precipitation "
                            "(needs network + free Earthdata Login)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--run", default="late",
                   help="early (~4 h latency) / late (~14 h, default) / "
                        "final (~3.5 months, gauge-adjusted)")
    s.add_argument("--accumulate", default="daily",
                   help="daily (default: 48 half-hourly mm/hr rates -> "
                        "mm/day totals) / native (half-hourly mm/hr rates)")
    s.add_argument("--stride-days", type=int, default=1,
                   help="keep every Nth day (default 1)")
    s.add_argument("--out", default="imerg_rain", help="output path prefix")
    s.set_defaults(func=cmd_fetch_imerg)

    s = sub.add_parser("rain-synthetic",
                       help="write a deterministic synthetic IMERG rain field (offline)")
    s.add_argument("--bbox", default="-125.0,25.0,-66.0,49.0",
                   help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--start", default="2024-01-01", help="start date YYYY-MM-DD")
    s.add_argument("--end", default="2024-01-05", help="end date YYYY-MM-DD")
    s.add_argument("--resolution", type=float, default=1.0,
                   help="synthetic grid resolution in degrees (default 1.0)")
    s.add_argument("--accumulate", default="daily",
                   help="daily (mm/day totals, default) / native (mm/hr rates)")
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="rain_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_rain_synthetic)

    s = sub.add_parser("fetch-blackmarble",
                       help="fetch NASA Black Marble daily night lights "
                            "(needs network + free Earthdata Login)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD "
                   "(record starts 2012-01-19)")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--product", default="daily",
                   help="daily (default: VNP46A2 V002 gap-filled NTL)")
    s.add_argument("--stride-days", type=int, default=1,
                   help="keep every Nth day (default 1)")
    s.add_argument("--resolution", type=float, default=0.05,
                   help="output grid spacing in degrees (default 0.05; "
                        "native 15 arc-second tiles are block-averaged)")
    s.add_argument("--out", default="blackmarble_lights", help="output path prefix")
    s.set_defaults(func=cmd_fetch_blackmarble)

    s = sub.add_parser("lights-synthetic",
                       help="write a deterministic synthetic night-lights field (offline)")
    s.add_argument("--bbox", default="-125.0,25.0,-66.0,49.0",
                   help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--start", default="2024-01-01", help="start date YYYY-MM-DD")
    s.add_argument("--end", default="2024-01-05", help="end date YYYY-MM-DD")
    s.add_argument("--resolution", type=float, default=1.0,
                   help="synthetic grid resolution in degrees (default 1.0)")
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="lights_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_lights_synthetic)

    s = sub.add_parser("fetch-grace",
                       help="fetch CSR GRACE/GRACE-FO RL06.3 terrestrial water "
                            "storage anomalies (needs network; keyless CSR "
                            "HTTPS; files cached)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD "
                   "(record starts 2002-04-01)")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--out", default="grace_water", help="output path prefix")
    s.set_defaults(func=cmd_fetch_grace)

    s = sub.add_parser("grace-synthetic",
                       help="write a deterministic synthetic water-storage "
                            "field (offline)")
    s.add_argument("--bbox", default="-125.0,30.0,-110.0,45.0",
                   help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--start", default="2020-01-01", help="start date YYYY-MM-DD")
    s.add_argument("--end", default="2020-12-01", help="end date YYYY-MM-DD")
    s.add_argument("--resolution", type=float, default=2.5,
                   help="synthetic grid resolution in degrees (default 2.5)")
    s.add_argument("--seed", type=int, default=11)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="grace_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_grace_synthetic)

    s = sub.add_parser("fetch-ibtracs",
                       help="fetch NOAA IBTrACS v4 tropical-cyclone best tracks "
                            "(needs network; keyless NCEI HTTPS; file cached)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD "
                   "(record starts 1980-01-01; 1842-01-01 with --full-archive)")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--min-wind", type=float, default=None,
                   help="keep storms whose lifetime max sustained wind "
                        "reaches this many kt (default: all)")
    s.add_argument("--storm-name", default=None,
                   help="select one named storm, e.g. katrina (case-insensitive)")
    s.add_argument("--full-archive", action="store_true",
                   help="use the 1842-present IBTrACS.ALL file "
                        "(default: 1980-present IBTrACS.since1980)")
    s.add_argument("--out", default="ibtracs_storms", help="output path prefix")
    s.set_defaults(func=cmd_fetch_ibtracs)

    s = sub.add_parser("storms-synthetic",
                       help="write a deterministic synthetic storm-track field (offline)")
    s.add_argument("--bbox", default="-100.0,10.0,-60.0,40.0",
                   help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--start", default="2024-08-01", help="start date YYYY-MM-DD")
    s.add_argument("--end", default="2024-09-15", help="end date YYYY-MM-DD")
    s.add_argument("--n-storms", type=int, default=3, help="storm count")
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="storms_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_storms_synthetic)

    s = sub.add_parser("sst-synthetic",
                       help="write a deterministic synthetic global SST field (offline)")
    s.add_argument("--nt", type=int, default=4)
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="sst_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_sst_synthetic)

    s = sub.add_parser("fetch-gebco",
                       help="fetch GEBCO 2024 topography/bathymetry subset "
                            "(needs network; tiles cached locally)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--resolution", default="0.25",
                   help="output grid spacing: '15s' native or degrees "
                        "(default 0.25; must be a multiple of 15 arc-seconds "
                        "that divides the 90° tiles, e.g. 0.25, 0.5, 1.0)")
    s.add_argument("--out", default="gebco_topo", help="output path prefix")
    s.set_defaults(func=cmd_fetch_gebco)

    s = sub.add_parser("fetch-naturalearth",
                       help="fetch Natural Earth coastline/countries vectors "
                            "(needs network; zips cached locally)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--scale", default="110m",
                   help="Natural Earth scale: 110m, 50m, or 10m (default 110m)")
    s.add_argument("--layers", default="coastline,countries",
                   help="comma-separated: coastline, countries (default both)")
    s.add_argument("--out", default="naturalearth", help="output path prefix")
    s.set_defaults(func=cmd_fetch_naturalearth)

    s = sub.add_parser("basemaps-synthetic",
                       help="write a deterministic synthetic topography field (offline)")
    s.add_argument("--bbox", default="-125.0,25.0,-66.0,49.0",
                   help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--resolution", type=float, default=1.0,
                   help="synthetic grid resolution in degrees (default 1.0)")
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="topo_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_basemaps_synthetic)

    s = sub.add_parser("fetch-era5",
                       help="fetch Copernicus ERA5 hourly reanalysis grids "
                            "(needs network + free CDS account)")
    s.add_argument("--variables", default="wind,msl",
                   help="comma-separated keys: wind,msl,t2m,tp (default wind,msl)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--stride-hours", type=int, default=6,
                   help="time stride in hours: 6 -> 00/06/12/18 UTC; "
                        "24 -> daily 12:00 UTC (default 6)")
    s.add_argument("--out", default="era5", help="output path prefix")
    s.set_defaults(func=cmd_fetch_era5)

    s = sub.add_parser("era5-synthetic",
                       help="write a deterministic synthetic ERA5 field (offline)")
    s.add_argument("--variables", default="wind,msl",
                   help="comma-separated keys: wind,msl,t2m,tp (default wind,msl)")
    s.add_argument("--nt", type=int, default=4)
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="era5_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_era5_synthetic)

    s = sub.add_parser("fetch-oscar",
                       help="fetch NASA PODAAC OSCAR v2.0 daily global surface "
                            "currents (needs network + Earthdata Login)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--stride-days", type=int, default=5,
                   help="granule sampling stride in days (default 5)")
    s.add_argument("--out", default="oscar_currents", help="output path prefix")
    s.set_defaults(func=cmd_fetch_oscar)

    s = sub.add_parser("fetch-cmems-currents",
                       help="fetch CMEMS global ocean physics daily currents "
                            "(needs network + free CMEMS account)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--stride-days", type=int, default=1,
                   help="keep every Nth daily timestep (default 1)")
    s.add_argument("--work-dir", default=None,
                   help="directory for the downloaded NetCDF "
                        "(default: a fresh temp dir)")
    s.add_argument("--out", default="cmems_currents", help="output path prefix")
    s.set_defaults(func=cmd_fetch_cmems_currents)

    s = sub.add_parser("currents-synthetic",
                       help="write deterministic synthetic OSCAR + CMEMS current "
                            "fields (offline)")
    s.add_argument("--nt", type=int, default=4)
    s.add_argument("--seed", type=int, default=13)
    s.add_argument("--out", default="currents_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_currents_synthetic)

    s = sub.add_parser("fetch-firms",
                       help="fetch NASA FIRMS active-fire detections "
                            "(needs network + free FIRMS MAP_KEY)")
    s.add_argument("--bbox", required=True, help="min_lon,min_lat,max_lon,max_lat "
                   "(conventional -180..180; antimeridian-crossing boxes wrap)")
    s.add_argument("--start", required=True, help="start date YYYY-MM-DD")
    s.add_argument("--end", required=True, help="end date YYYY-MM-DD")
    s.add_argument("--instruments", default="VIIRS_SNPP",
                   help="comma-separated: VIIRS_SNPP,VIIRS_NOAA20,"
                        "VIIRS_NOAA21,MODIS (default VIIRS_SNPP)")
    s.add_argument("--out", default="firms_fires", help="output path prefix")
    s.set_defaults(func=cmd_fetch_firms)

    s = sub.add_parser("fires-synthetic",
                       help="write deterministic synthetic FIRMS detections (offline)")
    s.add_argument("--bbox", default="-125.0,32.0,-114.0,42.0",
                   help="min_lon,min_lat,max_lon,max_lat")
    s.add_argument("--start", default="2024-08-01", help="start date YYYY-MM-DD")
    s.add_argument("--end", default="2024-08-07", help="end date YYYY-MM-DD")
    s.add_argument("--n", type=int, default=60, help="detection count")
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="fires_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_fires_synthetic)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (ImportError, RuntimeError, ValueError, KeyError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

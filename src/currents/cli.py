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


def cmd_sst_synthetic(a: argparse.Namespace) -> int:
    from .sst_global import SstField
    field = SstField.synthetic(nt=a.nt, seed=a.seed, source=a.source)
    path = f"{a.out.rstrip('/')}.json"
    field.to_json(path)
    print(f"wrote synthetic global SST field ({len(field.times)} steps) -> {path}")
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

    s = sub.add_parser("sst-synthetic",
                       help="write a deterministic synthetic global SST field (offline)")
    s.add_argument("--nt", type=int, default=4)
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--source", default="synthetic",
                   help="source label recorded on the field")
    s.add_argument("--out", default="sst_synthetic", help="output path prefix")
    s.set_defaults(func=cmd_sst_synthetic)

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

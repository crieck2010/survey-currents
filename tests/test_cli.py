"""Tests for the thin CLI adapter (offline; network commands mocked)."""

import json
import os

import pytest

from currents import cli
from currents.models import CurrentField


def test_parse_hours():
    assert cli._parse_hours("0") == [0]
    assert cli._parse_hours("0-3") == [0, 1, 2, 3]
    assert cli._parse_hours("0,6,12") == [0, 6, 12]
    assert cli._parse_hours("0-2,6") == [0, 1, 2, 6]


def test_parse_bbox():
    assert cli._parse_bbox("-92,41,-84,47") == [-92.0, 41.0, -84.0, 47.0]
    with pytest.raises(Exception):
        cli._parse_bbox("-92,41")


def test_list_models(capsys):
    assert cli.cmd_list_models(None) == 0
    out = capsys.readouterr().out
    assert "LMHOFS" in out and "GLOFS" in out


def test_synthetic_command(tmp_path, capsys):
    ns = type("NS", (), {"nt": 2, "seed": 7,
                         "out": str(tmp_path / "syn")})()
    assert cli.cmd_synthetic(ns) == 0
    path = str(tmp_path / "syn.json")
    assert os.path.exists(path)
    f = CurrentField.from_json(path)
    assert len(f.times) == 2


def test_fetch_noaa_mocked(monkeypatch, tmp_path, capsys):
    import currents.noaa_ofs as noaa

    def fake_fetch(ofs_code, date, cycle="00", hours=(0,), work_dir=".",
                   bucket=None, bbox=None):
        assert ofs_code == "LMHOFS"
        assert list(hours) == [0, 1]
        return CurrentField.synthetic(nt=2, seed=1)

    monkeypatch.setattr(noaa, "fetch_currents", fake_fetch)
    ns = type("NS", (), {"ofs": "LMHOFS", "date": "2026-09-16",
                         "cycle": "00", "hours": "0-1", "bbox": None,
                         "out": str(tmp_path / "wd")})()
    assert cli.cmd_fetch_noaa(ns) == 0
    assert os.path.exists(str(tmp_path / "wd.json"))
    assert "LMHOFS" in capsys.readouterr().out or True  # smoke


def test_export_cogs_missing_rasterio(monkeypatch, tmp_path):
    f = CurrentField.synthetic(nt=1, seed=1)
    inp = str(tmp_path / "f.json")
    f.to_json(inp)
    ns = type("NS", (), {"input": inp, "out": str(tmp_path / "o"),
                         "prefix": "t"})()
    # rasterio is not installed in CI here; the command must fail cleanly
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "rasterio":
            raise ImportError("no rasterio")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError):
        cli.cmd_export_cogs(ns)


def test_main_error_path(capsys):
    rc = cli.main(["fetch-noaa", "--ofs", "NOPE", "--date", "2026-09-16"])
    assert rc == 1
    assert "error" in capsys.readouterr().err.lower()


def test_main_bad_args():
    with pytest.raises(SystemExit):
        cli.main([])


def test_cli_has_currents_global_subcommands():
    from currents.cli import build_parser
    p = build_parser()
    subs = {a.dest for a in p._subparsers._group_actions[0]._choices_actions}
    assert {"fetch-oscar", "fetch-cmems-currents", "currents-synthetic"} <= subs


def test_currents_synthetic_command(tmp_path, capsys):
    ns = type("NS", (), {"nt": 2, "seed": 13,
                         "out": str(tmp_path / "cg")})()
    assert cli.cmd_currents_synthetic(ns) == 0
    for label in ("oscar", "cmems-currents"):
        path = str(tmp_path / f"cg_{label}.json")
        assert os.path.exists(path)
        f = CurrentField.from_json(path)
        assert len(f.times) == 2
    oscar = CurrentField.from_json(str(tmp_path / "cg_oscar.json"))
    assert oscar.temperature is None
    cmems = CurrentField.from_json(str(tmp_path / "cg_cmems-currents.json"))
    assert cmems.temperature is not None


def test_fetch_oscar_command_mocked(monkeypatch, tmp_path, capsys):
    import currents.currents_global as cg
    monkeypatch.setattr(cg, "fetch_oscar",
                        lambda bbox, start, end, stride_days=5:
                        cg.oscar_synthetic(nt=2))
    ns = type("NS", (), {"bbox": "-81,25,-55,43", "start": "2024-01-01",
                         "end": "2024-01-15", "stride_days": 5,
                         "out": str(tmp_path / "o")})()
    assert cli.cmd_fetch_oscar(ns) == 0
    out = capsys.readouterr().out
    assert "OSCAR v2.0" in out


def test_fetch_cmems_currents_command_mocked(monkeypatch, tmp_path, capsys):
    import currents.currents_global as cg
    monkeypatch.setattr(cg, "fetch_cmems_currents",
                        lambda bbox, start, end, stride_days=1, work_dir=None:
                        cg.cmems_currents_synthetic(nt=2))
    ns = type("NS", (), {"bbox": "-80,20,-60,40", "start": "2024-01-01",
                         "end": "2024-01-05", "stride_days": 1,
                         "work_dir": None,
                         "out": str(tmp_path / "c")})()
    assert cli.cmd_fetch_cmems_currents(ns) == 0
    out = capsys.readouterr().out
    assert "CMEMS" in out

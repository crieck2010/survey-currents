"""Tests for provenance sidecars."""

import json
import os

import pytest

from currents.provenance import (PROVENANCE_VERSION, read_provenance,
                                 sha256_file, sidecar_path,
                                 verify_provenance, write_provenance)


def test_write_and_read(tmp_path):
    data = tmp_path / "field.nc"
    data.write_bytes(b"payload")
    record = {"source": "s3://bucket/key.nc", "bbox": [-92, 41, -84, 47],
              "time_window": ["2026-09-16", "2026-09-17"]}
    sidecar = write_provenance(str(data), record)
    assert sidecar == sidecar_path(str(data))
    back = read_provenance(str(data))
    assert back["source"] == "s3://bucket/key.nc"
    assert back["provenance_version"] == PROVENANCE_VERSION
    assert back["tool"].startswith("survey-currents/")
    assert back["sha256"] == sha256_file(str(data))
    assert "downloaded_utc" in back


def test_write_does_not_clobber_explicit_fields(tmp_path):
    data = tmp_path / "f.nc"
    data.write_bytes(b"x")
    write_provenance(str(data), {"sha256": "pinned", "source": "s"})
    assert read_provenance(str(data))["sha256"] == "pinned"


def test_verify_ok(tmp_path):
    data = tmp_path / "f.nc"
    data.write_bytes(b"bytes")
    write_provenance(str(data), {"source": "s"})
    result = verify_provenance(str(data))
    assert result["ok"] is True


def test_verify_mismatch(tmp_path):
    data = tmp_path / "f.nc"
    data.write_bytes(b"bytes")
    write_provenance(str(data), {"source": "s"})
    data.write_bytes(b"tampered")
    result = verify_provenance(str(data))
    assert result["ok"] is False
    assert result["reason"] == "sha256 mismatch"


def test_verify_missing_sidecar(tmp_path):
    data = tmp_path / "f.nc"
    data.write_bytes(b"bytes")
    result = verify_provenance(str(data))
    assert result["ok"] is False
    assert "missing sidecar" in result["reason"]


def test_read_missing_sidecar(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_provenance(str(tmp_path / "nope.nc"))


def test_sidecar_is_json(tmp_path):
    data = tmp_path / "f.nc"
    data.write_bytes(b"bytes")
    sidecar = write_provenance(str(data), {"source": "s"})
    with open(sidecar) as fh:
        json.load(fh)  # must parse

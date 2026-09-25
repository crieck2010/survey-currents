"""Provenance sidecars: every downloaded file gets a JSON receipt.

A provenance record answers "exactly which bytes produced this field?"
six months later: source URL or dataset id, spatial/temporal window,
SHA-256 of the downloaded file, and the tool version that fetched it.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
from typing import Any, Dict, Optional

from . import __version__

PROVENANCE_VERSION = "1"
_SIDECAR_SUFFIX = ".provenance.json"


def sha256_file(path: str, chunk: int = 1 << 20) -> str:
    """Hex SHA-256 of a file on disk."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sidecar_path(data_path: str) -> str:
    """Provenance sidecar path for a downloaded data file."""
    return data_path + _SIDECAR_SUFFIX


def write_provenance(data_path: str, record: Dict[str, Any],
                     sidecar: Optional[str] = None) -> str:
    """Write the provenance sidecar for ``data_path``; return its path.

    ``record`` should carry at least ``source`` (URL or dataset id),
    ``bbox``, and ``time_window``. The writer stamps ``sha256``,
    ``downloaded_utc``, ``tool`` and ``provenance_version``.
    """
    record = dict(record)
    record.setdefault("sha256", sha256_file(data_path))
    record.setdefault("downloaded_utc",
                      _dt.datetime.now(_dt.timezone.utc).isoformat())
    record.setdefault("tool", f"survey-currents/{__version__}")
    record.setdefault("provenance_version", PROVENANCE_VERSION)
    out = sidecar or sidecar_path(data_path)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
    return out


def read_provenance(data_path: str,
                    sidecar: Optional[str] = None) -> Dict[str, Any]:
    """Read the provenance sidecar for ``data_path``."""
    path = sidecar or sidecar_path(data_path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"no provenance sidecar at {path!r}")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def verify_provenance(data_path: str,
                      sidecar: Optional[str] = None) -> Dict[str, Any]:
    """Re-hash ``data_path`` and compare with its sidecar.

    Returns ``{"ok": True, ...}`` on match, ``{"ok": False, "reason": ...}``
    on mismatch or a missing sidecar.
    """
    path = sidecar or sidecar_path(data_path)
    if not os.path.exists(path):
        return {"ok": False, "reason": f"missing sidecar {path!r}"}
    record = read_provenance(data_path, path)
    expected = record.get("sha256")
    if not expected:
        return {"ok": False, "reason": "sidecar has no sha256 field"}
    actual = sha256_file(data_path)
    if actual != expected:
        return {"ok": False, "reason": "sha256 mismatch",
                "expected": expected, "actual": actual}
    return {"ok": True, "sha256": actual, "record": record}

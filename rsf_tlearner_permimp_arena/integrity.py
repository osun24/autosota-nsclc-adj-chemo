"""Locked files, environment verification, hashes, and append-only ledger."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy
import pandas
import sklearn
import sksurv

from reactome_tlearner_arena import integrity as tlearner_integrity


ARENA_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = ARENA_DIR / "lock_manifest.json"
ENVIRONMENT_PATH = ARENA_DIR / "environment.json"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_sha256(value: object) -> str:
    return sha256_bytes(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def current_environment() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scikit-learn": sklearn.__version__,
        "scikit-survival": sksurv.__version__,
    }


def verify_environment() -> dict[str, str]:
    expected = json.loads(ENVIRONMENT_PATH.read_text(encoding="utf-8"))
    observed = current_environment()
    if observed != expected:
        raise RuntimeError(f"Environment lock mismatch: expected {expected}, observed {observed}")
    return observed


def verify_lock() -> None:
    tlearner_integrity.verify_lock()
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    body = {"schema_version": manifest["schema_version"], "files": manifest["files"]}
    if canonical_sha256(body) != manifest["manifest_sha256"]:
        raise RuntimeError("RSF T-learner PFI manifest is invalid")
    for relative, expected in manifest["files"].items():
        path = ARENA_DIR / relative
        if not path.exists() or sha256_file(path) != expected:
            raise RuntimeError(f"Locked RSF T-learner PFI file changed: {relative}")
    verify_environment()


def read_ledger(path: Path) -> list[dict]:
    if not path.exists() or not path.read_text(encoding="utf-8").strip():
        return []
    rows, previous = [], "GENESIS"
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stored = json.loads(line)
        supplied = stored.pop("row_sha256")
        if stored.get("previous_sha256") != previous:
            raise RuntimeError(f"Ledger chain broken at line {line_number}")
        actual = canonical_sha256(stored)
        if supplied != actual:
            raise RuntimeError(f"Ledger row hash mismatch at line {line_number}")
        stored["row_sha256"] = supplied
        rows.append(stored)
        previous = supplied
    return rows


def append_ledger(path: Path, payload: dict) -> dict:
    rows = read_ledger(path)
    row = {**payload, "previous_sha256": rows[-1]["row_sha256"] if rows else "GENESIS"}
    row["row_sha256"] = canonical_sha256(row)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(row, sort_keys=True) + "\n")
    return row

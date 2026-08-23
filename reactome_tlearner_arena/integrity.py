"""Integrity checks and append-only ledger for the T-learner arena."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from reactome_rsf_arena_v2 import integrity as v2_integrity


ARENA_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = ARENA_DIR / "lock_manifest.json"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def verify_lock() -> None:
    v2_integrity.verify_lock()
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    body = {"schema_version": manifest["schema_version"], "files": manifest["files"]}
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    if sha256_bytes(canonical) != manifest["manifest_sha256"]:
        raise RuntimeError("T-learner arena lock manifest is invalid")
    for relative, expected in manifest["files"].items():
        path = ARENA_DIR / relative
        if not path.exists() or sha256_file(path) != expected:
            raise RuntimeError(f"Locked T-learner arena file changed: {relative}")


def read_ledger(path: Path) -> list[dict]:
    if not path.exists() or not path.read_text(encoding="utf-8").strip():
        return []
    rows = []
    previous = "GENESIS"
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        row = json.loads(line)
        supplied = row.pop("row_sha256")
        if row.get("previous_sha256") != previous:
            raise RuntimeError(f"Ledger chain broken at line {line_number}")
        actual = sha256_bytes(json.dumps(row, sort_keys=True, separators=(",", ":")).encode())
        if supplied != actual:
            raise RuntimeError(f"Ledger row hash mismatch at line {line_number}")
        row["row_sha256"] = supplied
        rows.append(row)
        previous = supplied
    return rows


def append_ledger(path: Path, payload: dict) -> dict:
    rows = read_ledger(path)
    row = {**payload, "previous_sha256": rows[-1]["row_sha256"] if rows else "GENESIS"}
    canonical = json.dumps(row, sort_keys=True, separators=(",", ":")).encode()
    row["row_sha256"] = sha256_bytes(canonical)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(row, sort_keys=True) + "\n")
    return row


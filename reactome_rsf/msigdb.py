"""Download and verify the pinned MSigDB C2:CP:REACTOME symbol collection."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
from urllib.request import Request, urlopen

from .config import (
    DEFAULT_DATA_DIR,
    MSIGDB_COLLECTION,
    MSIGDB_FILENAME,
    MSIGDB_GENE_ID_TYPE,
    MSIGDB_RELEASE_CATALOG_URL,
    MSIGDB_RELEASE_NOTES_URL,
    MSIGDB_SHA256,
    MSIGDB_URL,
    MSIGDB_VERSION,
)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_gmt(path: str | Path) -> dict[str, tuple[str, ...]]:
    """Parse GMT into pathway -> ordered, de-duplicated gene symbols."""
    gene_sets: dict[str, tuple[str, ...]] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\n\r").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed GMT line {line_number}: expected >=3 fields")
            name = fields[0].strip()
            if not name:
                raise ValueError(f"Malformed GMT line {line_number}: empty set name")
            if name in gene_sets:
                raise ValueError(f"Duplicate GMT gene set name: {name}")
            genes = tuple(dict.fromkeys(g.strip() for g in fields[2:] if g.strip()))
            if not genes:
                raise ValueError(f"GMT gene set {name} has no gene members")
            gene_sets[name] = genes
    if not gene_sets:
        raise ValueError("GMT file contains no gene sets")
    return gene_sets


def _download_verified(destination: Path, force: bool = False) -> bool:
    """Return True after downloading, or False for a verified cache hit."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        actual = sha256_file(destination)
        if actual == MSIGDB_SHA256:
            return False
        if not force:
            raise RuntimeError(
                f"Cached GMT checksum mismatch at {destination}: {actual}; "
                "delete it or rerun with --force"
            )

    request = Request(MSIGDB_URL, headers={"User-Agent": "reactome-rsf/1.0"})
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", prefix=f".{MSIGDB_FILENAME}.", suffix=".part",
            dir=destination.parent, delete=False
        ) as temp_handle:
            temp_path = Path(temp_handle.name)
            with urlopen(request, timeout=120) as response:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    temp_handle.write(block)
        actual = sha256_file(temp_path)
        if actual != MSIGDB_SHA256:
            raise RuntimeError(
                f"Downloaded GMT checksum mismatch: expected {MSIGDB_SHA256}, got {actual}"
            )
        os.replace(temp_path, destination)
        temp_path = None
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
    return True


def retrieve_collection(data_dir: str | Path = DEFAULT_DATA_DIR, force: bool = False) -> dict:
    """Retrieve, verify, parse, and document the frozen MSigDB collection."""
    data_dir = Path(data_dir)
    gmt_path = data_dir / MSIGDB_FILENAME
    downloaded = _download_verified(gmt_path, force=force)
    gene_sets = parse_gmt(gmt_path)
    unique_genes = sorted({gene for members in gene_sets.values() for gene in members})

    genes_path = data_dir / f"reactome_unique_genes.v{MSIGDB_VERSION}.txt"
    genes_path.write_text("\n".join(unique_genes) + "\n", encoding="utf-8")
    verified_at = datetime.now(timezone.utc).isoformat()
    downloaded_at = datetime.fromtimestamp(
        gmt_path.stat().st_mtime, tz=timezone.utc
    ).isoformat()
    manifest = {
        "schema_version": 1,
        "source": "Molecular Signatures Database (MSigDB), Broad Institute",
        "collection": MSIGDB_COLLECTION,
        "msigdb_version": MSIGDB_VERSION,
        "gene_identifier_type": MSIGDB_GENE_ID_TYPE,
        "source_url": MSIGDB_URL,
        "release_catalog_url": MSIGDB_RELEASE_CATALOG_URL,
        "release_notes_url": MSIGDB_RELEASE_NOTES_URL,
        "gmt_filename": MSIGDB_FILENAME,
        "gmt_sha256": sha256_file(gmt_path),
        "gene_set_count": len(gene_sets),
        "unique_gene_count": len(unique_genes),
        "downloaded_at_utc": downloaded_at,
        "verified_at_utc": verified_at,
        "downloaded_this_run": downloaded,
        "license_note": "MSigDB content use is subject to the MSigDB license/terms.",
    }
    manifest_path = data_dir / "collection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {**manifest, "gmt_path": str(gmt_path), "genes_path": str(genes_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--force", action="store_true", help="replace a bad cached file")
    args = parser.parse_args()
    print(json.dumps(retrieve_collection(args.data_dir, force=args.force), indent=2))


if __name__ == "__main__":
    main()

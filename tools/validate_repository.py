#!/usr/bin/env python3
"""Fail fast on files that should never enter the public code repository."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAX_FILE_BYTES = 20 * 1024 * 1024
FORBIDDEN_SUFFIXES = {
    ".fastq",
    ".fq",
    ".rds",
    ".rdata",
    ".rda",
    ".biom",
    ".tif",
    ".tiff",
    ".gpkg",
    ".k2d",
}
REQUIRED = [
    ROOT / "README.md",
    ROOT / "LICENSE",
    ROOT / "CITATION.cff",
    ROOT / "DATA_AVAILABILITY.md",
    ROOT / "docs" / "provenance" / "SCRIPT_INDEX.tsv",
]


def tracked_files() -> list[Path]:
    return [
        path
        for path in ROOT.rglob("*")
        if path.is_file() and ".git" not in path.parts
    ]


def validate_script_index() -> list[str]:
    errors: list[str] = []
    index = ROOT / "docs" / "provenance" / "SCRIPT_INDEX.tsv"
    with index.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            canonical = row.get("canonical_file", "").strip()
            if not canonical:
                errors.append("SCRIPT_INDEX.tsv contains an empty canonical_file")
                continue
            if not (ROOT / canonical).is_file():
                errors.append(f"Indexed canonical file is missing: {canonical}")
    return errors


def main() -> int:
    errors: list[str] = []
    for required in REQUIRED:
        if not required.is_file():
            errors.append(f"Required file is missing: {required.relative_to(ROOT)}")

    for path in tracked_files():
        relative = path.relative_to(ROOT)
        lower_name = path.name.lower()
        if path.stat().st_size > MAX_FILE_BYTES:
            errors.append(f"File exceeds 20 MiB: {relative}")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES or lower_name.endswith(
            (".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz")
        ):
            errors.append(f"Forbidden data type in repository: {relative}")

    errors.extend(validate_script_index())

    with (ROOT / ".zenodo.json").open(encoding="utf-8") as handle:
        json.load(handle)

    if errors:
        print("REPOSITORY_VALIDATION_FAIL")
        for error in errors:
            print(f"- {error}")
        return 1

    print(f"REPOSITORY_VALIDATION_PASS: {len(tracked_files())} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


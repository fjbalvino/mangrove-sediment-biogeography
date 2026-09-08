#!/usr/bin/env python3
"""Freeze the 401-sample metadata cohort used by taxonomic analyses.

The four Giralia accessions are excluded by exact ``run_accession`` match.
No fuzzy project, locality, or substring filtering is permitted.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from metadata_pipeline_common import (
    atomic_write_csv,
    validate_coordinates,
    validate_primary_key,
    write_manifest,
)


DEFAULT_EXCLUDED = (
    "SRR29659508",
    "SRR29659599",
    "SRR29659629",
    "SRR29659854",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare the frozen 401-sample MAGENTA metadata cohort.")
    parser.add_argument("--base", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-out", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--lat-col", default="lat")
    parser.add_argument("--lon-col", default="lon")
    parser.add_argument("--expected-input-rows", type=int, default=405)
    parser.add_argument("--expected-output-rows", type=int, default=401)
    parser.add_argument("--exclude", nargs="*", default=list(DEFAULT_EXCLUDED))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base_path = Path(args.base)
    out_path = Path(args.out)
    audit_path = Path(args.audit_out)
    manifest_path = Path(args.manifest_out)
    if not base_path.is_file():
        raise FileNotFoundError(base_path)
    if out_path.resolve() == base_path.resolve():
        raise ValueError("--out debe ser distinto de --base")

    metadata = pd.read_csv(base_path, low_memory=False)
    validate_primary_key(metadata, args.id_col, "base")
    validate_coordinates(metadata, args.lat_col, args.lon_col, "base")
    if args.expected_input_rows > 0 and len(metadata) != args.expected_input_rows:
        raise ValueError(
            f"Filas de entrada observadas={len(metadata)}; esperadas={args.expected_input_rows}"
        )

    identifiers = metadata[args.id_col].astype("string").str.strip()
    excluded_set = set(args.exclude)
    observed = set(identifiers[identifiers.isin(excluded_set)].tolist())
    missing = sorted(excluded_set - observed)
    if missing:
        raise ValueError(f"Accesiones de exclusión ausentes en metadata: {missing}")

    excluded = metadata.loc[identifiers.isin(excluded_set)].copy()
    retained = metadata.loc[~identifiers.isin(excluded_set)].copy()
    retained = retained.sort_values(args.id_col, kind="stable").reset_index(drop=True)
    if args.expected_output_rows > 0 and len(retained) != args.expected_output_rows:
        raise ValueError(
            f"Filas retenidas observadas={len(retained)}; esperadas={args.expected_output_rows}"
        )
    validate_primary_key(retained, args.id_col, "set401")

    audit = pd.DataFrame(
        {
            args.id_col: identifiers,
            "excluded_from_set401": identifiers.isin(excluded_set),
            "filter_reason": [
                "exact_Giralia4_exclusion" if value else "retained"
                for value in identifiers.isin(excluded_set)
            ],
        }
    )
    atomic_write_csv(retained, out_path)
    atomic_write_csv(audit, audit_path)
    write_manifest(
        manifest_path,
        Path(__file__).name,
        args,
        [base_path],
        out_path,
        len(retained),
        len(metadata),
    )
    print(f"[OK] entrada: {len(metadata)}")
    print(f"[OK] excluidas: {len(excluded)}")
    print(f"[OK] set401: {len(retained)}")
    print(f"[OK] metadata: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

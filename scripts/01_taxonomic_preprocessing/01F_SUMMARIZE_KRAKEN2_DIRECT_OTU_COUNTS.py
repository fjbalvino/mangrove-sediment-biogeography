#!/usr/bin/env python3
"""Build a non-overlapping taxonomic OTU count table from Kraken2 reports.

Kraken2's third report column is the number of reads assigned directly to a
taxonomic node.  Summing these direct counts across non-root nodes counts each
classified read at most once, unlike clade counts, which are nested.  The
phyloseq builder subsequently uses the NCBI taxdump to retain Bacteria and
Archaea and reconstruct complete lineages.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Dict, Iterable, Iterator, Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize direct Kraken2 OTU counts.")
    parser.add_argument("--kraken-manifest", required=True)
    parser.add_argument("--out", required=True, help="Long TSV or TSV.GZ count table.")
    parser.add_argument("--audit-out", required=True)
    return parser.parse_args()


def text_reader(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def parse_report(path: Path) -> Tuple[Dict[str, int], Dict[str, int]]:
    direct_counts: Dict[str, int] = {}
    classified_from_root = 0
    unclassified = 0
    with text_reader(path) as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_number}: reporte Kraken2 con menos de seis columnas")
            try:
                clade_count = int(fields[1])
                direct_count = int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: conteo no entero") from error
            rank_code = fields[3].strip()
            taxid = fields[4].strip()
            if direct_count < 0 or clade_count < 0:
                raise ValueError(f"{path}:{line_number}: conteo negativo")
            if rank_code == "U" or taxid == "0":
                unclassified = max(unclassified, clade_count)
            if taxid == "1":
                classified_from_root = max(classified_from_root, clade_count)
            if taxid not in {"", "0", "1"} and direct_count > 0:
                direct_counts[taxid] = direct_counts.get(taxid, 0) + direct_count
    audit = {
        "classified_reads_report_root": classified_from_root,
        "unclassified_reads_report": unclassified,
        "direct_nonroot_reads": sum(direct_counts.values()),
        "n_direct_nonzero_taxa": len(direct_counts),
    }
    if sum(direct_counts.values()) > classified_from_root:
        raise ValueError(
            f"{path}: la suma de conteos directos excede el total clasificado del nodo raíz"
        )
    return direct_counts, audit


def output_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("wt", encoding="utf-8", newline="")


def main() -> int:
    args = parse_args()
    manifest_path = Path(args.kraken_manifest)
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = pd.read_csv(manifest_path, dtype="string")
    required = {"run_accession", "report"}
    if not required.issubset(manifest.columns):
        raise ValueError(f"kraken-manifest requiere {sorted(required)}")
    run_ids = manifest["run_accession"].str.strip()
    if run_ids.isna().any() or run_ids.eq("").any() or run_ids.duplicated().any():
        raise ValueError("run_accession debe ser completo y único")
    manifest = manifest.assign(run_accession=run_ids).sort_values("run_accession")

    out_path = Path(args.out)
    audit_rows = []
    total_rows = 0
    with output_writer(out_path) as output:
        output.write("run_accession\ttaxonomy_id\tcount\n")
        for row in manifest.itertuples(index=False):
            run = str(row.run_accession)
            report = Path(str(row.report))
            if not report.is_file():
                raise FileNotFoundError(f"{run}: reporte ausente: {report}")
            counts, audit = parse_report(report)
            for taxid in sorted(counts, key=lambda value: (len(value), value)):
                output.write(f"{run}\t{taxid}\t{counts[taxid]}\n")
                total_rows += 1
            audit_rows.append({"run_accession": run, "report": str(report), **audit})
            print(f"[PASS] {run}: {len(counts)} taxones directos", flush=True)

    audit_table = pd.DataFrame(audit_rows)
    audit_path = Path(args.audit_out)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_table.to_csv(audit_path, sep="\t", index=False)
    if len(audit_table) != len(manifest):
        raise RuntimeError("No se auditaron todas las muestras")
    print(f"[OK] muestras: {len(audit_table)}")
    print(f"[OK] filas taxón-muestra: {total_rows}")
    print(f"[OK] tabla: {out_path}")
    print(f"[OK] auditoría: {audit_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

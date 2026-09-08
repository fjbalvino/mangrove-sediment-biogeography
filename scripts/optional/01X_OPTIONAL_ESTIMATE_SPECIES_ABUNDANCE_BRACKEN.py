#!/usr/bin/env python3
"""Estimate species-level counts from Kraken2 reports with Bracken.

This optional stage is not called by the canonical workflow because the
current Methods text mentions only direct Kraken2 profiling. If it is used for
a separate sensitivity analysis, Bracken, its version, read length, database
distribution file, and abundance-estimation step must be documented.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Bracken species abundance estimation.")
    parser.add_argument("--kraken-manifest", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--threshold", type=int, default=0)
    parser.add_argument("--level", default="S", choices=["D", "P", "C", "O", "F", "G", "S"])
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def tool_version(binary: str) -> str:
    result = subprocess.run([binary, "-v"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    match = re.search(r"(\d+\.\d+(?:\.\d+)?)", result.stdout)
    return match.group(1) if match else "unknown"


def run_one(row: pd.Series, args: argparse.Namespace, binary: str) -> Dict[str, Any]:
    run = str(row["run_accession"])
    report = Path(str(row["report"]))
    output = Path(args.out_dir) / f"{run}.bracken.{args.level}.tsv"
    log = Path(args.out_dir) / f"{run}.bracken.log"
    if not report.is_file():
        raise FileNotFoundError(f"{run}: reporte Kraken2 ausente: {report}")
    if output.is_file() and output.stat().st_size > 0 and not args.overwrite:
        table = pd.read_csv(output, sep="\t")
        return {
            "run_accession": run,
            "bracken_table": str(output),
            "n_taxa": len(table),
            "estimated_reads": int(pd.to_numeric(table.get("new_est_reads"), errors="coerce").fillna(0).sum()),
            "status": "PASS_EXISTING",
        }
    command = [
        binary,
        "-d", str(args.db),
        "-i", str(report),
        "-o", str(output),
        "-r", str(args.read_length),
        "-l", args.level,
        "-t", str(args.threshold),
    ]
    with log.open("w", encoding="utf-8") as stream:
        stream.write("COMMAND\t" + " ".join(command) + "\n")
        stream.flush()
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
        return {"run_accession": run, "bracken_table": str(output), "n_taxa": 0, "estimated_reads": 0, "status": "FAIL"}
    table = pd.read_csv(output, sep="\t")
    required = {"name", "taxonomy_id", "taxonomy_lvl", "new_est_reads"}
    if not required.issubset(table.columns):
        raise ValueError(f"{run}: tabla Bracken sin columnas {sorted(required)}")
    return {
        "run_accession": run,
        "bracken_table": str(output),
        "n_taxa": len(table),
        "estimated_reads": int(pd.to_numeric(table["new_est_reads"], errors="coerce").fillna(0).sum()),
        "status": "PASS",
    }


def main() -> int:
    args = parse_args()
    if args.read_length < 1 or args.threshold < 0 or args.max_workers < 1:
        raise ValueError("read-length/workers deben ser positivos y threshold >= 0")
    binary = shutil.which("bracken")
    if not binary:
        raise RuntimeError("bracken no está en PATH")
    db = Path(args.db)
    distribution = db / f"database{args.read_length}mers.kmer_distrib"
    if not distribution.is_file():
        raise FileNotFoundError(
            f"La DB no contiene {distribution.name}. Ejecute bracken-build para esta longitud o cambie --read-length."
        )
    kraken = pd.read_csv(args.kraken_manifest, dtype="string")
    if not {"run_accession", "report"}.issubset(kraken.columns):
        raise ValueError("kraken-manifest requiere run_accession y report")
    if kraken["run_accession"].duplicated().any():
        raise ValueError("run_accession duplicado en kraken-manifest")
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    jobs = {}
    results: List[Dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        for _, row in kraken.iterrows():
            jobs[executor.submit(run_one, row, args, binary)] = str(row["run_accession"])
        for future in as_completed(jobs):
            result = future.result()
            results.append(result)
            print(f"[{result['status']}] {result['run_accession']}", flush=True)
    manifest = pd.DataFrame(results).sort_values("run_accession")
    output_path = Path(args.manifest_out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_path, index=False)
    if manifest["status"].str.startswith("PASS").sum() != len(kraken):
        raise RuntimeError(f"Bracken incompleto; revise {output_path}")
    print(f"[OK] Bracken {tool_version(binary)}: {len(manifest)} corridas")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

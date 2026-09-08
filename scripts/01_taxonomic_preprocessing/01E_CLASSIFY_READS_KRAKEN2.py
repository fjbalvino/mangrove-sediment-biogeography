#!/usr/bin/env python3
"""Classify trimmed MAGENTA reads with Kraken2 and the frozen PlusPF database."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Kraken2 taxonomic classification for MAGENTA reads.")
    parser.add_argument("--trimmed-manifest", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--threads-per-sample", type=int, default=40)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--confidence", type=float, default=0.0)
    parser.add_argument("--minimum-hit-groups", type=int, default=2)
    parser.add_argument("--require-version", default="2.1.2")
    parser.add_argument("--save-classifications", action="store_true")
    parser.add_argument("--memory-mapping", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def tool_version(binary: str) -> str:
    result = subprocess.run([binary, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    match = re.search(r"(?:version\s+)?(\d+\.\d+(?:\.\d+)?)", result.stdout, re.I)
    return match.group(1) if match else "unknown"


def sha256_small(path: Path, limit: int = 512 * 1024**2) -> str:
    if not path.is_file() or path.stat().st_size > limit:
        return "not_calculated"
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_database(db: Path) -> Dict[str, Any]:
    required = [db / "hash.k2d", db / "opts.k2d", db / "taxo.k2d"]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Kraken2 DB incompleta: " + ", ".join(missing))
    return {
        "kraken2_db": str(db),
        "hash_k2d_size": required[0].stat().st_size,
        "opts_k2d_sha256": sha256_small(required[1]),
        "taxo_k2d_sha256": sha256_small(required[2]),
    }


def report_counts(path: Path) -> Tuple[int, int, int]:
    classified = 0
    unclassified = 0
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            try:
                clade_reads = int(fields[1])
            except ValueError:
                continue
            rank = fields[3].strip()
            taxid = fields[4].strip()
            if rank == "U" or taxid == "0":
                unclassified = max(unclassified, clade_reads)
            elif taxid == "1":
                classified = max(classified, clade_reads)
    return classified, unclassified, classified + unclassified


def classify_one(row: pd.Series, args: argparse.Namespace, binary: str, db_audit: Dict[str, Any]) -> Dict[str, Any]:
    run = str(row["run_accession"])
    r1 = Path(str(row["r1"]))
    r2: Optional[Path] = None if pd.isna(row["r2"]) or not str(row["r2"]).strip() else Path(str(row["r2"]))
    if not r1.is_file() or (r2 is not None and not r2.is_file()):
        raise FileNotFoundError(f"{run}: FASTQ trimmed ausente")
    out_dir = Path(args.out_dir)
    report = out_dir / f"{run}.kraken.report.tsv"
    log = out_dir / f"{run}.kraken.log"
    classifications = out_dir / f"{run}.kraken.classifications.tsv" if args.save_classifications else None
    if report.is_file() and report.stat().st_size > 0 and not args.overwrite:
        classified, unclassified, total = report_counts(report)
        return {
            "run_accession": run,
            "layout": "PAIRED" if r2 is not None else "SINGLE",
            "report": str(report),
            "classifications": str(classifications) if classifications and classifications.is_file() else "",
            "classified_reads": classified,
            "unclassified_reads": unclassified,
            "total_reads": total,
            "classified_fraction": classified / total if total else float("nan"),
            "status": "PASS_EXISTING",
            **db_audit,
        }
    command: List[str] = [
        binary,
        "--db", str(args.db),
        "--threads", str(args.threads_per_sample),
        "--confidence", str(args.confidence),
        "--minimum-hit-groups", str(args.minimum_hit_groups),
        "--gzip-compressed",
        "--use-names",
        "--report", str(report),
        "--output", str(classifications) if classifications else "/dev/null",
    ]
    if args.memory_mapping:
        command.append("--memory-mapping")
    if r2 is not None:
        command.extend(["--paired", str(r1), str(r2)])
    else:
        command.append(str(r1))
    with log.open("w", encoding="utf-8") as stream:
        stream.write("COMMAND\t" + " ".join(command) + "\n")
        stream.flush()
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode != 0 or not report.is_file() or report.stat().st_size == 0:
        return {"run_accession": run, "status": "FAIL", "exit_code": result.returncode, **db_audit}
    classified, unclassified, total = report_counts(report)
    return {
        "run_accession": run,
        "layout": "PAIRED" if r2 is not None else "SINGLE",
        "report": str(report),
        "classifications": str(classifications) if classifications else "",
        "classified_reads": classified,
        "unclassified_reads": unclassified,
        "total_reads": total,
        "classified_fraction": classified / total if total else float("nan"),
        "status": "PASS",
        **db_audit,
    }


def main() -> int:
    args = parse_args()
    if args.threads_per_sample < 1 or args.max_workers < 1:
        raise ValueError("threads/workers deben ser >= 1")
    if not 0 <= args.confidence <= 1:
        raise ValueError("--confidence debe estar entre 0 y 1")
    binary = shutil.which("kraken2")
    if not binary:
        raise RuntimeError("kraken2 no está en PATH")
    observed_version = tool_version(binary)
    if args.require_version and observed_version != args.require_version:
        raise RuntimeError(f"Kraken2 {observed_version}; el manuscrito requiere {args.require_version}")
    db = Path(args.db)
    db_audit = audit_database(db)
    trimmed = pd.read_csv(args.trimmed_manifest, dtype="string")
    required = {"run_accession", "r1", "r2"}
    if not required.issubset(trimmed.columns):
        raise ValueError(f"trimmed-manifest requiere columnas {sorted(required)}")
    if trimmed["run_accession"].duplicated().any():
        raise ValueError("run_accession duplicado en trimmed-manifest")
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    jobs = {}
    results: List[Dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        for _, row in trimmed.iterrows():
            jobs[executor.submit(classify_one, row, args, binary, db_audit)] = str(row["run_accession"])
        for future in as_completed(jobs):
            result = future.result()
            results.append(result)
            print(f"[{result['status']}] {result['run_accession']}", flush=True)
    manifest = pd.DataFrame(results).sort_values("run_accession")
    manifest_path = Path(args.manifest_out)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, index=False)
    if manifest["status"].str.startswith("PASS").sum() != len(trimmed):
        raise RuntimeError(f"Kraken2 incompleto; revise {manifest_path}")
    print(f"[OK] Kraken2 {observed_version}: {len(manifest)} corridas")
    print(f"[OK] manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

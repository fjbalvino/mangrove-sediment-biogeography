#!/usr/bin/env python3
"""Trim MAGENTA reads with an explicit, auditable fastp parameterization.

The executable command uses Q15 as the definition of an unqualified base,
requires at least 50 nt after filtering, permits no more than five N bases per
read, and performs automatic adapter detection for paired-end libraries.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Trim MAGENTA reads using fastp.")
    parser.add_argument("--raw-manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--report-dir", required=True)
    parser.add_argument("--threads-per-sample", type=int, default=8)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--qualified-quality-phred", type=int, default=15)
    parser.add_argument("--unqualified-percent-limit", type=int, default=40)
    parser.add_argument("--length-required", type=int, default=50)
    parser.add_argument("--n-base-limit", type=int, default=5)
    parser.add_argument("--require-version", default="0.23.4")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def fastp_version(binary: str) -> str:
    result = subprocess.run([binary, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    match = re.search(r"(\d+\.\d+(?:\.\d+)?)", result.stdout)
    return match.group(1) if match else "unknown"


def nonempty(path: Optional[Path]) -> bool:
    return path is not None and path.is_file() and path.stat().st_size > 0


def trim_one(row: pd.Series, args: argparse.Namespace, binary: str) -> Dict[str, Any]:
    run = str(row["run_accession"])
    r1 = Path(str(row["r1"]))
    r2 = None if pd.isna(row["r2"]) or not str(row["r2"]).strip() else Path(str(row["r2"]))
    out_dir = Path(args.out_dir)
    report_dir = Path(args.report_dir)
    r1_out = out_dir / f"{run}_1.trimmed.fastq.gz"
    r2_out = out_dir / f"{run}_2.trimmed.fastq.gz" if r2 is not None else None
    json_out = report_dir / f"{run}.fastp.json"
    html_out = report_dir / f"{run}.fastp.html"
    log_out = report_dir / f"{run}.fastp.log"
    if not r1.is_file() or (r2 is not None and not r2.is_file()):
        raise FileNotFoundError(f"{run}: faltan FASTQ crudos")
    if not args.overwrite and nonempty(r1_out) and (r2_out is None or nonempty(r2_out)) and nonempty(json_out):
        return {
            "run_accession": run,
            "layout": "PAIRED" if r2 is not None else "SINGLE",
            "r1": str(r1_out),
            "r2": str(r2_out) if r2_out is not None else "",
            "fastp_json": str(json_out),
            "fastp_html": str(html_out),
            "status": "PASS_EXISTING",
        }
    command: List[str] = [
        binary,
        "--thread", str(args.threads_per_sample),
        "--in1", str(r1),
        "--out1", str(r1_out),
        "--qualified_quality_phred", str(args.qualified_quality_phred),
        "--unqualified_percent_limit", str(args.unqualified_percent_limit),
        "--length_required", str(args.length_required),
        "--n_base_limit", str(args.n_base_limit),
        "--json", str(json_out),
        "--html", str(html_out),
    ]
    if r2 is not None:
        command.extend(["--in2", str(r2), "--out2", str(r2_out), "--detect_adapter_for_pe"])
    with log_out.open("w", encoding="utf-8") as log:
        log.write("COMMAND\t" + " ".join(command) + "\n")
        log.flush()
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    ok = result.returncode == 0 and nonempty(r1_out) and (r2_out is None or nonempty(r2_out)) and nonempty(json_out)
    return {
        "run_accession": run,
        "layout": "PAIRED" if r2 is not None else "SINGLE",
        "r1": str(r1_out) if nonempty(r1_out) else "",
        "r2": str(r2_out) if nonempty(r2_out) else "",
        "fastp_json": str(json_out) if nonempty(json_out) else "",
        "fastp_html": str(html_out) if nonempty(html_out) else "",
        "fastp_exit_code": result.returncode,
        "status": "PASS" if ok else "FAIL",
    }


def main() -> int:
    args = parse_args()
    if args.threads_per_sample < 1 or args.max_workers < 1:
        raise ValueError("threads/workers deben ser >= 1")
    binary = shutil.which("fastp")
    if not binary:
        raise RuntimeError("fastp no está en PATH")
    observed_version = fastp_version(binary)
    if args.require_version and observed_version != args.require_version:
        raise RuntimeError(f"fastp {observed_version}; el manuscrito requiere {args.require_version}")
    raw = pd.read_csv(args.raw_manifest, dtype="string")
    required = {"run_accession", "r1", "r2"}
    if not required.issubset(raw.columns):
        raise ValueError(f"raw-manifest requiere columnas {sorted(required)}")
    if raw["run_accession"].duplicated().any():
        raise ValueError("run_accession duplicado en raw-manifest")
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    Path(args.report_dir).mkdir(parents=True, exist_ok=True)
    jobs = {}
    results: List[Dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        for _, row in raw.iterrows():
            jobs[executor.submit(trim_one, row, args, binary)] = str(row["run_accession"])
        for future in as_completed(jobs):
            result = future.result()
            results.append(result)
            print(f"[{result['status']}] {result['run_accession']}", flush=True)
    manifest = pd.DataFrame(results).sort_values("run_accession")
    output_path = Path(args.manifest_out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_path, index=False)
    if manifest["status"].str.startswith("PASS").sum() != len(raw):
        raise RuntimeError(f"fastp incompleto; revise {output_path}")
    print(f"[OK] fastp {observed_version}: {len(manifest)} corridas")
    print(f"[OK] manifiesto trimmed: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

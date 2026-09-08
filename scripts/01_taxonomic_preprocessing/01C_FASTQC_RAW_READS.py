#!/usr/bin/env python3
"""Run FastQC v0.11.9-compatible QC on the frozen raw-read manifest."""

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
    parser = argparse.ArgumentParser(description="Run FastQC and optional MultiQC on raw MAGENTA reads.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--status-out", required=True)
    parser.add_argument("--threads-per-job", type=int, default=2)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--require-version", default="0.11.9")
    parser.add_argument("--skip-multiqc", action="store_true")
    return parser.parse_args()


def version(binary: str) -> str:
    result = subprocess.run([binary, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    match = re.search(r"(\d+\.\d+(?:\.\d+)?)", result.stdout)
    return match.group(1) if match else "unknown"


def run_fastqc(run: str, files: List[str], out_dir: Path, threads: int, binary: str) -> Dict[str, Any]:
    command = [binary, "--threads", str(threads), "--outdir", str(out_dir), *files]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    expected = [out_dir / (Path(path).name.removesuffix(".gz").removesuffix(".fastq") + "_fastqc.zip") for path in files]
    ok = result.returncode == 0 and all(path.is_file() for path in expected)
    return {
        "run_accession": run,
        "n_fastq": len(files),
        "fastqc_exit_code": result.returncode,
        "outputs_present": ok,
        "status": "PASS" if ok else "FAIL",
        "stderr_tail": result.stderr[-1000:].replace("\n", " | "),
    }


def main() -> int:
    args = parse_args()
    if args.threads_per_job < 1 or args.max_workers < 1:
        raise ValueError("threads/workers deben ser >= 1")
    binary = shutil.which("fastqc")
    if not binary:
        raise RuntimeError("fastqc no está en PATH")
    observed_version = version(binary)
    if args.require_version and observed_version != args.require_version:
        raise RuntimeError(f"FastQC {observed_version}; el manuscrito requiere {args.require_version}")

    manifest = pd.read_csv(args.manifest, dtype="string")
    required = {"run_accession", "r1", "r2"}
    if not required.issubset(manifest.columns):
        raise ValueError(f"El manifiesto requiere columnas {sorted(required)}")
    if manifest["run_accession"].duplicated().any():
        raise ValueError("run_accession duplicado en el manifiesto")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = {}
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        for _, row in manifest.iterrows():
            files = [str(row["r1"])]
            if pd.notna(row["r2"]) and str(row["r2"]).strip():
                files.append(str(row["r2"]))
            missing = [path for path in files if not Path(path).is_file()]
            if missing:
                raise FileNotFoundError(f"{row['run_accession']}: FASTQ ausentes: {missing}")
            future = executor.submit(run_fastqc, str(row["run_accession"]), files, out_dir, args.threads_per_job, binary)
            jobs[future] = str(row["run_accession"])
        results: List[Dict[str, Any]] = []
        for future in as_completed(jobs):
            result = future.result()
            results.append(result)
            print(f"[{result['status']}] {result['run_accession']}", flush=True)
    status = pd.DataFrame(results).sort_values("run_accession")
    status_path = Path(args.status_out)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status.to_csv(status_path, index=False)
    if not (status["status"] == "PASS").all():
        raise RuntimeError(f"FastQC falló; revise {status_path}")
    if not args.skip_multiqc:
        multiqc = shutil.which("multiqc")
        if multiqc:
            result = subprocess.run([multiqc, str(out_dir), "--outdir", str(out_dir), "--force"])
            if result.returncode != 0:
                raise RuntimeError("MultiQC falló")
        else:
            print("[WARN] multiqc no está en PATH; FastQC individual está completo")
    print(f"[OK] FastQC {observed_version}: {len(status)} corridas")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

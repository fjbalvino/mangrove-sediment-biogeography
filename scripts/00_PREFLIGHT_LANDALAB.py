#!/usr/bin/env python3
"""Read-only preflight for the MAGENTA taxonomic/environmental package."""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


EXPECTED_VERSIONS = {"fastqc": "0.11.9", "fastp": "0.23.4", "kraken2": "2.1.2"}


def parse_env(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip("\"").strip("'")
    return values


def version_of(binary: str) -> str:
    for flag in ("--version", "-v"):
        result = subprocess.run([binary, flag], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        match = re.search(r"(\d+\.\d+(?:\.\d+)?)", result.stdout)
        if match:
            return match.group(1)
    return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(description="MAGENTA Landa Lab preflight")
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--scope", choices=["taxonomic", "environmental", "all"], default="all")
    parser.add_argument("--allow-missing-environmental-assets", action="store_true")
    args = parser.parse_args()
    config_path = Path(args.config)
    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    config = parse_env(config_path)
    checks: List[Dict[str, Any]] = []
    taxonomic_scope = args.scope in {"taxonomic", "all"}
    environmental_scope = args.scope in {"environmental", "all"}

    def add(check: str, passed: bool, observed: Any, expected: Any, severity: str = "error") -> None:
        checks.append(
            {
                "check": check,
                "pass": bool(passed),
                "observed": str(observed),
                "expected": str(expected),
                "severity": severity,
            }
        )

    common_keys = ["BIOGEOGRAPHY_BASE", "MAGENTA_METADATA_CURATED", "PYTHON_BIN"]
    taxonomic_keys = ["MAGENTA_ROOT", "KRAKEN2_DB", "NCBI_TAXDUMP_DIR", "R_BIN"] if taxonomic_scope else []
    for key in common_keys + taxonomic_keys:
        value = config.get(key, "")
        add(f"config_{key}", bool(value), value or "EMPTY", "non-empty")

    directory_keys = ["BIOGEOGRAPHY_BASE"]
    if taxonomic_scope:
        directory_keys.extend(["MAGENTA_ROOT", "KRAKEN2_DB", "NCBI_TAXDUMP_DIR"])
    for key in directory_keys:
        path = Path(config.get(key, ""))
        add(f"directory_{key}", path.is_dir(), path, "existing directory")
    metadata = Path(config.get("MAGENTA_METADATA_CURATED", ""))
    add("curated_metadata", metadata.is_file(), metadata, "existing file")

    if taxonomic_scope:
        db = Path(config.get("KRAKEN2_DB", ""))
        for filename in ["hash.k2d", "opts.k2d", "taxo.k2d"]:
            add(f"kraken_db_{filename}", (db / filename).is_file(), db / filename, "existing file")
        taxdump = Path(config.get("NCBI_TAXDUMP_DIR", ""))
        for filename in ["names.dmp", "nodes.dmp"]:
            candidates = [taxdump / filename, taxdump / "taxonomy" / filename, taxdump / "taxdump" / filename]
            matches = [path for path in candidates if path.is_file()]
            add(f"taxdump_{filename}", len(matches) == 1, matches or candidates, "one existing file")

        for binary_name, expected in EXPECTED_VERSIONS.items():
            binary = shutil.which(binary_name)
            observed = version_of(binary) if binary else "missing"
            add(f"binary_{binary_name}", binary is not None, binary or "missing", "available in PATH")
            add(f"version_{binary_name}", observed == expected, observed, expected)
        binary = shutil.which("gzip")
        add("binary_gzip", binary is not None, binary or "missing", "available in PATH")
        download_clients = [name for name in ["aria2c", "curl"] if shutil.which(name)]
        add(
            "download_client",
            bool(download_clients),
            ",".join(download_clients) if download_clients else "missing",
            "aria2c or curl available in PATH",
        )

    python_bin = config.get("PYTHON_BIN", "python3")
    python_imports = "pandas,numpy,rasterio,geopandas,tqdm" if environmental_scope else "pandas,numpy"
    python_command = [
        python_bin,
        "-c",
        f"import {python_imports}; print('PASS')",
    ]
    python_result = subprocess.run(python_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) if Path(python_bin).is_file() else None
    add(
        "python_packages",
        python_result is not None and python_result.returncode == 0,
        python_result.stdout.strip() if python_result else "python missing",
        python_imports + " importable",
    )
    if taxonomic_scope:
        r_bin = config.get("R_BIN", "Rscript")
        r_code = "pkgs<-c('phyloseq','data.table'); quit(status=ifelse(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)),0,1))"
        r_result = subprocess.run([r_bin, "-e", r_code], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) if Path(r_bin).is_file() else None
        add(
            "R_packages",
            r_result is not None and r_result.returncode == 0,
            "PASS" if r_result is not None and r_result.returncode == 0 else (r_result.stdout[-1000:] if r_result else "R missing"),
            "phyloseq,data.table importable",
        )

    if environmental_scope:
        for key in ["HUMAN_FOOTPRINT_2020_TIF", "WORLDCLIM_BIO_DIR", "MEOW_VECTOR"]:
            value = config.get(key, "")
            path = Path(value) if value else None
            exists = path is not None and path.exists()
            severity = "warning" if args.allow_missing_environmental_assets else "error"
            add(f"environmental_asset_{key}", exists, path or "EMPTY", "existing path", severity)

    errors = [row for row in checks if not row["pass"] and row["severity"] == "error"]
    warnings = [row for row in checks if not row["pass"] and row["severity"] == "warning"]
    payload = {
        "status": "PASS" if not errors else "FAILED",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(config_path),
        "scope": args.scope,
        "n_checks": len(checks),
        "n_errors": len(errors),
        "n_warnings": len(warnings),
        "checks": checks,
    }
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Download the frozen MAGENTA run cohort and verify ENA checksums.

ENA FASTQ links are preferred.  ``fasterq-dump`` is used only when public
FASTQ links are absent.  Internal runs such as TX1_3 can be supplied through a
local mapping table with columns run_accession,r1,r2.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib.parse import urlparse

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download and validate MAGENTA FASTQ files.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--local-fastq-map", default=None)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--fasterq-threads", type=int, default=8)
    parser.add_argument("--run-accession", action="append", default=[])
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def split_field(value: Any) -> List[str]:
    if value is None or pd.isna(value):
        return []
    return [token.strip() for token in str(value).split(";") if token.strip()]


def normalize_url(value: str) -> str:
    if "://" in value:
        return value
    return "https://" + value.lstrip("/")


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_gzip(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"FASTQ ausente o vacío: {path}")
    with gzip.open(path, "rb") as stream:
        stream.read(1)
    result = subprocess.run(["gzip", "-t", str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(f"gzip -t falló para {path}: {result.stderr.decode(errors='replace')}")


def run_download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    aria2 = shutil.which("aria2c")
    curl = shutil.which("curl")
    if aria2:
        command = [
            aria2,
            "--continue=true",
            "--max-connection-per-server=8",
            "--split=8",
            "--min-split-size=16M",
            "--max-tries=10",
            "--retry-wait=5",
            "--file-allocation=none",
            "--dir", str(destination.parent),
            "--out", destination.name,
            url,
        ]
    elif curl:
        command = [curl, "--fail", "--location", "--retry", "10", "--continue-at", "-", "--output", str(destination), url]
    else:
        raise RuntimeError("Se requiere aria2c o curl")
    result = subprocess.run(command)
    if result.returncode != 0:
        raise RuntimeError(f"Descarga fallida ({result.returncode}): {url}")


def copy_local(source: Path, destination: Path, overwrite: bool) -> None:
    if not source.is_file():
        raise FileNotFoundError(f"FASTQ local no existe: {source}")
    if destination.exists() and not overwrite:
        validate_gzip(destination)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    shutil.copy2(source, temporary)
    validate_gzip(temporary)
    os.replace(temporary, destination)


def compress_fastq(source: Path, destination: Path) -> None:
    pigz = shutil.which("pigz")
    if pigz:
        with destination.open("wb") as stream:
            result = subprocess.run([pigz, "-c", str(source)], stdout=stream)
        if result.returncode != 0:
            raise RuntimeError(f"pigz falló para {source}")
    else:
        with source.open("rb") as incoming, gzip.open(destination, "wb", compresslevel=6) as outgoing:
            shutil.copyfileobj(incoming, outgoing, length=8 * 1024**2)


def fasterq_fallback(run: str, layout: str, out_dir: Path, threads: int, overwrite: bool) -> Tuple[Path, Optional[Path]]:
    binary = shutil.which("fasterq-dump")
    if not binary:
        raise RuntimeError(f"{run}: sin enlaces FASTQ y fasterq-dump no está en PATH")
    r1_final = out_dir / f"{run}_1.fastq.gz"
    r2_final = out_dir / f"{run}_2.fastq.gz" if layout == "PAIRED" else None
    if r1_final.exists() and (r2_final is None or r2_final.exists()) and not overwrite:
        validate_gzip(r1_final)
        if r2_final:
            validate_gzip(r2_final)
        return r1_final, r2_final
    with tempfile.TemporaryDirectory(prefix=f"{run}.", dir=out_dir) as temporary_name:
        temporary = Path(temporary_name)
        command = [binary, run, "--split-files", "--threads", str(threads), "--outdir", str(temporary), "--temp", str(temporary)]
        result = subprocess.run(command)
        if result.returncode != 0:
            raise RuntimeError(f"{run}: fasterq-dump falló ({result.returncode})")
        r1_source = temporary / f"{run}_1.fastq"
        single_source = temporary / f"{run}.fastq"
        if not r1_source.exists() and single_source.exists():
            r1_source = single_source
        if not r1_source.exists():
            raise RuntimeError(f"{run}: fasterq-dump no generó R1")
        r1_partial = r1_final.with_suffix(r1_final.suffix + ".partial")
        compress_fastq(r1_source, r1_partial)
        validate_gzip(r1_partial)
        os.replace(r1_partial, r1_final)
        if r2_final is not None:
            r2_source = temporary / f"{run}_2.fastq"
            if not r2_source.exists():
                raise RuntimeError(f"{run}: metadata PAIRED pero fasterq-dump no generó R2")
            r2_partial = r2_final.with_suffix(r2_final.suffix + ".partial")
            compress_fastq(r2_source, r2_partial)
            validate_gzip(r2_partial)
            os.replace(r2_partial, r2_final)
    return r1_final, r2_final


def download_run(row: pd.Series, args: argparse.Namespace, local_map: Dict[str, Tuple[Path, Optional[Path]]]) -> Dict[str, Any]:
    run = str(row[args.id_col]).strip()
    layout = str(row.get("library_layout", "")).strip().upper() or "UNKNOWN"
    out_dir = Path(args.out_dir)
    r1_final = out_dir / f"{run}_1.fastq.gz"
    r2_final = out_dir / f"{run}_2.fastq.gz" if layout == "PAIRED" else None
    source = ""
    expected_md5 = split_field(row.get("fastq_md5"))

    if run in local_map:
        local_r1, local_r2 = local_map[run]
        copy_local(local_r1, r1_final, args.overwrite)
        if local_r2 is not None:
            r2_final = out_dir / f"{run}_2.fastq.gz"
            copy_local(local_r2, r2_final, args.overwrite)
            layout = "PAIRED"
        else:
            r2_final = None
            layout = "SINGLE"
        source = "local_map"
    else:
        urls = [normalize_url(value) for value in split_field(row.get("fastq_ftp"))]
        if urls:
            expected_files = 2 if layout == "PAIRED" else 1
            if len(urls) != expected_files:
                raise RuntimeError(f"{run}: layout={layout}, pero ENA reporta {len(urls)} FASTQ")
            destinations: List[Path] = [r1_final]
            if len(urls) == 2:
                r2_final = out_dir / f"{run}_2.fastq.gz"
                destinations.append(r2_final)
            for index, (url, destination) in enumerate(zip(urls, destinations)):
                if destination.exists() and not args.overwrite:
                    validate_gzip(destination)
                else:
                    partial = destination.with_suffix(destination.suffix + ".partial")
                    run_download(url, partial)
                    validate_gzip(partial)
                    if index < len(expected_md5) and md5sum(partial).lower() != expected_md5[index].lower():
                        partial.unlink(missing_ok=True)
                        raise RuntimeError(f"{run}: MD5 incorrecto para {url}")
                    os.replace(partial, destination)
            source = "ENA_fastq"
        else:
            r1_final, r2_final = fasterq_fallback(run, layout, out_dir, args.fasterq_threads, args.overwrite)
            source = "fasterq_dump"

    validate_gzip(r1_final)
    if r2_final is not None:
        validate_gzip(r2_final)
    return {
        "run_accession": run,
        "layout": "PAIRED" if r2_final is not None else "SINGLE",
        "r1": str(r1_final),
        "r2": str(r2_final) if r2_final is not None else "",
        "r1_size_bytes": r1_final.stat().st_size,
        "r2_size_bytes": r2_final.stat().st_size if r2_final is not None else 0,
        "r1_md5": md5sum(r1_final),
        "r2_md5": md5sum(r2_final) if r2_final is not None else "",
        "download_source": source,
        "status": "PASS",
    }


def load_local_map(path: Optional[str]) -> Dict[str, Tuple[Path, Optional[Path]]]:
    if not path:
        return {}
    table = pd.read_csv(path, sep="\t" if str(path).endswith(".tsv") else ",", dtype="string")
    required = {"run_accession", "r1"}
    if not required.issubset(table.columns):
        raise ValueError(f"local-fastq-map requiere columnas {sorted(required)}")
    if table["run_accession"].duplicated().any():
        raise ValueError("local-fastq-map contiene run_accession duplicados")
    mapping: Dict[str, Tuple[Path, Optional[Path]]] = {}
    for _, row in table.iterrows():
        r2 = None if "r2" not in table.columns or pd.isna(row.get("r2")) or not str(row.get("r2")).strip() else Path(str(row["r2"]))
        mapping[str(row["run_accession"]).strip()] = (Path(str(row["r1"])), r2)
    return mapping


def main() -> int:
    args = parse_args()
    if args.max_workers < 1 or args.fasterq_threads < 1:
        raise ValueError("Los hilos/workers deben ser >= 1")
    metadata = pd.read_csv(args.metadata, low_memory=False)
    if args.id_col not in metadata.columns:
        raise ValueError(f"Falta {args.id_col} en metadata")
    metadata[args.id_col] = metadata[args.id_col].astype("string").str.strip()
    if metadata[args.id_col].isna().any() or metadata[args.id_col].duplicated().any():
        raise ValueError(f"{args.id_col} debe ser completo y único")
    selected = set(args.run_accession)
    if selected:
        missing = selected - set(metadata[args.id_col].astype(str))
        if missing:
            raise ValueError(f"run_accession solicitados ausentes: {sorted(missing)}")
        metadata = metadata.loc[metadata[args.id_col].isin(selected)].copy()
    if args.limit > 0:
        metadata = metadata.head(args.limit).copy()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    local_map = load_local_map(args.local_fastq_map)
    results: List[Dict[str, Any]] = []
    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {executor.submit(download_run, row, args, local_map): str(row[args.id_col]) for _, row in metadata.iterrows()}
        for future in as_completed(futures):
            run = futures[future]
            try:
                result = future.result()
                results.append(result)
                print(f"[PASS] {run}", flush=True)
            except Exception as error:
                errors.append(f"{run}: {type(error).__name__}: {error}")
                print(f"[FAIL] {errors[-1]}", flush=True)
    manifest = pd.DataFrame(results).sort_values("run_accession") if results else pd.DataFrame()
    manifest_path = Path(args.manifest_out)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, index=False)
    if errors:
        error_path = manifest_path.with_name(manifest_path.stem + ".errors.txt")
        error_path.write_text("\n".join(errors) + "\n", encoding="utf-8")
        raise RuntimeError(f"Fallaron {len(errors)} corridas; revise {error_path}")
    if len(manifest) != len(metadata):
        raise RuntimeError(f"Manifiesto incompleto: {len(manifest)} de {len(metadata)} corridas")
    print(f"[OK] FASTQ manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

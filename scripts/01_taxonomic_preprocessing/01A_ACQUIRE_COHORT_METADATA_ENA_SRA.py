#!/usr/bin/env python3
"""Freeze and enrich the curated MAGENTA cohort with ENA/SRA run metadata.

The manually curated table defines the cohort.  Public APIs are queried only
for those ``run_accession`` values, avoiding a time-dependent broad search that
would silently change the study population.  Curated values are preserved;
public metadata fill missing fields and conflicts are exported for review.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FIELDS = [
    "run_accession",
    "sample_accession",
    "study_accession",
    "scientific_name",
    "tax_id",
    "library_strategy",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "collection_date",
    "country",
    "location",
    "lat",
    "lon",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "read_count",
    "base_count",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Enrich a frozen run_accession cohort from ENA/SRA.")
    parser.add_argument("--curated-metadata", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--public-out", required=True)
    parser.add_argument("--conflicts-out", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    parser.add_argument("--max-attempts", type=int, default=5)
    parser.add_argument("--allow-missing-public", action="store_true")
    parser.add_argument("--allow-non-wgs-illumina", action="store_true")
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator, dtype="string", low_memory=False)


def clean_text(value: Any) -> Any:
    if value is None or pd.isna(value):
        return pd.NA
    text = str(value).strip()
    return pd.NA if text.lower() in {"", "na", "n/a", "nan", "none", "null"} else text


def collapse_curated(table: pd.DataFrame, id_col: str) -> pd.DataFrame:
    if id_col not in table.columns:
        aliases = [column for column in table.columns if column.lower() in {"run", "run_accession"}]
        if len(aliases) != 1:
            raise ValueError(f"No se encontró una llave inequívoca {id_col}")
        table = table.rename(columns={aliases[0]: id_col})
    work = table.copy()
    work[id_col] = work[id_col].map(clean_text).astype("string")
    if work[id_col].isna().any():
        raise ValueError(f"{id_col} contiene valores vacíos")
    rows: List[Dict[str, Any]] = []
    conflicts: List[str] = []
    for run, group in work.groupby(id_col, sort=False, dropna=False):
        record: Dict[str, Any] = {id_col: run}
        for column in work.columns:
            if column == id_col:
                continue
            values = [clean_text(value) for value in group[column]]
            unique = sorted({str(value) for value in values if not pd.isna(value)})
            if len(unique) > 1:
                conflicts.append(f"{run}:{column}={'|'.join(unique[:5])}")
            record[column] = unique[0] if unique else pd.NA
        rows.append(record)
    if conflicts:
        raise ValueError("Metadata curada duplicada con conflictos: " + "; ".join(conflicts[:20]))
    output = pd.DataFrame(rows, columns=work.columns)
    if output[id_col].duplicated().any():
        raise RuntimeError("La metadata colapsada conserva run_accession duplicados")
    return output


def batched(values: Sequence[str], size: int) -> Iterable[Sequence[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def request_ena(run_ids: Sequence[str], args: argparse.Namespace) -> pd.DataFrame:
    query = " OR ".join(f'run_accession=\"{run}\"' for run in run_ids)
    parameters = {
        "result": "read_run",
        "query": query,
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "limit": 0,
    }
    request = Request(
        ENA_URL + "?" + urlencode(parameters),
        headers={"User-Agent": "MAGENTA-taxonomic-workflow/1.0", "Accept": "text/tab-separated-values"},
    )
    last_error = "unknown"
    for attempt in range(1, args.max_attempts + 1):
        try:
            with urlopen(request, timeout=args.timeout_seconds) as response:
                payload = response.read().decode("utf-8", errors="replace")
            if not payload.strip():
                return pd.DataFrame(columns=ENA_FIELDS)
            return pd.read_csv(io.StringIO(payload), sep="\t", dtype="string")
        except Exception as error:
            last_error = f"{type(error).__name__}: {error}"
            if attempt < args.max_attempts:
                time.sleep(min(2 ** (attempt - 1), 30))
    raise RuntimeError(f"ENA falló después de {args.max_attempts} intentos: {last_error}")


def atomic_csv(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(handle)
    try:
        table.to_csv(temporary, index=False)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    if args.batch_size < 1 or args.batch_size > 200:
        raise ValueError("--batch-size debe estar entre 1 y 200")
    curated_path = Path(args.curated_metadata)
    output_path = Path(args.out)
    public_path = Path(args.public_out)
    conflicts_path = Path(args.conflicts_out)
    manifest_path = Path(args.manifest_out)
    if not curated_path.is_file():
        raise FileNotFoundError(f"No existe metadata curada: {curated_path}")

    curated = collapse_curated(read_table(curated_path), args.id_col)
    run_ids = curated[args.id_col].astype(str).tolist()
    public_parts: List[pd.DataFrame] = []
    for index, chunk in enumerate(batched(run_ids, args.batch_size), start=1):
        print(f"[ENA] lote {index}: {len(chunk)} accesiones", flush=True)
        public_parts.append(request_ena(chunk, args))
    public = pd.concat(public_parts, ignore_index=True) if public_parts else pd.DataFrame(columns=ENA_FIELDS)
    if not public.empty:
        public = public.drop_duplicates("run_accession", keep="first")
    public_ids = set(public.get("run_accession", pd.Series(dtype="string")).dropna().astype(str))
    missing_public = [run for run in run_ids if run not in public_ids]
    if missing_public and not args.allow_missing_public:
        raise RuntimeError(
            f"{len(missing_public)} accesiones sin registro público ENA/SRA: {missing_public[:20]}. "
            "Use --allow-missing-public únicamente para muestras internas documentadas."
        )

    if not public.empty and not args.allow_non_wgs_illumina:
        strategy = public["library_strategy"].fillna("").str.upper()
        platform = public["instrument_platform"].fillna("").str.upper()
        invalid = ~(strategy.eq("WGS") & platform.eq("ILLUMINA"))
        if invalid.any():
            examples = public.loc[invalid, ["run_accession", "library_strategy", "instrument_platform"]]
            raise RuntimeError("Corridas públicas no WGS/ILLUMINA:\n" + examples.head(20).to_string(index=False))

    public_lookup = public.set_index("run_accession") if not public.empty else pd.DataFrame()
    enriched = curated.copy()
    conflicts: List[Dict[str, Any]] = []
    for column in ENA_FIELDS:
        if column == "run_accession":
            continue
        incoming = enriched[args.id_col].map(public_lookup[column]) if not public.empty and column in public_lookup else pd.Series(pd.NA, index=enriched.index, dtype="string")
        if column in enriched.columns:
            existing = enriched[column].map(clean_text)
            incoming_clean = incoming.map(clean_text)
            conflict = existing.notna() & incoming_clean.notna() & existing.astype(str).ne(incoming_clean.astype(str))
            for row_index in enriched.index[conflict]:
                conflicts.append(
                    {
                        "run_accession": enriched.at[row_index, args.id_col],
                        "field": column,
                        "curated_value": existing.at[row_index],
                        "public_value": incoming_clean.at[row_index],
                    }
                )
            enriched[column] = existing.fillna(incoming_clean)
        else:
            enriched[column] = incoming
    enriched["metadata_public_status"] = enriched[args.id_col].map(
        lambda run: "ENA_SRA_matched" if str(run) in public_ids else "curated_internal_or_unmatched"
    )
    enriched["metadata_frozen_timestamp_utc"] = datetime.now(timezone.utc).isoformat()

    atomic_csv(public, public_path)
    atomic_csv(pd.DataFrame(conflicts, columns=["run_accession", "field", "curated_value", "public_value"]), conflicts_path)
    atomic_csv(enriched, output_path)
    manifest = {
        "schema_version": "1.0",
        "script": Path(__file__).name,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "curated_metadata": str(curated_path),
        "curated_sha256": sha256(curated_path),
        "output": str(output_path),
        "output_sha256": sha256(output_path),
        "n_curated_runs": len(curated),
        "n_public_matches": len(public_ids & set(run_ids)),
        "n_missing_public": len(missing_public),
        "missing_public_run_accession": missing_public,
        "n_curated_public_conflicts": len(conflicts),
        "primary_key": args.id_col,
        "cohort_policy": "curated_metadata_defines_frozen_cohort",
        "ena_fields": ENA_FIELDS,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[OK] metadata enriquecida: {output_path}")
    print(f"[OK] corridas: {len(enriched)}; públicas: {len(public_ids & set(run_ids))}; no públicas: {len(missing_public)}")
    print(f"[OK] conflictos auditados: {len(conflicts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

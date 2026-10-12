#!/usr/bin/env python3
"""Validate the final environmental metadata table and export missingness/QC."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd


GAM_COLUMNS = [
    "human_footprint_2020",
    "bio1", "bio2", "bio3", "bio13", "bio15",
    "gbif_shannon_H_2km", "gbif_richness_species_2km",
    "soc_mean", "ph_mean", "sand_mean", "silt_mean", "clay_mean", "bdod_mean",
]
MEOW_COLUMNS = ["REALM", "PROVINCE", "ECOREGION"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate MAGENTA environmental metadata.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--lat-col", default="lat")
    parser.add_argument("--lon-col", default="lon")
    parser.add_argument("--expected-rows", type=int, default=341)
    parser.add_argument("--allow-missing-meow", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    path = Path(args.metadata)
    if not path.is_file():
        raise FileNotFoundError(path)
    metadata = pd.read_csv(path, low_memory=False)
    if "ph_mean" not in metadata and "phh2o_mean" in metadata:
        metadata["ph_mean"] = metadata["phh2o_mean"]
    required = [args.id_col, args.lat_col, args.lon_col, *GAM_COLUMNS, *MEOW_COLUMNS]
    missing_columns = [column for column in required if column not in metadata.columns]
    if missing_columns:
        raise ValueError(f"Faltan columnas: {missing_columns}")
    ids = metadata[args.id_col].astype("string").str.strip()
    if ids.isna().any() or ids.eq("").any() or ids.duplicated().any():
        raise ValueError("run_accession debe ser completo y único")
    lat = pd.to_numeric(metadata[args.lat_col], errors="coerce")
    lon = pd.to_numeric(metadata[args.lon_col], errors="coerce")
    missing_coordinates = lat.isna() & lon.isna()
    invalid_coordinates = (lat.isna() ^ lon.isna()) | (lat.notna() & lon.notna() & (~lat.between(-90, 90) | ~lon.between(-180, 180)))
    if invalid_coordinates.any():
        raise ValueError(f"{int(invalid_coordinates.sum())} coordenadas WGS84 inválidas")
    if args.expected_rows > 0 and len(metadata) != args.expected_rows:
        raise ValueError(f"Filas observadas={len(metadata)}; esperadas={args.expected_rows}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    missingness: List[Dict[str, Any]] = []
    for column in [*GAM_COLUMNS, *MEOW_COLUMNS]:
        values = metadata[column]
        if column in MEOW_COLUMNS:
            missing = values.isna() | values.astype("string").str.strip().isin(["", "NA", "Unassigned"])
        else:
            numeric = pd.to_numeric(values, errors="coerce")
            missing = ~np.isfinite(numeric)
        missingness.append(
            {
                "variable": column,
                "n_total": len(metadata),
                "n_missing": int(missing.sum()),
                "n_complete": int((~missing).sum()),
                "complete_fraction": float((~missing).mean()),
            }
        )
    missingness_table = pd.DataFrame(missingness)
    missingness_table.to_csv(out_dir / "02G_environmental_variable_missingness.tsv", sep="\t", index=False)
    complete_gam = metadata[GAM_COLUMNS].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
    sample_audit = pd.DataFrame(
        {
            "run_accession": ids,
            "complete_14_predictor_GAM": complete_gam,
            "n_missing_GAM_predictors": metadata[GAM_COLUMNS].apply(pd.to_numeric, errors="coerce").isna().sum(axis=1),
            "MEOW_complete": metadata[MEOW_COLUMNS].notna().all(axis=1),
        }
    )
    sample_audit.to_csv(out_dir / "02G_environmental_sample_completeness.tsv", sep="\t", index=False)
    if not args.allow_missing_meow and not sample_audit["MEOW_complete"].all():
        raise RuntimeError("Hay muestras sin clasificación MEOW; revise la auditoría")
    summary = {
        "status": "PASS",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "metadata": str(path),
        "metadata_sha256": sha256(path),
        "n_rows": len(metadata),
        "n_unique_run_accession": ids.nunique(),
        "n_complete_14_predictor_GAM": int(complete_gam.sum()),
        "n_incomplete_14_predictor_GAM": int((~complete_gam).sum()),
        "n_realms": int(metadata["REALM"].nunique(dropna=True)),
        "primary_key": args.id_col,
        "coordinate_crs": "EPSG:4326",
        "n_missing_coordinate_pairs": int(missing_coordinates.sum()),
    }
    (out_dir / "02G_validation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

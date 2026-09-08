#!/usr/bin/env python3
"""Extract the configured 2020 Human Footprint raster at MAGENTA coordinates.

The raster is supplied explicitly because the manuscript must cite and archive
the exact dataset version.  Coastal NoData values are resolved through the
nearest valid raster cell within a declared distance, and the actual anchor
cell is exported for auditing.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import rasterio

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable=None, **_: Any):
        return iterable

from metadata_pipeline_common import (
    assert_row_contract,
    atomic_write_csv,
    derived_path,
    fill_from_diagnostics,
    nearest_valid_raster_cell,
    select_coordinate_work,
    validate_coordinates,
    validate_primary_key,
    write_manifest,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract the archived 2020 Human Footprint raster.")
    parser.add_argument("--base", required=True, help="Input metadata CSV.")
    parser.add_argument("--human-footprint-tif", required=True, help="Archived GeoTIFF used in the manuscript.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--diag-out", default=None)
    parser.add_argument("--manifest-out", default=None)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--lat-col", default="lat")
    parser.add_argument("--lon-col", default="lon")
    parser.add_argument("--value-col", default="human_footprint_2020")
    parser.add_argument("--dataset-label", default="Human Footprint 2020")
    parser.add_argument("--max-anchor-distance-m", type=float, default=2000.0)
    parser.add_argument("--allow-unresolved", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base_path = Path(args.base)
    raster_path = Path(args.human_footprint_tif)
    output_path = Path(args.out)
    if not base_path.is_file():
        raise FileNotFoundError(f"No existe base: {base_path}")
    if not raster_path.is_file():
        raise FileNotFoundError(f"No existe raster Human Footprint: {raster_path}")
    if output_path.resolve() == base_path.resolve():
        raise ValueError("--out debe ser distinto de --base")
    if args.max_anchor_distance_m < 0:
        raise ValueError("--max-anchor-distance-m debe ser >= 0")

    diag_path = Path(args.diag_out) if args.diag_out else derived_path(output_path, ".diagnostics.csv")
    manifest_path = Path(args.manifest_out) if args.manifest_out else derived_path(output_path, ".manifest.json")

    metadata = pd.read_csv(base_path, low_memory=False)
    validate_primary_key(metadata, args.id_col, "base")
    validate_coordinates(metadata, args.lat_col, args.lon_col, "base")
    if args.value_col not in metadata.columns:
        metadata[args.value_col] = np.nan
    coordinate_work = select_coordinate_work(
        metadata,
        [args.value_col],
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )

    diagnostics: List[Dict[str, Any]] = []
    with rasterio.open(raster_path) as raster:
        if raster.count != 1:
            raise ValueError(f"Se esperaba un raster de una banda; se observaron {raster.count}")
        if raster.crs is None:
            raise ValueError("El raster Human Footprint no tiene CRS")
        for _, point in tqdm(coordinate_work.iterrows(), total=len(coordinate_work), desc="Human Footprint"):
            lat = float(point[args.lat_col])
            lon = float(point[args.lon_col])
            result = nearest_valid_raster_cell(raster, lat, lon, args.max_anchor_distance_m)
            record: Dict[str, Any] = {
                "latlon_key": str(point["latlon_key"]),
                "sample_lat": lat,
                "sample_lon": lon,
                "human_footprint_dataset_label": args.dataset_label,
                "human_footprint_source_file": raster_path.name,
                "human_footprint_raster_crs": str(raster.crs),
            }
            if result is None:
                record.update(
                    {
                        args.value_col: np.nan,
                        "human_footprint_anchor_method": "not_found",
                        "human_footprint_anchor_lat": np.nan,
                        "human_footprint_anchor_lon": np.nan,
                        "human_footprint_anchor_distance_m": np.nan,
                        "human_footprint_raster_row": np.nan,
                        "human_footprint_raster_col": np.nan,
                        "human_footprint_error": "no_valid_cell_within_radius",
                    }
                )
            else:
                record.update(
                    {
                        args.value_col: result["value"],
                        "human_footprint_anchor_method": result["method"],
                        "human_footprint_anchor_lat": result["anchor_lat"],
                        "human_footprint_anchor_lon": result["anchor_lon"],
                        "human_footprint_anchor_distance_m": result["distance_m"],
                        "human_footprint_raster_row": result["row"],
                        "human_footprint_raster_col": result["col"],
                        "human_footprint_error": "",
                    }
                )
            diagnostics.append(record)

    diag = pd.DataFrame(diagnostics)
    atomic_write_csv(diag, diag_path)
    unresolved = diag[args.value_col].isna() if args.value_col in diag.columns else pd.Series(dtype=bool)
    if unresolved.any() and not args.allow_unresolved:
        raise RuntimeError(
            f"{int(unresolved.sum())} coordenadas sin Human Footprint. Revise {diag_path} "
            "o use --allow-unresolved de forma explícita."
        )

    before = metadata.copy()
    enriched = fill_from_diagnostics(
        metadata,
        diag,
        [args.value_col],
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )
    assert_row_contract(before, enriched, args.id_col)
    atomic_write_csv(enriched, output_path)
    write_manifest(
        manifest_path,
        Path(__file__).name,
        args,
        [base_path, raster_path],
        output_path,
        len(enriched),
        len(coordinate_work),
    )
    print(f"[OK] metadata: {output_path}")
    print(f"[OK] diagnóstico: {diag_path}")
    print(f"[OK] manifiesto: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

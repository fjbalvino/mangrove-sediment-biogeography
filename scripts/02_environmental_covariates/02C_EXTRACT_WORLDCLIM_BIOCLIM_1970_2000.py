#!/usr/bin/env python3
"""Extract WorldClim v2.1 BIO variables from one common raster cell.

WorldClim temperature BIO layers are distributed as values scaled by ten in
the standard GeoTIFF archive.  The default ``--temperature-scale 0.1`` exports
degrees Celsius for BIO1, BIO2, BIO5-BIO11 while preserving the conventional
units of BIO3, BIO4 and the precipitation variables.
"""

from __future__ import annotations

import argparse
import glob
from contextlib import ExitStack
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import rasterio

try:
    from tqdm import tqdm
except ImportError:  # progress display is optional
    def tqdm(iterable=None, **_: Any):
        return iterable

from metadata_pipeline_common import (
    assert_row_contract,
    atomic_write_csv,
    derived_path,
    ensure_same_raster_grid,
    fill_from_diagnostics,
    nearest_valid_raster_cell,
    raster_value_at_cell,
    select_coordinate_work,
    validate_coordinates,
    validate_primary_key,
    write_manifest,
)


def parse_bios(specification: str) -> List[int]:
    values: List[int] = []
    for token in specification.replace(" ", "").split(","):
        if not token:
            continue
        if "-" in token:
            start, end = token.split("-", 1)
            values.extend(range(int(start), int(end) + 1))
        else:
            values.append(int(token))
    result = sorted(set(values))
    if not result or any(value < 1 or value > 19 for value in result):
        raise ValueError("--bios debe contener números entre 1 y 19")
    return result


def detect_prefix(climate_dir: Path) -> str:
    matches = sorted(glob.glob(str(climate_dir / "*_bio_1.tif")))
    if matches:
        return Path(matches[0]).name[: -len("_1.tif")]
    matches = sorted(glob.glob(str(climate_dir / "*_1.tif")))
    if matches:
        return Path(matches[0]).name[: -len("_1.tif")]
    raise FileNotFoundError(f"No se encontró el raster BIO1 en {climate_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract WorldClim v2.1 BIO1-BIO19, normals 1970-2000."
    )
    parser.add_argument("--base", required=True)
    parser.add_argument("--clim-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--diag-out", default=None)
    parser.add_argument("--manifest-out", default=None)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--lat-col", default="lat")
    parser.add_argument("--lon-col", default="lon")
    parser.add_argument("--bios", default="1-19")
    parser.add_argument("--ref-bio", type=int, default=1)
    parser.add_argument("--prefix", default=None)
    parser.add_argument(
        "--temperature-scale",
        type=float,
        default=0.1,
        help="Multiplier for temperature BIO layers; standard WorldClim GeoTIFFs use 0.1.",
    )
    parser.add_argument("--max-anchor-distance-m", type=float, default=2000.0)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base_path = Path(args.base)
    climate_dir = Path(args.clim_dir)
    output_path = Path(args.out)
    if not base_path.is_file():
        raise FileNotFoundError(f"No existe base: {base_path}")
    if not climate_dir.is_dir():
        raise FileNotFoundError(f"No existe clim-dir: {climate_dir}")
    if output_path.resolve() == base_path.resolve():
        raise ValueError("--out debe ser distinto de --base")
    bios = parse_bios(args.bios)
    if not np.isfinite(args.temperature_scale) or args.temperature_scale <= 0:
        raise ValueError("--temperature-scale debe ser positivo")
    if args.ref_bio not in range(1, 20):
        raise ValueError("--ref-bio debe estar entre 1 y 19")
    prefix = args.prefix or detect_prefix(climate_dir)
    raster_paths = {
        bio: climate_dir / f"{prefix}_{bio}.tif"
        for bio in sorted(set(bios + [args.ref_bio]))
    }
    missing = [str(path) for path in raster_paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Faltan rasters WorldClim: " + ", ".join(missing))

    diag_path = Path(args.diag_out) if args.diag_out else derived_path(output_path, ".diagnostics.csv")
    manifest_path = (
        Path(args.manifest_out)
        if args.manifest_out
        else derived_path(output_path, ".manifest.json")
    )

    metadata = pd.read_csv(base_path, low_memory=False)
    validate_primary_key(metadata, args.id_col, "base")
    validate_coordinates(metadata, args.lat_col, args.lon_col, "base")
    bio_columns = [f"bio{bio}" for bio in bios]
    for column in bio_columns:
        if column not in metadata.columns:
            metadata[column] = np.nan
    coordinate_work = select_coordinate_work(
        metadata,
        bio_columns,
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )

    diagnostics: List[Dict[str, Any]] = []
    with ExitStack() as stack:
        sources = {
            bio: stack.enter_context(rasterio.open(path))
            for bio, path in raster_paths.items()
        }
        reference = sources[args.ref_bio]
        if reference.crs is None:
            raise ValueError("El raster WorldClim de referencia no tiene CRS")
        for bio, source in sources.items():
            ensure_same_raster_grid(reference, source, f"BIO{bio}")

        temperature_bios = {1, 2, 5, 6, 7, 8, 9, 10, 11}
        for _, point in tqdm(
            coordinate_work.iterrows(),
            total=len(coordinate_work),
            desc="WorldClim",
        ):
            lat = float(point[args.lat_col])
            lon = float(point[args.lon_col])
            anchor = nearest_valid_raster_cell(
                reference,
                lat,
                lon,
                max_distance_m=args.max_anchor_distance_m,
            )
            record: Dict[str, Any] = {
                "latlon_key": str(point["latlon_key"]),
                "sample_lat": lat,
                "sample_lon": lon,
                "worldclim_version": "2.1",
                "worldclim_period": "1970-2000",
                "worldclim_raster_prefix": prefix,
                "worldclim_x_resolution": abs(float(reference.transform.a)),
                "worldclim_y_resolution": abs(float(reference.transform.e)),
                "worldclim_reference_bio": args.ref_bio,
                "worldclim_temperature_scale": args.temperature_scale,
                "worldclim_crs": str(reference.crs),
            }
            if anchor is None:
                record.update(
                    {
                        "worldclim_anchor_method": "not_found",
                        "worldclim_anchor_lat": np.nan,
                        "worldclim_anchor_lon": np.nan,
                        "worldclim_anchor_distance_m": np.nan,
                        "worldclim_raster_row": np.nan,
                        "worldclim_raster_col": np.nan,
                        "worldclim_error": "no_valid_reference_cell_within_radius",
                    }
                )
                record.update({column: np.nan for column in bio_columns})
            else:
                row = int(anchor["row"])
                col = int(anchor["col"])
                missing_at_anchor: List[str] = []
                for bio in bios:
                    value = raster_value_at_cell(sources[bio], row, col)
                    if np.isfinite(value) and bio in temperature_bios:
                        value *= args.temperature_scale
                    record[f"bio{bio}"] = value
                    if not np.isfinite(value):
                        missing_at_anchor.append(f"bio{bio}")
                record.update(
                    {
                        "worldclim_anchor_method": anchor["method"],
                        "worldclim_anchor_lat": anchor["anchor_lat"],
                        "worldclim_anchor_lon": anchor["anchor_lon"],
                        "worldclim_anchor_distance_m": anchor["distance_m"],
                        "worldclim_raster_row": row,
                        "worldclim_raster_col": col,
                        "worldclim_error": (
                            ""
                            if not missing_at_anchor
                            else "missing_at_common_cell:" + ",".join(missing_at_anchor)
                        ),
                    }
                )
            diagnostics.append(record)

    diag = pd.DataFrame(diagnostics)
    before = metadata.copy()
    enriched = fill_from_diagnostics(
        metadata,
        diag,
        bio_columns,
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )
    assert_row_contract(before, enriched, args.id_col)
    atomic_write_csv(enriched, output_path)
    atomic_write_csv(diag, diag_path)
    write_manifest(
        manifest_path,
        Path(__file__).name,
        args,
        [base_path] + [raster_paths[bio] for bio in sorted(raster_paths)],
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

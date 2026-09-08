#!/usr/bin/env python3
"""Assign MEOW realm, province, and ecoregion from the official vector layer.

The official MEOW file is mandatory. Boundary points are handled with
``covers`` and points outside polygons are assigned only to the actual nearest
polygon within a declared maximum distance.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # progress display is optional
    def tqdm(iterable=None, **_: Any):
        return iterable

from metadata_pipeline_common import (
    assert_row_contract,
    atomic_write_csv,
    derived_path,
    fill_from_diagnostics,
    offset_point,
    select_coordinate_work,
    validate_coordinates,
    validate_primary_key,
    write_manifest,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assign official Marine Ecoregions of the World (MEOW) regions."
    )
    parser.add_argument("--base", required=True)
    parser.add_argument("--meow", required=True, help="Official MEOW vector file (.shp or .gpkg).")
    parser.add_argument("--out", required=True)
    parser.add_argument("--diag-out", default=None)
    parser.add_argument("--manifest-out", default=None)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--lat-col", default="lat")
    parser.add_argument("--lon-col", default="lon")
    parser.add_argument("--max-anchor-distance-m", type=float, default=2000.0)
    parser.add_argument("--allow-unassigned", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def find_column(columns: Sequence[str], candidates: Sequence[str], label: str) -> str:
    lookup = {str(column).lower(): str(column) for column in columns}
    for candidate in candidates:
        found = lookup.get(candidate.lower())
        if found is not None:
            return found
    raise ValueError(f"MEOW: no se encontró la columna {label}; columnas={list(columns)}")


def clean_classification(value: Any) -> Optional[str]:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "desconocida", "unassigned"}:
        return None
    return text


def classification(row: pd.Series, mapping: Dict[str, str]) -> Dict[str, str]:
    output: Dict[str, str] = {}
    for canonical, source in mapping.items():
        value = clean_classification(row[source])
        if value is None:
            raise ValueError(f"MEOW: polígono con {canonical} vacío")
        output[canonical] = value
    return output


def resolve_hits(
    candidates: Any,
    point: Any,
    mapping: Dict[str, str],
) -> Optional[Tuple[pd.Series, str]]:
    covers = candidates.loc[candidates.geometry.apply(lambda geometry: geometry.covers(point))]
    if covers.empty:
        return None
    strict = covers.loc[covers.geometry.apply(lambda geometry: geometry.contains(point))]
    selected = strict if not strict.empty else covers
    classes = {
        tuple(classification(row, mapping).items())
        for _, row in selected.iterrows()
    }
    if len(classes) > 1:
        raise ValueError("Punto MEOW ambiguo: varios polígonos con clasificaciones distintas")
    return selected.iloc[0], "contains" if not strict.empty else "covers_boundary"


def spatial_candidates(gdf: Any, point: Any, bounds: Tuple[float, float, float, float]) -> Any:
    try:
        indices = list(gdf.sindex.intersection(bounds))
        return gdf.iloc[indices]
    except Exception:
        min_lon, min_lat, max_lon, max_lat = bounds
        return gdf.loc[
            gdf.geometry.bounds.minx.le(max_lon)
            & gdf.geometry.bounds.maxx.ge(min_lon)
            & gdf.geometry.bounds.miny.le(max_lat)
            & gdf.geometry.bounds.maxy.ge(min_lat)
        ]


def local_utm_epsg(lat: float, lon: float) -> int:
    zone = max(1, min(60, int(math.floor((lon + 180.0) / 6.0)) + 1))
    return (32600 if lat >= 0 else 32700) + zone


def nearest_polygon(
    gdf: Any,
    point: Any,
    lat: float,
    lon: float,
    max_distance_m: float,
    mapping: Dict[str, str],
) -> Optional[Dict[str, Any]]:
    import geopandas as gpd
    from pyproj import Transformer
    from shapely.geometry import Point
    from shapely.ops import nearest_points, transform

    north = offset_point(lat, lon, max_distance_m, 0.0)
    east = offset_point(lat, lon, max_distance_m, 90.0)
    south = offset_point(lat, lon, max_distance_m, 180.0)
    west = offset_point(lat, lon, max_distance_m, 270.0)
    bounds = (west[1], south[0], east[1], north[0])
    candidates = spatial_candidates(gdf, point, bounds)
    if candidates.empty:
        return None

    epsg = local_utm_epsg(lat, lon)
    projected = candidates.to_crs(epsg=epsg)
    point_series = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326").to_crs(epsg=epsg)
    point_projected = point_series.iloc[0]
    distances = projected.geometry.distance(point_projected)
    within = distances[distances <= max_distance_m + 1e-6]
    if within.empty:
        return None
    minimum = float(within.min())
    tied_indices = within.index[np.isclose(within.to_numpy(), minimum, atol=1e-6)]
    tied_rows = candidates.loc[tied_indices]
    classes = {
        tuple(classification(row, mapping).items())
        for _, row in tied_rows.iterrows()
    }
    if len(classes) > 1:
        raise ValueError("Ancla MEOW ambigua: polígonos equidistantes con clasificaciones distintas")
    selected_index = tied_indices[0]
    selected_original = candidates.loc[selected_index]
    selected_projected = projected.loc[selected_index].geometry
    _, anchor_projected = nearest_points(point_projected, selected_projected)
    transformer = Transformer.from_crs(epsg, 4326, always_xy=True)
    anchor_lon, anchor_lat = transformer.transform(anchor_projected.x, anchor_projected.y)
    return {
        **classification(selected_original, mapping),
        "method": "nearest_polygon",
        "anchor_lat": float(anchor_lat),
        "anchor_lon": float(anchor_lon),
        "distance_m": minimum,
        "utm_epsg": epsg,
    }


def lookup_meow(
    gdf: Any,
    lat: float,
    lon: float,
    max_distance_m: float,
    mapping: Dict[str, str],
) -> Optional[Dict[str, Any]]:
    from shapely.geometry import Point

    point = Point(lon, lat)
    direct_candidates = spatial_candidates(gdf, point, point.bounds)
    direct = resolve_hits(direct_candidates, point, mapping)
    if direct is not None:
        row, method = direct
        return {
            **classification(row, mapping),
            "method": method,
            "anchor_lat": lat,
            "anchor_lon": lon,
            "distance_m": 0.0,
            "utm_epsg": np.nan,
        }
    return nearest_polygon(gdf, point, lat, lon, max_distance_m, mapping)


def main() -> int:
    args = parse_args()
    try:
        import geopandas as gpd
    except ImportError as error:
        raise RuntimeError("Se requiere geopandas para la clasificación MEOW oficial") from error

    base_path = Path(args.base)
    meow_path = Path(args.meow)
    output_path = Path(args.out)
    if not base_path.is_file():
        raise FileNotFoundError(f"No existe base: {base_path}")
    if not meow_path.is_file():
        raise FileNotFoundError(f"No existe MEOW: {meow_path}")
    if output_path.resolve() == base_path.resolve():
        raise ValueError("--out debe ser distinto de --base")
    diag_path = Path(args.diag_out) if args.diag_out else derived_path(output_path, ".diagnostics.csv")
    manifest_path = (
        Path(args.manifest_out)
        if args.manifest_out
        else derived_path(output_path, ".manifest.json")
    )

    metadata = pd.read_csv(base_path, low_memory=False)
    validate_primary_key(metadata, args.id_col, "base")
    validate_coordinates(metadata, args.lat_col, args.lon_col, "base")
    region_columns = ["REALM", "PROVINCE", "ECOREGION"]
    for column in region_columns:
        if column not in metadata.columns:
            metadata[column] = pd.Series(pd.NA, index=metadata.index, dtype="string")
        else:
            metadata[column] = metadata[column].map(clean_classification).astype("string")

    gdf = gpd.read_file(meow_path)
    if gdf.empty:
        raise ValueError("El archivo MEOW no contiene polígonos")
    if gdf.crs is None:
        raise ValueError("El archivo MEOW no tiene CRS; no se asumirá EPSG:4326")
    gdf = gdf.to_crs(4326)
    if hasattr(gdf.geometry, "make_valid"):
        gdf.geometry = gdf.geometry.make_valid()
    else:
        gdf.geometry = gdf.geometry.buffer(0)
    invalid = gdf.geometry.isna() | gdf.geometry.is_empty | ~gdf.geometry.is_valid
    if invalid.any():
        raise ValueError(f"MEOW contiene {int(invalid.sum())} geometrías inválidas")
    mapping = {
        "REALM": find_column(gdf.columns, ["REALM", "realm"], "REALM"),
        "PROVINCE": find_column(gdf.columns, ["PROVINCE", "province"], "PROVINCE"),
        "ECOREGION": find_column(gdf.columns, ["ECOREGION", "ecoregion"], "ECOREGION"),
    }

    coordinate_work = select_coordinate_work(
        metadata,
        region_columns,
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )
    diagnostics: List[Dict[str, Any]] = []
    for _, point in tqdm(
        coordinate_work.iterrows(),
        total=len(coordinate_work),
        desc="MEOW",
    ):
        lat = float(point[args.lat_col])
        lon = float(point[args.lon_col])
        result = lookup_meow(
            gdf,
            lat,
            lon,
            args.max_anchor_distance_m,
            mapping,
        )
        record: Dict[str, Any] = {
            "latlon_key": str(point["latlon_key"]),
            "sample_lat": lat,
            "sample_lon": lon,
            "meow_source_file": meow_path.name,
        }
        if result is None:
            record.update(
                {
                    "REALM": pd.NA,
                    "PROVINCE": pd.NA,
                    "ECOREGION": pd.NA,
                    "meow_anchor_method": "not_found",
                    "meow_anchor_lat": np.nan,
                    "meow_anchor_lon": np.nan,
                    "meow_anchor_distance_m": np.nan,
                    "meow_local_utm_epsg": np.nan,
                    "meow_error": "no_polygon_within_radius",
                }
            )
        else:
            record.update(
                {
                    "REALM": result["REALM"],
                    "PROVINCE": result["PROVINCE"],
                    "ECOREGION": result["ECOREGION"],
                    "meow_anchor_method": result["method"],
                    "meow_anchor_lat": result["anchor_lat"],
                    "meow_anchor_lon": result["anchor_lon"],
                    "meow_anchor_distance_m": result["distance_m"],
                    "meow_local_utm_epsg": result["utm_epsg"],
                    "meow_error": "",
                }
            )
        diagnostics.append(record)

    diag = pd.DataFrame(diagnostics)
    unresolved = diag.loc[diag.get("meow_error", pd.Series(dtype=str)).eq("no_polygon_within_radius")]
    atomic_write_csv(diag, diag_path)
    if not unresolved.empty and not args.allow_unassigned:
        raise RuntimeError(
            f"{len(unresolved)} coordenadas no pudieron asignarse a MEOW dentro de "
            f"{args.max_anchor_distance_m:g} m. Revise {diag_path} o use --allow-unassigned."
        )

    before = metadata.copy()
    enriched = fill_from_diagnostics(
        metadata,
        diag,
        region_columns,
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )
    assert_row_contract(before, enriched, args.id_col)
    atomic_write_csv(enriched, output_path)
    meow_inputs = sorted(
        path for path in meow_path.parent.glob(meow_path.stem + ".*") if path.is_file()
    )
    write_manifest(
        manifest_path,
        Path(__file__).name,
        args,
        [base_path] + meow_inputs,
        output_path,
        len(enriched),
        len(coordinate_work),
    )
    print(f"[OK] metadata final: {output_path}")
    print(f"[OK] diagnóstico: {diag_path}")
    print(f"[OK] manifiesto: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Shared, dependency-light utilities for the MAGENTA metadata pipeline.

All tables are keyed by ``run_accession`` and coordinate-level lookups are
cached/joined through a deterministic WGS84 coordinate key.  The functions in
this module never silently delete or reorder metadata rows.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


MISSING_ID_TOKENS = {"", "na", "n/a", "nan", "none", "null"}


def require_columns(table: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{label}: faltan columnas requeridas: {missing}")


def clean_ids(values: Iterable[Any]) -> pd.Series:
    series = pd.Series(values, copy=False, dtype="string").str.strip()
    return series.mask(series.str.lower().isin(MISSING_ID_TOKENS))


def validate_primary_key(table: pd.DataFrame, id_col: str, label: str) -> None:
    require_columns(table, [id_col], label)
    ids = clean_ids(table[id_col])
    missing = ids.isna()
    if missing.any():
        rows = table.index[missing].tolist()[:20]
        raise ValueError(f"{label}: {id_col} vacío en filas {rows}")
    duplicated = ids.duplicated(keep=False)
    if duplicated.any():
        examples = sorted(ids[duplicated].astype(str).unique().tolist())[:20]
        raise ValueError(f"{label}: {id_col} duplicado: {examples}")
    table[id_col] = ids.astype(str)


def validate_coordinates(
    table: pd.DataFrame,
    lat_col: str,
    lon_col: str,
    label: str,
    require_complete: bool = False,
) -> None:
    require_columns(table, [lat_col, lon_col], label)
    lat = pd.to_numeric(table[lat_col], errors="coerce")
    lon = pd.to_numeric(table[lon_col], errors="coerce")
    partial = lat.isna() ^ lon.isna()
    if partial.any():
        raise ValueError(f"{label}: {int(partial.sum())} filas tienen solo una coordenada")
    if require_complete and (lat.isna() | lon.isna()).any():
        raise ValueError(f"{label}: se requieren coordenadas completas")
    invalid = lat.notna() & lon.notna() & (~lat.between(-90, 90) | ~lon.between(-180, 180))
    if invalid.any():
        raise ValueError(f"{label}: {int(invalid.sum())} coordenadas fuera de rango WGS84")
    table[lat_col] = lat
    table[lon_col] = lon


def coordinate_key(lat: float, lon: float, decimals: int = 6) -> str:
    return f"{float(lat):.{decimals}f}|{float(lon):.{decimals}f}"


def select_coordinate_work(
    metadata: pd.DataFrame,
    target_columns: Sequence[str],
    overwrite: bool,
    lat_col: str,
    lon_col: str,
) -> pd.DataFrame:
    validate_coordinates(metadata, lat_col, lon_col, "metadata")
    valid = metadata[lat_col].notna() & metadata[lon_col].notna()
    if overwrite:
        needed = valid
    else:
        for column in target_columns:
            if column not in metadata.columns:
                metadata[column] = np.nan
        needed = valid & metadata[list(target_columns)].isna().any(axis=1)
    work = metadata.loc[needed, [lat_col, lon_col]].copy()
    work["latlon_key"] = [
        coordinate_key(lat, lon) for lat, lon in zip(work[lat_col], work[lon_col])
    ]
    duplicate_coordinates = work.groupby("latlon_key")[[lat_col, lon_col]].nunique()
    if ((duplicate_coordinates > 1).any(axis=1)).any():
        raise ValueError("Una clave de coordenada representa valores incompatibles")
    return work.drop_duplicates("latlon_key", keep="first").reset_index(drop=True)


def fill_from_diagnostics(
    metadata: pd.DataFrame,
    diagnostics: pd.DataFrame,
    value_columns: Sequence[str],
    overwrite: bool,
    lat_col: str,
    lon_col: str,
) -> pd.DataFrame:
    if diagnostics.empty:
        return metadata.copy()
    require_columns(diagnostics, ["latlon_key", *value_columns], "diagnostics")
    if diagnostics["latlon_key"].duplicated().any():
        raise ValueError("diagnostics: latlon_key duplicado")
    output = metadata.copy()
    keys = pd.Series(pd.NA, index=output.index, dtype="string")
    valid = output[lat_col].notna() & output[lon_col].notna()
    keys.loc[valid] = [
        coordinate_key(lat, lon)
        for lat, lon in zip(output.loc[valid, lat_col], output.loc[valid, lon_col])
    ]
    lookup = diagnostics.set_index("latlon_key")
    for column in value_columns:
        if column not in output.columns:
            output[column] = np.nan
        incoming = keys.map(lookup[column])
        if overwrite:
            mask = incoming.notna()
        else:
            mask = output[column].isna() & incoming.notna()
        output.loc[mask, column] = incoming.loc[mask].to_numpy()
    return output


def assert_row_contract(before: pd.DataFrame, after: pd.DataFrame, id_col: str) -> None:
    if len(before) != len(after):
        raise ValueError(f"Contrato de filas roto: {len(before)} -> {len(after)}")
    require_columns(before, [id_col], "before")
    require_columns(after, [id_col], "after")
    if clean_ids(before[id_col]).tolist() != clean_ids(after[id_col]).tolist():
        raise ValueError("El orden o conjunto de run_accession cambió durante el enriquecimiento")


def derived_path(path: Path, suffix: str) -> Path:
    return path.with_name(path.stem + suffix)


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if pd.isna(value):
        return None
    raise TypeError(f"No serializable: {type(value).__name__}")


def atomic_write_json(payload: Mapping[str, Any], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True, default=_json_default)
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_write_csv(table: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(handle)
    try:
        table.to_csv(temporary, index=False)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def sha256_file(path: Path, maximum_bytes: int = 2 * 1024**3) -> Optional[str]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size > maximum_bytes:
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> Dict[str, Any]:
    path = Path(path)
    record: Dict[str, Any] = {"path": str(path), "exists": path.is_file()}
    if path.is_file():
        stat = path.stat()
        record.update(
            {
                "size_bytes": stat.st_size,
                "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                "sha256": sha256_file(path),
            }
        )
    return record


def write_manifest(
    path: Path,
    script_name: str,
    args: Any,
    inputs: Sequence[Path],
    output: Path,
    n_rows: int,
    n_coordinates_processed: int,
) -> None:
    payload = {
        "schema_version": "1.0",
        "script": script_name,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "arguments": vars(args),
        "inputs": [file_record(Path(item)) for item in inputs],
        "output": file_record(Path(output)),
        "n_rows": int(n_rows),
        "n_coordinates_processed": int(n_coordinates_processed),
    }
    atomic_write_json(payload, Path(path))


def cache_file(cache_dir: Path, namespace: str, key: str, parameters: Mapping[str, Any]) -> Path:
    canonical = json.dumps(parameters, sort_keys=True, separators=(",", ":"), default=_json_default)
    token = hashlib.sha256(f"{namespace}|{key}|{canonical}".encode("utf-8")).hexdigest()
    return Path(cache_dir) / namespace / f"{token}.json"


def read_json_cache(path: Path) -> Optional[Dict[str, Any]]:
    path = Path(path)
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as stream:
            payload = json.load(stream)
        return payload if isinstance(payload, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def offset_point(lat: float, lon: float, distance_m: float, bearing_deg: float) -> Tuple[float, float]:
    radius_m = 6_371_008.8
    angular = float(distance_m) / radius_m
    bearing = math.radians(float(bearing_deg))
    lat1 = math.radians(float(lat))
    lon1 = math.radians(float(lon))
    lat2 = math.asin(
        math.sin(lat1) * math.cos(angular)
        + math.cos(lat1) * math.sin(angular) * math.cos(bearing)
    )
    lon2 = lon1 + math.atan2(
        math.sin(bearing) * math.sin(angular) * math.cos(lat1),
        math.cos(angular) - math.sin(lat1) * math.sin(lat2),
    )
    return math.degrees(lat2), ((math.degrees(lon2) + 180.0) % 360.0) - 180.0


def haversine_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    radius_m = 6_371_008.8
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    value = math.sin(dp / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    return 2.0 * radius_m * math.asin(min(1.0, math.sqrt(value)))


def raster_value_at_cell(source: Any, row: int, col: int) -> float:
    if row < 0 or col < 0 or row >= source.height or col >= source.width:
        return np.nan
    value = source.read(1, window=((row, row + 1), (col, col + 1)), masked=True)
    if value.size != 1 or np.ma.is_masked(value[0, 0]):
        return np.nan
    numeric = float(value[0, 0])
    return numeric if np.isfinite(numeric) else np.nan


def ensure_same_raster_grid(reference: Any, candidate: Any, label: str) -> None:
    if reference.crs != candidate.crs:
        raise ValueError(f"{label}: CRS distinto al raster de referencia")
    if reference.width != candidate.width or reference.height != candidate.height:
        raise ValueError(f"{label}: dimensiones distintas al raster de referencia")
    if not np.allclose(tuple(reference.transform), tuple(candidate.transform), rtol=0, atol=1e-12):
        raise ValueError(f"{label}: malla/transformación distinta al raster de referencia")


def _pixel_size_m(source: Any, lat: float) -> float:
    xres = abs(float(source.transform.a))
    yres = abs(float(source.transform.e))
    if getattr(source.crs, "is_geographic", False):
        x_m = xres * 111_320.0 * max(0.05, math.cos(math.radians(lat)))
        y_m = yres * 110_574.0
        return max(1.0, min(x_m, y_m))
    return max(1.0, min(xres, yres))


def nearest_valid_raster_cell(
    source: Any,
    lat: float,
    lon: float,
    max_distance_m: float,
) -> Optional[Dict[str, Any]]:
    from rasterio.warp import transform as transform_coordinates
    from rasterio.windows import Window

    if source.crs is None:
        raise ValueError("El raster no tiene CRS")
    x, y = transform_coordinates("EPSG:4326", source.crs, [float(lon)], [float(lat)])
    row0, col0 = source.index(x[0], y[0])
    direct = raster_value_at_cell(source, row0, col0)
    if np.isfinite(direct):
        anchor_x, anchor_y = source.xy(row0, col0)
        anchor_lon, anchor_lat = transform_coordinates(source.crs, "EPSG:4326", [anchor_x], [anchor_y])
        distance = haversine_m(lat, lon, anchor_lat[0], anchor_lon[0])
        if distance <= max_distance_m + 1e-6:
            return {
                "value": direct,
                "row": int(row0),
                "col": int(col0),
                "anchor_lat": float(anchor_lat[0]),
                "anchor_lon": float(anchor_lon[0]),
                "distance_m": float(distance),
                "method": "sample_cell",
            }

    radius = int(math.ceil(max_distance_m / _pixel_size_m(source, lat))) + 2
    radius = min(max(radius, 1), 4096)
    row_start = max(0, row0 - radius)
    col_start = max(0, col0 - radius)
    row_stop = min(source.height, row0 + radius + 1)
    col_stop = min(source.width, col0 + radius + 1)
    window = Window(col_start, row_start, col_stop - col_start, row_stop - row_start)
    data = source.read(1, window=window, masked=True)
    best: Optional[Dict[str, Any]] = None
    for local_row, local_col in np.argwhere(~np.ma.getmaskarray(data)):
        value = float(data[local_row, local_col])
        if not np.isfinite(value):
            continue
        row = row_start + int(local_row)
        col = col_start + int(local_col)
        anchor_x, anchor_y = source.xy(row, col)
        anchor_lon, anchor_lat = transform_coordinates(source.crs, "EPSG:4326", [anchor_x], [anchor_y])
        distance = haversine_m(lat, lon, anchor_lat[0], anchor_lon[0])
        if distance > max_distance_m + 1e-6:
            continue
        if best is None or distance < best["distance_m"]:
            best = {
                "value": value,
                "row": row,
                "col": col,
                "anchor_lat": float(anchor_lat[0]),
                "anchor_lon": float(anchor_lon[0]),
                "distance_m": float(distance),
                "method": "nearest_valid_cell",
            }
    return best

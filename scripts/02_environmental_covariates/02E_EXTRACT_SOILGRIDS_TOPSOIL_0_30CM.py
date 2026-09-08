#!/usr/bin/env python3
"""Extract SoilGrids v2.0 topsoil properties for 0-30 cm.

Only the 0-5, 5-15, and 15-30 cm layers are accepted. Values are converted
from SoilGrids integer encodings and combined with thickness weights 5:10:15.
No shallow-depth fallback is permitted.
"""

from __future__ import annotations

import argparse
import json
import math
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

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
    atomic_write_json,
    cache_file,
    coordinate_key,
    derived_path,
    fill_from_diagnostics,
    offset_point,
    read_json_cache,
    select_coordinate_work,
    validate_coordinates,
    validate_primary_key,
    write_manifest,
)


SOILGRIDS_URL = "https://rest.isric.org/soilgrids/v2.0/properties/query"
USER_AGENT = "MAGENTA-metadata-pipeline/2.0"
DEPTHS = ("0-5cm", "5-15cm", "15-30cm")
THICKNESS_CM = {"0-5cm": 5.0, "5-15cm": 10.0, "15-30cm": 15.0}
DEFAULT_PROPERTIES = ("clay", "soc", "sand", "silt", "bdod", "cec", "cfvo", "phh2o")

# Scale factors convert official SoilGrids integer encodings to analysis units.
PROPERTY_METADATA: Dict[str, Dict[str, Any]] = {
    "clay": {"factor": 0.1, "unit": "%", "range": (0.0, 100.0)},
    "sand": {"factor": 0.1, "unit": "%", "range": (0.0, 100.0)},
    "silt": {"factor": 0.1, "unit": "%", "range": (0.0, 100.0)},
    "soc": {"factor": 0.1, "unit": "g/kg", "range": (0.0, 1000.0)},
    "bdod": {"factor": 0.01, "unit": "kg/dm3", "range": (0.0, 3.0)},
    "cec": {"factor": 0.1, "unit": "cmol(c)/kg", "range": (0.0, 200.0)},
    "cfvo": {"factor": 0.1, "unit": "cm3/dm3", "range": (0.0, 1000.0)},
    "phh2o": {"factor": 0.1, "unit": "pH", "range": (0.0, 14.0)},
}


class RateLimiter:
    def __init__(self, rate_per_second: float):
        self.rate = float(rate_per_second)
        self.capacity = max(1.0, self.rate * 2.0)
        self.tokens = self.capacity
        self.updated = time.monotonic()
        self.lock = threading.Lock()

    def acquire(self) -> None:
        while True:
            with self.lock:
                now = time.monotonic()
                self.tokens = min(
                    self.capacity,
                    self.tokens + (now - self.updated) * self.rate,
                )
                self.updated = now
                if self.tokens >= 1.0:
                    self.tokens -= 1.0
                    return
            time.sleep(0.02)


RATE_LIMITER: RateLimiter


def request_json(parameters: Mapping[str, Any], args: argparse.Namespace) -> Dict[str, Any]:
    backoff = 1.0
    last_error = "unknown"
    for attempt in range(1, args.max_attempts + 1):
        RATE_LIMITER.acquire()
        request = Request(
            SOILGRIDS_URL + "?" + urlencode(parameters, doseq=True),
            headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        )
        try:
            with urlopen(request, timeout=args.timeout_seconds) as response:
                payload = json.loads(response.read().decode("utf-8", errors="replace"))
            if not isinstance(payload, dict):
                raise RuntimeError("SoilGrids devolvió una respuesta no estructurada")
            return payload
        except HTTPError as error:
            last_error = f"HTTP {error.code}"
            if error.code not in {429, 500, 502, 503, 504}:
                break
            retry_after = error.headers.get("Retry-After") if error.headers else None
            wait = float(retry_after) if retry_after and retry_after.isdigit() else backoff
        except (URLError, TimeoutError, json.JSONDecodeError) as error:
            last_error = f"{type(error).__name__}: {error}"
            wait = backoff
        except Exception as error:
            last_error = f"{type(error).__name__}: {error}"
            wait = backoff
        if attempt < args.max_attempts:
            time.sleep(min(wait, args.backoff_cap_seconds))
            backoff = min(backoff * 2.0, args.backoff_cap_seconds)
    raise RuntimeError(
        f"Consulta SoilGrids fallida después de {args.max_attempts} intentos: {last_error}"
    )


def depth_label(depth: Mapping[str, Any]) -> Optional[str]:
    label = depth.get("label")
    if label:
        return str(label).replace(" ", "")
    depth_range = depth.get("range") or {}
    top = depth_range.get("top_depth")
    bottom = depth_range.get("bottom_depth")
    unit = depth_range.get("unit_depth", "cm")
    if top is None or bottom is None:
        return None
    return f"{top}-{bottom}{unit}".replace(".0", "")


def parse_layer_values(payload: Mapping[str, Any], prop: str) -> Dict[str, float]:
    layers = ((payload.get("properties") or {}).get("layers") or [])
    matching = [layer for layer in layers if str(layer.get("name", prop)) == prop]
    if not matching and len(layers) == 1:
        matching = layers
    if len(matching) != 1:
        raise ValueError(f"Respuesta SoilGrids ambigua para {prop}: {len(matching)} capas")
    output: Dict[str, float] = {}
    for depth in matching[0].get("depths") or []:
        label = depth_label(depth)
        value = (depth.get("values") or {}).get("mean")
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if label in DEPTHS and np.isfinite(numeric):
            output[label] = numeric
    return output


def converted_topsoil_values(raw: Mapping[str, float], prop: str) -> Dict[str, Any]:
    missing = [depth for depth in DEPTHS if depth not in raw]
    if missing:
        return {"value": np.nan, "depth_values": {}, "error": "missing_depths:" + ",".join(missing)}
    metadata = PROPERTY_METADATA[prop]
    converted = {depth: float(raw[depth]) * metadata["factor"] for depth in DEPTHS}
    lower, upper = metadata["range"]
    outside = [depth for depth, value in converted.items() if value < lower or value > upper]
    if outside:
        return {
            "value": np.nan,
            "depth_values": converted,
            "error": "converted_value_out_of_range:" + ",".join(outside),
        }
    weighted = sum(converted[depth] * THICKNESS_CM[depth] for depth in DEPTHS) / 30.0
    return {"value": float(weighted), "depth_values": converted, "error": ""}


def query_property(
    lat: float,
    lon: float,
    prop: str,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    key = coordinate_key(lat, lon)
    parameters_for_cache = {
        "endpoint": "soilgrids_v2_properties",
        "property": prop,
        "depths": DEPTHS,
        "value": "mean",
        "schema": 2,
    }
    path = cache_file(Path(args.cache_dir), "soilgrids_topsoil", key, parameters_for_cache)
    cached = None if args.refresh_cache else read_json_cache(path)
    if cached is not None:
        cached["cache_status"] = "hit"
        return cached
    parameters = {
        "lat": f"{lat:.6f}",
        "lon": f"{lon:.6f}",
        "property": prop,
        "depth": list(DEPTHS),
        "value": "mean",
    }
    payload = request_json(parameters, args)
    raw = parse_layer_values(payload, prop)
    result = converted_topsoil_values(raw, prop)
    result.update(
        {
            "property": prop,
            "unit": PROPERTY_METADATA[prop]["unit"],
            "raw_depth_values": raw,
            "query_timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "cache_status": "miss",
        }
    )
    atomic_write_json(result, path)
    return result


def candidate_points(
    lat: float,
    lon: float,
    max_distance_m: float,
    radial_steps_m: Sequence[float],
    bearings: int,
) -> List[Tuple[float, float, float, float]]:
    output = [(lat, lon, 0.0, np.nan)]
    for distance in sorted(set(float(value) for value in radial_steps_m)):
        if distance <= 0 or distance > max_distance_m:
            continue
        for bearing in np.linspace(0.0, 360.0, bearings, endpoint=False):
            candidate_lat, candidate_lon = offset_point(lat, lon, distance, float(bearing))
            output.append((candidate_lat, candidate_lon, distance, float(bearing)))
    return output


def find_anchor(
    lat: float,
    lon: float,
    anchor_prop: str,
    args: argparse.Namespace,
) -> Tuple[Optional[Dict[str, Any]], List[str]]:
    attempts: List[str] = []
    steps = [float(value) for value in args.radial_steps_m.split(",") if value.strip()]
    for candidate_lat, candidate_lon, distance, bearing in candidate_points(
        lat,
        lon,
        args.max_anchor_distance_m,
        steps,
        args.search_bearings,
    ):
        result = query_property(candidate_lat, candidate_lon, anchor_prop, args)
        if result["error"]:
            attempts.append(f"{distance:g}m@{bearing:g}:{result['error']}")
            continue
        return (
            {
                "lat": candidate_lat,
                "lon": candidate_lon,
                "distance_m": distance,
                "bearing_deg": bearing,
                "method": "query_point" if distance == 0 else "nearest_valid_query_point",
                "anchor_result": result,
            },
            attempts,
        )
    return None, attempts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="SoilGrids v2.0 properties, thickness-weighted 0-30 cm."
    )
    parser.add_argument("--base", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--diag-out", default=None)
    parser.add_argument("--manifest-out", default=None)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--lat-col", default="lat")
    parser.add_argument("--lon-col", default="lon")
    parser.add_argument("--properties", default=",".join(DEFAULT_PROPERTIES))
    parser.add_argument("--anchor-property", default="clay")
    parser.add_argument("--max-anchor-distance-m", type=float, default=2000.0)
    parser.add_argument("--radial-steps-m", default="250,500,1000,2000")
    parser.add_argument("--search-bearings", type=int, default=8)
    parser.add_argument("--global-rps", type=float, default=1.0)
    parser.add_argument("--timeout-seconds", type=int, default=45)
    parser.add_argument("--max-attempts", type=int, default=6)
    parser.add_argument("--backoff-cap-seconds", type=float, default=60.0)
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    properties = [value.strip() for value in args.properties.split(",") if value.strip()]
    unknown = sorted(set(properties + [args.anchor_property]) - set(PROPERTY_METADATA))
    if unknown:
        raise ValueError(f"Propiedades SoilGrids sin conversión definida: {unknown}")
    if args.search_bearings < 1:
        raise ValueError("--search-bearings debe ser >= 1")
    global RATE_LIMITER
    RATE_LIMITER = RateLimiter(args.global_rps)

    base_path = Path(args.base)
    output_path = Path(args.out)
    if not base_path.is_file():
        raise FileNotFoundError(f"No existe base: {base_path}")
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
    property_columns = [f"{prop}_mean" for prop in properties]
    if "phh2o" in properties:
        property_columns.append("ph_mean")
    for column in property_columns:
        if column not in metadata.columns:
            metadata[column] = np.nan
    coordinate_work = select_coordinate_work(
        metadata,
        [f"{prop}_mean" for prop in properties],
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )

    diagnostics: List[Dict[str, Any]] = []
    for _, point in tqdm(
        coordinate_work.iterrows(),
        total=len(coordinate_work),
        desc="SoilGrids 0-30 cm",
    ):
        lat = float(point[args.lat_col])
        lon = float(point[args.lon_col])
        key = str(point["latlon_key"])
        record: Dict[str, Any] = {
            "latlon_key": key,
            "sample_lat": lat,
            "sample_lon": lon,
            "soilgrids_version": "2.0",
            "soilgrids_depth_basis": "0-30cm_thickness_weighted_5_10_15cm",
            "soilgrids_anchor_property": args.anchor_property,
        }
        anchor, attempts = find_anchor(lat, lon, args.anchor_property, args)
        if anchor is None:
            record.update(
                {
                    "soilgrids_anchor_method": "not_found",
                    "soilgrids_anchor_lat": np.nan,
                    "soilgrids_anchor_lon": np.nan,
                    "soilgrids_anchor_distance_m": np.nan,
                    "soilgrids_anchor_bearing_deg": np.nan,
                    "soilgrids_error": "no_complete_topsoil_anchor_within_radius",
                    "soilgrids_anchor_attempts": "|".join(attempts),
                }
            )
            record.update({f"{prop}_mean": np.nan for prop in properties})
        else:
            record.update(
                {
                    "soilgrids_anchor_method": anchor["method"],
                    "soilgrids_anchor_lat": anchor["lat"],
                    "soilgrids_anchor_lon": anchor["lon"],
                    "soilgrids_anchor_distance_m": anchor["distance_m"],
                    "soilgrids_anchor_bearing_deg": anchor["bearing_deg"],
                    "soilgrids_anchor_attempts": "|".join(attempts),
                }
            )
            errors: List[str] = []
            for prop in properties:
                result = (
                    anchor["anchor_result"]
                    if prop == args.anchor_property
                    else query_property(anchor["lat"], anchor["lon"], prop, args)
                )
                record[f"{prop}_mean"] = result["value"]
                record[f"{prop}_unit"] = PROPERTY_METADATA[prop]["unit"]
                for depth, value in result.get("depth_values", {}).items():
                    safe_depth = depth.replace("-", "_").replace("cm", "cm")
                    record[f"{prop}_{safe_depth}"] = value
                if result["error"]:
                    errors.append(f"{prop}:{result['error']}")
            record["soilgrids_error"] = "|".join(errors)
        if "phh2o_mean" in record:
            record["ph_mean"] = record["phh2o_mean"]
        diagnostics.append(record)

    diag = pd.DataFrame(diagnostics)
    before = metadata.copy()
    enriched = fill_from_diagnostics(
        metadata,
        diag,
        property_columns,
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
        [base_path],
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

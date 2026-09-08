#!/usr/bin/env python3
"""Extract reproducible GBIF plant-diversity metrics in fixed 2-km buffers.

Occurrence records are restricted to accepted ``speciesKey`` values, duplicate
occurrences are removed, records are spatially thinned within species, and
diversity is rarefied to a common record count with a deterministic seed.
Raw query summaries are cached by coordinate and parameter set.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple
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
    derived_path,
    fill_from_diagnostics,
    read_json_cache,
    select_coordinate_work,
    validate_coordinates,
    validate_primary_key,
    write_manifest,
)


GBIF_OCCURRENCE_URL = "https://api.gbif.org/v1/occurrence/search"
GBIF_RADIUS_M = 2000
GBIF_KINGDOM_KEY = 6
PAGE_LIMIT = 300
USER_AGENT = "MAGENTA-metadata-pipeline/2.0"


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


def request_json(
    parameters: Mapping[str, Any],
    timeout_seconds: int,
    max_attempts: int,
    backoff_cap_seconds: float,
) -> Dict[str, Any]:
    backoff = 1.0
    last_error = "unknown"
    for attempt in range(1, max_attempts + 1):
        RATE_LIMITER.acquire()
        request = Request(
            GBIF_OCCURRENCE_URL + "?" + urlencode(parameters, doseq=True),
            headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        )
        try:
            with urlopen(request, timeout=timeout_seconds) as response:
                payload = json.loads(response.read().decode("utf-8", errors="replace"))
            if not isinstance(payload, dict):
                raise RuntimeError("GBIF devolvió una respuesta no estructurada")
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
        if attempt < max_attempts:
            time.sleep(min(wait, backoff_cap_seconds))
            backoff = min(backoff * 2.0, backoff_cap_seconds)
    raise RuntimeError(f"Consulta GBIF fallida después de {max_attempts} intentos: {last_error}")


def query_parameters(lat: float, lon: float, args: argparse.Namespace) -> Dict[str, Any]:
    parameters: Dict[str, Any] = {
        "hasCoordinate": "true",
        "hasGeospatialIssue": "false",
        "occurrenceStatus": "PRESENT",
        "kingdomKey": GBIF_KINGDOM_KEY,
        "geoDistance": f"{lat:.6f},{lon:.6f},{GBIF_RADIUS_M}m",
    }
    if args.year_from is not None:
        year_to = args.year_to if args.year_to is not None else datetime.now().year
        parameters["year"] = f"{args.year_from},{year_to}"
    return parameters


def spatial_cell(
    record_lat: float,
    record_lon: float,
    center_lat: float,
    center_lon: float,
    grid_m: float,
) -> Tuple[int, int]:
    north_m = (record_lat - center_lat) * 111_320.0
    east_m = (
        (record_lon - center_lon)
        * 111_320.0
        * max(0.01, math.cos(math.radians(center_lat)))
    )
    return math.floor(east_m / grid_m), math.floor(north_m / grid_m)


def fetch_coordinate(
    key: str,
    lat: float,
    lon: float,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    cache_parameters = {
        "api": "gbif_occurrence_search_v1",
        "radius_m": GBIF_RADIUS_M,
        "kingdom_key": GBIF_KINGDOM_KEY,
        "year_from": args.year_from,
        "year_to": args.year_to,
        "thinning_grid_m": args.thinning_grid_m,
        "max_records": args.max_records,
        "schema": 2,
    }
    path = cache_file(Path(args.cache_dir), "gbif_2km", key, cache_parameters)
    cached = None if args.refresh_cache else read_json_cache(path)
    if cached is not None:
        cached["gbif_cache_status"] = "hit"
        return cached

    common = query_parameters(lat, lon, args)
    first_parameters = dict(common, limit=0)
    first = request_json(
        first_parameters,
        args.timeout_seconds,
        args.max_attempts,
        args.backoff_cap_seconds,
    )
    total_reported = int(first.get("count", 0) or 0)
    if total_reported > args.max_records and not args.allow_truncated:
        raise RuntimeError(
            f"{key}: GBIF reporta {total_reported} registros; excede --max-records "
            f"{args.max_records}. Use un GBIF download o --allow-truncated."
        )
    fetch_target = min(total_reported, args.max_records)
    offset = 0
    seen_occurrences: set[str] = set()
    seen_species_cells: set[Tuple[int, int, int]] = set()
    species_counter: Counter[int] = Counter()
    fetched_records = 0
    records_with_species = 0

    while fetched_records < fetch_target:
        parameters = dict(common)
        parameters["limit"] = min(PAGE_LIMIT, fetch_target - fetched_records)
        parameters["offset"] = offset
        page = request_json(
            parameters,
            args.timeout_seconds,
            args.max_attempts,
            args.backoff_cap_seconds,
        )
        results = page.get("results") or []
        if not results:
            break
        for record in results:
            occurrence_id = record.get("key") or record.get("gbifID")
            if occurrence_id is None:
                occurrence_id = hashlib.sha256(
                    json.dumps(record, sort_keys=True, default=str).encode("utf-8")
                ).hexdigest()
            occurrence_id = str(occurrence_id)
            if occurrence_id in seen_occurrences:
                continue
            seen_occurrences.add(occurrence_id)
            species_key = record.get("speciesKey")
            if not isinstance(species_key, int):
                continue
            records_with_species += 1
            record_lat = record.get("decimalLatitude")
            record_lon = record.get("decimalLongitude")
            try:
                record_lat = float(record_lat)
                record_lon = float(record_lon)
            except (TypeError, ValueError):
                continue
            cell_x, cell_y = spatial_cell(
                record_lat,
                record_lon,
                lat,
                lon,
                args.thinning_grid_m,
            )
            species_cell = (species_key, cell_x, cell_y)
            if species_cell in seen_species_cells:
                continue
            seen_species_cells.add(species_cell)
            species_counter[species_key] += 1
        page_size = len(results)
        fetched_records += page_size
        offset += page_size
        if page_size < parameters["limit"]:
            break

    truncated = total_reported > fetched_records
    payload: Dict[str, Any] = {
        "latlon_key": key,
        "sample_lat": lat,
        "sample_lon": lon,
        "gbif_query_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "gbif_radius_m": GBIF_RADIUS_M,
        "gbif_kingdom_key": GBIF_KINGDOM_KEY,
        "gbif_occurrences_reported_2km": total_reported,
        "gbif_occurrences_fetched_2km": fetched_records,
        "gbif_occurrences_unique_2km": len(seen_occurrences),
        "gbif_occurrences_with_species_2km": records_with_species,
        "gbif_occurrences_thinned_2km": int(sum(species_counter.values())),
        "gbif_species_before_rarefaction_2km": len(species_counter),
        "gbif_truncated": bool(truncated),
        "gbif_thinning_grid_m": float(args.thinning_grid_m),
        "species_counts": {str(key_): int(value) for key_, value in sorted(species_counter.items())},
        "gbif_cache_status": "miss",
    }
    atomic_write_json(payload, path)
    return payload


def diversity_metrics(counts: np.ndarray) -> Dict[str, float]:
    counts = np.asarray(counts, dtype=int)
    counts = counts[counts > 0]
    richness = int(counts.size)
    total = int(counts.sum())
    if richness == 0 or total == 0:
        return {
            "richness": 0.0,
            "shannon": np.nan,
            "pielou": np.nan,
            "simpson": np.nan,
            "chao1": np.nan,
        }
    proportions = counts / total
    shannon = float(-np.sum(proportions * np.log(proportions)))
    singletons = int(np.sum(counts == 1))
    doubletons = int(np.sum(counts == 2))
    chao1 = (
        richness + (singletons * singletons) / (2.0 * doubletons)
        if doubletons > 0
        else richness + singletons * (singletons - 1) / 2.0
    )
    return {
        "richness": float(richness),
        "shannon": shannon,
        "pielou": shannon / math.log(richness) if richness > 1 else np.nan,
        "simpson": float(1.0 - np.sum(proportions * proportions)),
        "chao1": float(chao1),
    }


def rarefy_counts(counts: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    counts = np.asarray(counts, dtype=int)
    if n < 0 or n > int(counts.sum()):
        raise ValueError("Profundidad de rarefacción inválida")
    remaining_draws = int(n)
    remaining_total = int(counts.sum())
    output = np.zeros_like(counts)
    for index in range(max(0, len(counts) - 1)):
        good = int(counts[index])
        bad = remaining_total - good
        draw = int(rng.hypergeometric(good, bad, remaining_draws))
        output[index] = draw
        remaining_draws -= draw
        remaining_total -= good
    if len(counts):
        output[-1] = remaining_draws
    return output


def coordinate_seed(base_seed: int, key: str) -> int:
    digest = hashlib.sha256(f"{base_seed}|{key}".encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def rarefied_metrics(
    counts: np.ndarray,
    n: int,
    repetitions: int,
    seed: int,
) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    metrics = [
        diversity_metrics(rarefy_counts(counts, n, rng))
        for _ in range(repetitions)
    ]
    output: Dict[str, float] = {}
    for name in ["richness", "shannon", "pielou", "simpson", "chao1"]:
        values = np.array([row[name] for row in metrics], dtype=float)
        output[name] = float(np.nanmean(values)) if np.isfinite(values).any() else np.nan
    return output


def parse_rarefaction_depth(specification: str, raw_results: List[Dict[str, Any]]) -> int:
    if specification.lower() == "min":
        positive = [
            int(result["gbif_occurrences_thinned_2km"])
            for result in raw_results
            if int(result["gbif_occurrences_thinned_2km"]) > 0
        ]
        if not positive:
            raise ValueError("No hay registros GBIF válidos después del filtrado")
        return min(positive)
    value = int(specification)
    if value < 1:
        raise ValueError("--rarefaction-n debe ser 'min' o un entero positivo")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GBIF Plantae diversity in fixed 2-km buffers with spatial thinning and rarefaction."
    )
    parser.add_argument("--base", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--diag-out", default=None)
    parser.add_argument("--manifest-out", default=None)
    parser.add_argument("--id-col", default="run_accession")
    parser.add_argument("--lat-col", default="lat")
    parser.add_argument("--lon-col", default="lon")
    parser.add_argument("--year-from", type=int, default=None)
    parser.add_argument("--year-to", type=int, default=None)
    parser.add_argument("--thinning-grid-m", type=float, default=100.0)
    parser.add_argument(
        "--rarefaction-n",
        default="min",
        help="Common record count or 'min'. 'min' requires --overwrite.",
    )
    parser.add_argument("--rarefaction-reps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--max-records", type=int, default=100000)
    parser.add_argument("--allow-truncated", action="store_true")
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--global-rps", type=float, default=2.0)
    parser.add_argument("--timeout-seconds", type=int, default=45)
    parser.add_argument("--max-attempts", type=int, default=6)
    parser.add_argument("--backoff-cap-seconds", type=float, default=60.0)
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.year_from is not None and args.year_to is not None and args.year_from > args.year_to:
        raise ValueError("--year-from no puede ser mayor que --year-to")
    if args.thinning_grid_m <= 0:
        raise ValueError("--thinning-grid-m debe ser positivo")
    if args.rarefaction_reps < 1:
        raise ValueError("--rarefaction-reps debe ser >= 1")
    if args.rarefaction_n.lower() == "min" and not args.overwrite:
        raise ValueError(
            "--rarefaction-n min requiere --overwrite para calcular una profundidad común "
            "sobre todas las coordenadas. En una actualización incremental indique el mismo "
            "valor entero usado en la corrida original."
        )

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
    metric_columns = [
        "gbif_shannon_H_2km",
        "gbif_richness_species_2km",
        "gbif_pielou_J_2km",
        "gbif_simpson_1minD_2km",
        "gbif_chao1_2km",
        "gbif_occurrences_effective",
        "gbif_occurrences_raw_2km",
        "gbif_occurrences_thinned_2km",
        "gbif_rarefaction_n_2km",
    ]
    for column in metric_columns:
        if column not in metadata.columns:
            metadata[column] = np.nan
    coordinate_work = select_coordinate_work(
        metadata,
        ["gbif_shannon_H_2km", "gbif_richness_species_2km"],
        overwrite=args.overwrite,
        lat_col=args.lat_col,
        lon_col=args.lon_col,
    )

    raw_by_key: Dict[str, Dict[str, Any]] = {}
    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {}
        for _, point in coordinate_work.iterrows():
            key = str(point["latlon_key"])
            lat = float(point[args.lat_col])
            lon = float(point[args.lon_col])
            futures[executor.submit(fetch_coordinate, key, lat, lon, args)] = key
        for future in tqdm(as_completed(futures), total=len(futures), desc="GBIF 2 km"):
            key = futures[future]
            try:
                raw_by_key[key] = future.result()
            except Exception as error:
                errors.append(f"{key}: {error}")
    if errors:
        raise RuntimeError("Fallaron consultas GBIF:\n" + "\n".join(errors[:20]))

    raw_results = [raw_by_key[str(key)] for key in coordinate_work["latlon_key"]]
    rarefaction_n = parse_rarefaction_depth(args.rarefaction_n, raw_results)
    diagnostics: List[Dict[str, Any]] = []
    for raw in raw_results:
        key = str(raw["latlon_key"])
        counts = np.array(
            [int(value) for _, value in sorted(raw["species_counts"].items())],
            dtype=int,
        )
        available = int(counts.sum())
        record = {name: value for name, value in raw.items() if name != "species_counts"}
        record["gbif_rarefaction_n_2km"] = rarefaction_n
        record["gbif_rarefaction_reps"] = args.rarefaction_reps
        record["gbif_rarefaction_seed"] = coordinate_seed(args.seed, key)
        record["gbif_occurrences_raw_2km"] = raw["gbif_occurrences_reported_2km"]
        record["gbif_occurrences_effective"] = rarefaction_n if available >= rarefaction_n else np.nan
        if available == 0:
            metrics = diversity_metrics(counts)
            record["gbif_error"] = "no_species_records_after_filtering"
        elif available < rarefaction_n:
            metrics = {name: np.nan for name in ["richness", "shannon", "pielou", "simpson", "chao1"]}
            record["gbif_error"] = "insufficient_records_for_common_rarefaction"
        else:
            metrics = rarefied_metrics(
                counts,
                rarefaction_n,
                args.rarefaction_reps,
                record["gbif_rarefaction_seed"],
            )
            record["gbif_error"] = ""
        record.update(
            {
                "gbif_richness_species_2km": metrics["richness"],
                "gbif_shannon_H_2km": metrics["shannon"],
                "gbif_pielou_J_2km": metrics["pielou"],
                "gbif_simpson_1minD_2km": metrics["simpson"],
                "gbif_chao1_2km": metrics["chao1"],
            }
        )
        diagnostics.append(record)

    diag = pd.DataFrame(diagnostics)
    before = metadata.copy()
    enriched = fill_from_diagnostics(
        metadata,
        diag,
        metric_columns,
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
    print(f"[OK] rarefacción común: {rarefaction_n} registros")
    print(f"[OK] diagnóstico: {diag_path}")
    print(f"[OK] manifiesto: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

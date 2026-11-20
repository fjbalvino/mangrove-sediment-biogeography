#!/usr/bin/env python3
"""Fill only missing marine salinity pairs using progressive N/E/S/W queries.

Reuses the geodesic offsets of SoilGrids (02E / metadata_pipeline_common).
Existing salinity values, all other fields, row order and IDs stay unchanged.
This is a regional raster extraction, not a coordinate correction or a new
hydrological validation. Prior warnings remain in the external audit.
"""
from __future__ import annotations

import argparse
import collections
import csv
import importlib.util
import json
import math
import platform
import re
import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from metadata_pipeline_common import offset_point

SPEC = importlib.util.spec_from_file_location("biooracle_source", Path(__file__).with_name("02H_EXTRACT_BIOORACLE_SALINITY.py"))
bio = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bio)
DIRECTIONS = (("N", 0), ("E", 90), ("S", 180), ("W", 270))
DEFAULT_STEPS = "250,500,1000,2000,5000,10000,20000,30000,50000"


def missing_pair(row):
    missing = [row[c].strip().lower() in bio.MISSING for c in bio.COLUMNS]
    if any(missing) and not all(missing):
        raise ValueError("Partial salinity pair: cannot silently mix cells or periods")
    if not any(missing) and any(not math.isfinite(float(row[c])) or float(row[c]) < 0 for c in bio.COLUMNS):
        raise ValueError("Invalid existing salinity value")
    return all(missing)


def point_url(index):
    iy, ix = index
    bounds = (iy, iy, ix, ix)
    selection = f"[0:1][{iy}:1:{iy}][{ix}:1:{ix}]"
    return bio.BASE_URL + "/griddap/" + bio.DATASET + ".csv?" + ",".join(v + selection for v in bio.VARIABLES), bounds


class CellSource:
    """Read validated old patches, then request only previously unknown cells."""
    def __init__(self, cache_dir, offline=False, workers=2):
        self.cache_dir, self.offline, self.workers = Path(cache_dir), offline, workers
        self.cells, self.records = {}, {}
        prefix = bio.BASE_URL + "/griddap/" + bio.DATASET + ".csv?"
        pattern = re.compile(r"^so_mean\[0:1\]\[(\d+):1:(\d+)\]\[(\d+):1:(\d+)\],so_range\[0:1\]\[\1:1:\2\]\[\3:1:\4\]$")
        for path in sorted(self.cache_dir.glob("*.json")):
            record = json.loads(path.read_text())
            url = record.get("url", "")
            if not url.startswith(prefix):
                continue
            match = pattern.fullmatch(urllib.parse.unquote(url[len(prefix):]))
            if not match:
                continue
            data, record = bio.fetch(url, self.cache_dir, offline=True)
            self.add(bio.parse_patch(data, tuple(map(int, match.groups()))), record)

    def add(self, cells, record):
        for index, periods in cells.items():
            if index in self.cells and self.cells[index] != periods:
                raise ValueError("Cached responses disagree for the same source cell")
            self.cells[index] = periods
        self.records[record["url"]] = record

    def ensure(self, indices):
        needed = sorted(set(indices) - self.cells.keys())
        def retrieve(index):
            url, bounds = point_url(index)
            data, record = bio.fetch(url, self.cache_dir, self.offline)
            return bio.parse_patch(data, bounds), record
        with ThreadPoolExecutor(max_workers=self.workers) as executor:
            for cells, record in executor.map(retrieve, needed):
                self.add(cells, record)


def cardinal_search(points, source, steps, max_distance_m):
    """Origin, then increasing radii; first complete N/E/S/W cell wins.

    Repeated cells are logged but not re-queried. Both the query-point radius
    and actual source-cell-center distance must meet the stated bound.
    """
    remaining = set(points)
    selected, attempts = {}, []
    seen = {point: set() for point in remaining}
    for radius in [0.0] + sorted(set(steps)):
        if radius > max_distance_m or not remaining:
            continue
        queries = {}
        for point in sorted(remaining):
            positions = [("origin", None, point)] if radius == 0 else [(name, bearing, offset_point(*point, radius, bearing)) for name, bearing in DIRECTIONS]
            queries[point] = [(name, bearing, query, bio.cell_index(*query)) for name, bearing, query in positions]
        source.ensure(index for qs in queries.values() for _, _, _, index in qs)
        for point, positions in queries.items():
            for direction, bearing, query, index in positions:
                center = bio.cell_center(*index)
                distance = bio.distance_km(point, center)
                periods = source.cells[index]
                values = bio.aggregate(periods)
                result = "repeated_cell" if index in seen[point] else "no_complete_decadal_pair"
                if index not in seen[point] and values is not None:
                    result = "outside_cell_center_distance_bound" if distance * 1000 > max_distance_m + 1e-6 else "selected"
                attempts.append({"lat": point[0], "lon": point[1], "radius_m": radius, "direction": direction, "bearing_deg": bearing if bearing is not None else "", "query_lat": query[0], "query_lon": query[1], "cell_lat": center[0], "cell_lon": center[1], "cell_distance_km": round(distance, 6), "result": result})
                seen[point].add(index)
                if result == "selected":
                    selected[point] = {**attempts[-1], "values": values, "periods": periods}
                    remaining.remove(point)
                    break
        print(f"CARDINAL radius_m={radius:g} recovered_coordinates={len(selected)} unresolved_coordinates={len(remaining)}", flush=True)
    return selected, attempts


def run(args):
    fields, rows = bio.read_metadata(args.base)
    if not set(bio.COLUMNS).issubset(fields):
        raise ValueError("Run 02H first: both salinity columns are required")
    if Path(args.base).resolve() == Path(args.out).resolve():
        raise ValueError("Input and output must be different paths")
    _, prior = bio.read_metadata(args.previous_audit)
    if [r["run_accession"] for r in rows] != [r["run_accession"] for r in prior]:
        raise ValueError("Previous audit and metadata accession order differs")
    for row, diag in zip(rows, prior):
        if any(row[c] != diag[c] for c in ("lat", "lon", *bio.COLUMNS)):
            raise ValueError("Stale previous audit: coordinates or salinity differs")
    missing = [missing_pair(r) for r in rows]
    steps = sorted(set(float(v) for v in args.radial_steps_m.split(",")))
    if not steps or any(not math.isfinite(v) or v <= 0 for v in steps):
        raise ValueError("Radial steps must be finite positive distances")
    if not math.isfinite(args.max_anchor_distance_m) or args.max_anchor_distance_m <= 0:
        raise ValueError("Maximum anchor distance must be finite and positive")
    if max(steps) < args.max_anchor_distance_m:
        steps.append(args.max_anchor_distance_m)
    steps = [v for v in steps if v <= args.max_anchor_distance_m]
    source_data, source_record = bio.fetch(bio.BASE_URL + "/info/" + bio.DATASET + "/index.json", args.cache_dir, args.offline)
    bio.validate_source(source_data)
    source = CellSource(args.cache_dir, args.offline, args.workers)
    points = {bio.coordinates(r) for r, needed in zip(rows, missing) if needed and bio.coordinates(r) is not None}
    chosen, attempts = cardinal_search(points, source, steps, args.max_anchor_distance_m)
    output, diagnostics, changes = [], [], []
    for row, old, needed in zip(rows, prior, missing):
        result, diag = dict(row), dict(old)
        diag.update(previous_status=old["status"], previous_reason=old.get("reason", ""), recovery_query_radius_m="", recovery_direction="", recovery_query_lat="", recovery_query_lon="", salinity_changed="false", recovery_note="")
        if needed:
            point = bio.coordinates(row)
            if point is None:
                diag["status"] = "missing_or_invalid_coordinates"
                diag["reason"] = "Cardinal search requires a valid coordinate pair; row retained"
            elif point not in chosen:
                diag["status"] = "no_cardinal_coverage_within_limit"
                diag["reason"] = f"No complete native/cardinal cell within {args.max_anchor_distance_m:g} m; no coordinate correction performed"
            else:
                selection = chosen[point]
                for i, column in enumerate(bio.COLUMNS):
                    result[column] = diag[column] = format(selection["values"][i], ".10g")
                    for j, decade in enumerate(("2000_2009", "2010_2019")):
                        diag[bio.VARIABLES[i] + "_" + decade] = selection["periods"][bio.TIMES[j]][i]
                diag.update(status="assigned_recovered_native" if selection["direction"] == "origin" else "assigned_cardinal", reason="Regional value recovered by origin then increasing-radius N/E/S/W search; prior warnings retained", cell_lat=selection["cell_lat"], cell_lon=selection["cell_lon"], distance_km=selection["cell_distance_km"], recovery_query_radius_m=selection["radius_m"], recovery_direction=selection["direction"], recovery_query_lat=selection["query_lat"], recovery_query_lon=selection["query_lon"], salinity_changed="true", recovery_note="Numerical regional extraction; not approval of coordinates or hydrological equivalence. See previous_status and previous_reason.")
                changes.append({"run_accession": row["run_accession"]} | {"old_" + c: row[c] for c in bio.COLUMNS} | {"new_" + c: result[c] for c in bio.COLUMNS})
        output.append(result)
        diagnostics.append(diag)
    for before, after, was_missing in zip(rows, output, missing):
        if any(before[c] != after[c] for c in fields if not (was_missing and c in bio.COLUMNS)):
            raise AssertionError("A previously available value or original field changed")
    audit_dir = Path(args.audit_dir)
    audit_dir.mkdir(parents=True, exist_ok=True)
    bio.write_csv(args.out, fields, output)
    bio.write_csv(audit_dir / "sample_salinity_audit.csv", list(diagnostics[0]), diagnostics)
    bio.write_csv(audit_dir / "cardinal_search_attempts.csv", ["lat", "lon", "radius_m", "direction", "bearing_deg", "query_lat", "query_lon", "cell_lat", "cell_lon", "cell_distance_km", "result"], attempts)
    bio.write_csv(audit_dir / "recovered_accessions.csv", ["run_accession"] + ["old_" + c for c in bio.COLUMNS] + ["new_" + c for c in bio.COLUMNS], changes)
    bio.write_csv(audit_dir / "registros_sin_salinidad.csv", list(diagnostics[0]), [d for d in diagnostics if not d["status"].startswith("assigned_")])
    manifest = {"status": "PASS", "dataset": bio.DATASET, "reference_period": "2000-2019", "aggregation": "same complete two-decadal means as 02H", "mode": "fill_missing_only_cardinal", "geometry_source": "metadata_pipeline_common.offset_point (also used by 02E SoilGrids)", "direction_order": [name for name, _ in DIRECTIONS], "radial_steps_m": steps, "max_query_and_cell_center_distance_m": args.max_anchor_distance_m, "coordinate_warnings": "preserved, not used to block requested regional extraction", "n_rows": len(rows), "n_previously_available": len(rows)-sum(missing), "n_recovered": len(changes), "n_assigned": len(rows)-sum(missing)+len(changes), "n_missing": sum(missing)-len(changes), "status_counts": dict(collections.Counter(d["status"] for d in diagnostics)), "all_previous_values_and_other_fields_preserved": True, "input_sha256": bio.sha(Path(args.base).read_bytes()), "output_sha256": bio.sha(Path(args.out).read_bytes()), "previous_audit_sha256": bio.sha(Path(args.previous_audit).read_bytes()), "script_sha256": bio.sha(Path(__file__).read_bytes()), "geometry_module_sha256": bio.sha(Path(__file__).with_name("metadata_pipeline_common.py").read_bytes()), "source_module_sha256": bio.sha(Path(SPEC.origin).read_bytes()), "timestamp_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(), "command": sys.argv, "retrievals": [source_record] + [source.records[url] for url in sorted(source.records)]}
    bio.write_json(audit_dir / "cardinal_recovery_manifest.json", manifest)
    print(json.dumps({k: manifest[k] for k in ("status", "n_rows", "n_previously_available", "n_recovered", "n_assigned", "n_missing", "status_counts")}, indent=2))
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True)
    parser.add_argument("--previous-audit", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-dir", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--radial-steps-m", default=DEFAULT_STEPS)
    parser.add_argument("--max-anchor-distance-m", type=float, default=50000)
    parser.add_argument("--workers", type=int, default=2, choices=range(1, 5))
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    try:
        rc = run(args)
    except Exception as exc:
        bio.write_json(Path(args.audit_dir) / "STATUS.json", {"status": "FAIL", "error": str(exc)})
        raise
    bio.write_json(Path(args.audit_dir) / "STATUS.json", {"status": "PASS"})
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

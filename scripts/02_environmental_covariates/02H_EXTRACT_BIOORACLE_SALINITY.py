#!/usr/bin/env python3
"""Append two Bio-ORACLE v3 marine-reference covariates without filtering rows.

Only native grid cells are assigned automatically. Alternative coastal cells need
an explicit, coordinate-bound review with supporting evidence. All diagnostics
remain outside the main metadata table. Requires Python 3.9+, standard library.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import math
import platform
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

DATASET = "so_baseline_2000_2019_depthsurf"
BASE_URL = "https://erddap.bio-oracle.org/erddap"
TIMES = ("2000-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
VARIABLES = ("so_mean", "so_range")
COLUMNS = ("marine_salinity_mean", "marine_salinity_annual_range")
MISSING = {"", "na", "nan", "null", "none"}
BLOCKED_COORDINATES = ("conflicto", "por_verificar", "por_resolver", "pendientes_validacion")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_csv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_metadata(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if not fields or len(fields) != len(set(fields)):
        raise ValueError("Missing or duplicate column names")
    if not {"run_accession", "lat", "lon"}.issubset(fields):
        raise ValueError("Required columns: run_accession, lat, lon")
    ids = [r["run_accession"] for r in rows]
    if not rows or any(not x.strip() or x != x.strip() for x in ids) or len(set(ids)) != len(ids):
        raise ValueError("run_accession must be nonempty, unique and without surrounding whitespace")
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("Malformed CSV row")
    return fields, rows


def coordinates(row):
    try:
        lat, lon = float(row["lat"]), float(row["lon"])
    except (ValueError, TypeError):
        return None
    return (lat, lon) if math.isfinite(lat) and math.isfinite(lon) and -90 <= lat <= 90 and -180 <= lon <= 180 else None


def cell_index(lat, lon):
    # Cell containing the point; exact boundaries use north/east, poles clamp.
    return min(3599, max(0, math.floor(round((lat + 90) / .05, 10)))), min(7199, max(0, math.floor(round((lon + 180) / .05, 10))))


def cell_center(iy, ix):
    return round(-89.975 + iy * .05, 6), round(-179.975 + ix * .05, 6)


def distance_km(a, b):
    lat1, lat2 = math.radians(a[0]), math.radians(b[0])
    dlat, dlon = lat2 - lat1, math.radians(b[1] - a[1])
    h = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    return 6371.0088 * 2 * math.asin(math.sqrt(min(1, max(0, h))))


def fetch(url, cache_dir, offline=False):
    path = Path(cache_dir) / (sha(url.encode()) + ".payload")
    record = path.with_suffix(".json")
    if path.exists() and record.exists():
        data, provenance = path.read_bytes(), json.loads(record.read_text())
        if provenance["url"] != url or provenance["sha256"] != sha(data):
            raise ValueError(f"Cache checksum/URL mismatch: {path.name}")
        return data, provenance
    if offline:
        raise FileNotFoundError(f"Missing verified offline cache: {url}")
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "MAGENTA-BioORACLE/1.0 scientific reproducible extraction"})
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
            break
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_bytes(data)
    temporary.replace(path)
    provenance = {"url": url, "sha256": sha(data), "bytes": len(data), "retrieved_utc": datetime.now(timezone.utc).isoformat()}
    write_json(record, provenance)
    return data, provenance


def validate_source(data):
    rows = json.loads(data)["table"]["rows"]
    attrs = {(r[1], r[2]): r[4] for r in rows if r[0] == "attribute"}
    dims = {r[1]: r[4] for r in rows if r[0] == "dimension"}
    if "V3.0" not in attrs.get(("NC_GLOBAL", "source"), ""):
        raise ValueError("Unexpected Bio-ORACLE version")
    if not dims.get("time", "").startswith("nValues=2,"):
        raise ValueError("Expected exactly two decadal layers")
    for axis, count in (("latitude", 3600), ("longitude", 7200)):
        if not dims.get(axis, "").startswith(f"nValues={count},"):
            raise ValueError("Unexpected grid dimensions")
    for name, expected in (("geospatial_lat_min", -89.975), ("geospatial_lat_max", 89.975), ("geospatial_lon_min", -179.975), ("geospatial_lon_max", 179.975)):
        if not math.isclose(float(attrs[("NC_GLOBAL", name)]), expected, abs_tol=1e-5):
            raise ValueError("Unexpected grid extent")
    if attrs.get(("NC_GLOBAL", "time_coverage_start")) != TIMES[0] or attrs.get(("NC_GLOBAL", "time_coverage_end")) != TIMES[1]:
        raise ValueError("Unexpected baseline time coordinates")
    for axis, expected in (("lat", .05), ("lon", .05)):
        if not math.isclose(float(attrs[("NC_GLOBAL", f"geospatial_{axis}_resolution")]), expected):
            raise ValueError("Unexpected grid resolution")
    for var in VARIABLES:
        if attrs.get((var, "units")) != "PSU":
            raise ValueError("Unexpected salinity units")
    return attrs


def patch_url(iy, ix):
    # Seven by seven cells, approximately 20 km to each side. No edge wrapping:
    # candidates at the dateline are conservative; native cells are still valid.
    y0, y1, x0, x1 = max(0, iy - 3), min(3599, iy + 3), max(0, ix - 3), min(7199, ix + 3)
    selection = f"[0:1][{y0}:1:{y1}][{x0}:1:{x1}]"
    query = ",".join(v + selection for v in VARIABLES)
    return BASE_URL + "/griddap/" + DATASET + ".csv?" + urllib.parse.quote(query, safe=",:[]"), (y0, y1, x0, x1)


def parse_patch(data, bounds):
    reader = csv.reader(io.StringIO(data.decode("utf-8")))
    if next(reader) != ["time", "latitude", "longitude", *VARIABLES]:
        raise ValueError("Unexpected ERDDAP CSV header")
    if next(reader) != ["UTC", "degrees_north", "degrees_east", "PSU", "PSU"]:
        raise ValueError("Unexpected ERDDAP CSV units")
    cells = collections.defaultdict(dict)
    y0, y1, x0, x1 = bounds
    for row in reader:
        if len(row) != 5 or row[0] not in TIMES:
            raise ValueError("Unexpected ERDDAP row/time")
        lat, lon = float(row[1]), float(row[2])
        index = cell_index(lat, lon)
        if not (y0 <= index[0] <= y1 and x0 <= index[1] <= x1):
            raise ValueError("Source cell outside request")
        if any(abs(a-b) > 1e-4 for a, b in zip((lat, lon), cell_center(*index))):
            raise ValueError("Source coordinates do not match grid centers")
        if row[0] in cells[index]:
            raise ValueError("Duplicate source cell/time")
        values = tuple(float(v) for v in row[3:])
        if any(math.isfinite(v) and v < 0 and not math.isclose(v, -9999.9) for v in values):
            raise ValueError("Unexpected negative salinity statistic")
        cells[index][row[0]] = tuple(v if math.isfinite(v) and v >= 0 else None for v in values)
    if len(cells) != (y1-y0+1)*(x1-x0+1) or any(set(v) != set(TIMES) for v in cells.values()):
        raise ValueError("Incomplete source grid/decadal response")
    return dict(cells)


def aggregate(periods):
    if set(periods) != set(TIMES) or any(v is None for t in TIMES for v in periods[t]):
        return None
    return tuple(sum(periods[t][i] for t in TIMES) / 2 for i in range(2))


def load_reviews(path, rows):
    if not path:
        return {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        decisions = list(csv.DictReader(handle))
    by_id = {r["run_accession"]: r for r in rows}
    out = {}
    for d in decisions:
        key = d["run_accession"]
        if key in out or key not in by_id or d["decision"] not in {"hold", "approve"}:
            raise ValueError("Review must have unique existing IDs and hold/approve decisions")
        if (d["lat"], d["lon"]) != (by_id[key]["lat"], by_id[key]["lon"]):
            raise ValueError(f"Stale review coordinates: {key}")
        if not d.get("reason", "").strip() or not d.get("evidence", "").strip():
            raise ValueError(f"Review needs a reason and evidence: {key}")
        out[key] = d
    return out


def run(args):
    fields, rows = read_metadata(args.base)
    if Path(args.base).resolve() == Path(args.out).resolve():
        raise ValueError("Write a separate output; do not overwrite the source snapshot")
    if any(c in fields for c in COLUMNS) and not args.replace_salinity:
        raise ValueError("Salinity columns already exist; explicitly use --replace-salinity")
    audits = Path(args.audit_dir)
    audits.mkdir(parents=True, exist_ok=True)
    reviews = load_reviews(args.review, rows)
    source, source_record = fetch(BASE_URL + "/info/" + DATASET + "/index.json", args.cache_dir, args.offline)
    validate_source(source)
    (audits / "source_metadata.json").write_bytes(source)
    groups = sorted({cell_index(*coordinates(r)) for r in rows if coordinates(r) is not None})
    def retrieve(index):
        url, bounds = patch_url(*index)
        data, record = fetch(url, args.cache_dir, args.offline)
        return index, parse_patch(data, bounds), record
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        fetched = list(executor.map(retrieve, groups))
    patches = {index: cells for index, cells, _ in fetched}
    candidates, sample_audit, output = [], [], []
    coord_seen = set()
    for original in rows:
        row = dict(original)
        coord = coordinates(row)
        diag = {k: row.get(k, "") for k in ("run_accession", "lat", "lon", "Country", "Locality", "coordinate_review_status", "review_flags")}
        diag.update(status="missing_or_invalid_coordinates", reason="", cell_lat="", cell_lon="", distance_km="", native_cell_lat="", native_cell_lon="", nearest_candidate_km="", review_evidence="")
        values = None
        selected_periods = None
        if coord is not None:
            native = cell_index(*coord)
            cells = patches[native]
            diag["native_cell_lat"], diag["native_cell_lon"] = cell_center(*native)
            ranked = sorted((distance_km(coord, cell_center(*idx)), idx, aggregate(v)) for idx, v in cells.items() if aggregate(v) is not None)
            if ranked:
                diag["nearest_candidate_km"] = round(ranked[0][0], 6)
            if coord not in coord_seen:
                for distance, idx, val in ranked:
                    candidates.append({"lat": row["lat"], "lon": row["lon"], "cell_lat": cell_center(*idx)[0], "cell_lon": cell_center(*idx)[1], "distance_km": round(distance, 6), "is_native": idx == native, COLUMNS[0]: val[0], COLUMNS[1]: val[1]})
                coord_seen.add(coord)
            choice = native
            decision = reviews.get(row["run_accession"])
            blocked = any(flag in row.get("coordinate_review_status", "").lower() for flag in BLOCKED_COORDINATES)
            if decision:
                diag["reason"], diag["review_evidence"] = decision["reason"], decision["evidence"]
                if decision["decision"] == "hold":
                    diag["status"] = "held_by_spatial_review"
                    choice = None
                else:
                    target = (float(decision["cell_lat"]), float(decision["cell_lon"]))
                    choice = cell_index(*target)
                    if any(abs(a-b) > 1e-6 for a, b in zip(target, cell_center(*choice))):
                        raise ValueError("Reviewed coordinates must be exact cell centers")
                    if distance_km(coord, target) > args.max_review_distance_km:
                        raise ValueError("Reviewed cell exceeds explicit maximum distance")
                    if choice not in cells or aggregate(cells[choice]) is None:
                        raise ValueError("Reviewed cell missing from complete cached patch")
                    diag["status"] = "assigned_reviewed_native" if choice == native else "assigned_reviewed_coastal"
            elif blocked:
                diag["status"] = "held_existing_coordinate_warning"
                diag["reason"] = row["coordinate_review_status"]
                choice = None
            elif aggregate(cells[native]) is None:
                diag["status"] = "no_native_coverage"
                diag["reason"] = "No complete native cell; nearby values are candidates only, pending waterbody review"
                choice = None
            else:
                diag["status"] = "assigned_native"
                diag["reason"] = "Regional marine reference at supplied coordinates; local hydrological representativeness unvalidated"
            if choice is not None:
                selected_periods = cells[choice]
                values = aggregate(selected_periods)
                diag["cell_lat"], diag["cell_lon"] = cell_center(*choice)
                diag["distance_km"] = round(distance_km(coord, cell_center(*choice)), 6)
        for i, column in enumerate(COLUMNS):
            row[column] = format(values[i], ".10g") if values else "NA"
            diag[column] = row[column]
            for j, decade in enumerate(("2000_2009", "2010_2019")):
                diag[VARIABLES[i] + "_" + decade] = selected_periods[TIMES[j]][i] if selected_periods else ""
        sample_audit.append(diag)
        output.append(row)
    # Every old field is compared as a string; no NA inference or number coercion.
    old_fields = [f for f in fields if f not in COLUMNS]
    assert len(rows) == len(output)
    assert all(all(a[f] == b[f] for f in old_fields) for a, b in zip(rows, output))
    final_fields = list(fields)
    for c in COLUMNS:
        if c not in final_fields:
            final_fields.insert(final_fields.index("bio19") + 1 + list(COLUMNS).index(c) if "bio19" in final_fields else len(final_fields), c)
    write_csv(args.out, final_fields, output)
    write_csv(audits / "sample_salinity_audit.csv", list(sample_audit[0]), sample_audit)
    candidate_fields = ["lat", "lon", "cell_lat", "cell_lon", "distance_km", "is_native", *COLUMNS]
    write_csv(audits / "candidate_cells.csv", candidate_fields, candidates)
    summary = {"status": "PASS", "dataset": DATASET, "source_variables": list(VARIABLES), "output_columns": list(COLUMNS), "reference_period": "2000-2019", "decade_labels": list(TIMES), "aggregation": "equal mean of the two complete decadal statistics; range is mean annual max-minus-min", "surface_depth_m": [0, .49], "resolution_degrees": .05, "units_as_published": "PSU", "interpretation": "Marine surface climatological context, not measured porewater salinity", "assignment_policy": "native cell only unless explicit coordinate-bound review; known coordinate warnings held", "max_review_distance_km": args.max_review_distance_km, "n_rows": len(rows), "n_unique_valid_coordinates": len(coord_seen), "n_assigned": sum(d["status"].startswith("assigned_") for d in sample_audit), "status_counts": dict(collections.Counter(d["status"] for d in sample_audit)), "all_original_fields_preserved": True, "input_sha256": sha(Path(args.base).read_bytes()), "output_sha256": sha(Path(args.out).read_bytes()), "review_sha256": sha(Path(args.review).read_bytes()) if args.review else None, "script_sha256": sha(Path(__file__).read_bytes()), "python": platform.python_version(), "timestamp_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv, "retrievals": [source_record] + [r for _, _, r in fetched]}
    write_json(audits / "salinity_manifest.json", summary)
    write_json(audits / "salinity_dictionary.json", {COLUMNS[0]: {"source": "so_mean", "definition_es": "Salinidad marina superficial media de referencia 2000–2019"}, COLUMNS[1]: {"source": "so_range", "definition_es": "Media del rango anual de salinidad superficial 2000–2019; no es el rango extremo de todo el periodo"}, "shared": {"dataset": DATASET, "units": "PSU (convención del proveedor)", "missing": "NA; motivo en sample_salinity_audit.csv", "source_documentation": "https://www.bio-oracle.org/documentation.php", "doi": "10.1111/geb.13813"}})
    print(json.dumps({k: summary[k] for k in ("status", "n_rows", "n_unique_valid_coordinates", "n_assigned", "status_counts")}, indent=2))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--audit-dir", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--review", help="External CSV: run_accession,lat,lon,decision,cell_lat,cell_lon,reason,evidence")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--replace-salinity", action="store_true")
    parser.add_argument("--workers", type=int, default=2, choices=range(1, 5))
    parser.add_argument("--max-review-distance-km", type=float, default=10)
    args = parser.parse_args(argv)
    if not math.isfinite(args.max_review_distance_km) or not 0 < args.max_review_distance_km <= 20:
        parser.error("--max-review-distance-km must be >0 and <=20")
    try:
        result = run(args)
    except Exception as exc:
        write_json(Path(args.audit_dir) / "STATUS.json", {"status": "FAIL", "error": str(exc)})
        raise
    write_json(Path(args.audit_dir) / "STATUS.json", {"status": "PASS"})
    return result


if __name__ == "__main__":
    raise SystemExit(main())

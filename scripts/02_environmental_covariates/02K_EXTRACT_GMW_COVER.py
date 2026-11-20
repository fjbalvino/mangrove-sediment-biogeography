#!/usr/bin/env python3
"""Mapped mangrove area as percent of a fixed geodesic 2 km landscape buffer.

Pinned Global Mangrove Watch 4.1.12, year 2020. No coordinate displacement,
cohort filtering or substitution of zero by NA. Audit stays outside metadata.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import platform
import shutil
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyogrio
from pyproj import CRS, Geod, Transformer
from shapely import make_valid, transform
from shapely.geometry import Polygon
from shapely.ops import unary_union

SPEC = importlib.util.spec_from_file_location("metadata_io", Path(__file__).with_name("02H_EXTRACT_BIOORACLE_SALINITY.py"))
bio = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bio)
RECORD_ID = "21346457"
VERSION = "4.1.12"
YEAR = 2020
RADIUS_M = 2000
FIELD = "mangrove_cover_2020_2km_pct"
FILENAME = "gmw_v4112_2020_mng_ext_cntry_info_vec.gpkg.gz"
SOURCE_URL = f"https://zenodo.org/api/records/{RECORD_ID}/files/{FILENAME}/content"
SOURCE_MD5 = "5d75f466c31f6a48e54d616ee1b75c91"
SOURCE_BYTES = 167847900
GEOD = Geod(ellps="WGS84")


def file_hash(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acquire_source(cache_dir, offline=False):
    """Download only the pinned year and verify the publisher's checksum."""
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    archive = cache / FILENAME
    if not archive.exists():
        if offline:
            raise FileNotFoundError(f"Missing source archive: {archive}")
        temporary = archive.with_suffix(archive.suffix + ".part")
        with urllib.request.urlopen(SOURCE_URL, timeout=90) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
        if temporary.stat().st_size != SOURCE_BYTES or file_hash(temporary, "md5") != SOURCE_MD5:
            raise ValueError("GMW archive does not match the pinned publisher checksum")
        temporary.replace(archive)
    if archive.stat().st_size != SOURCE_BYTES or file_hash(archive, "md5") != SOURCE_MD5:
        raise ValueError("Cached GMW archive checksum mismatch")
    gpkg = archive.with_suffix("")
    receipt_path = cache / "gmw_source_verification.json"
    receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    if not gpkg.exists() or receipt.get("gpkg_sha256") != file_hash(gpkg):
        temporary = gpkg.with_suffix(".gpkg.part")
        with gzip.open(archive, "rb") as source, temporary.open("wb") as output:
            shutil.copyfileobj(source, output, length=1024 * 1024)
        temporary.replace(gpkg)
    info = pyogrio.read_info(gpkg)
    if info["crs"] != "EPSG:4326" or info["layer_name"] != FILENAME.removesuffix(".gpkg.gz"):
        raise ValueError("Unexpected GMW layer/CRS")
    receipt = {"source_url": SOURCE_URL, "doi": f"10.5281/zenodo.{RECORD_ID}",
               "version": VERSION, "reference_year": YEAR, "license": "CC-BY-4.0",
               "archive_md5": SOURCE_MD5, "archive_sha256": file_hash(archive),
               "gpkg_sha256": file_hash(gpkg), "bytes": SOURCE_BYTES,
               "verified_utc": datetime.now(timezone.utc).isoformat(),
               "layer": info["layer_name"], "features": int(info["features"])}
    bio.write_json(receipt_path, receipt)
    return gpkg, receipt


def buffer_geometry(lat, lon, radius_m=RADIUS_M, vertices=1440):
    """Ellipsoidal geodesic boundary, measured in a local equal-area projection."""
    if not (math.isfinite(lat) and math.isfinite(lon) and -90 <= lat <= 90 and -180 <= lon <= 180):
        raise ValueError("Invalid coordinate")
    if not math.isfinite(radius_m) or radius_m <= 0 or vertices < 360:
        raise ValueError("Invalid geodesic buffer parameters")
    bearings = np.linspace(0, 360, vertices, endpoint=False)
    lons, lats, _ = GEOD.fwd(np.full(vertices, lon), np.full(vertices, lat), bearings, np.full(vertices, radius_m))
    projection = CRS.from_proj4(f"+proj=laea +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs")
    project = Transformer.from_crs("EPSG:4326", projection, always_xy=True).transform
    x, y = project(lons, lats)
    buffer = Polygon(zip(x, y))
    ymin, ymax = float(min(lats)), float(max(lats))
    if max(lons) - min(lons) > 180:
        boxes = [(float(min(lons[lons >= 0])), ymin, 180, ymax),
                 (-180, ymin, float(max(lons[lons < 0])), ymax)]
    else:
        boxes = [(float(min(lons)), ymin, float(max(lons)), ymax)]
    return buffer, project, boxes


def measure_geometries(geometries, buffer, project):
    """Union clipped polygons so borders/overlap cannot double-count area."""
    clipped, repaired = [], 0
    for geometry in geometries:
        if geometry is None or geometry.is_empty:
            continue
        if not geometry.is_valid:
            geometry = make_valid(geometry)
            repaired += 1
        local = transform(geometry, lambda xy: np.column_stack(project(xy[:, 0], xy[:, 1])))
        if not local.is_valid:
            local = make_valid(local)
        part = local.intersection(buffer)
        if not part.is_empty:
            clipped.append(part)
    area = float(unary_union(clipped).area) if clipped else 0.0
    denominator = float(buffer.area)
    if not math.isfinite(area) or not 0 <= area <= denominator * (1 + 1e-9):
        raise ValueError("Impossible mangrove area")
    pct = min(100.0, 100 * area / denominator)
    return {FIELD: pct, "mangrove_area_m2": area, "buffer_area_m2": denominator,
            "invalid_source_geometries_repaired": repaired}


def measure_point(gpkg, point):
    buffer, project, boxes = buffer_geometry(*point)
    features = {}
    for bbox in boxes:
        frame = pyogrio.read_dataframe(gpkg, bbox=bbox, columns=[], fid_as_index=True)
        features.update(zip(frame.index, frame.geometry))
    return {**measure_geometries(features.values(), buffer, project), "source_features_in_bbox": len(features)}


def run(args):
    fields, rows = bio.read_metadata(args.base)
    if FIELD in fields:
        raise ValueError(f"{FIELD} already exists; use the archived input before this stage")
    if Path(args.base).resolve() == Path(args.out).resolve():
        raise ValueError("Input and output paths must differ")
    source, source_record = acquire_source(args.cache_dir, args.offline)
    points = sorted({bio.coordinates(row) for row in rows if bio.coordinates(row) is not None})
    measurements = {}
    for i, point in enumerate(points):
        measurements[point] = measure_point(source, point)
        if (i + 1) % 20 == 0 or i + 1 == len(points):
            print(f"GMW coordinates processed: {i+1}/{len(points)}", flush=True)
    output, audit = [], []
    for row in rows:
        point = bio.coordinates(row)
        value = measurements.get(point)
        result = dict(row)
        result[FIELD] = format(value[FIELD], ".8g") if value is not None else "NA"
        diag = {k: row.get(k, "") for k in ("run_accession", "lat", "lon", "Country", "Locality", "coordinate_review_status", "review_flags")}
        diag.update({FIELD: result[FIELD], "status": "missing_or_invalid_coordinates" if value is None else ("assigned_zero_mapped" if value[FIELD] == 0 else "assigned_mapped_cover"),
                     "reference_year": YEAR, "gmw_version": VERSION, "buffer_radius_m": RADIUS_M,
                     "mangrove_area_m2": "", "buffer_area_m2": "", "source_features_in_bbox": "", "invalid_source_geometries_repaired": "",
                     "interpretation": "Mapped mangrove extent / full landscape buffer (land and water); not canopy density. Zero is no mapped mangrove, not field-confirmed absence.",
                     "coordinate_note": "Original coordinates and inherited warnings retained; no displacement or new coordinate validation."})
        if value is not None:
            diag.update({k: v for k, v in value.items() if k != FIELD})
        output.append(result)
        audit.append(diag)
    insert_at = fields.index("human_footprint_2020") + 1 if "human_footprint_2020" in fields else len(fields)
    out_fields = fields[:insert_at] + [FIELD] + fields[insert_at:]
    assert all(all(before[k] == after[k] for k in fields) for before, after in zip(rows, output))
    audit_dir = Path(args.audit_dir)
    bio.write_csv(args.out, out_fields, output)
    bio.write_csv(audit_dir / "sample_mangrove_cover_audit.csv", list(audit[0]), audit)
    bio.write_csv(audit_dir / "zero_mapped_mangrove.csv", list(audit[0]), [r for r in audit if r["status"] == "assigned_zero_mapped"])
    manifest = {"status": "PASS", "field": FIELD, "n_rows": len(rows), "n_unique_coordinates": len(points),
                "status_counts": dict(collections.Counter(r["status"] for r in audit)),
                "radius_m": RADIUS_M, "buffer_vertices": 1440, "buffer_geodesy": "WGS84 ellipsoid",
                "area_method": "Local Lambert azimuthal equal-area; union of clipped mapped polygons; denominator includes land and water",
                "national_accuracy_correction_applied": False, "all_original_fields_and_order_preserved": True,
                "input_sha256": file_hash(args.base), "output_sha256": file_hash(args.out), "source": source_record,
                "script_sha256": file_hash(__file__), "python": platform.python_version(),
                "packages": {p: importlib.metadata.version(p) for p in ("numpy", "pyproj", "shapely", "pyogrio", "geopandas")},
                "timestamp_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv}
    bio.write_json(audit_dir / "mangrove_cover_manifest.json", manifest)
    bio.write_json(audit_dir / "STATUS.json", {"status": "PASS"})
    print(json.dumps({k: manifest[k] for k in ("status", "n_rows", "n_unique_coordinates", "status_counts")}, indent=2))
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--audit-dir", required=True)
    p.add_argument("--cache-dir", required=True)
    p.add_argument("--offline", action="store_true")
    args = p.parse_args()
    try:
        return run(args)
    except Exception as exc:
        bio.write_json(Path(args.audit_dir) / "STATUS.json", {"status": "FAIL", "error": str(exc)})
        raise


if __name__ == "__main__":
    raise SystemExit(main())

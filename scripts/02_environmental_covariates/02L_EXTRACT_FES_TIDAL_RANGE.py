#!/usr/bin/env python3
"""Prepare or compute regional mean tidal range from an external FES2022b atlas.

AVISO access is required separately. Missing atlas files are an input blocker,
never evidence of missing tides. --preflight writes a per-accession pending audit
without adding an entirely empty variable to the master metadata.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

SPEC = importlib.util.spec_from_file_location("metadata_io", Path(__file__).with_name("02H_EXTRACT_BIOORACLE_SALINITY.py"))
bio = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bio)
FIELD = "tidal_range_mean_m"
CONSTITUENTS = ("2N2", "Eps2", "J1", "K1", "K2", "L2", "Lambda2", "M2", "M3", "M4", "M6", "M8", "MKS2", "MN4", "MS4", "MSf", "Mf", "Mm", "Msqm", "Mtm", "Mu2", "N2", "N4", "Nu2", "O1", "P1", "Q1", "R2", "S1", "S2", "S4", "Sa", "Ssa", "T2")
MODEL_DOI = "10.24400/527896/A01-2024.004"
EXPECTED_PYFES = "2026.5.2"
START, END = "2000-01-01", "2020-01-01"
STEP_MINUTES = 15


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atlas_paths(directory):
    return {c: Path(directory) / f"{c.lower()}_fes2022.nc" for c in CONSTITUENTS}


def check_atlas(paths):
    """Require standard non-extrapolated 1/30 degree Cartesian cm/degree files."""
    from netCDF4 import Dataset
    if any("extrapolat" in str(p.parent).lower() for p in paths.values()):
        raise ValueError("Use the non-extrapolated ocean_tide_20241025 atlas")
    hashes = {}
    for name, path in paths.items():
        with Dataset(path) as ds:
            if not {"lat", "lon", "amplitude", "phase"}.issubset(ds.variables):
                raise ValueError(f"Missing FES variables in {path.name}")
            for axis in ("lat", "lon"):
                v = np.asarray(ds[axis][:])
                if v.ndim != 1 or len(v) < 2 or not np.allclose(np.abs(np.diff(v)), 1 / 30, rtol=1e-3, atol=1e-6):
                    raise ValueError(f"Expected the FES2022b 1/30 degree grid: {path.name}/{axis}")
            if getattr(ds["amplitude"], "units", "").lower() not in ("cm", "centimeter", "centimeters", "centimetres"):
                raise ValueError(f"FES amplitude must be in cm: {path.name}")
            if getattr(ds["phase"], "units", "").lower() not in ("degree", "degrees", "deg"):
                raise ValueError(f"FES phase must be in degrees: {path.name}")
        hashes[name] = file_hash(path)
    return hashes


def range_summary(heights_m):
    """Mean |high-low| over successive interior extrema, with parabolic refinement.

Each rising/falling half-cycle has equal weight. Endpoints do not become extrema.
No day boundaries or monthly maxima are used. Any nonfinite height fails.
"""
    heights = np.asarray(heights_m, dtype=float)
    if heights.ndim != 1 or len(heights) < 5 or not np.isfinite(heights).all():
        raise ValueError("Tide series must be finite, one-dimensional and complete")
    if np.ptp(heights) <= 1e-12:
        return {FIELD: 0.0, "n_half_cycles": 0, "range_p05_m": 0.0, "range_p95_m": 0.0}
    maxima, _ = find_peaks(heights)
    minima, _ = find_peaks(-heights)
    extrema = sorted([(int(i), 1) for i in maxima] + [(int(i), -1) for i in minima])
    if len(extrema) < 3:
        raise ValueError("Insufficient turning points for mean tidal range")
    values, kinds = [], []
    for i, kind in extrema:
        before, center, after = heights[i - 1:i + 2]
        denominator = before - 2 * center + after
        offset = 0.5 * (before - after) / denominator if denominator else 0.0
        refined = center - (before - after) ** 2 / (8 * denominator) if denominator and abs(offset) <= 1 else center
        values.append(refined)
        kinds.append(kind)
    if any(a == b for a, b in zip(kinds, kinds[1:])):
        raise ValueError("Tide extrema do not alternate")
    ranges = np.abs(np.diff(values))
    return {FIELD: float(np.mean(ranges)), "n_half_cycles": len(ranges),
            "range_p05_m": float(np.quantile(ranges, .05)), "range_p95_m": float(np.quantile(ranges, .95))}


def read_review(path, rows, max_distance_km=5):
    if path is None:
        return {}
    _, review = bio.read_metadata(path)
    metadata = {r["run_accession"]: r for r in rows}
    result = {}
    for record in review:
        accession = record["run_accession"]
        original = metadata.get(accession)
        if original is None or any(record[c] != original[c] for c in ("lat", "lon")):
            raise ValueError("Review does not match original accession/coordinate strings")
        if record.get("decision") not in ("approve", "hold") or not record.get("reason", "").strip() or not record.get("evidence", "").strip():
            raise ValueError("Each review needs approve/hold, reason and evidence")
        if record["decision"] == "approve":
            point = bio.coordinates({"lat": record.get("anchor_lat"), "lon": record.get("anchor_lon")})
            original_point = bio.coordinates(original)
            if point is None or original_point is None or bio.distance_km(point, original_point) > max_distance_km:
                raise ValueError("Approved anchor missing, invalid, or beyond the 5 km bound")
            record["point"] = point
        result[accession] = record
    return result


class Predictor:
    def __init__(self, paths, audit_dir):
        import pyfes
        import yaml
        if pyfes.__version__ != EXPECTED_PYFES:
            raise ValueError(f"Use validated pyfes=={EXPECTED_PYFES}; found {pyfes.__version__}")
        self.pyfes = pyfes
        self.config_path = Path(audit_dir) / "fes_runtime_config.yaml"
        self.config_path.parent.mkdir(parents=True, exist_ok=True)
        self.config_path.write_text(yaml.safe_dump({"engine": "darwin", "tide": {"cartesian": {
            "amplitude": "amplitude", "phase": "phase", "latitude": "lat", "longitude": "lon", "dynamic": ["A5"],
            "paths": {c: str(p.resolve()) for c, p in paths.items()}}}}))
        self.dates = np.arange(np.datetime64(START), np.datetime64(END), np.timedelta64(STEP_MINUTES, "m")).astype("datetime64[us]")
        self.cache = {}

    def predict(self, point):
        if point in self.cache:
            return self.cache[point]
        lat, lon = point
        # Small per-point windows avoid loading all 34 global arrays into RAM.
        # The primary dataset has no sample at the date line. Refuse to silently
        # clip a crossing window; a reviewed anchor or separate implementation is needed.
        if abs(lon) > 179.8 or abs(lat) > 89.8:
            result = {"status": "grid_boundary_requires_review"}
            self.cache[point] = result
            return result
        config = self.pyfes.config.load(self.config_path, bbox=(lon - .2, lat - .2, lon + .2, lat + .2))
        config.settings.with_num_threads(2).with_time_tolerance(0)
        # Reject missing/extrapolated locations before computing twenty years.
        preview = self.pyfes.evaluate_tide(config.models["tide"], self.dates[:2], np.full(2, lon), np.full(2, lat), settings=config.settings)
        if not (np.isfinite(preview[0]).all() and np.isfinite(preview[1]).all()) or np.any(preview[2] <= 0):
            result = {"status": "no_nonextrapolated_model_prediction"}
        else:
            heights, flags = [], []
            for start in range(0, len(self.dates), 100000):
                dates = self.dates[start:start + 100000]
                tide, lp, flag = self.pyfes.evaluate_tide(config.models["tide"], dates, np.full(len(dates), lon), np.full(len(dates), lat), settings=config.settings)
                heights.append((tide + lp) / 100)
                flags.append(flag)
            heights, flags = np.concatenate(heights), np.concatenate(flags)
            if not np.isfinite(heights).all() or np.any(flags <= 0):
                raise ValueError("Incomplete or extrapolated tidal prediction")
            result = {"status": "computed", **range_summary(heights), "minimum_interpolation_points": int(min(flags)), "n_time_steps": len(heights)}
        self.cache[point] = result
        return result


def run(args):
    fields, rows = bio.read_metadata(args.base)
    audit_dir = Path(args.audit_dir)
    paths = atlas_paths(args.atlas_dir)
    missing = [str(p) for p in paths.values() if not p.is_file()]
    base_manifest = {"model": "FES2022b ocean_tide_20241025 (non-extrapolated)", "model_doi": MODEL_DOI,
                     "reference_period": "2000-2019", "step_minutes": STEP_MINUTES, "field": FIELD,
                     "n_rows": len(rows), "missing_atlas_files": missing, "input_sha256": file_hash(args.base),
                     "timestamp_utc": datetime.now(timezone.utc).isoformat(), "script_sha256": file_hash(__file__), "command": sys.argv}
    if args.preflight:
        audit = [{"run_accession": r["run_accession"], "lat": r["lat"], "lon": r["lon"], FIELD: "NA",
                  "status": "missing_or_invalid_coordinates" if bio.coordinates(r) is None else ("awaiting_fes2022_atlas" if missing else "ready_not_extracted")}
                 for r in rows]
        status = "AWAITING_INPUT" if missing else "READY_NOT_EXTRACTED"
        bio.write_csv(audit_dir / "sample_tidal_range_audit.csv", list(audit[0]), audit)
        bio.write_json(audit_dir / "tidal_range_manifest.json", {**base_manifest, "status": status, "master_metadata_modified": False})
        bio.write_json(audit_dir / "STATUS.json", {"status": status})
        print(json.dumps({"status": status, "n_rows": len(rows), "missing_atlas_files": len(missing)}))
        return 2 if missing else 0
    if missing:
        raise FileNotFoundError("Authorized FES2022b atlas required; run --preflight for the exact missing-file list")
    if not args.out or Path(args.out).resolve() == Path(args.base).resolve() or FIELD in fields:
        raise ValueError("Use a different output path and input without an existing tidal-range column")
    hashes = check_atlas(paths)
    reviews = read_review(args.review, rows)
    predictor = Predictor(paths, audit_dir)
    output, audit = [], []
    for row in rows:
        point = bio.coordinates(row)
        review = reviews.get(row["run_accession"])
        anchor = review["point"] if review and review["decision"] == "approve" else point
        diag = {k: row.get(k, "") for k in ("run_accession", "lat", "lon", "Country", "Locality", "coordinate_review_status", "review_flags")}
        diag.update({FIELD: "NA", "status": "missing_or_invalid_coordinates", "anchor_lat": "", "anchor_lon": "", "distance_km": "",
                     "minimum_interpolation_points": "", "n_time_steps": "", "n_half_cycles": "", "range_p05_m": "", "range_p95_m": "",
                     "reason": review["reason"] if review else "", "evidence": review["evidence"] if review else ""})
        if review and review["decision"] == "hold":
            diag["status"] = "held_by_spatial_review"
        elif anchor is not None:
            prediction = predictor.predict(anchor)
            diag.update(anchor_lat=anchor[0], anchor_lon=anchor[1], distance_km=bio.distance_km(point, anchor), status=prediction["status"])
            if prediction["status"] == "computed":
                diag.update({k: v for k, v in prediction.items() if k not in ("status", FIELD)})
                if prediction["minimum_interpolation_points"] < 4 and not review:
                    diag["status"] = "partial_coastal_interpolation_requires_review"
                else:
                    diag["status"] = "assigned_reviewed" if review else "assigned_native_interpolation"
                    diag[FIELD] = format(prediction[FIELD], ".8g")
        output.append({**row, FIELD: diag[FIELD]})
        audit.append(diag)
    insertion = fields.index("marine_salinity_annual_range") + 1 if "marine_salinity_annual_range" in fields else len(fields)
    out_fields = fields[:insertion] + [FIELD] + fields[insertion:]
    assert all(all(old[k] == new[k] for k in fields) for old, new in zip(rows, output))
    bio.write_csv(args.out, out_fields, output)
    bio.write_csv(audit_dir / "sample_tidal_range_audit.csv", list(audit[0]), audit)
    bio.write_json(audit_dir / "tidal_range_manifest.json", {**base_manifest, "status": "PASS", "atlas_sha256": hashes,
                   "pyfes": EXPECTED_PYFES, "output_sha256": file_hash(args.out), "all_original_fields_preserved": True,
                   "review_sha256": file_hash(args.review) if args.review else None,
                   "status_counts": dict(collections.Counter(r["status"] for r in audit)),
                   "interpretation": "Regional astronomical tide range; not plot hydroperiod, current speed, surge or measured sediment conditions"})
    bio.write_json(audit_dir / "STATUS.json", {"status": "PASS"})
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True)
    parser.add_argument("--out")
    parser.add_argument("--atlas-dir", required=True)
    parser.add_argument("--audit-dir", required=True)
    parser.add_argument("--review")
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    try:
        return run(args)
    except Exception as exc:
        bio.write_json(Path(args.audit_dir) / "STATUS.json", {"status": "FAIL", "error": str(exc)})
        raise


if __name__ == "__main__":
    raise SystemExit(main())

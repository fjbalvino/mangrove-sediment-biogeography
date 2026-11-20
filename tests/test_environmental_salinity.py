"""Offline scientific and row-preservation contracts for marine covariates."""
import csv
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/02_environmental_covariates/02H_EXTRACT_BIOORACLE_SALINITY.py"
SPEC = importlib.util.spec_from_file_location("biooracle_salinity", SCRIPT)
bio = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bio)


def fixture_patch(bounds, missing_native=False):
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(["time", "latitude", "longitude", *bio.VARIABLES])
    writer.writerow(["UTC", "degrees_north", "degrees_east", "PSU", "PSU"])
    for decade, timestamp in enumerate(bio.TIMES):
        for iy in range(bounds[0], bounds[1]+1):
            for ix in range(bounds[2], bounds[3]+1):
                vals = (0, 2) if decade == 0 else (4, 6)
                if missing_native and (iy, ix) == (1800, 3600):
                    vals = ("NaN", "-9999.9")
                writer.writerow([timestamp, *bio.cell_center(iy, ix), *vals])
    return out.getvalue().encode()


class SalinityContracts(unittest.TestCase):
    def test_two_complete_decades_and_mean_annual_range(self):
        self.assertEqual(bio.aggregate(dict(zip(bio.TIMES, [(0, 2), (4, 6)]))), (2, 4))
        self.assertIsNone(bio.aggregate({bio.TIMES[0]: (0, 2)}))
        self.assertIsNone(bio.aggregate(dict(zip(bio.TIMES, [(0, 2), (4, None)]))))

    def test_boundaries_and_dateline_distance(self):
        self.assertEqual(bio.cell_index(0, 0), (1800, 3600))
        self.assertEqual(bio.cell_index(-90, -180), (0, 0))
        self.assertEqual(bio.cell_index(90, 180), (3599, 7199))
        self.assertLess(bio.distance_km((0, 179.99), (0, -179.99)), 3)
        self.assertIsNone(bio.coordinates({"lat": "NaN", "lon": "0"}))
        self.assertIsNone(bio.coordinates({"lat": "91", "lon": "0"}))

    def test_source_completeness_units_fill_and_zero(self):
        bounds = (1800, 1800, 3600, 3600)
        data = fixture_patch(bounds)
        self.assertEqual(bio.aggregate(bio.parse_patch(data, bounds)[(1800, 3600)]), (2, 4))
        self.assertIsNone(bio.aggregate(bio.parse_patch(fixture_patch(bounds, True), bounds)[(1800, 3600)]))
        with self.assertRaises(ValueError):
            bio.parse_patch(data.replace(b"PSU", b"ppt"), bounds)
        with self.assertRaises(ValueError):
            bio.parse_patch(b"\n".join(data.splitlines()[:-1]) + b"\n", bounds)
        with self.assertRaises(ValueError):
            bio.parse_patch(data + data.splitlines()[-1] + b"\n", bounds)

    def test_no_row_loss_no_coercion_and_no_automatic_nearest_fill(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fields = ["run_accession", "lat", "lon", "bio19", "untouched", "coordinate_review_status"]
            rows = [dict(zip(fields, x)) for x in [
                ("S1", "0.0100", "0.0100", "001", "NA", ""),
                ("S2", "0.0100", "0.0100", "2.000", "comma,quote\"", ""),
                ("S3", "", "", "", "", ""),
                ("S4", "0.11", "0.11", "5", "keep", "conflicto_por_resolver"),
                ("S5", "0.11", "0.11", "5", "keep", ""),
            ]]
            bio.write_csv(root / "in.csv", fields, rows)
            args = SimpleNamespace(base=root / "in.csv", out=root / "out.csv", audit_dir=root / "audit", cache_dir=root / "cache", review=None, offline=True, replace_salinity=False, workers=1, max_review_distance_km=10)
            def fake_fetch(url, cache, offline):
                if "/info/" in url:
                    return b"{}", {"url": url}
                for index in [(1800, 3600), (1802, 3602)]:
                    target, bounds = bio.patch_url(*index)
                    if target == url:
                        return fixture_patch(bounds, True), {"url": url}
                raise AssertionError(url)
            with patch.object(bio, "fetch", side_effect=fake_fetch), patch.object(bio, "validate_source"):
                bio.run(args)
            out_fields, out = bio.read_metadata(args.out)
            self.assertEqual(len(out), len(rows))
            self.assertEqual(len(out_fields), len(fields)+2)
            self.assertTrue(all(a[k] == b[k] for a, b in zip(rows, out) for k in fields))
            self.assertEqual([r[bio.COLUMNS[0]] for r in out], ["NA", "NA", "NA", "NA", "2"])
            self.assertEqual([r[bio.COLUMNS[1]] for r in out], ["NA", "NA", "NA", "NA", "4"])
            with (root / "audit/sample_salinity_audit.csv").open() as handle:
                audit = list(csv.DictReader(handle))
            self.assertEqual(audit[0]["status"], "no_native_coverage")
            self.assertTrue(audit[0]["nearest_candidate_km"])
            first = args.out.read_bytes()
            with patch.object(bio, "fetch", side_effect=fake_fetch), patch.object(bio, "validate_source"):
                bio.run(args)
            self.assertEqual(first, args.out.read_bytes())
            decisions = [
                {"run_accession": "S1", "lat": "0.0100", "lon": "0.0100", "decision": "approve", "cell_lat": ".075", "cell_lon": ".025", "reason": "same coastal waterbody", "evidence": "reviewed source"},
                {"run_accession": "S5", "lat": "0.11", "lon": "0.11", "decision": "hold", "cell_lat": "", "cell_lon": "", "reason": "barrier separates native cell", "evidence": "reviewed source"},
            ]
            args.review = root / "reviews.csv"
            bio.write_csv(args.review, list(decisions[0]), decisions)
            with patch.object(bio, "fetch", side_effect=fake_fetch), patch.object(bio, "validate_source"):
                bio.run(args)
            _, reviewed = bio.read_metadata(args.out)
            self.assertEqual([r[bio.COLUMNS[0]] for r in reviewed], ["2", "NA", "NA", "NA", "NA"])
            args.max_review_distance_km = 1
            with patch.object(bio, "fetch", side_effect=fake_fetch), patch.object(bio, "validate_source"), self.assertRaisesRegex(ValueError, "maximum distance"):
                bio.run(args)

    def test_duplicate_ids_and_stale_reviews_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "test.csv"
            p.write_text("run_accession,lat,lon\nA,1,2\nA,3,4\n")
            with self.assertRaises(ValueError):
                bio.read_metadata(p)
            p.write_text("run_accession,lat,lon,decision,reason,evidence\nA,1,3,hold,test,source\n")
            with self.assertRaises(ValueError):
                bio.load_reviews(p, [{"run_accession": "A", "lat": "1", "lon": "2"}])

    def test_cache_tampering_or_missing_offline_cache_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                bio.fetch("https://example.org/test", tmp, True)
            p = Path(tmp) / (bio.sha(b"https://example.org/test") + ".payload")
            p.write_bytes(b"changed")
            bio.write_json(p.with_suffix(".json"), {"url": "https://example.org/test", "sha256": bio.sha(b"original")})
            with self.assertRaises(ValueError):
                bio.fetch("https://example.org/test", tmp, True)


if __name__ == "__main__":
    unittest.main()

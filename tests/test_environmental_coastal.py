"""Geometry, zero/missing distinction, tidal ranges and metadata preservation."""
import csv
import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import geopandas as gpd
import numpy as np
from shapely.geometry import box

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts/02_environmental_covariates"


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gmw = load("gmw", "02K_EXTRACT_GMW_COVER.py")
tide = load("tidal_range", "02L_EXTRACT_FES_TIDAL_RANGE.py")


class MangroveGeometry(unittest.TestCase):
    def test_buffer_area_and_latitude_independence(self):
        areas = [gmw.buffer_geometry(lat, -88)[0].area for lat in (0, 20, 60)]
        for area in areas:
            self.assertAlmostEqual(area / (math.pi * 2000 ** 2), 1, places=5)
        self.assertLess(max(areas) - min(areas), 1)

    def test_union_prevents_double_count_and_clips_edges(self):
        buffer = box(-1, -1, 1, 1)
        left_half = box(-2, -2, 0, 2)
        result = gmw.measure_geometries([left_half, left_half], buffer, lambda x, y: (x, y))
        self.assertEqual(result[gmw.FIELD], 50)

    def test_empty_extent_is_zero(self):
        result = gmw.measure_geometries([], box(-1, -1, 1, 1), lambda x, y: (x, y))
        self.assertEqual(result[gmw.FIELD], 0)

    def test_dateline_window_is_split(self):
        buffer, _, boxes = gmw.buffer_geometry(0, 179.999)
        self.assertEqual(len(boxes), 2)
        self.assertLess(sum(b[2] - b[0] for b in boxes), .1)
        self.assertAlmostEqual(buffer.area / (math.pi * 2000 ** 2), 1, places=5)

    def test_full_run_preserves_rows_strings_order_and_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            source = p / "tiny.gpkg"
            gpd.GeoDataFrame(geometry=[box(-.005, -.005, .005, .005)], crs="EPSG:4326").to_file(source, driver="GPKG")
            fields = ["run_accession", "lat", "lon", "human_footprint_2020", "marine_salinity_mean", "review_flags"]
            rows = [dict(zip(fields, values)) for values in [
                ("A", "0.000", "0", "001", "34.230", "warning"),
                ("B", "1", "1", "NA", "0", ""),
                ("C", "", "", "9", "NA", "missing")]]
            gmw.bio.write_csv(p / "in.csv", fields, rows)
            args = SimpleNamespace(base=p / "in.csv", out=p / "out.csv", audit_dir=p / "audit", cache_dir=p / "cache", offline=True)
            with patch.object(gmw, "acquire_source", return_value=(source, {"fixture": True})):
                self.assertEqual(gmw.run(args), 0)
            out_fields, out = gmw.bio.read_metadata(p / "out.csv")
            self.assertEqual([c for c in out_fields if c != gmw.FIELD], fields)
            self.assertEqual([{k: r[k] for k in fields} for r in out], rows)
            self.assertGreater(float(out[0][gmw.FIELD]), 0)
            self.assertEqual(out[1][gmw.FIELD], "0")
            self.assertEqual(out[2][gmw.FIELD], "NA")
            with self.assertRaises(ValueError):
                gmw.run(SimpleNamespace(**{**vars(args), "base": p / "out.csv", "out": p / "other.csv"}))


class TidalRangeContracts(unittest.TestCase):
    def test_range_is_twice_amplitude_not_daily_extreme(self):
        hours = np.arange(0, 30 * 24, .25)
        for period in (12, 24):
            summary = tide.range_summary(1.5 * np.cos(2 * np.pi * hours / period + .37) + 9)
            self.assertAlmostEqual(summary[tide.FIELD], 3, places=4)

    def test_nonfinite_values_are_not_silently_dropped(self):
        with self.assertRaises(ValueError):
            tide.range_summary([0, 1, np.nan, -1, 0, 1, 0])

    def test_constant_tide_has_zero_range(self):
        self.assertEqual(tide.range_summary(np.full(20, 3))[tide.FIELD], 0)

    def test_preflight_reports_blocker_without_modifying_master(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            rows = [{"run_accession": "A", "lat": "20", "lon": "-88"}, {"run_accession": "B", "lat": "", "lon": ""}]
            tide.bio.write_csv(p / "in.csv", list(rows[0]), rows)
            original = (p / "in.csv").read_bytes()
            args = SimpleNamespace(base=p / "in.csv", out=p / "out.csv", atlas_dir=p / "absent", audit_dir=p / "audit", preflight=True, review=None)
            self.assertEqual(tide.run(args), 2)
            self.assertFalse((p / "out.csv").exists())
            self.assertEqual((p / "in.csv").read_bytes(), original)
            report = json.loads((p / "audit/tidal_range_manifest.json").read_text())
            self.assertEqual(report["status"], "AWAITING_INPUT")
            self.assertEqual(len(report["missing_atlas_files"]), 34)
            args.preflight = False
            with self.assertRaises(FileNotFoundError):
                tide.run(args)

    def test_review_is_bound_to_original_coordinates_and_distance(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "review.csv"
            rows = [{"run_accession": "A", "lat": "20.0", "lon": "-88.0"}]
            review = {**rows[0], "decision": "approve", "anchor_lat": "20.01", "anchor_lon": "-88", "reason": "same estuary", "evidence": "test review"}
            tide.bio.write_csv(p, list(review), [review])
            self.assertIn("A", tide.read_review(p, rows))
            review["anchor_lat"] = "21"
            tide.bio.write_csv(p, list(review), [review])
            with self.assertRaises(ValueError):
                tide.read_review(p, rows)
            review["lat"] = "20"
            tide.bio.write_csv(p, list(review), [review])
            with self.assertRaises(ValueError):
                tide.read_review(p, rows)

    def test_real_pyfes_engine_with_synthetic_constituents(self):
        import pyfes
        dates = np.arange(np.datetime64("2001-01-01"), np.datetime64("2001-02-01"), np.timedelta64(15, "m")).astype("datetime64[us]")
        constants = {c: (0., 0.) for c in tide.CONSTITUENTS}
        constants["S2"] = (100., 0.)
        settings = pyfes.FESSettings().with_num_threads(1).with_time_tolerance(0).with_compute_long_period_equilibrium(False)
        short, long = pyfes.evaluate_tide_from_constituents(constants, dates, latitude=20, settings=settings)
        result = tide.range_summary((short + long) / 100)
        self.assertAlmostEqual(result[tide.FIELD], 2, places=3)

    def test_full_tidal_pipeline_on_synthetic_atlas(self):
        from netCDF4 import Dataset
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            atlas = p / "synthetic_atlas"
            atlas.mkdir()
            for name, path in tide.atlas_paths(atlas).items():
                with Dataset(path, "w") as ds:
                    for axis, size, origin in (("lat", 16, 10), ("lon", 19, 270)):
                        ds.createDimension(axis, size)
                        ds.createVariable(axis, "f8", (axis,))[:] = origin + np.arange(size) / 30
                    amplitude = ds.createVariable("amplitude", "f4", ("lat", "lon"))
                    amplitude.units = "cm"
                    amplitude[:] = 100 if name == "S2" else 0
                    phase = ds.createVariable("phase", "f4", ("lat", "lon"))
                    phase.units = "degrees"
                    phase[:] = 0
            rows = [{"run_accession": "A", "lat": "10.250", "lon": "-89.75", "unchanged": "001"},
                    {"run_accession": "B", "lat": "", "lon": "", "unchanged": "NA"}]
            tide.bio.write_csv(p / "input.csv", list(rows[0]), rows)
            args = SimpleNamespace(base=p / "input.csv", out=p / "output.csv", atlas_dir=atlas, audit_dir=p / "audit", preflight=False, review=None)
            with patch.object(tide, "START", "2001-01-01"), patch.object(tide, "END", "2001-02-01"):
                self.assertEqual(tide.run(args), 0)
            _, output = tide.bio.read_metadata(p / "output.csv")
            self.assertEqual([{k: r[k] for k in rows[0]} for r in output], rows)
            self.assertAlmostEqual(float(output[0][tide.FIELD]), 2, delta=.02)
            self.assertEqual(output[1][tide.FIELD], "NA")


if __name__ == "__main__":
    unittest.main()

"""Check image-coordinate geometry against independently known synthetic controls."""
import importlib.util
import math
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/'scripts/12_gene_environment/12C_RECOMPUTE_FIGURE_GEOREFERENCES.py'
SPEC=importlib.util.spec_from_file_location('figure_georef',PATH)
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GeoreferenceTests(unittest.TestCase):
    def setUp(self):
        p=[[0,0],[10,0],[0,10],[10,10],[5,2],[3,8],[7,9]]
        q=[[2*x+50,3*y+100] for x,y in p]
        self.evidence={'panels':{'A':{'source_points':p,'reference_points':q,
            'reference_export':{'width':1000,'height':1000,'extent':
                {'xmin':0,'xmax':10000,'ymin':0,'ymax':10000,'spatialReference':{'wkid':3857}}}}},
            'sites':[{'label':'synthetic','panel':'A','source_center_px':[4,6],'source_marker_radius_px':2}]}

    def test_actual_extent_pixel_center_and_axis_direction(self):
        result=MODULE.recompute(self.evidence)
        site=result['sites'][0]
        self.assertAlmostEqual(site['lon_unrounded'],math.degrees(585/6378137),places=10)
        expected_lat=math.degrees(2*math.atan(math.exp(8815/6378137))-math.pi/2)
        self.assertAlmostEqual(site['lat_unrounded'],expected_lat,places=10)
        self.assertLess(result['panels']['A']['loo_max_pixels'],1e-8)

    def test_geographic_reference_is_not_misread_as_mercator(self):
        self.evidence['panels']['A']['reference_export']['extent']['spatialReference']['wkid']=4326
        with self.assertRaisesRegex(ValueError,'Web Mercator'):
            MODULE.recompute(self.evidence)


if __name__=='__main__':
    unittest.main()

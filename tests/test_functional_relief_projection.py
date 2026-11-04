"""Map projection checks independent of external tile services."""
from pathlib import Path
import importlib.util
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('relief',ROOT/'scripts/12_gene_environment/template/atlas/scripts/relief_projection.py')
relief=importlib.util.module_from_spec(spec);spec.loader.exec_module(relief)


class ReliefProjection(unittest.TestCase):
    def test_mercator_pixel_centres_round_trip(self):
        for height in [1,32,2160]:
            lat=relief.row_latitudes(height)
            y=(1-np.arcsinh(np.tan(np.radians(lat)))/np.pi)/2
            np.testing.assert_allclose(y,(np.arange(height)+0.5)/height,atol=1e-14)
            np.testing.assert_allclose(lat,-lat[::-1],atol=1e-12)
            self.assertTrue((np.abs(lat)<=relief.MAX_LATITUDE).all())
        self.assertAlmostEqual(relief.MAX_LATITUDE,85.0511287798066)

    def test_land_colour_and_longitude_are_not_flattened(self):
        rgb=np.zeros((90,180,3),dtype='uint8')
        rgb[:,:,0]=np.arange(180)[None,:]
        rgb[:,:,1]=np.arange(90)[:,None]*2
        rgb[:,:,2]=127
        projected=relief.reproject_rgb(rgb)
        self.assertEqual(projected.shape,(180,180,3))
        np.testing.assert_array_equal(projected[:,:,0],np.tile(np.arange(180),(180,1)))
        self.assertTrue((projected[:,:,2]==127).all())
        expected=np.rint((90-relief.row_latitudes(180))-1).astype('uint8')
        np.testing.assert_array_equal(projected[:,45,1],expected)
        self.assertGreater(np.unique(projected[:,:,1]).size,70)
        with self.assertRaises(ValueError):relief.reproject_rgb(np.zeros((40,40,3)))


if __name__=='__main__':unittest.main()

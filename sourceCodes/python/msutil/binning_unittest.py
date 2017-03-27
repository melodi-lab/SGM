#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Usage: source <base>/env.sh
#        python2.5 <base>/code/msutil/binning_unittest.py

"""
Unit test for binning.py.
"""

__authors__ = [ 'Adam Gustafson <amg81@ee.washington.edu>' ]

import os
import unittest

from msutil.binning import uniform, quantile, simple_uniform, histogram_spectra
from msutil.ms2util import MS2_iterator

class uniformTestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')
        self.ms2_file = os.path.join(self.testdata_dir, 'spectra.ms2')
        self.num_bins = 4

    def tearDown(self):
        pass

    def test_bad_input(self):
        axis = 'mz'
        use_ints = False
        self.assertRaises(ValueError, uniform, self.ms2_file, 'm/z', 
                          self.num_bins, use_ints)  
        self.assertRaises(TypeError, uniform, self.ms2_file, axis, 
                          self.num_bins, 'whoops')
        self.assertRaises(TypeError, uniform, self.ms2_file, axis, 
                          'whoops', use_ints)                          
        self.assertRaises(TypeError, uniform, self.ms2_file, 1, 
                          self.num_bins, use_ints)                          
        self.assertRaises(TypeError, uniform, 1, axis,
                          self.num_bins, use_ints)                
                              
    def test_mz(self):
        axis = 'mz'
        use_ints = True
        tick_points = uniform(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertEqual(tick_points,
                         [(244, 252), (252, 260), (260, 268), (268, 276)])
        use_ints = False
        tick_points = uniform(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertAlmostEqual(tick_points[0][0], 244.0)
        self.assertAlmostEqual(tick_points[0][1], 252.0)
        self.assertAlmostEqual(tick_points[1][0], 252.0)
        self.assertAlmostEqual(tick_points[1][1], 260.0)
        self.assertAlmostEqual(tick_points[2][0], 260.0)
        self.assertAlmostEqual(tick_points[2][1], 268.0)
        self.assertAlmostEqual(tick_points[3][0], 268.0)
        self.assertAlmostEqual(tick_points[3][1], 276.0)
        
    def test_intensity(self):
        axis = 'intensity'
        use_ints = True
        tick_points = uniform(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertEqual(tick_points,
                         [(27, 186), (186, 345), (345, 504), (504, 667)])
        use_ints = False
        tick_points = uniform(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertAlmostEqual(tick_points[0][0], 27.002161026)
        self.assertAlmostEqual(tick_points[0][1], 186.88978910450001)
        self.assertAlmostEqual(tick_points[1][0], 186.88978910450001)
        self.assertAlmostEqual(tick_points[1][1], 346.77741718300001)
        self.assertAlmostEqual(tick_points[2][0], 346.77741718300001)
        self.assertAlmostEqual(tick_points[2][1], 506.66504526150004)
        self.assertAlmostEqual(tick_points[3][0], 506.66504526150004)
        self.assertAlmostEqual(tick_points[3][1], 666.55267334000007)

    def test_intensity_float(self):
        axis = "intensity"
        use_ints = False
        tick_points = uniform(self.ms2_file, axis, self.num_bins, use_ints)
        
class quantileTestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')
        self.ms2_file = os.path.join(self.testdata_dir, 'spectra.ms2')
        self.num_bins = 4

    def tearDown(self):
        pass

    def test_bad_input(self):
        axis = 'mz'
        use_ints = False
        self.assertRaises(ValueError, quantile, self.ms2_file, 'm/z',
                          self.num_bins, use_ints)  
        self.assertRaises(TypeError, quantile, self.ms2_file, axis, 
                          self.num_bins, 'whoops')
        self.assertRaises(TypeError, quantile, self.ms2_file, axis, 
                          'whoops', use_ints)                          
        self.assertRaises(TypeError, quantile, self.ms2_file, 1, 
                          self.num_bins, use_ints)                          
        self.assertRaises(TypeError, quantile, 1, axis,
                          self.num_bins, use_ints) 
                              
    def test_mz(self):
        axis = 'mz'
        use_ints = True
        tick_points = quantile(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertEqual(tick_points,
                         [(244,257), (257, 266), (266, 271), (271, 276)])
        use_ints = False
        tick_points = quantile(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertAlmostEqual(tick_points[0][0], 244.0)
        self.assertAlmostEqual(tick_points[0][1], 256.8)
        self.assertAlmostEqual(tick_points[1][0], 256.8)
        self.assertAlmostEqual(tick_points[1][1], 266.3)         
        self.assertAlmostEqual(tick_points[2][0], 266.3)
        self.assertAlmostEqual(tick_points[2][1], 271.3)
        self.assertAlmostEqual(tick_points[3][0], 271.3)
        self.assertAlmostEqual(tick_points[3][1], 276.0)

    def test_intensity(self):
        axis = 'intensity'
        use_ints = True
        tick_points = quantile(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertEqual(tick_points, [(27, 34), (34, 45), (45, 92), (92, 667)])
        use_ints = False
        tick_points = quantile(self.ms2_file, axis, self.num_bins, use_ints)
        self.assertAlmostEqual(tick_points[0][0], 27.002161026)
        self.assertAlmostEqual(tick_points[0][1], 34.2786712646)
        self.assertAlmostEqual(tick_points[1][0], 34.2786712646)
        self.assertAlmostEqual(tick_points[1][1], 44.694499969500001)
        self.assertAlmostEqual(tick_points[2][0], 44.694499969500001)
        self.assertAlmostEqual(tick_points[2][1], 92.4735870361)
        self.assertAlmostEqual(tick_points[3][0], 92.4735870361)
        self.assertAlmostEqual(tick_points[3][1], 666.55267333999996)       
    
    
class simple_uniformTestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')
        self.ms2_file = os.path.join(self.testdata_dir, 'spectra.ms2')
        self.num_bins = 4

    def tearDown(self):
        pass

    def test_bad_input(self):
        axis = 'm/z'
        use_ints = False
        self.assertRaises(TypeError, simple_uniform, 'whoops', 10, 2)
        self.assertRaises(TypeError, simple_uniform, 0, 'whoops', 2)
        self.assertRaises(TypeError, simple_uniform, 0, 10, 'whoops')
        self.assertRaises(ValueError, simple_uniform, 11, 10, 2)
        self.assertRaises(ValueError, simple_uniform, 0, 0, 0)

    def test_result(self):
        tick_points = simple_uniform(0,10,2)
        self.assertEqual(tick_points[0][0], 0)
        self.assertEqual(tick_points[0][1], 5)
        self.assertEqual(tick_points[1][0], 5)
        self.assertEqual(tick_points[1][1], 10)

    def test_too_many_bins(self):
        self.assertRaises(AssertionError, simple_uniform, 0, 100, 101)
        self.assertRaises(AssertionError, simple_uniform, 0, 100, 1000)
        self.assertRaises(AssertionError, simple_uniform, 0, 100, 250)

    def test_tight(self):
        ranges = simple_uniform(0, 100, 100)
        for r in ranges:
            self.assertEqual(r[0]+1, r[1])

    def test_histogram_spectra_coarse(self):
        spectra = list(MS2_iterator(self.ms2_file))
        spectrum = spectra[0]

        ranges = simple_uniform(0, 300, 2)
        bins = histogram_spectra(spectrum, ranges, max, use_mz = True)

        self.assertTrue(min(spectrum.mz) > 150)
        self.assertEqual(ranges[0], (0, 150))
        self.assertEqual(bins[0], 0)

        self.assertEqual(ranges[1], (150, 300))
        self.assertEqual(bins[1], max(spectrum.intensity))

    def test_histogram_spectra_fine(self):
        spectra = list(MS2_iterator(self.ms2_file))
        spectrum = spectra[0]

        ranges = simple_uniform(0, 300, 300)
        bins = histogram_spectra(spectrum, ranges, max, use_mz = True)

        mz = [ 244, 245, 248, 255, 256, 259, 262, 264, 265, 266, 268, 269,
               270, 271, 272, 274, 275 ]
        height = [ 39.5469398499, 117.0206604, 666.55267334,
                   126.886894226, 27.002161026, 34.2933921814,
                   32.0194511414, 76.8686294556, 39.6276931763,
                   130.83052063, 55.613243103, 92.4735870361,
                   34.2786712646, 53.1507644653, 123.959289551,
                   44.6944999695, 28.8996353149 ]

        for x, y in zip(mz, height):
            self.assertAlmostEqual(bins[x], y)

if __name__ == '__main__':
    unittest.main()

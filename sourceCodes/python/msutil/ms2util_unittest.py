#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Usage: source <base>/env.sh
#        python2.5 <base>/code/msutil/ms2util_unittest.py

"""Unit test for ms2util.py.
"""

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import copy
import os
import string
import unittest

from msutil.ms2util import MS2_iterator

class MS2IteratorTestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')

    def tearDown(self):
        pass

    def validate_spectrum(self, spectrum):
        self.assertEqual(spectrum.spectrum_id, 5453)
        self.assertAlmostEqual(spectrum.precursor_mass, 1718.80296)
        self.assertAlmostEqual(spectrum.retention_time, 43.1522)
        self.assertEqual(spectrum.charge, 2)
        self.assertEqual(len(spectrum.mz), len(spectrum.intensity))
        self.assertEqual(len(spectrum.mz), 17)
        self.assertAlmostEqual(spectrum.mz[0], 244.0)
        self.assertAlmostEqual(spectrum.mz[4], 256.8)
        self.assertAlmostEqual(spectrum.mz[10], 268.2)
        self.assertAlmostEqual(spectrum.mz[16], 275.5)
        self.assertAlmostEqual(spectrum.m[0], 244.0*2)
        self.assertAlmostEqual(spectrum.m[4], 256.8*2)
        self.assertAlmostEqual(spectrum.m[10], 268.2*2)
        self.assertAlmostEqual(spectrum.m[16], 275.5*2)
        self.assertAlmostEqual(spectrum.intensity[0], 39.5469398499)
        self.assertAlmostEqual(spectrum.intensity[4], 27.002161026)
        self.assertAlmostEqual(spectrum.intensity[10], 55.613243103)
        self.assertAlmostEqual(spectrum.intensity[16], 28.8996353149)

    def find_in_path(self, progname):
        paths = string.split(os.environ.get('PATH', ''), os.pathsep)
        locs = [ p for p in paths if os.path.exists(os.path.join(p, progname)) ]
        return len(locs) > 0

    def test_parse_text_spectrum(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        num_spectra = 0
        for spectrum in MS2_iterator(testfile):
            self.validate_spectrum(spectrum)
            num_spectra = num_spectra + 1
        self.assertEqual(num_spectra, 1)

    def test_parse_compressed_spectrum(self):
        testfile = os.path.join(self.testdata_dir, 'compressed_spectrum.ms2.gz')
        num_spectra = 0
        for spectrum in MS2_iterator(testfile):
            self.validate_spectrum(spectrum)
            num_spectra = num_spectra + 1
        self.assertEqual(num_spectra, 1)

    def test_parse_compressed_spectrum_using_gzcat(self):
        """Test the use of gzcat and pipes to decompress the gzipped spectra."""
        if self.find_in_path('gzcat'):
            testfile = os.path.join(self.testdata_dir,
                                    'compressed_spectrum.ms2.gz')
            num_spectra = 0
            for spectrum in MS2_iterator(testfile, True):
                self.validate_spectrum(spectrum)
                num_spectra = num_spectra + 1
            self.assertEqual(num_spectra, 1)

    def test_parse_pickled_spectrum(self):
        testfile = os.path.join(self.testdata_dir, 'spectra.pickle.gz')
        spectrum_id = [5453, 7811, 32]
        charge = [2, 3, 2]
        mass = [1718.80296, 718.9, 818.802960]
        spectrum_length = [18, 20, 21]
        num_spectra = 0
        for spectrum in MS2_iterator(testfile):
            spectrum.validate()
            self.assertEqual(spectrum.spectrum_id, spectrum_id[num_spectra])
            self.assertAlmostEqual(spectrum.precursor_mass, mass[num_spectra])
            self.assertEqual(spectrum.charge, charge[num_spectra])
            self.assertEqual(len(spectrum.mz), spectrum_length[num_spectra])
            num_spectra = num_spectra + 1
        self.assertEqual(num_spectra, 3)

    def test_iterate_over_spectra(self):
        testfile = os.path.join(self.testdata_dir, 'spectra.ms2')
        num_spectra = 0
        spectrum_id = [5453, 7811, 32]
        charge = [2, 3, 2]
        mass = [1718.80296, 718.9, 818.802960]
        spectrum_length = [18, 20, 21]
        for spectrum in MS2_iterator(testfile):
            spectrum.validate()
            self.assertEqual(spectrum.spectrum_id, spectrum_id[num_spectra])
            self.assertAlmostEqual(spectrum.precursor_mass, mass[num_spectra])
            self.assertEqual(spectrum.charge, charge[num_spectra])
            self.assertEqual(len(spectrum.mz), spectrum_length[num_spectra])
            num_spectra = num_spectra + 1
        self.assertEqual(num_spectra, 3)

    def test_mz_field(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2_iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]

        self.assertAlmostEqual(244.0, s.mz[0])
        self.assertAlmostEqual(39.5469398499, s.intensity[0])
        self.assertAlmostEqual(275.5, s.mz[-1])
        self.assertAlmostEqual(28.8996353149, s.intensity[-1])

        testfile = os.path.join(self.testdata_dir, 'spectra.ms2')
        spectra = list(MS2_iterator(testfile))
        self.assertEqual(len(spectra), 3)

        self.assertAlmostEqual(244.0, spectra[0].mz[0])
        self.assertAlmostEqual(39.5469398499, spectra[0].intensity[0])
        self.assertAlmostEqual(275.5, spectra[0].mz[-1])
        self.assertAlmostEqual(28.8996353149, spectra[0].intensity[-1])

        self.assertAlmostEqual(244.0, spectra[2].mz[0])
        self.assertAlmostEqual(39.5469398499, spectra[2].intensity[0])
        self.assertAlmostEqual(276.0, spectra[2].mz[-1])
        self.assertAlmostEqual(28.8996353149, spectra[2].intensity[-1])

    def test_log_normalization(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2_iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]
        s.log_normalize()
        for v in s.intensity:
            self.assertTrue(v >= 0.0)

    def test_rank_normalization(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2_iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]
        r = copy.deepcopy(s)

        s.rank_normalize()
        s.validate()
        self.assertAlmostEqual(max(s.intensity), 1.0)
        self.assertAlmostEqual(min(s.intensity), 1.0/float(len(s.intensity)))
        self.assertEqual(len(s.mz), len(r.mz))
        self.assertEqual(len(s.m), len(r.m))
        self.assertEqual(len(s.intensity), len(r.intensity))

        for a, b in zip(s.mz, r.mz):
            self.assertAlmostEqual(a, b)

        for a, b in zip(s.m, r.m):
            self.assertAlmostEqual(a, b)

    def test_remove_low_peaks(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2_iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]
        r = copy.deepcopy(s)

        frac = 0.10
        intensities = copy.deepcopy(s.intensity)
        intensities.sort()
        s.remove_low_peaks(frac)
        s.validate()
        self.assertTrue(len(s.intensity) < len(r.intensity))

        # Check that we removed intensities that are smaller than anything
        # remaining in the transformed spectrum.
        for h in intensities[0:int(frac*len(intensities))]:
            self.assertTrue(h < min(s.intensity))

    def test_hash_cmp_len(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2_iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s1 = spectra[0]
        s2 = copy.deepcopy(s1)
        s2.spectrum_id = s2.spectrum_id + 1 # artificially change sid

        # length(spectrum)
        self.assertEqual(s1.length, len(s1.mz))

        # __hash__
        dic = { }
        dic[s1] = 1
        dic[s2] = 2
        self.assertEqual(len(dic), 2)

        # __cmp__
        self.assertTrue(s1 < s2)
        self.assertFalse(s1 > s2)
        self.assertFalse(s1 == s2)
        self.assertTrue(s1 == s1)
        self.assertTrue(s2 == s2)
        

if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python

import os
import copy
import unittest

from msutil.spectrum import MS2Iterator

class MS2TestCase(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')

    def tearDown(self):
        pass

    def validate_spectrum(self, spectrum):
        self.assertEqual(spectrum.spectrum_id, 5453)
        self.assertEqual(len(spectrum.mz), len(spectrum.intensity))
        self.assertEqual(len(spectrum.mz), 17)
        self.assertAlmostEqual(spectrum.mz[0], 244.0)
        self.assertAlmostEqual(spectrum.mz[4], 256.8)
        self.assertAlmostEqual(spectrum.mz[10], 268.2)
        self.assertAlmostEqual(spectrum.mz[16], 275.5)
        self.assertAlmostEqual(spectrum.intensity[0], 39.5469398499)
        self.assertAlmostEqual(spectrum.intensity[4], 27.002161026)
        self.assertAlmostEqual(spectrum.intensity[10], 55.613243103)
        self.assertAlmostEqual(spectrum.intensity[16], 28.8996353149)

    def test_parse_text_spectrum(self):
        testfile = os.path.join(self.testdata_dir, 'multiz.ms2')
        spectra = list(MS2Iterator(testfile))
        self.assertEqual(len(spectra), 2)
        self.assertEqual(spectra[0].spectrum_id, 5453)
        self.assertEqual(spectra[1].spectrum_id, 9453)
        self.assertEqual(spectra[0].precursor_mz, 859.90545)
        self.assertEqual(spectra[1].precursor_mz, 859.90545)
        self.assertEqual(spectra[0].charge_list, [(2,1718.80296),(3,718.80296)])
        self.assertEqual(spectra[1].charge_list, [(2,1718.80296)])
        for s in spectra:
            self.assertEqual(len(s.mz), len(s.intensity))
            self.assertEqual(min(s.mz), 244.0)
            self.assertEqual(max(s.mz), 275.5)
            self.assertEqual(min(s.intensity), 27.002161026)
            self.assertEqual(max(s.intensity), 666.55267334)

    def test_parse_compressed_spectrum(self):
        testfile = os.path.join(self.testdata_dir, 'compressed_spectrum.ms2.gz')
        spectra = list(MS2Iterator(testfile))
        self.validate_spectrum(spectra[0])
        self.assertEqual(len(spectra), 1)

        spectra = list(MS2Iterator(testfile, True))
        self.validate_spectrum(spectra[0])
        self.assertEqual(len(spectra), 1)

        #spectra = list(MS2Iterator('../../data/strict-orbitrap.ms2.gz'))
        #print len(spectra)

    #def test_big_file(self):
        #spectra = list(MS2Iterator('/g/ssli/projects/msms/data/spectra/yeast-02.ms2'))
        #spectra = list(MS2Iterator('/s0/ajit/spectra/yeast-01.ms2'))
        #spectra = list(MS2Iterator(os.path.join(self.testdata_dir, 'debug.ms2')))

    def test_mz_field(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2Iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]

        self.assertAlmostEqual(244.0, s.mz[0])
        self.assertAlmostEqual(39.5469398499, s.intensity[0])
        self.assertAlmostEqual(275.5, s.mz[-1])
        self.assertAlmostEqual(28.8996353149, s.intensity[-1])

        testfile = os.path.join(self.testdata_dir, 'spectra.ms2')
        spectra = list(MS2Iterator(testfile))
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
        spectra = list(MS2Iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]
        s.log_normalize()
        for v in s.intensity:
            self.assertTrue(v >= 0.0)

    def test_rank_normalization(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2Iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]
        r = copy.deepcopy(s)

        s.rank_normalize()
        self.assertAlmostEqual(max(s.intensity), 1.0)
        self.assertAlmostEqual(min(s.intensity), 1.0/float(len(s.intensity)))
        self.assertEqual(len(s.mz), len(r.mz))
        self.assertEqual(len(s.intensity), len(r.intensity))

        for a, b in zip(s.mz, r.mz):
            self.assertAlmostEqual(a, b)

    def test_remove_low_peaks(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2Iterator(testfile))
        self.assertEqual(len(spectra), 1)
        s = spectra[0]
        r = copy.deepcopy(s)

        frac = 0.10
        intensities = copy.deepcopy(s.intensity)
        intensities.sort()
        s.remove_low_peaks(frac)
        self.assertTrue(len(s.intensity) < len(r.intensity))

        # Check that we removed intensities that are smaller than anything
        # remaining in the transformed spectrum.
        for h in intensities[0:int(frac*len(intensities))]:
            self.assertTrue(h < min(s.intensity))

    def test_hash_cmp_len(self):
        testfile = os.path.join(self.testdata_dir, 'spectrum.ms2')
        spectra = list(MS2Iterator(testfile))
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

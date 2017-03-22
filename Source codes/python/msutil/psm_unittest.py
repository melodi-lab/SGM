#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Usage: source <base>/env.sh
#        python2.5 <base>/code/msutil/psm_unittest.py

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os
import unittest

from msutil.psm import PeptideSpectrumMatches
from protein.peptide import Peptide

class PeptideSpectrumMatchesTestCase(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')

    def tearDown(self):
        pass

    def data(self, filename):
        return os.path.join(self.testdata_dir, filename)

    def test_init(self):
        m = PeptideSpectrumMatches(self.data('test.psm'))

    def test_init_failure(self):
        m = PeptideSpectrumMatches( )
        self.assertRaises(RuntimeError,
                          PeptideSpectrumMatches.parse_psm_file,
                          m, self.data('bad_id.psm'), False)

    def test_duplicate_psm(self):
        m = PeptideSpectrumMatches( )
        self.assertRaises(KeyError,
                          PeptideSpectrumMatches.parse_psm_file,
                          m, self.data('dup_spectra.psm'), False)
        self.assertRaises(KeyError,
                          PeptideSpectrumMatches.parse_psm_file,
                          m, self.data('dup_peptide.psm'), False)

    def test_set_get(self):
        """Test the bijection invariant between peptide and spectrum ids."""
        m = PeptideSpectrumMatches( )
        p = Peptide('LANK')
        s = 5463
        m[s] = p
        self.assertRaises(KeyError,
                          PeptideSpectrumMatches.__setitem__, m, p, s)

        self.assertEquals(len(m), 1)
        self.assertEquals(m[p], s)
        self.assertEquals(m[s], p)

        p2 = Peptide('LANKR')
        s2 = s + 1
        m[p2] = s2
        self.assertEquals(len(m), 2)
        self.assertEquals(m[p2], s2)
        self.assertEquals(m[s2], p2)

        self.assertRaises(KeyError,
                          PeptideSpectrumMatches.__setitem__,
                          m, p2, s)

    def test_iteration(self):
        # Test __iter__ routine.
        m = PeptideSpectrumMatches(self.data('test.psm'))
        c = 0
        for s, p in m:
            c = c + 1
        self.assertEquals(c, 10)

    def test_peptides(self):
        m = PeptideSpectrumMatches(self.data('test.psm'))
        c = 0
        for p in m.peptides():
            c = c + 1
        self.assertEqual(c, 10)

if __name__ == '__main__':
    unittest.main()

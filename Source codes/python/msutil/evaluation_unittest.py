#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Usage: source <base>/env.sh
#        python2.5 <base>/code/msutil/evaluation_unittest.py

"""
Unit test for evaluation.py.
"""

__authors__ = [ 'Adam Gustafson <amg81@ee.washington.edu>' ]

import os
import unittest
from protein.peptide import Peptide
from msutil.evaluation import rel_rank, abs_rank
from msutil.psm import PeptideSpectrumMatches


class rel_rankTestCase(unittest.TestCase):
    def setUp(self):
        # make test true database
        self.pep_list = \
            [Peptide("arkvc"), Peptide("peptide"), Peptide("hahargcvkr")]
        self.pep_scores = [-600, -400, -250]
        # pick the true peptide used as evidence
        self.true_pep = self.pep_list[0]
        # make test decoy database
        self.d_pep_list = []
        for pep in self.pep_list:
            self.d_pep_list.append(Peptide.shuffle(pep,False))        
        self.d_pep_scores = [-601, -401, -251]
        self.K = 3        
        self.norm = True

        self.test_dir = os.path.dirname(os.path.realpath(__file__))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')

    def tearDown(self):
        pass

    def test_bad_input(self):
        for K in [-1, 0, 4, 5, 6]:    
            self.assertRaises(ValueError, rel_rank, 
                              self.pep_list, self.pep_scores,
                              self.d_pep_list, self.d_pep_scores,
                              self.true_pep, K, self.norm)
        self.assertRaises(TypeError, rel_rank,              # bad pep_list
                          'whoops', self.pep_scores,  
                          self.d_pep_list, self.d_pep_scores,
                          self.true_pep, self.K, self.norm)
        self.assertRaises(TypeError, rel_rank,              # bad pep_scores
                          self.pep_list, 'whoops',
                          self.d_pep_list, self.d_pep_scores,
                          self.true_pep, self.K, self.norm)
        self.assertRaises(TypeError, rel_rank,              # bad d_pep_list
                          self.pep_list, self.pep_scores,   
                          'whoops', self.d_pep_scores,
                          self.true_pep, self.K, self.norm)
        self.assertRaises(TypeError, rel_rank,              # bad d_pep_scores
                          self.pep_list, self.pep_scores,   
                          self.d_pep_list, 'whoops',
                          self.true_pep, self.K, self.norm)
        self.assertRaises(TypeError, rel_rank,              # bad true_pep
                          self.pep_list, self.pep_scores,   
                          self.d_pep_list, self.d_pep_scores,
                          'whoops', self.K)                          
        self.assertRaises(TypeError, rel_rank,              # bad K
                          self.pep_list, self.pep_scores,   
                          self.d_pep_list, self.d_pep_scores,
                          self.true_pep, 'whoops', self.norm)
        self.assertRaises(TypeError, rel_rank,              # bad norm
                          self.pep_list, self.pep_scores,   
                          self.d_pep_list, self.d_pep_scores,
                          self.true_pep, K, 'whoops')
  
    def test_results(self):
        performance, true_pep_pos = rel_rank(self.pep_list, self.pep_scores,
                                            self.d_pep_list, self.d_pep_scores,
                                            self.true_pep, self.K, self.norm)
        self.assertEqual(performance, float(2)/3)
        self.assertEqual(true_pep_pos, 5)


class abs_rankTestCase(unittest.TestCase):
    def setUp(self):
        # form the psms from test.psm (ignoring any duplicates)
        self.psm = PeptideSpectrumMatches()
        self.test_dir = os.path.dirname(os.path.realpath(__file__))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')
        testfile = os.path.join(self.testdata_dir, 'test.psm')
        self.psm.parse_psm_file(testfile, True)
        # form the lists of the spectrum_IDs and dummy inferred spectra
        self.spectrum_IDs = []
        self.guessed_peps = []
        test_file = open(os.path.join(self.testdata_dir, 'test.psm'), 'r')
        for line in test_file:
            tokens = line.split()
            ID = int(tokens[0])
            sequence = tokens[1]
            self.spectrum_IDs.append(ID)
            pep = Peptide(sequence)
            if ID % 2 == 0:  # wrong inference for even-numberd spectrum IDs
                pep = pep.shuffle(False)
            self.guessed_peps.append(pep)

            
    def tearDown(self):
        pass

    def test_bad_input(self):
        dummy_IDs = [5453, 4653, 6127, 13483, 2501, 5273, 6817]
        self.assertRaises(RuntimeError, abs_rank, 
                          dummy_IDs, self.guessed_peps, self.psm)
        self.assertRaises(TypeError, abs_rank, 
                          'whoops', self.guessed_peps, self.psm)
        self.assertRaises(TypeError, abs_rank, 
                          self.spectrum_IDs, 'whoops', self.psm)       
        self.assertRaises(TypeError, abs_rank, 
                          self.spectrum_IDs, self.guessed_peps, 'whoops')
                          
    def test_results(self):
        performance, correct_IDs = \
                abs_rank(self.spectrum_IDs, self.guessed_peps, self.psm)
        self.assertEqual(performance, 0.7)
        self.assertEqual(correct_IDs, 
                        [5453, 4653, 6127, 13483, 2501, 5273, 6817])


if __name__ == '__main__':
    unittest.main()

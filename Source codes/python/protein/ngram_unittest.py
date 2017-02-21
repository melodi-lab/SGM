#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Usage: source <base>/env.sh
#        python2.5 <base>/code/msutil/binning_unittest.py

"""
Unit test for ngram.py.
"""

__authors__ = [ 'Adam Gustafson <amg81@ee.washington.edu>' ]

import os
import unittest

from random import uniform
from protein.ngram import load_peptides, prior, transition_matrix

class priorTestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')
        self.d = 20

    def tearDown(self):
        pass

    def test_prior_uniform(self):
        filename = os.path.join(self.testdata_dir, 'test-prior-uniform.pep')
        pep_list = load_peptides(filename)
        for prob in prior(pep_list): 
            self.assertAlmostEqual(prob, 1./self.d)
    
    def test_prior_uniform_smoothing(self):
        alpha = uniform(0,1)
        N = 20
        filename = os.path.join(self.testdata_dir, 'test-prior-uniform.pep')
        for i, prob in enumerate(prior(load_peptides(filename), alpha)):
                self.assertAlmostEqual(prob, (1.0 + alpha)/(N + self.d*alpha))

    def test_prior_single(self):
        filename = os.path.join(self.testdata_dir, 'test-prior-single.pep')
        for i, prob in enumerate(prior(load_peptides(filename))):
            if i == 0: # 'A'
                self.assertAlmostEqual(prob, 1.0)
            else:
                self.assertAlmostEqual(prob, 0.0)
    
    def test_prior_single_smoothing(self):
        alpha = uniform(0,1)
        N = 1
        filename = os.path.join(self.testdata_dir, 'test-prior-single.pep')
        pep_list = load_peptides(filename)
        for i, prob in enumerate(prior(pep_list, alpha)):
            if i == 0: # 'A'
                self.assertAlmostEqual(prob, (1.0 + alpha)/(N + self.d*alpha))
            else: 
                self.assertAlmostEqual(prob, (0.0 + alpha)/(N + self.d*alpha))


class transition_matrixTestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(os.path.realpath( __file__ ))
        self.testdata_dir = os.path.join(self.test_dir, 'testdata')

    def tearDown(self):
        pass

    def test_ngram_uniform(self):
        f = []
        f.append(os.path.join(self.testdata_dir, 'test-bigram-uniform.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-trigram-uniform.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-quadgram-uniform.pep'))
        for i in range(3):
            n = i + 2          # bigram, trigram, etc.
            d = pow(20, i + 1) # cardinality of conditional distribution
            Q = transition_matrix(load_peptides(f[i]), n)
            for i, p in enumerate(Q): # ith row of transition matrix
                for j, prob in enumerate(p): 
                    self.assertAlmostEqual(prob, 1.0/d)
    
    def test_ngram_uniform_smoothing(self):
        alpha = uniform(0,1)
        f = []
        f.append(os.path.join(self.testdata_dir, 'test-bigram-uniform.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-trigram-uniform.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-quadgram-uniform.pep'))
        for n in range(2,5):    # bigram, trigram, etc. 
            d = pow(20, n - 1)  # cardinality of conditional distribution
            N = d               # instances of transitions to any amino acid
            Q = transition_matrix(load_peptides(f[n-2]), n)
            for i, p in enumerate(Q): # ith row of transition matrix
                for j, prob in enumerate(p): 
                    self.assertAlmostEqual(prob, (1.0 + alpha)/(N + alpha*d))
    
    def test_ngram_single(self):
        f = []
        f.append(os.path.join(self.testdata_dir, 'test-bigram-single.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-trigram-single.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-quadgram-single.pep'))
        for n in range(2,5):    # bigram, trigram, etc.
            Q = transition_matrix(load_peptides(f[n-2]), n)
            for i, p in enumerate(Q): # ith row of transition matrix
                for j, prob in enumerate(p):
                    if j == 0: # transitions from 'A', 'AA', etc.
                        self.assertAlmostEqual(prob, 1.0)
                    else:
                        self.assertAlmostEqual(prob, 0.0)
    
    def test_ngram_single_smoothing(self):
        alpha = uniform(0,1)
        f = []
        f.append(os.path.join(self.testdata_dir, 'test-bigram-single.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-trigram-single.pep'))
        f.append(os.path.join(self.testdata_dir, 'test-quadgram-single.pep'))
        for n in range(2,5):    # bigram, trigram, etc.
            d = pow(20, n - 1)  # cardinality of conditional distribution
            N = 1
            Q = transition_matrix(load_peptides(f[n-2]), n, alpha)
            for i, p in enumerate(Q): # ith row of transition matrix
                for j, prob in enumerate(p):
                    if j == 0: # transitions from 'A', 'AA', etc.
                        self.assertAlmostEqual(prob, (1 + alpha)/(N + d*alpha))
                    else:
                        self.assertAlmostEqual(prob, (0 + alpha)/(N + d*alpha))


if __name__ == '__main__':
    unittest.main()
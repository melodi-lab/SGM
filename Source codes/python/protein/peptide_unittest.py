#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Usage: source <base>/env.sh
#        python2.5 <base>/code/protein/peptide_unittest.py

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import unittest
import pickle
import tempfile

from protein import peptide
from protein.peptide import Peptide

class PeptideTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        p = Peptide()
        self.assertEqual(len(p.seq), 0)
        p = Peptide('ranlrag')
        self.assertEqual(str(p), 'RANLRAG')
        self.assertRaises(ValueError, Peptide, 'ZAN')
        self.assertRaises(ValueError, Peptide, '!an')

    def test_length(self):
        p = Peptide('RANQG')
        self.assertEqual(p.length, 5)

    def test_printer(self):
        seq = 'RANQG'
        p = Peptide(seq)
        self.assertEqual(str(p), seq)

    def test_copy(self):
        p = Peptide('RANQGNR')
        p = Peptide()
        self.assertEqual(p.length, 0)

    def test_shuffle(self):
        """Test Peptide.shuffle( ), assuming non-tryptic shuffle."""
        seq = 'ANQGNR'
        shuffled_peptides = [ ]
        for _ in range(100):
            shuffled_peptides.append(Peptide(seq).shuffle())
        is_different = [str(p) != seq for p in shuffled_peptides]
        self.assertTrue(any(is_different))

    def test_shuffle_tryptic(self):
        """Test Peptide.shuffle(True), a tryptic preserving shuffle."""
        seq = 'ANQGNR'
        for _ in range(100):
            peptide = Peptide(seq).shuffle(True)
            self.assertEqual(str(peptide)[-1], 'R')

    def test_reverse(self):
        self.assertTrue(str(Peptide('').reverse()) == '')
        self.assertTrue(str(Peptide('').reverse(False)) == '')

        self.assertTrue(str(Peptide('A').reverse()) == 'A')
        self.assertTrue(str(Peptide('A').reverse(False)) == 'A')

        self.assertTrue(str(Peptide('AR').reverse()) == 'AR')
        self.assertTrue(str(Peptide('AR').reverse(False)) == 'RA')

        self.assertTrue(str(Peptide('ANQGNR').reverse()) == 'NGQNAR')
        self.assertTrue(str(Peptide('ANQGNR').reverse(False)) == 'RNGQNA')

    def test_random(self):
        self.assertTrue(str(Peptide('').random()) == '')
        self.assertTrue(Peptide('A').random().length == 1)

        p = Peptide('ANQGNR')
        shuffled_peptides = [ ]
        for _ in xrange(100):
            shuffled_peptides.append(p.random())
        is_different = [ str(q) != p.seq for q in shuffled_peptides ]
        self.assertTrue(any(is_different))

    def test_hashability(self):
        """Test that Peptides hash by their IUPAC string representation;
        not their address in memory, which is the default for objects.
        """
        seqs = ['ANE', 'QGH', 'EGH', 'LKM']
        peps = [Peptide(s) for s in seqs]
        # __cmp__ is required for hashability
        for s1,p1 in zip(seqs,peps):
            for s2,p2 in zip(seqs,peps):
                self.assertEqual(s1 < s2, p1 < p2)
                self.assertEqual(s1 == s2, p1 == p2)
                self.assertEqual(s1 > s2, p1 > p2)
        hash_map = { }
        for p in peps:
            hash_map[p] = None
        for s in seqs:
            hash_map[Peptide(s)] = s
        self.assertEqual(len(hash_map), len(peps))
        for p, s in list(hash_map.items()):
            self.assertEqual(p.seq, s)

    def test_mass_table(self):
        """Test that the output of peptide.mass_table() is safely mutable."""
        self.assertRaises(ValueError, peptide.mass_table, 'askdaksd')
        self.assertRaises(ValueError, peptide.mass_table, 'mono')
        self.assertRaises(ValueError, peptide.mass_table, 'avg')
        mono = peptide.mass_table('monoisotopic')
        avg = peptide.mass_table('average')
        self.assertEqual(len(mono), 20)
        self.assertEqual(len(avg), 20)
        mono = { }
        avg = { }
        self.assertEqual(len(peptide.mass_table('monoisotopic')), 20)
        self.assertEqual(len(peptide.mass_table('average')), 20)

    def test_iupac_table(self):
        """Test that the IUPAC naming tables were typed in correctly.

        The test is a partial (weak) test, probing only a few entries
        in the naming tables.
        """
        pairs = {'M' : 'Met', 'P' : 'Pro', 'G' : 'Gly'}
        one_to_three = peptide.IUPAC_table()
        three_to_one = peptide.IUPAC_table('3-to-1')
        self.assertRaises(ValueError, peptide.IUPAC_table, '1to3')
        self.assertRaises(ValueError, peptide.IUPAC_table, '3to1')
        for k, v in list(pairs.items()):
            self.assertEqual(one_to_three[k], v)
            self.assertEqual(three_to_one[v], k)

    def test_amino_acids(self):
        """Test that the amino acid list has been computed."""
        self.assertEqual(len(peptide.amino_acids()), 20)

    def test_mass_dictionaries(self):
        """The monoisotopic and average mass dictionaries should be the same
        length, and have almost the same masses for each amino acid. The average
        and monoisotopic dictionaries should be interchangeable in code.
        """
        mono = peptide.mass_table('monoisotopic')
        avg = peptide.mass_table('average')
        self.assertEqual(len(mono), len(avg))
        for (key_mono, key_avg) in zip(mono,avg):
            self.assertEqual(key_mono, key_avg)
        values = list(zip(iter(list(mono.values())),
                          iter(list(avg.values()))))
        for (a, b) in values:
            self.assertTrue(abs(a - b) < 1) # 1 Dalton

    def test_mass(self):
        """Test computed properties, Peptide.mass and Peptide.average_mass."""
        sequences = ['G', 'GAN', 'WYV']
        true_mass = [57.02146, 242.1015, 448.21105]
        true_avg_mass = [57.0519, 242.2345, 448.5218]

        peps = [Peptide(s) for s in sequences]
        mass = [p.mass for p in peps]
        avg_mass = [p.average_mass for p in peps]
        for (m, t) in zip(mass, true_mass):
            self.assertAlmostEqual(m, t)
        for (m, t) in zip(avg_mass, true_avg_mass):
            self.assertAlmostEqual(m, t)

    def test_serialize(self):
        """Test that Peptide objects are serializable."""
        peptide = Peptide('GANWYVRKA')
        temp = tempfile.NamedTemporaryFile()
        pickle.dump(peptide, temp)
        temp.seek(0)
        new_peptide = pickle.load(temp)
        self.assertEqual(peptide, new_peptide)
        temp.close()

    def test_ideal_fragment_masses(self):
        """Test the synthetic fragmentation routine."""
        p = Peptide('')
        (left, right) = p.ideal_fragment_masses( )
        self.assertTrue(1 == len(left) and 1 == len(right))
        self.assertEqual(left[0], 0)
        self.assertEqual(right[0], 0)

        p = Peptide('A')
        (left, right) = p.ideal_fragment_masses( )
        self.assertTrue(2 == len(left) and 2 == len(right))
        self.assertEqual(left[0], 0)
        self.assertEqual(right[0], p.mass)
        self.assertEqual(left[1], p.mass)
        self.assertEqual(right[1], 0)

        (left, right) = p.ideal_fragment_masses('average')
        self.assertEqual(left[0], 0)
        self.assertEqual(right[0], p.average_mass)

        p = Peptide('GAN')
        s = p.seq
        (left, right) = p.ideal_fragment_masses( )
        self.assertTrue(4 == len(left) and 4 == len(right))
        for i in range(len(s)+1):
            self.assertEqual(s[0:i] + s[i:3], s)
            self.assertAlmostEqual(left[i], Peptide(s[0:i]).mass)
            self.assertAlmostEqual(right[i], Peptide(s[i:3]).mass)
            self.assertAlmostEqual(left[i] + right[i], p.mass)

        p = Peptide('GAN')
        s = p.seq
        (left, right) = p.ideal_fragment_masses('average')
        self.assertTrue(4 == len(left) and 4 == len(right))
        for i in range(len(s)+1):
            self.assertEqual(s[0:i] + s[i:3], s)
            self.assertAlmostEqual(left[i], Peptide(s[0:i]).average_mass)
            self.assertAlmostEqual(right[i], Peptide(s[i:3]).average_mass)
            self.assertAlmostEqual(left[i] + right[i], p.average_mass)

    # TODO(ajit): Add unit-test for mass_op mode.

if __name__ == '__main__':
    unittest.main()


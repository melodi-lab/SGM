#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

"""Peptides and operations on peptides.

A peptide is a short sequence of amino acids, usually a piece of a protein.
Amino acid tables are taken from [1].

[1] http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
"""

import random
import string

# Monoisotopic: elements are assumed to be the most common isotopes.
__amino_acid_mass = { 'A' : 71.03711,
                      'R' : 156.10111,
                      'N' : 114.04293,
                      'D' : 115.02694,
                      'C' : 103.00919 + 57.021464,
                      'E' : 129.04259,
                      'Q' : 128.05858,
                      'G' : 57.02146,
                      'H' : 137.05891,
                      'I' : 113.08406,
                      'L' : 113.08406,
                      'K' : 128.09496,
                      'M' : 131.04049,
                      'F' : 147.06841,
                      'P' : 97.05276,
                      'S' : 87.03203,
                      'T' : 101.04768,
                      'W' : 186.07931,
                      'Y' : 163.06333,
                      'V' : 99.06841 }

__amino_acid_mass_pure = { 'A' : 71.03711,
                           'R' : 156.10111,
                           'N' : 114.04293,
                           'D' : 115.02694,
                           'C' : 103.00919,
                           'E' : 129.04259,
                           'Q' : 128.05858,
                           'G' : 57.02146,
                           'H' : 137.05891,
                           'I' : 113.08406,
                           'L' : 113.08406,
                           'K' : 128.09496,
                           'M' : 131.04049,
                           'F' : 147.06841,
                           'P' : 97.05276,
                           'S' : 87.03203,
                           'T' : 101.04768,
                           'W' : 186.07931,
                           'Y' : 163.06333,
                           'V' : 99.06841 }

# Average: weighted average of isotopic masses, weighted by isotope abundance.
__amino_acid_average_mass = { 'A' : 71.0788,
                              'R' : 156.1875,
                              'N' : 114.1038,
                              'D' : 115.0886,
                              'C' : 103.1388 + 57.021464,
                              'E' : 129.1155,
                              'Q' : 128.1307,
                              'G' : 57.0519,
                              'H' : 137.1411,
                              'I' : 113.1594,
                              'L' : 113.1594,
                              'K' : 128.1741,
                              'M' : 131.1926,
                              'F' : 147.1766,
                              'P' : 97.1167,
                              'S' : 87.0782,
                              'T' : 101.1051,
                              'W' : 186.2132,
                              'Y' : 163.1760,
                              'V' : 99.1326 }

__amino_acid_average_mass_pure = { 'A' : 71.0788,
                                   'R' : 156.1875,
                                   'N' : 114.1038,
                                   'D' : 115.0886,
                                   'C' : 103.1388,
                                   'E' : 129.1155,
                                   'Q' : 128.1307,
                                   'G' : 57.0519,
                                   'H' : 137.1411,
                                   'I' : 113.1594,
                                   'L' : 113.1594,
                                   'K' : 128.1741,
                                   'M' : 131.1926,
                                   'F' : 147.1766,
                                   'P' : 97.1167,
                                   'S' : 87.0782,
                                   'T' : 101.1051,
                                   'W' : 186.2132,
                                   'Y' : 163.1760,
                                   'V' : 99.1326 }

__amino_acid_letter_codes = {'A' : 'Ala',
                             'R' : 'Arg',
                             'N' : 'Asn',
                             'D' : 'Asp',
                             'C' : 'Cys',
                             'E' : 'Glu',
                             'Q' : 'Gln',
                             'G' : 'Gly',
                             'H' : 'His',
                             'I' : 'Ile',
                             'L' : 'Leu',
                             'K' : 'Lys',
                             'M' : 'Met',
                             'F' : 'Phe',
                             'P' : 'Pro',
                             'S' : 'Ser',
                             'T' : 'Thr',
                             'W' : 'Trp',
                             'Y' : 'Tyr',
                             'V' : 'Val' }

__amino_acids = __amino_acid_mass.keys()
__amino_acids.sort()

def mass_table(tbl = 'monoisotopic', pure = False):
    """Return a dictionary of IUPAC amino acid codes to masses (in Daltons).
    """
    if tbl == 'monoisotopic':
        if pure: return __amino_acid_mass_pure
        else: return __amino_acid_mass
    elif tbl == 'average':
        if pure: return __amino_acid_average_mass_pure
        else: return __amino_acid_average_mass
    else:
        raise ValueError("Argument must be one of ['monoisotopic', 'average'].")

def IUPAC_table(tbl = '1-to-3'):
    """Return a dictionary of IUPAC 1-letter to 3-letter codes, or IUPAC
    3-letter to 1-letter codes: e.g., 'M' -> 'Met' or 'Met' -> 'M'
    """
    if tbl == '1-to-3':
        return __amino_acid_letter_codes
    elif tbl == '3-to-1':
        return dict([(v, k) for k,v in list(__amino_acid_letter_codes.items())])
    else:
        raise ValueError("Argument must be one of ['1-to-3', '3-to-1'].")

def amino_acids():
    return __amino_acids

def amino_acids_to_indices():
    """Return a dictionary that maps amino acid letter codes to its index.

    The returned dictionary can be used to map letter like 'A' to their integer
    index, which is extensively used in CPTs which condition on amino acid
    identity.

    """
    dic = { }
    for index, aa in enumerate(amino_acids()):
        dic[aa] = index
    return dic

def transition_probabilities(peptides, dirichlet_scale = 1e-10):
    """Compute the probability of transitions between amino acids.

    Args:
        peptides: Sequence of peptides.
        dirichlet_scale: Strength of uniform Dirichlet smoothing.

    Returns:
        A square matrix of size n, where n is the number of amino
        acids specified by protein.peptide.amino_acids(). The i^th
        row of the matrix represents the multinomial distribution
        P(AA_{t+1} | AA_{t} = i), where AA is a discrete variable
        representing an amino acid.

        The MLE is smoothed.

    """
    n = len(amino_acids())
    count = [ [ 0 for col in range(n) ] for row in range(n)]

    def pairs(sequence):
        for i in xrange(0, len(sequence)-1):
            yield (sequence[i], sequence[i+1])

    indexer = amino_acids_to_indices()
    for p in peptides:
        for a, b in pairs(p.seq):
            i = indexer[a]
            j = indexer[b]
            count[i][j] = count[i][j] + 1

    for i in range(n):
        total = float(sum(count[i]))
        count[i] = [ v/total for v in count[i] ]
        count[i] = [ v + dirichlet_scale for v in count[i] ] #smoothing
    return count


class Peptide(object):
    """Short sequence of amino acids.

    Peptides represent a sequence of amino acids, a string where each
    letter corresponds to the IUPAC code for an amino acid (e.g., A = Arginine).
    This class encapsulates the string sequence, and common operations on
    peptides. Peptides hash by their amino acid sequence; not their object id
    in memory.
    """

    def __init__(self, sequence = ''):
        """Test validity of peptide sequence string."""
        self.seq = sequence.upper() # String representation: e.g., 'KVRN'
        if not set(self.seq) <= set(amino_acids()):
            diff = list(set(self.seq) - set(amino_acids()))
            raise ValueError("Argument 'sequence' contains letters not in "
                             "IUPAC amino acid tables: " + str(diff))

    def __cmp__(self,other):
        if self.seq < other.seq:
            return -1
        elif self.seq > other.seq:
            return 1
        else:
            return 0

    def __hash__(self):
        return hash(self.seq)

    def __str__(self):
        return self.seq

    @property
    def length(self):
        """Number of amino acids in the peptide. Read-only computed property."""
        return len(self.seq)

    @property
    def mass(self):
        """Mass of the peptide. Read-only computed property."""
        return self.__compute_mass(mass_table('monoisotopic', pure = False)) + 18.010564684

    @property
    def average_mass(self):
        """Isotope weighted mass of the peptide. Read-only computed property."""
        return self.__compute_mass(mass_table('average', pure = False)) + 18.0153

    @property
    def pure_average_mass(self):
        """Isotope weighted mass of the peptide. Read-only computed property."""
        return self.__compute_mass(mass_table('average', pure = True)) + 18.0153

    def __compute_mass(self, mass):
        """Compute the mass of a peptide using a selected mass table."""
        #return sum(map(lambda aa: mass[aa], list(self.seq)))
        return sum([mass[aa] for aa in list(self.seq)])

    def shuffle(self, tryptic = True):
        """Return a shuffled copy of the peptide sequence, as a Peptide object.

        Randomly permutes the amino acids. If tryptic is set to true, then
        we guarantee that the last position (which is usually
        'R' or 'K' in a trypic peptide) is not shuffled.

        """
        residue_list = list(self.seq)
        if tryptic:
            last_char = residue_list[-1]
            residue_list = residue_list[:-1]
            random.shuffle(residue_list)
            residue_list.append(last_char)
        else:
            random.shuffle(residue_list)
        return Peptide(string.join(residue_list, sep = ''))

    def reverse(self, tryptic = True):
        """Return a reverse copy of the peptide sequence, as a Peptide object.

        If tryptic is set to true, then we guarantee that the last position,
        which is 'R' or 'K' in a tryptic peptide, is not reversed.

        """
        if len(self.seq) == 0:
            return Peptide(self.seq)

        residue_list = list(self.seq)
        if tryptic:
            last_char = residue_list.pop()
            residue_list.reverse()
            residue_list.append(last_char)
        else:
            residue_list.reverse()
        return Peptide(string.join(residue_list, sep = ''))

    def random(self):
        """Return a peptide of the same length, generated at random.

        Each residue occurs uniformly at random.

        """
        if self.length == 0:
            return Peptide(self.seq)
        l = amino_acids()
        residue_list = list(l[i] for i in (random.randint(0, len(l)-1) for _ in
                                           xrange(self.length)))
        return Peptide(string.join(residue_list, sep = ''))

    def ideal_fragment_masses(self, tbl = 'monoisotopic',
                              mass_op = lambda x: x):
        """Return the masses of the synthetic fragmentation products of self.

        Consider a peptide of length n as a character sequence: p[0..n-1].
        We return two arrays of length n+1, (left_mass, right_mass), where

        for each i = 0...n
        left_mass[i] = mass(p[0..i-1])  [ i.e., left_mass[0] = 0 ]
        right_mass[i] = mass(p[i..n-1])  [ i.e., right_mass[n] = 0 ]

        and where mass( ) is the mass of a peptide sequence. We do
        *not* consider the effect of mobile proton losses, or neutral
        losses.

        Args:
            tbl: string indicating the type of mass table to use. See
                protein.peptide.mass_table() for all possible values.
            mass_op: transformation of each of the mass table elements,
                i.e., a scalar function which takes a float and returns
                a numeric value.

        Returns:
            A tuple containing two vectors.
            left_mass:
            right_mass:

        """
        residue_mass = mass_table(tbl)
        mass = lambda aa: mass_op(residue_mass[aa])
        left_mass = [ 0 ]
        for i, aa in enumerate(self.seq):
            left_mass.append(left_mass[i] + mass(aa))

        right_mass = [ 0 ]
        for i, aa in enumerate(self.seq[::-1]): # reversed string
            right_mass.append(right_mass[i] + mass(aa))
        right_mass.reverse()
        return (left_mass, right_mass)

    def n_masses(self, tbl = 'monoisotopic'):
        """ Return the n-term masses at each amide bond location.
        Args:
            tbl:    (optional) string indicating the type of mass table to use.
                    Defaults to 'monoisotopic'; otherwise, use 'average.'
        Returns:
            A list of the n-term masses for the peptide.
        """
        residue_mass = mass_table(tbl)
        mass = lambda aa: residue_mass[aa]
        n_masses = []
        # omit the last amino acid -- just want n-masses at amide bonds
        for i, aa in enumerate(self.seq[:-1]):
            if i == 0:
                n_masses.append(mass(aa))
            else:
                n_masses.append(n_masses[i-1] + mass(aa))
        return n_masses

    def c_masses(self, tbl = 'monoisotopic'):
        """ Return the c-term masses at each amide bond location.
        Args:
            tbl:    (optional) string indicating the type of mass table to use.
                    Defaults to 'monoisotopic'; otherwise, use 'average.'
        Returns:
            A list of the c-term masses for the peptide.
        """
        residue_mass = mass_table(tbl)
        mass = lambda aa: residue_mass[aa]
        c_masses = []
        # omit the first amino acid -- just want c-masses at amide bonds
        for i, aa in enumerate(self.seq[:0:-1]):
            if i == 0:
                c_masses.append(mass(aa))
            else:
                c_masses.append(c_masses[i-1] + mass(aa))
        return c_masses

    def predict_ions(self, ion_type, charge = 1, tbl = 'monoisotopic'):
        """Return the m/z series of a given ion type, assuming monoisotopic
        masses.
        Args:
            ion_type: a, b, or y
            charge: 1 (default) or 2
            tbl: (optional) string indicating the type of mass table to use.
        Returns:
            a list of ion fragmentation m/z locations
        """
        ###### John note: Below is a mix of the average and monoisotopic masses of elements,
        ###### ie mass(elem == 'H') is the average mass of H, mass(elem == 'C') is the monoisotopic mass of C,
        ###### mass(elem == 'N') is the average mass of N, and mass(eleme == 'O') is both the average and
        ###### monoisotopic mass of O to 4 significant digits
        # elem = {'H': 1.00794, 'C': 12, 'N': 14.00674, 'O': 15.9994}

        # Fixed elem:
        elem = {'H': 1.0078246, 'C': 12, 'N': 14.003074, 'O': 15.9994}
        if charge != 1 and charge != 2:
            raise ValueError("charge must be 1 or 2")
        if ion_type.lower() == 'a':
            ions = self.n_masses(tbl)
            for i in range(0,len(ions)):
                ions[i] += charge*elem['H'] - elem['C'] - elem['O']
                ions[i] = ions[i]/charge
        elif ion_type.lower() == 'b':
            ions = self.n_masses(tbl)
            for i in range(0,len(ions)):
                ions[i] = (ions[i] + charge*elem['H'])/charge
        elif ion_type.lower() == 'y':
            ions = self.c_masses(tbl)
            for i in range(0, len(ions)):
                ions[i] = (ions[i] + (charge + 2)*elem['H'] + elem['O'])/charge
        else:
            raise ValueError("ion_type must be one of a, b, or y")
        return ions


def load_peptide_list(filename):
    """Load a file of peptide sequences into a list of peptides.

    Args:
        filename: A file with peptide sequences, one per line.

    Return:
        A list of Peptide objects, one for each line in filename.

    """
    try:
        lst = [ ]
        f = open(filename)
        for line in f:
            p = Peptide(line.strip())
            lst.append(p)
    except IOError:
        pass
    return lst

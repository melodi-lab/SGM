#!/usr/bin/env python
#
# Copyright 2011 <fill in later>

"""Peptide database populated by a list of peptides, instead of the results
from crux search-for-matches, which is what crux.peptide_db does."""

from __future__ import with_statement

import os
import protein.peptide
import util.file
import random
import types

class PeptideDB(object):

    def __init__(self, source, tbl = 'monoisotopic'):
        """Load a peptide database from a file.

        There is no specific assignment of peptides to spectra.

        Arguments:
           source: If a filename, load the peptides. If an iterable, use
               the iterable to populate the list of peptides.

        """
        self.peptides = [ ]
        self.masses = [ ]
        self.filename = None
        self.tbl = tbl

        if type(source) == types.StringType:
            source = os.path.realpath(os.path.expanduser(source))
            assert os.path.exists(source), 'Bad path %s' % source
            self.peptides = list(self._parser(source))
            self.filename = source
        else:
            self.peptides = list(source)
            self.filename = '<iterable>'

        if not self.peptides:
            raise RuntimeError('No peptides loaded.')

        # These must be M+H+ masses, since get_peptides uses the
        # MS2Spectrum.precursor_mass field to select peptides, and that is
        # an M+H+ value. For a peptide, this is the sum of the residue masses
        # + water
        if tbl == 'monoisotopic':
            self.masses = [ p.mass for p in self.peptides ]
        elif tbl == 'average':
            self.masses = [ p.average_mass for p in self.peptides ]

    def __str__(self):
        out = '%s: from file %s. Contains %d peptides.' % (
            self.__class__.__name__, self.filename, len(self.peptides))
        return out

    def __len__(self):
        return len(self.peptides)

    def _parser(self, filename):
        log = util.file.CommentedFile(open(filename))
        for line in log:
            yield protein.peptide.Peptide(line.strip())

    def make_decoys(self):
        decoys = [ p.shuffle() for p in self.peptides ]
        tl = set(self.peptides)
        for i in xrange(len(decoys)):
            j = 0
            while decoys[i] in tl and j < 10:
                decoys[i] = random.choice(self.peptides).shuffle()
                j = j + 1
        db = PeptideDB(decoys, self.tbl)
        assert len(self) == len(db)
        return db

    def make_filtered(self, spectra, mass_tol = 3.0):
        return PeptideDB(self._filter(spectra, mass_tol), self.tbl)

    def _filter(self, spectra, mass_tol):
        for p, m in zip(self.peptides, self.masses):
            if any(abs(m - s.precursor_mass) <= mass_tol for s in spectra):
                yield p

    def get_peptides(self, spectrum, mass_tol = 3.0):
        """Return the peptides that are within the given mass_tolerance.

        Arguments:
           spectrum: Instance of MS2Spectrum.
           mass_tol: Return all peptides within +/- mass_tol, in Da.

        Returns:
           peptides: returns an iterable of peptides which statisfy the
           mass tolerance criterion.

        """
        for p, m in zip(self.peptides, self.masses):
            if abs(m - spectrum.precursor_mass) <= mass_tol:
                yield p

#!/usr/bin/env python
#
# Copyright 2010 <fill in later>

"""Tools for creating a peptide database, usually using the output generated
by crux search-for-matches: i.e., search.{target,decoy}.txt.

"""

from __future__ import with_statement

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

from bisect import bisect_left, bisect_right
from operator import itemgetter
from protein.peptide import Peptide
import sys
import csv

class CruxPeptideDB(object):

    def __init__(self, filename, relevant_ids = None):
        """Load the Crux peptide database from a file.

        Only the peptides matched to spectra in relevant_ids are stored,
        which saves a substantial amount of memory.

        Uses the 'diverted' file format, where the output of crux
        search-for-matches is preprocessed into a file with the following
        record format:

        <sid> <peptide> <Sp score> <Xcorr score>

        The script used to convert search.{target,decoy}.txt to the diverted
        format is in scripts/crux/divertPeptides.py.

        Args:
            filename: Name of the peptide database (e.g.,
                search.target.diverted.txt)
            relevant_ids: Set of spectrum ids we care about.

        """
        self.db = [ ]
        self.keys = None
        self.filename = filename
        self._parser(filename, relevant_ids)

    def __str__(self):
        nspectra = len(set(self.keys))
        npeptides = len(self.db)
        out = '%s: from file %s. Contains %d spectra and %d peptides.' % (
            self.__class__.__name__, self.filename, nspectra, npeptides)
        return out

    def _parser(self, filename, relevant_ids):
        print >> sys.stderr, 'Loading matches from %s' % filename
        with open(filename) as log:
            for line in log:
                tokens = line.split()
                spectrum_id = int(tokens[0])
                if relevant_ids == None or (spectrum_id in relevant_ids):
                    self.db.append( (spectrum_id, tokens[1], float(tokens[2])) )
            self.db.sort(key = itemgetter(0))
            self.keys = [ r[0] for r in self.db ]

    def peptides(self, spectrum_id, k = None, paranoid = False):
        """Return the peptides that were listed as matches for a spectrum
        in the peptide database.

        Args:
            spectrum_id: Spectrum identifier in the Crux peptide database.
            k: Use only the k peptides with highest SP score.

        Returns:
            An iterable over peptides matched to the spectrum in the database,
            i.e., an iterable over protein.peptide.Peptide instances.

        Raises:
            ValueError: If no matches can be found for the spectrum_id.

        """
        low = bisect_left(self.keys, spectrum_id)
        high = bisect_right(self.keys, spectrum_id)

        if paranoid:
            for record in self.db[low:high]:
                assert(record[0] == spectrum_id)

        if k:
            lst = self.db[low:high]
            lst.sort(key = itemgetter(2), reverse = True)
            res = ( Peptide(p[1]) for p in lst[0:k] )
        else:
            res = ( Peptide(p[1]) for p in self.db[low:high] )
        return res

    def records(self, spectrum_id):
        """Yield the records for a particular spectrum."""

        low = bisect_left(self.keys, spectrum_id)
        high = bisect_right(self.keys, spectrum_id)
        return self.db[low:high]

class SimplePeptideDB(object):
    """A simpler version of a peptide database that takes the peptides.txt
    file produces by crux create-index --peptide-list T <fasta-file> <index-dir>

    The peptide database is a text file where each line is of the form
    <peptide> <neutral mass>

    Where peptide is a IUPAC sequence, and neutral mass it the average mass
    of a peptide (c.f.,
    http://noble.gs.washington.edu/proj/crux/crux-search-for-matches.html)

    """
    def __init__(self, filename):
        self.filename = filename
        self.peptides = [ ]
        self.masses = None
        self._parser(filename)

    def _parser(self, filename):
        records = [ ]
        # first detect what is the used delimiter; crux switched from a space to a tab when shifting to 2.0
        with open(filename, "r") as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(), delimiters=' \t')
            csvfile.seek(0)
            for line in csv.reader(csvfile, dialect):
                if line: records.append( (line[0], float(line[1])) )
            records.sort(key = lambda r: r[1])
        self.peptides, self.masses = zip(*records)

    def filter(self, precursor_mass, window = 3, ppm = False):
        """Return the sequence of the peptides that are within 'window'
        Daltons of the mass given: i.e., where the peptide's mass is within
        [precursor_mass - window, precursor_mass + window]

        Note: avgmass is what is reported on the neutral mass line, and
        that should be a more reasonable mass calculation since the precursor
        mass comes from MS1, which reflects the average mass.

        """
        if not ppm:
            l = bisect_left(self.masses, precursor_mass - window)
            h = bisect_right(self.masses, precursor_mass + window, lo = l)
        else:
            wl = window*1e-6
            # print "precursor_mass=%f" % (precursor_mass)
            # print "l=%f" % (precursor_mass / (1.0 + window*0.000001))
            # print "h=%f" % (precursor_mass / (1.0 - window*0.000001))
            l = bisect_left(self.masses, precursor_mass / (1.0 + window*0.000001))
            h = bisect_right(self.masses, precursor_mass / (1.0 - window*0.000001), lo = l)
            # print "l=%f" % (precursor_mass * (1.0 - wl))
            # print "h=%f" % (precursor_mass * (1.0 + wl))
            # l = bisect_left(self.masses, precursor_mass * (1.0 - wl))
            # h = bisect_right(self.masses, precursor_mass * (1.0 + wl), lo = l)
            # print "l=%d" % (l)
            # print "h=%d" % (h)
            # print self.peptides[l:h]
        return self.peptides[l:h]
        # if l <= h:
        #     sorted_by_closeness = sorted(zip(self.peptides[l:h], self.masses[l:h]),
        #                                  key = lambda r: abs(r[1] - precursor_mass))
        #     return list(p for p, _ in sorted_by_closeness)
        # else: # empty set of candidates
        #     return []

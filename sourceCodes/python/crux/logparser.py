#!/usr/bin/env python
#
# Copyright 2010 <fill in later>

from __future__ import with_statement

"""Tools for parsing logs output by Crux and Tide.

References:
[1] http://noble.gs.washington.edu/proj/crux/

"""

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import csv
import os
import sys

from operator import itemgetter
from protein.peptide import Peptide
from itertools import islice

from progress import ProgressBar, Percentage, Bar, ETA
from buzhug import Base

class ParseError(Exception):
    def __init__(self, prefix, filename, lineno, line):
        self.prefix = prefix
        self.filename = filename
        self.lineno = lineno
        self.line = line
    def __str__(self):
        return repr('%s%s(%d): %s' %
                    (self.prefix, self.filename, self.lineno, self.line))

def parse_crux_search_txt(filename):
    """Iterate over records in a search.{target,decoy}.txt.

    Crux txt format files are tab-delimited with 30 fields*, described
    in the online documentation [1]. This function returns an iterator
    which yields a dictionary with the fields and their values.

    * 'decoy q-value (p-value)' is not output by Crux, at least as of v1.33.

    [1] http://noble.gs.washington.edu/proj/crux/txt-format.html

    Arguments:
       filename: Name of the crux search-for-matches output.

    Returns:
       Dictionary that maps field names to values. Only fields that
       are non-empty in the input exist in the returned dictionary.
       Many of the fields are not usually set in the output of crux
       search-for-matches, and will not be available.

    """
    fields = ['scan', # int
              'charge', # int
              'spectrum precursor m/z', # float
              'spectrum neutral mass', # float
              'peptide mass', # float
              'delta_cn', # float
              'sp score', # float
              'sp rank', # float
              'xcorr score', # float
              'xcorr rank', # int
              'p-value', # float
              'Weibull est. q-value', # float
              'decoy q-value (xcorr)', # float
              'percolator score', # float
              'percolator rank', # int
              'percolator q-value', # float
              'q-ranker score', # float
              'q-ranker q-value', # float
              'b/y ions matched', # int
              'b/y ions total', # int
              'matches/spectrum', # int
              'sequence', # string
              'cleavage type', # string
              'protein id', # string
              'flanking aa', # string
              'unshuffled sequence', # string
              'eta', # float
              'beta', # float
              'shift', # float
              'corr'] # float
    casts = [ int, int, float, float, float, float, float, float, float, int,
              float, float, float, float, int, float, float, float, int, int,
              int, str, str, str, str, str, float, float, float, float ]
    assert(len(fields) == len(casts))

    _mandatories = [ 'scan', 'charge', 'spectrum precursor m/z',
                     'spectrum neutral mass', 'xcorr score',
                     'xcorr rank', 'sequence' ]

    def conv(f, value):
        value = value.strip()
        if len(value):
            return f(value)

    def validate(record):
        return all(record.has_key(m) for m in _mandatories)

    widgets = [ Percentage(), Bar(), ETA() ]
    progress = ProgressBar(widgets = widgets,
                           maxval = os.path.getsize(filename)).start()

    with open(filename) as f:
        reader = csv.reader(f, delimiter='\t')
        # Header
        row = reader.next()
        if row != fields:
            raise ParseError('Header: ', filename, 1, ' '.join(row))
        # Body
        for row in reader:
            progress.update(f.tell())
            if len(row) != len(fields):
                raise ParseError('Line: ', filename, reader.line_num,
                                 ' '.join(row))

            r = dict((k, conv(f,x)) for k, f, x in zip(fields, casts, row))
            if r:
                if not validate(r):
                    raise ParseError('Missing: ', filename, reader.line_num,
                                     ' '.join(row))
                yield r

    progress.finish()
    sys.stdout.write('\n')

def parse_crux_search_txt_dict(filename, scorer = 'xcorr', sortlists = True,
                               nrecords = None):
    """Extract scored peptide-spectrum matches from crux search-for-matches.

    It is the caller's responsibility to make sure that the crux log file
    has the required type of score. XCorr is always returned, but other scores
    may or may not be defined in 'filename'.

    Arguments:
       filename: Text output generated by Crux, usu. search.{target,decoy}.txt.
       scorer: Score type to use, one of {'xcorr', 'percolator', 'q-ranker',
               'delta_cn', 'sp'}
       sortlists: If true, the list of candidates is sorted in order of
                  decreasing score.
       nrecords: Number of records to read in. If None, read all the records.

    Returns:
       A dictionary which maps scan id to a list of (peptide, score) tuples,
       which are the peptides tested against the spectrum.

    """
    mapper = { 'xcorr' : 'xcorr score',
               'percolator' : 'percolator score',
               'q-ranker' : 'q-ranker score',
               'delta_cn' : 'delta_cn',
               'sp' : 'sp score' }
    if scorer not in mapper:
        raise ValueError('"%s" is not a valid scorer' % scorer)

    dic = { }
    for record in islice(parse_crux_search_txt(filename), nrecords):
        sid = record['scan']
        if sid not in dic:
            dic[sid] = [ ]

        score = record[mapper[scorer]]
        peptide = Peptide(record['sequence'])
        dic[sid].append((peptide, score))

    if sortlists:
        for sid, lst in dic.iteritems():
            dic[sid] = sorted(lst, key = itemgetter(1), reverse = True)

    return dic


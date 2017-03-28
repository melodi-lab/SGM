#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Get the m/z range for a set of spectra or PSM file
#
# Usage: ./mzrange.py -s <base>/data/foo.ms2.gz -p <base>/data/foo.psm

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os
import platform
import sys

from optparse import OptionParser

from msutil.spectrum import MS2Iterator as MS2_iterator
from protein.peptide import Peptide
from msutil.psm import PeptideSpectrumMatches

def mz_range_spectra(filename):
    # TODO(ajit): Figure out how to get DBL_MAX in python.
    mz_min = 10e10
    mz_max = 0
    fn = os.path.abspath(filename)
    for spectrum in MS2_iterator(fn):
        mz_min = min(mz_min, min(spectrum.mz))
        mz_max = max(mz_max, max(spectrum.mz))
    return (mz_min, mz_max)

def mz_range_psm(filename):
    matches = PeptideSpectrumMatches(filename, True)
    masses = [ p.mass for _, p in matches ]
    masses.extend( [ p.average_mass for _, p in matches ] )
    return (min(masses), max(masses))

if __name__ == "__main__":
    parser = OptionParser(usage = "usage: %prog -i SPECTRA_FILE")
    parser.add_option("-s", "--spectra",
                      action = "store", type = "string", dest = "spectra",
                      default = None,
                      help = "Name of MS2 formatted spectra file to read.")
    parser.add_option("-p", "--psm",
                      action = "store", type = "string", dest = "psm",
                      default = None,
                      help = "Nam eof PSM file to read.")
    (options, args) = parser.parse_args()

    if not options.spectra and not options.psm:
        parser.print_help()
        exit(-1)

    if options.spectra != None:
        (low, high) = mz_range_spectra(options.spectra)
        print('Spectra m/z range is [%f, %f]' % (low, high))

    if options.psm:
        (low, high) = mz_range_psm(options.psm)
        print('PSM m/z range is [%f, %f]' % (low, high))


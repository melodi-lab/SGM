#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Compute intensity quantiles of a set of spectra.
#
# Usage:

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os
import sys

from optparse import OptionParser

from msutil.ms2util import MS2_iterator
from util.statistics import mean, variance
from recipes.quantile import quantile

if __name__ == "__main__":
    parser = OptionParser(usage = "usage: %prog [options]")
    parser.add_option("-i", "--input",
                      action = "store", type = "string", dest = "filename",
                      help = "Name of MS2 formatted spectra file to read.")
    parser.add_option("-z", "--use_gzcat",
                      action = "store_true", dest = "gzcat",
                      default = False, help = "Use gzcat to decompress.")
    (options, args) = parser.parse_args()

    if options.filename == '':
        parser.print_help()
        exit(-1)

    intensities = [ ]
    fn = os.path.abspath(options.filename)
    for s in MS2_iterator(fn, options.gzcat):
        intensities.extend(s.intensity)
    intensities = sorted(intensities)

    qtype = 8
    low = quantile(intensities, 1.0/3.0, qtype = qtype, issorted = True)
    med = quantile(intensities, 0.5, qtype = qtype, issorted = True)
    high = quantile(intensities, 2.0/3.0, qtype = qtype, issorted = True)

    print (low, med, high)


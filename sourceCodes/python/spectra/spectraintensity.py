#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Count spectra in an MS2 formatted file (.ms2, .ms2.gz)
#
# Usage: ./spectraintensity -i <base>/data/spectrum.ms2.gz [-z]

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os
import sys

from optparse import OptionParser

from msutil.ms2util import MS2_iterator
from util.statistics import mean, variance

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

    fn = os.path.abspath(options.filename)
    for s in MS2_iterator(fn, options.gzcat):
        print("Spectrum %d: mean intensity %.4f, intensity variance %.4f" %
              (s.spectrum_id, mean(s.intensity), variance(s.intensity)))


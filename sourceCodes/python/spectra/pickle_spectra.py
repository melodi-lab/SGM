#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Convert an .ms2 or .ms2.gz file to Python's binary pickle format,
# which can be faster to load. Pickled files should be treated as
# temporary, cache, files. Unlike .ms2/.ms2.gz files, there is no
# guarantee that a pickled collection of spectra will be readable by
# future versions of MS2_iterator.
#
# Usage: ./pickle_spectra.py -i SPECTRUM_FILE -o OUTPUT_FILE [-z]
#
# Example: ./pickle_spectra.py -i ../data/spectrum.ms2 -o /tmp/spectrum.gzpickle

from __future__ import with_statement

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import cPickle
import gzip
import os
import sys

import util.file
from optparse import OptionParser

from msutil.ms2util import MS2_iterator

if __name__ == "__main__":
    parser = OptionParser(usage = "usage: %prog [options]")
    parser.add_option("-i", "--input",
                      action = "store", type = "string", dest = "input",
                      help = "Name of MS2 formatted spectra file to read.")
    parser.add_option("-o", "--output",
                      action = "store", type = "string", dest = "output",
                      help = "Name of file to pickle to.")
    parser.add_option("-z", "--use_gzcat",
                      action = "store_true", dest = "gzcat",
                      default = False, help = "Use gzcat to decompress.")
    (options, args) = parser.parse_args()

    if not options.input or not options.output:
        parser.print_help()
        exit(-1)

    if os.path.exists(options.output):
        print('File already exists: %s' % options.output)
        exit(-2)

    with util.file.GzipFile(options.output, 'wb') as outfile:
        fn = os.path.abspath(options.input)
        for s in MS2_iterator(fn, options.gzcat):
            cPickle.dump(s, outfile, cPickle.HIGHEST_PROTOCOL)


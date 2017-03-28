#!/usr/bin/env python
#
# Copyright 2010 <fill in later>
#
# Print a single spectrum using Matplotlib.

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os
import sys

try:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
except ImportError:
    print >> sys.stderr, \
    'Module "pylab" not available. You must install matplotlib to use this.'
    exit(-1)
from optparse import OptionParser

from msutil.ms2util import MS2_iterator


def print_spectrum(spectrum, filename):
    x = spectrum.mz
    y = spectrum.intensity
    pylab.stem(x, y, markerfmt='b,')
    pylab.xlabel('m/z')
    pylab.ylabel('Transformed intensity')
    pylab.grid()
    pylab.xlim(0.95*min(x), 1.05*max(x))
    pylab.savefig(filename)

if __name__ == "__main__":
    parser = OptionParser(usage = "usage: %prog [options]")
    parser.add_option("-i", "--input",
                      action = "store", type = "string", dest = "filename",
                      help = "Name of MS2 formatted spectra file to read.")
    parser.add_option("-z", "--use_gzcat",
                      action = "store_true", dest = "gzcat",
                      default = False, help = "Use gzcat to decompress.")
    parser.add_option("-o", "--output", default = None,
                      action = "store", type = "string", dest = "output",
                      help = "Name of file to create.")
    parser.add_option("--id",
                      action = "store", type = "int", dest = "spectrum_id",
                      default = None, help = "ID of spectrum to print.")
    (options, args) = parser.parse_args()

    if options.filename == '':
        parser.print_help()
        exit(-1)

    if not options.spectrum_id:
        parser.print_help()
        exit(-1)

    if not options.output:
        option.output = '%s.pdf' % options.spectrum_id

    fn = os.path.abspath(options.filename)
    for spectrum in MS2_iterator(fn, options.gzcat):
        if spectrum.spectrum_id == options.spectrum_id:
            print_spectrum(spectrum, options.output)
            exit(0)


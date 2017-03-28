#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

from __future__ import with_statement

import optparse
import os
import sys

#from msutil.ms2util import MS2_iterator_sample
from msutil.spectrum import MS2IteratorSample
from util.args import read_file_callback, write_file_callback

if __name__ == "__main__":
    usage = 'Usage: %prog --spectra MS2FILE --output MS2FILE --num_spectra INT'
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('--spectra', default = None, action = 'callback',
                      type = 'string', callback = read_file_callback,
                      dest = 'spectra_file',
                      help = 'File containing spectra.')
    parser.add_option('--output', default = None, action = 'callback',
                      type = 'string', callback = write_file_callback,
                      dest = 'output_file',
                      help = 'New spectrum file to create.')
    parser.add_option('--num_spectra', default = None, type = 'int',
                      action = 'store', dest = 'num_spectra',
                      help = 'Number of spectra to sample.')
    parser.add_option("-z", "--use_gzcat",
                      action = "store_true", dest = "gzcat",
                      default = False, help = "Use gzcat to decompress.")
    (options, args) = parser.parse_args()

    mandatories = ['spectra_file', 'output_file', 'num_spectra']
    for m in mandatories:
        if not options.__dict__[m]:
            print('Mandatory argument "%s" is missing.' % m)
            parser.print_help()
            exit(-1)

    with open(options.output_file, 'w') as f:
        for spectrum in MS2IteratorSample(options.spectra_file,
                                          options.num_spectra,
                                          options.gzcat):
            spectrum.write(f)



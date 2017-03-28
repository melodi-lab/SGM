#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

from __future__ import with_statement

import optparse
import os
import sys

from util.file import extension
from itertools import ifilter
from msutil.ms2util import MS2_iterator
from util.args import read_file_callback, write_file_callback
from util.iterables import grouper

if __name__ == "__main__":
    usage = 'Usage: %prog --spectra MS2FILE --num_spectra INT [-z]'
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('--spectra', default = None, action = 'callback',
                      type = 'string', callback = read_file_callback,
                      dest = 'spectra_file',
                      help = 'File containing spectra.')
    parser.add_option('--num_spectra', default = None, type = 'int',
                      action = 'store', dest = 'num_spectra',
                      help = 'Number of spectra per-chunk.')
    parser.add_option("-z", "--use_gzcat",
                      action = "store_true", dest = "gzcat",
                      default = False, help = "Use gzcat to decompress.")
    (options, args) = parser.parse_args()

    mandatories = ['spectra_file', 'num_spectra']
    for m in mandatories:
        if not options.__dict__[m]:
            print('Mandatory argument "%s" is missing.' % m)
            parser.print_help()
            exit(-1)

    _, basename = extension(options.spectra_file)
    spectra = MS2_iterator(options.spectra_file, options.gzcat)

    for part, group in enumerate(grouper(options.num_spectra, spectra)):
        with open('%s-part%d.ms2' % (basename, part), 'w') as f:
            print >> sys.stderr, 'Writing %s' % f.name
            for spectrum in ifilter(lambda x: x, group):
                spectrum.write(f)



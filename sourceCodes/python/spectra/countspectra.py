#!/usr/bin/env python
#
# Copyright 2010 <fill in later>
#
# Count spectra in an MS2 formatted file (.ms2, .ms2.gz), or all .ms2 files under
# the directory given as input.
#
# Usage: ./countSpectra.py -i <base>/data/spectrum.ms2.gz [-z]
# Usage: ./countSpectra.py -i <dir> -z

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import fnmatch
import os
import sys

from optparse import OptionParser

from msutil.spectrum import MS2Iterator


def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.abspath(os.path.join(root, basename))
        yield filename

def count_in_file(filename, gzcat = False):
    fn = os.path.abspath(filename)
    num_spectra = sum(1 for _ in MS2Iterator(fn, gzcat))
    return num_spectra


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

    num_spectra = 0    
    if os.path.isdir(options.filename):
      for fn in find_files(options.filename, '*.ms2*'):
        count = count_in_file(fn, True)
        num_spectra += count
        print("Found %d spectra in %s" % (count, os.path.split(fn)[1]))
      print("There are %d spectra in directory %s" % (num_spectra, options.filename))
    else:
      num_spectra = count_in_file(options.filename, options.gzcat)
      print("There are %d spectra in %s" % (num_spectra, os.path.split(options.filename)[1]))

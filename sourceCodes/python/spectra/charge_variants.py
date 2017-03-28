#!/usr/bin/env python
#
# Copyright 2011 <fill in later>
#
# Count charge variations among spectra in an MS2 file

import collections
import optparse
import os

from msutil.spectrum import MS2Iterator

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = "usage: %prog [options]")
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

    charges = collections.defaultdict(int)
    for s in MS2Iterator(fn, options.gzcat):
        c = tuple(r[0] for r in s.charge_lines)
        charges[c] += 1

    for c, count in charges.iteritems():
        print '%s occurs in %d spectra.' % (str(c), count)


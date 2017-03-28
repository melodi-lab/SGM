#!/usr/bin/env python
#
# Copyright 2011 <fill in later>
#
# Filter spectra by charge.
#
# Usage: To retain only spectra with a reported charge of 2 or 3.
#        filter_spectra.py -z -i <ms2file> -c 2,3 > foo.ms2

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os
import sys

from optparse import OptionParser

from msutil.spectrum import MS2Iterator

def prune_chargelines(s, chargelist):
    """Ensure that there is exactly one charge line for each spectra.

    Some spectra have multiple Z-lines with the same charge, but different
    m/z values (I don't know why). This routine picks the first such chargeline,
    and remove all extra chargelines.

    The resulting files are pretty useful when we only want to test charges +2 and +3,
    but it's obviously a bad idea to do this from an end-point sort of view.

    Arguments:
      charges: set of charges
      s: instance of MS2Spectrum

    """
    newcl = [ ]
    for charge, mass in s.charge_lines:
        if charge in [ c for c, m in newcl ]:
            continue
        if charge not in chargelist:
            continue
        newcl.append((charge, mass))
    s.charge_lines = newcl       


if __name__ == "__main__":
    parser = OptionParser(usage = "usage: %prog [options]")
    parser.add_option("-i", "--input",
                      action = "store", type = "string", dest = "filename",
                      help = "Name of MS2 formatted spectra file to read.")
    parser.add_option("-z", "--use_gzcat",
                      action = "store_true", dest = "gzcat",
                      default = False, help = "Use gzcat to decompress.")
    parser.add_option("-c", "--spectrum-charges",
                      action = "store", type = "string", dest = "charges",
                      help = "Spectrum charges to retain.")
    (options, args) = parser.parse_args()

    if options.filename == '':
        parser.print_help()
        exit(-1)

    charges = set(int(c) for c in options.charges.split(','))
    print >> sys.stderr, 'INFO: Retaining spectra with charges %s' % (
        str(charges))

    count = 0
    fn = os.path.abspath(options.filename)
    for s in MS2Iterator(fn,options.gzcat):
        if set(s.charges) & charges:
            prune_chargelines(s, charges)
            s.write(sys.stdout)
            count = count + 1

    print >> sys.stderr, 'INFO: Wrote %d spectra with charges %s' % (
        count, str(charges))

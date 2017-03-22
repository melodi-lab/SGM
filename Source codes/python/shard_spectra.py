#!/usr/bin/env python
#
# Copyright 2011
# authors: John Halloran, Ajit Singh

from __future__ import with_statement

import cPickle as pickle
import collections
import csv
import itertools
import math
import time
import optparse
import os
import random
import sys

from crux.peptide_db import SimplePeptideDB
from msutil.spectrum import MS2Iterator, SampleSpectra
#from progress import ETA, ProgressBar
import util.args
import util.iterables
from visualize.parsers import load_ident

def shard_spectra_file(options):
    """Split a spectra file into pieces, or shards.

    Arguments:
        options.spectra: File with the spectra.
        options.gzcat: Use gzcat to decompress, if file is gzipped.
        options.num_spectra: Number of spectra to select, at random.

    Returns:
        (spectra, n) - An iterator of MS2Spectrum objects, and the
        number of spectra contained in the iterator.

    """

    # Consider pre-specified charges
    validcharges = set(int(c) for c in options.charges.split(','))
    spectra = list(s for s in MS2Iterator(options.spectra, options.gzcat) if
                   set(s.charges) & validcharges) 
    print >> sys.stderr, '%s has %d spectra with charges %s' % (
        options.spectra, len(spectra),
        '{' + ','.join('+%d' % c for c in validcharges) + '}')

    # # Only use spectra where 2+ is a possible charge state.
    # spectra = list(s for s in MS2Iterator(options.spectra, options.gzcat) if
    #                set(s.charges) & set([2]))

    # # Only use spectra where 3+ is a possible charge state.
    # spectra = list(s for s in MS2Iterator(options.spectra, options.gzcat) if
    #                set(s.charges) & set([3]))

    minMz = 2000.0
    maxMz = -1
    if options.filter_ident:
        print "Filtering spectra from ident %s" % (options.ident)
        targets, decoys = load_ident(options.ident)
        sids = [ t[0] for t in targets ]

        s = []
        for spec in spectra:
            if spec.mz[0] < minMz:
                minMz = spec.mz[0]
            if spec.mz[-1] > maxMz:
                maxMz = spec.mz[-1]
            if(spec.length <= options.max_spectrum_length): 
                if(spec.spectrum_id not in sids):
                # if(spec.spectrum_id in sids):
                    s.append(spec)
    else:
        s = []
        for spec in spectra:
            if spec.mz[0] < minMz:
                minMz = spec.mz[0]
            if spec.mz[-1] > maxMz:
                maxMz = spec.mz[-1]
            if(spec.length >= options.min_spectrum_length): 
            # if(spec.length <= options.max_spectrum_length): 
                s.append(spec)
                
    spectra = s
    print "minMz=%f" % minMz
    print "maxMz=%f" % maxMz

    if len(spectra) == 0:
        print >> sys.stderr, 'There are no spectra with charge_line (+2).'
        exit(-1)
    elif options.num_spectra and options.num_spectra < len(spectra):
        spectra = SampleSpectra(spectra, options.num_spectra)
        return spectra, options.num_spectra
    else:
        random.shuffle(spectra) # comment this out later
        return spectra, len(spectra)

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--min_spectrum_length', type = 'int', action = 'store',
                      dest = 'min_spectrum_length', default = 1)
    parser.add_option('--max_spectrum_length', type = 'int', action = 'store',
                      dest = 'max_spectrum_length', default = 10000)
    parser.add_option('--ident', type = 'string', action= 'store')
    parser.add_option('--spectra', type = 'string', action = 'store')
    parser.add_option('--targets', type = 'string', action = 'store')
    parser.add_option('--decoys', type = 'string', action = 'store')
    parser.add_option('--output_dir', type = 'string', action = 'store')
    parser.add_option('--shards', type = 'int', action = 'store')
    parser.add_option('--num_spectra', type = 'int', action = 'store')
    parser.add_option('--tol', type = 'float', action = 'store', default = 3.0)
    parser.add_option('--ppm', action = 'store_true',
                      default = False,
                      help = "Is mass tolerance measured in part per million.")
    parser.add_option('--filter_ident', action = 'store_true',
                      default = False,
                      help = "Filter sids by ident")
    parser.add_option('--charges', action = "store", type = "string",
                      default = "2",
                      help = "Spectrum charges to retain.")
    parser.add_option('-z', '--use_gzcat', action = 'store_true',
                      dest = 'gzcat', default = False,
                      help = "Use gzcat to decompress .ms2.gz files, if avail.")

    (options, args) = parser.parse_args()
    req = [ 'spectra', 'shards', 'targets', 'output_dir' ]
    if not util.args.mandatories_defined(req, options):
        parser.print_help()
        exit(-1)
    assert(options.shards > 0)

    # Compute the number of spectra in each slice.
    spectra, n = shard_spectra_file(options)
    if n < options.shards:
        print('More shards than spectra: %d vs. %d, defaulting to %d shards' % (
        options.shards, n, max(int(n/4), 1)))
        options.shards = max(int(n/4), 1)
    # assert n >= options.shards, 'More shards than spectra: %d vs. %d' % (
    #     options.shards, n)
    sz = int(math.floor(float(n)/options.shards))
    print options.shards
    print >> sys.stderr, 'Each shard has at most %d spectra' % sz

    targets = SimplePeptideDB(options.targets)
    decoys = SimplePeptideDB(options.decoys)

    mass_h = 1.00727646677
    # mass_h = 1.00782503207
    # mass_h = 1.01
    validcharges = set(int(c) for c in options.charges.split(','))
    for part, group in enumerate(util.iterables.grouper(sz, spectra)):
        #print >> sys.stderr, '%d' % part
        spectra_app = []
        spectra_list = list(s for s in group if s)
        emptySpectra = 0
        # Use neutral mass to select peptides (Z-lines report M+H+ mass)
        t = collections.defaultdict(list)
        d = collections.defaultdict(list)
        for s in spectra_list:
            spec_added = 0
            repSpec = 0
            for c, m in s.charge_lines:
                if c in validcharges:
                    tc = targets.filter(m - mass_h, options.tol, options.ppm)
                    dc = decoys.filter(m - mass_h, options.tol, options.ppm)
                    if len(tc) > 0 and len(dc) > 0:
                        print "sid=%d, charge=%d, num targets=%d, num decoys=%d" % (s.spectrum_id, c, len(tc), len(dc))
                        if (s.spectrum_id,c) in t: # high res MS1 may have multiple precursor charge estimates
                            t[(s.spectrum_id,c)] += tc
                            d[(s.spectrum_id,c)] += dc
                            repSpec = 1
                        else:
                            t[(s.spectrum_id,c)] = tc
                            d[(s.spectrum_id,c)] = dc
                        if not spec_added:
                            spectra_app.append(s)
                            spec_added = 1
                    else:
                        emptySpectra += 1
            if repSpec: # prune any multiply added peptide candidates
                for c in validcharges:
                    t[(s.spectrum_id, c)] = list(set(t[(s.spectrum_id, c)]))
                    d[(s.spectrum_id, c)] = list(set(d[(s.spectrum_id, c)]))
                

        print "%d spectra with empty candidate sets" % emptySpectra
        data = { 'spectra' : spectra_app,
                 'target' : t,
                 'decoy' : d,
                 'shardname' : 'shard-%d' % part }
        outfn = os.path.join(options.output_dir, 'shard-%d.pickle' % part)
        with open(outfn, 'w') as out:
            pickle.dump(data, out, pickle.HIGHEST_PROTOCOL)

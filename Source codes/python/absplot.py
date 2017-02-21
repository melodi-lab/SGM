#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2011 <fill in later>
#
# Usage: absplot.py --output foo.png Crux:crux-60cm.txt \"DBN A\":dbn-60cm.txt

import optparse
import sys
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    print >> sys.stderr, \
        'Module "matplotlib" not available.'
    exit(-1)

import util.args
import visualize.absolute
import visualize.parsers
import visualize.util


if __name__ == '__main__':

    usage = ('Usage: %prog [options] label1:IDENT1 ... labeln:IDENTn\n\n'
             'Example: %prog --output foo.png Crux:crux-60cm-ident.txt\n'
             'Example %prog --output foo.png "DBN Model":dbn-60cm-ident.txt')
    desc = ('The positional arguments IDENT1, IDENT2, ..., IDENTn are the '
            'names of spectrum identification files: tab-separated files with '
            'Kind, Sid, Peptide, and Score fields. It is assumed that the '
            'files are ordered in terms of specificity: the spectra matched '
            'in IDENT{n} are a subset of the spectra identified in IDENT{n+1}. '
            'Only the spectra identified in IDENT1 are used to generate the '
            'absolute ranking curve. The typical usage is that each of the '
            'ident files represent a particular algorithm, and we want to '
            'compare their performance under absolute ranking.\n\n '
            'The labeli prefixes are the name used to describe that method '
            'in the resulting plot. e.g., "Crux 1":crux-60cm-ident.txt.')

    parser = optparse.OptionParser(usage=usage, description=desc)

    help_output = 'Name of the file where the plot is stored.'
    parser.add_option('--output', type='string', help=help_output,
                      default=None)

    help_q = 'Maximum q-value to plot to: 0 < q <= 1.0'
    parser.add_option('--maxq', type='float', default=1.0, help=help_q)
    parser.add_option('--q', type='float', default=-10)

    help_publish = 'Generate plot for b/w publication.'
    parser.add_option('--publish', action='store_true', default=False,
                      help=help_publish)
    (options, args) = parser.parse_args()

    assert len(args) >= 1, 'No identification files listed.'

    labels = []
    scorelists = []
    sids = None

    if options.q < 0.0 and not options.output:
        print >> sys.stderr, 'You need to define one of --q, --output'
        exit(-2)

    def process(arg, sids=None, silent=False):
        desc, fn = util.args.parse_positional_pair(arg)
        labels.append(desc)
        targets, decoys = visualize.parsers.load_ident(fn, sids)
        scorelists.append((targets, decoys))
        if not silent:
            print 'Loaded %d spectrum identifications from %s' % (
                len(targets), fn)
        return set(r[0] for r in targets)

    sids = process(args[0])
    for argument in args[1:]:
        process(argument, sids)

    scorelists = visualize.util.refine(scorelists)

    if options.output:
        qrange = (0.0, options.maxq)
        visualize.absolute.plot(scorelists, options.output, qrange=qrange,
                                labels=labels, paranoid=True,
                                publish=options.publish)

    print "hi"

    if options.q > 0:
        if options.q < 0.0 or options.q > 1.0:
            print >> sys.stderr, '--q <float> must be a q-value in [0,1]'
            exit(-2)
        num_accepted = visualize.absolute.identifications(
            scorelists, options.q)
        for lbl, n in zip(labels, num_accepted):
            print '%s accepted %d spectra at q = %.3g' % (lbl, n, options.q)

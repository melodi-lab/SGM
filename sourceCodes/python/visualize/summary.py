#!/usr/bin/env python
#
# Copyright 2011 <fill in later>

"""Summary visualization, which are especially useful which analyzing behaviour:

- Histogram: plot the density of the score of target vs. decoy PSMs.
- Scatterplot : plot the correlation between two peptide identifications

"""

import os
try:
    import sys
    import matplotlib
    import pylab
    if not sys.modules.has_key('matplotlib'):
        matplotlib.use('Agg')
except ImportError:
    print >> sys.stderr, \
    'Module "pylab" not available. You must install matplotlib to use this.'
    exit(-1)


def histogram(targets, decoys, fn, bins = 40):
    """Histogram of the score distribution between target and decoy PSMs.

    Arguments:
        targets: Iterable of floats, each the score of a target PSM.
        decoys: Iterable of floats, each the score of a decoy PSM.
        fn: Name of the output file. The format is inferred from the
            extension: e.g., foo.png -> PNG, foo.pdf -> PDF. The image
            formats allowed are those supported by matplotlib: png,
            pdf, svg, ps, eps, tiff.
        bins: Number of bins in the histogram [default: 40].

    Effects:
        Outputs the image to the file specified in 'fn'.

    """
    pylab.clf()
    pylab.hold(True)
    pylab.xlabel('Score')
    pylab.ylabel('Pr(Score)')

    l = min(min(decoys), min(targets))
    h = max(max(decoys), max(targets))
    l = -0.2
    h = 0.0
    _, _, h1 = pylab.hist(targets, bins = bins, range = (l,h), normed = 1,
                          color = 'b', alpha = 0.25)
    _, _, h2 = pylab.hist(decoys, bins = bins, range = (l,h), normed = 1,
                          color = 'm', alpha = 0.25)
    pylab.legend((h1[0], h2[0]), ('Target PSMs', 'Decoy PSMs'), loc = 2)
    pylab.savefig('%s' % fn)

def scatterplot(first_method, second_method, fn, labels = None):
    """Scatterplot of the PSM scores for two methods.

    Assumes that the (target, decoy) tuples are ordered, so that
    target[i] and decoy[i] refer to the same spectrum.

    Arguments:
        first_method: (targets, decoys) tuple, where each entry in
           targets, and in decoys, is a triplet of the form (s, p, f)
           where s is a spectrum id, p is a peptide, and f is the score
           of the match. See the outputs of load_ident.
        second_method: Same form as first_method, but represents the
           identifications produced by another algorithm
        fn: Name of the file to plot to.
        labels: (first_method_name, second_method_name), the labels for
           the two methods, used to labelled in the x and y axes of the
           scatterplot.

    Effects:
        Creates a plot in file 'fn'. The position of each red point encodes
        the score of the top decoy PSM, for a spectrum. The position of
        each blue point encodes the score of the top target PSM, for a
        spectrum.

    """
    if len(first_method[0]) != len(second_method[0]):
        raise Exception('Methods evaluated over different spectra.')

    t1 = list(r[2] for r in first_method[0])
    d1 = list(r[2] for r in first_method[1])
    t2 = list(r[2] for r in second_method[0])
    d2 = list(r[2] for r in second_method[1])

    pylab.clf()
    pylab.hold(True)
    pylab.scatter(t1, t2, color = 'b', alpha = 0.20, s = 2)
    pylab.scatter(d1, d2, color = 'r', alpha = 0.10, s = 1)
    pylab.xlim( (min(min(t1), min(d1)), max(max(t1), max(d1))) )
    if labels:
        pylab.xlabel(labels[0])
        pylab.ylabel(labels[1])

    pylab.savefig(fn)

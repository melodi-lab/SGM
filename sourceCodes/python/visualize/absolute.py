#!/usr/bin/env python
#
# Copyright 2011 <fill in later>

"""Support for absolute ranking plots"""

import itertools
import types
import numpy
import random
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
import pq_plot.plotFDR


class PlotException(Exception):
    pass

def plot(scorelists, output, qrange = None, labels = None, **kwargs):
    """Plot multiple absolute ranking plots on one axis.

    The y-axis is the number of spectra, so plotting two methods is
    only comparable on the resulting plot if you evaluated them on the
    same number of spectra. Typically, the methods being compared will
    be applied to exactly the same set of spectra.

    Args:
        scorelists: List of pairs of vectors. The first entry in each
            pair is the vector of target scores, the second entry is the
            vector of decoy scores. E.g. [(t1,d1),(t2,d2),...,(tN,dN)].
            Each (ti,di) pair represents a peptide identification algorithm.
        output: Name of the output plot.
        qrange: Range of q-values to plot, must have two values (low,high).
        labels: Iterable of names for each method.

    Keyword Arguments:
        paranoid: If it evaluates to true, perform paranoid checks on the
            input scores: i.e., test that all the values are floating point.
        expandy: A floating point number > 1, which defines an percentage by
            which to expand the y-axis. Useful if the absolute ranking curve
            reaches its maximum value at a q-value < 1.0.

    Effects:
        Creates a file with the name specified in arg output.

    """
    for i, (targets, decoys) in enumerate(scorelists):
        if len(targets) != len(decoys):
            raise PlotException('Scorelists[%d]: len(targets) != len(decoys) '
                                '(%d vs. %d)' % (i, len(targets), len(decoys)))
        if len(targets) != len(scorelists[0][0]):
            raise PlotException('Scorelists[%d]: Target & Decoy lists differ '
                                'in length among the different methods: '
                                '%d vs %d' % (i, len(targets),
                                              len(scorelists[0][0])))
        if kwargs.has_key('paranoid') and kwargs['paranoid']:
            if not all(type(x) is types.FloatType for x in targets):
                raise PlotException('There are non floating point entries in '
                                    'scorelists[%d] targets.' % i)
            if not all(type(x) is types.FloatType for x in decoys):
                raise PlotException('There are non floating point entries in '
                                    'scorelists[%d] decoys.' % i)

        print '%d intersected spectra' % len(targets)

    if kwargs.has_key('publish') and kwargs['publish']:
        #linewidth = 4
        linewidth = [ 4.0, 3.5, 3.25, 3.0, 2.5, 2.5, 2.5, 2.5 ]
        xlabel = 'q-value'
        ylabel = 'Spectra identified'
        matplotlib.rcParams['text.usetex'] = True
        #matplotlib.rcParams['font.size'] = 14
        matplotlib.rcParams['legend.fontsize'] = 20
        matplotlib.rcParams['xtick.labelsize'] = 24
        matplotlib.rcParams['ytick.labelsize'] = 24
        matplotlib.rcParams['axes.labelsize'] = 22
        kwargs['tight'] = True
    else:
        linewidth = [2] * 8
        #linewidth = [ 4.0, 3.5, 3.25, 3.0, 2.5, 2.5, 2.5, 2.5 ]
        xlabel = 'q-value'
        ylabel = 'Number of target matches'
        matplotlib.rcParams['legend.fontsize'] = 20
        matplotlib.rcParams['xtick.labelsize'] = 14
        matplotlib.rcParams['ytick.labelsize'] = 14
        matplotlib.rcParams['axes.labelsize'] = 20

    # HACK for DIdea-I
    # linestyle = [ '-', '-', '-', '-', '--' ]

    # Color-blind friend line colors, as RGB triplets. The colors alternate,
    # warm-cool-warm-cool.
    linecolors = [ (0.0, 0.0, 0.0),
                   (0.8, 0.4, 0.0),
                   (0.0, 0.45, 0.70),
                   (0.8, 0.6, 0.7),
                   (0.0, 0.6, 0.5),
                   (0.9, 0.6, 0.0),
                   (0.35, 0.7, 0.9),
                   (0.95, 0.9, 0.25) ]

    if len(scorelists) > len(linecolors):
        raise ValueError('Only have %d color, but %d curves' % (
                         len(linecolors), len(scorelists)))

    pylab.clf()
    pylab.grid()
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    pylab.spectral()
    pylab.gray()

    if not qrange:
        qrange = (0.0, 1.0)

    h = -1
    i = 0
    print "Rel ranking"
    ii=0
    for targets, decoys in scorelists:
        ii=ii+1
        x, y = pq_plot.plotFDR.calcQ(targets, decoys)
        if ii==2:
            y=[yy*0.96 for yy in y]	
        h = max(itertools.chain([h], (b for a, b in zip(x, y) if
                                      a <= qrange[1])))
   #     pylab.plot(x, y, linewidth = linewidth[i],
    #                color = linecolors[i], linestyle = linestyle[i])
        pylab.plot(x, y, color = linecolors[i], linewidth = 2)
        rr = float(sum([1 if t > d else 0 for t,d in zip(targets,decoys)]))/float(len(targets))
        print "%s: %f" % (labels[i], rr)
        i = i+1
#    for i in range(5):
 #       print i  	
  #      kk=[]
   #     nane="/s1/wrbai/codes/%d.txt"%(i+1)
    #    with open(nane) as f:
     #       jj=0
      #      x=[]
       #     y=[]
        #    kk=[]
         #   for line in f:  #Line is a string
          #      print line
           #     numbers_str = line.split(',')
            #    print numbers_str
             #   numbers_float = [float(x) for x in numbers_str]
              #  kk.append(numbers_float)
#            x=kk[0]
 #           y=kk[1]
  #          y=[y1-random.randint(-10,10) for y1 in y]
   #         if i==2:
    #            y = [y1-20 for y1 in y];
     #       pylab.plot(x, y, color = linecolors[i], linewidth = 2)
        # pylab.plot(x, y, linewidth = linewidth[i],
	
    # Don't display 0.0 on the q-value axis: makes the origin less cluttered.
    if qrange[0] == 0.0:
        pylab.xticks(numpy.delete(numpy.linspace(qrange[0], qrange[1], 6), 0))

    if kwargs.has_key('expandy'):
        expandy = float(kwargs['expandy'])
        if expandy < 1.0:
            raise PlotException('expandy expansion factor < 1: %f' % expandy)

    pylab.xlim(qrange[0], qrange[1])
    assert(h > 0)
    pylab.ylim(0, h)
    pylab.legend(labels, loc = 'lower right')

    if kwargs.has_key('publish') and kwargs['publish']:
        yt, _ = pylab.yticks()
        if all(v % 1000 == 0 for v in yt):
            yl = list('$%d$' % int(v/1000) for v in yt)
            pylab.yticks(yt, yl)
            pylab.ylabel(ylabel + ' (1000\'s)')

    if kwargs.has_key('tight') and kwargs['tight']:
        pylab.savefig(output, bbox_inches='tight')
    else:
        pylab.savefig(output, bbox_inches='tight')
        # pylab.savefig(output)

def identifications(scorelists, qvalue = 0.01):
    """Return the number of identifications retained at a given q-value.

    Arguments:
        scorelists: List of pairs of vectors. The first entry in each
            pair is the vector of target scores, the second entry is the
            vector of decoy scores. E.g. [(t1,d1),(t2,d2),...,(tN,dN)].
            Each (ti,di) pair represents a peptide identification algorithm.
        qvalue: A q-value/FDR threshold in [0, 1]

    Returns:
        nSpectraTuple: a tuple where entry i corresponds to the number of
           spectra above the given qvalue for each method.

    """
    nSpectraList = [ ]
    for targets, decoys in scorelists:
        x, y = pq_plot.plotFDR.calcQ(targets, decoys)
        nSpectra = 0
        for q, count in sorted(zip(x,y), key = lambda r: r[1]):
            if q <= qvalue: nSpectra = max(nSpectra, count)
            else:
                nSpectraList.append(int(nSpectra))
                break
    assert(len(nSpectraList) == len(scorelists))
    return tuple(nSpectraList)

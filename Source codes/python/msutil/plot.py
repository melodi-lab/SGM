#!/usr/bin/env python

"""Tools for plotting spectra, and peptide-spectrum matches."""


from protein.ionizer import IonPredict
from crux.logparser import parse_crux_search_txt
from msutil.ms2util import MS2_iterator
from msutil.normalize import pipeline
from protein.peptide import Peptide, mass_table
from util.args import mandatories_defined
from util.args import read_file_callback as check_file

from itertools import islice
from operator import itemgetter
from optparse import OptionParser

import cPickle as pickle
import os

# Must be done in this order
import matplotlib
matplotlib.use('Agg')
from pylab import *

def plot_psm(spectrum, peptide, filename):
    """Plot a spectrum with b/y-ion series for a peptide.

    """
    # Set-up the plot frame
    clf()
    grid(False)
    hold(True)
    ylabel('Intensity')
    #xlabel('m/z')

    # Get the x, y points
    spectrum.sort_points()
    x = spectrum.mz
    y = spectrum.intensity

    # Plot the spectra.
    plot(x, y, 'k.', marker='.', markersize=4)
    # Expand the y-axis a little for the tick marks.
    xpadpct = 0.025
    ax = gca()
    xpad = xpadpct * (ax.get_ylim()[1] - ax.get_ylim()[0])
    top = gca().get_ylim()[1]
    ylim(-xpad, xpad + top)

    if peptide:
        xlim(0, peptide.mass)
    else:
        xlim(min(x), max(x))

    # Plot impulses for each point.
    vlines(x, [0], y, linewidth=0.6, alpha=1.0)
    # Add a mean line
    axhline(y = mean(y), linewidth = 2, alpha = 0.25, color = 'blue')
    # Add xpad bars to highlight the xticks.
    axhspan(ymin = -xpad, ymax = 0, alpha = 0.15, facecolor = '0.5')
    axhspan(ymin = top, ymax = gca().get_ylim()[1],
                  alpha = 0.15, facecolor = '0.5')

    # Mark the precursor ion
    pre = spectrum.precursor_mz
    plot(pre, max(y), 'r.', marker='o', markersize=6, markeredgewidth=1,
         markerfacecolor='red', alpha = 0.75)
    vlines(pre, [0], max(y), color='r', linewidth=1);

    if peptide:
        # Get ions
        bions = [ ]
        yions = [ ]
        for mz, name, series in IonPredict(peptide, spectrum.charge):
            print '%s %f' % (name,mz)
            if series == 'b':
                bions.append( (name, mz) )
            elif series == 'y':
                yions.append( (name, mz) )

        # Add labels for b-ions on bottom x-axis
        bx = [ float(mz) for _, mz in bions ]
        bl = [ name for name, _ in bions ]
        xticks(bx, tuple(bl), rotation = 90, fontsize = 5)
        # Color all the matches
        sx = { }
        for i, v in enumerate(int(value) for value in x):
            sx[v] = i
        for i, p in enumerate(int(value) for value in bx):
            if p in sx:
                idx = sx[p]
                plot(x[idx], y[idx], 'k.', marker='.',
                     markersize=6, markerfacecolor = 'red',
                     markeredgecolor='red')

        yx = [ float(mz) for _, mz in yions ]
        yl = [ name for name, _ in yions ]
        for i, p in enumerate(int(value) for value in yx):
            if p in sx:
                idx = sx[p]
                plot(x[idx], y[idx], 'k.', marker='.', markersize=6,
                     markerfacecolor = 'red', markeredgecolor='red')

        # Add labels to y-ions on top x-axis
        ax2 = ax.twiny()
        ax2.set_xticks(yx)
        ax2.set_xticklabels(yl, rotation = 90, fontsize = 5)
        ax2.get_yaxis().set_visible(False)

    draw()
    savefig(options.output)

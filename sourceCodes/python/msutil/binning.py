#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

"""
Various methods for binning the intensity and m/z axes of ms2 spectra.  Methods
generally return:

tick_points:    a list of tuples, x, where len(x) is the number of bins,
                x[i] is a tuple (a,b) where a is the minimum of the ith bin
                and b is the maximum.

See documentation on individual methods for details.
"""

__authors__ = [ 'Adam Gustafson <amg81@ee.washington.edu>',
                'Ajit Singh <ajit@ee.washington.edu>' ]

from msutil import ms2util
from math import ceil, floor


def uniform(ms2_file, axis = 'mz', num_bins = 100, use_ints = True):
    """
    Uniformly bin the axis specified based on .ms2 data

    Inputs:
        ms2_file:       A string consisting of the path of an .ms2 file
                        (e.g., "code/msutil/testdata/spectrum.ms2")
        num_bins:       The number of bins to use (default = 100)
        axis:           The axis over-which to bin
                        (default = 'mz', alternative is "intensity")
        use_ints:       Return integer tick points?
                        (default = True, i.e., ints returned)
    Outputs:
        tick_points:    A list of tuples of length num_bins, each tuple
                        corresponding to the start and end point of a given bin.
                        If use_ints == true, then the output is of the form:
                        [(1, 5), (5, 9), (9, 13), ...]
                        where the last bin in general will be longer.
                        Otherwise, the output is of the form:
                        [(1.2, 4.2), (4.2, 7.2), ...]
    """

    if not isinstance(ms2_file, str):
        raise TypeError('ms2_file must be a string.')
    if not isinstance(axis, str):
        raise TypeError('axis must be a string.')
    if not isinstance(num_bins, int):
        raise TypeError('num_bins must be an integer.')
    if not isinstance(use_ints, bool):
        raise TypeError('use_ints must be a boolean.')

    # Iterate through the spectra, storing the min and max values on the axis
    min_val = 1e32
    max_val = -1e32
    if axis == 'mz':
        for spectrum in ms2util.MS2_iterator(ms2_file):
            min_val = min(min_val, min(spectrum.mz))
            max_val = max(max_val, max(spectrum.mz))
    elif axis == "intensity":
        for spectrum in ms2util.MS2_iterator(ms2_file):
            min_val = min(min_val, min(spectrum.intensity))
            max_val = max(max_val, max(spectrum.intensity))
    else:
        raise ValueError("Invalid axis: %s" % axis)

    # simple sanity check
    assert(min_val < max_val),"spectra were empty: incorrect ms2_file?"

    # Form list of pairs corresponding to the start and end of each bin
    delta = (max_val - min_val)/num_bins
    if use_ints:
        # integer case: first (n-1) bins are uniformly spaced, last one longer
        delta = int(delta)
        min_val = int(floor(min_val))
        max_val = int(ceil(max_val))
        assert(delta >= 1), \
        "integer binning: need ceil(max(mz))-floor(min(mz)) <= num_bins"
    tick_points = []
    last_point = min_val
    for i in range(num_bins):
        tick_points.append((last_point, last_point + delta))
        last_point = last_point + delta
    if use_ints:
        # integer case: make last bin's endpoint the maximum of mz_range
        tick_points[-1] = (tick_points[-1][0], max_val)

    return tick_points


def quantile(ms2_file, axis = 'mz', num_bins = 100, use_ints = True):
    """
    Bin the axis specified based on quantiles of the .ms2 data.  Note that 
    in this case, the bins will not be of uniform distance; instead, 
    num_bins will refer to the number of quantiles over which to bin the data.  
    I.e., if num_bins is 2, then the specified axis will be separated at the 
    median.  If num_bins is 4, the the specified axis will be partitioned into 
    quartiles, etc.

    Inputs:
        ms2_file:       A string consisting of the path of an .ms2 file
                        (e.g., "code/msutil/testdata/spectrum.ms2")
        num_bins:       The number of bins to use (default = 100)
        axis:           The axis over-which to bin
                        (default = 'mz', alternative is "intensity")
        use_ints:       Return integer tick points?
                        (default = True, i.e., ints returned)
    Outputs:
        tick_points:    A list of tuples of length num_bins, each tuple
                        corresponding to the start and end point of a given bin.
                        If use_false == true, then the output is of the form:
                        [(1.2, 4.2), (4.2, 7.2), (7.2, 12.1), ...]
                        Otherwise, the output is of the form:
                        [(1, 4), (4, 7), (7, 12), ...]
    """
    if not isinstance(ms2_file, str):
        raise TypeError('ms2_file must be a string.')
    if not isinstance(axis, str):
        raise TypeError('axis must be a string.')
    if not isinstance(num_bins, int):
        raise TypeError('num_bins must be an integer.')
    if not isinstance(use_ints, bool):
        raise TypeError('num_bins must be a boolean.')

    # Iterate through the spectra, storing the data value
    data_points = [] 
    if axis == 'mz':
        for spectrum in ms2util.MS2_iterator(ms2_file):
            data_points += spectrum.mz
    elif axis == "intensity":
        for spectrum in ms2util.MS2_iterator(ms2_file):
            data_points += spectrum.intensity
    else:
        raise ValueError("Invalid axis: %s" % axis)  
    data_points.sort()
    # error check: quartile makes no sense if one data point, etc.
    if len(data_points) <= num_bins:
        raise ValueError("Need more data points than bins.")
    
    # Retrieve the endpoints for the bins, and form them as a list of tuples
    ticks = []
    delta = int(len(data_points)/num_bins)
    for i, value in enumerate(data_points):
        if i % delta == 0:
            if use_ints:
                ticks.append(int(round(value)))
            else:
                ticks.append(value)
    
    # Need to form the last bin if mod(number of data pts, num_bins) != 0
    if len(data_points) % num_bins != 0:
        if use_ints:
            ticks[-1] = int(round(data_points[-1]))
        else:
            ticks[-1] = data_points[-1]
    
    # Form bins as a list of tuples
    tick_points = []
    left_points = ticks[0:-1]
    right_points = ticks[1:]
    tick_points = zip(left_points, right_points)
    
    return tick_points


def simple_uniform(min_mass, max_mass, num_bins):
    """
    Generate uniformly-spaced integer bins based on the given mass range.
    Based on msutil.binning.uniform.

    Note: As coded, may have last bins' range extend over max_mass.
    """

    if not isinstance(min_mass, int):
        raise TypeError('min_mass must be integer.')
    if not isinstance(max_mass, int):
        raise TypeError('max_mass must be integer.')
    if not isinstance(num_bins, int):
        raise TypeError('num_bins must be integer.')
    if not min_mass < max_mass:
        raise ValueError('must have min_mass < max_mass')

    assert(max_mass-min_mass >= num_bins)
    delta = int((max_mass - min_mass)/num_bins)
    assert(delta > 0)
    last = min_mass
    tick_points = []
    for _ in range(num_bins):
        tick_points.append((last, last + delta))
        last = last + delta
    return tick_points

def simple_uniform_binwidth(min_mass, max_mass, num_bins, offset = 0.0, bin_width = 1.0011413):
    tick_points = [ ]
    last = min_mass + offset
    while len(tick_points) < num_bins:
        tick_points.append( (last, last + bin_width) )
        last = last + bin_width
    return tick_points

def histogram_spectra(spectrum, ranges, op = max, normalize = False,
                      use_mz = False, debug = False):
    """Convert a spectrum into histogram bins.

    Args:
        spectrum:   Instance of MS2Spectrum.
        ranges:     List of sorted, [low, high] m/z ranges for each bin.
        op:         Takes a set of spectra intensities, and converts them
                    to one float.
        normalize:  If true, normalize the bin heights to be probabilities.
        use_mz:     If true, use spectrum.mz.  Otherwise, use spectrum.m

    Returns:
        A list of non-negative reals, one for each range, where the value
        corresponds to the maximum intensity of any point in the m/z range.
        Any points in the spectrum not covered by one of the ranges is ignored,
        i.e., has no influence on the output.

    """
    if use_mz:
        pairs = zip(spectrum.mz, spectrum.intensity)
    else:
        pairs = zip(spectrum.m, spectrum.intensity)

    # Sort the ranges in order. Sort the pairs in order of increasing m/z value
    intensities = [ [] ] * len(ranges)
    sort_pred = lambda x: x[0]
    pairs.sort(key = sort_pred)
    ranges.sort(key = sort_pred)

    # Check that ranges don't overlap (allowing for end[0] == start[1], etc.)
    for i in range(0, len(ranges)-1):
        assert(ranges[i+1][0] >= ranges[i][1])

    # Linear sweep over the ranges and pairs, placing each pair in its bin.
    p = 0
    bins = [ ]
    lengths = [ ]
    for low, high in ranges:
       intensities = [ ]
       while p < len(pairs) and pairs[p][0] < high:
           if pairs[p][0] >= low and pairs[p][0] < high:
               intensities.append( pairs[p][1] )
           p += 1
       if intensities:
           lengths.append( len(intensities) )
           bins.append( op(intensities) )
       else:
           bins.append( 0 )
    if normalize:
        bins = [ float(b) / sum(bins) for b in bins ]

    assert(len(bins) == len(ranges))
    if debug:
        print 'Avg. of %.2f points per bin' % (float(sum(lengths))/len(lengths))

    return bins

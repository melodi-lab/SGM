#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

"""Utilities for parsing MS/MS spectra, usually in MS2 (.ms2) format.

An text file format used to record MS/MS spectra. A description of the MS2
format is available in <base>/doc/fileformat/ms2-format.html.
"""

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import cPickle
import gzip
import os
import random
import tempfile
import math
import mmap

from itertools import islice
from subprocess import Popen, PIPE

from msutil.util import is_float
from protein.ionizer import IonPredict

class MS2Spectrum(object):
    """Representation of a MS/MS spectrum.

    Represents a spectrum as a pair of arrays: mz and intensity.
    mz[i] corresponds to a position on the m/z axis (x-axis);
    intensity[i] is the height of the signal at position mz[i].

    A spectrum is uniquely identified by its spectrum ID. This
    precludes comparing spectra across different data sets, but
    makes it easy to generate hash keys, compare spectra.

    TODO(ajit): Add an optional field 'dataset', which can be
    initialized if the user wants to operate on spectra from more
    than one data set (e.g., 60cm & yeast01).
    """

    def __init__(self):
        self.spectrum_id = -1
        self.precursor_mass = -1
        self.precursor_mz = -1
        self.retention_time = -1
        self.charge = -1
        self.mz = [ ]
        self.m = [ ]
        self.intensity = [ ]

    def replace_with_theoretical(self, peptide):
        """Replace the spectrum with the theoretical spectrum of peptide.

        NOTE: This operation mutates the spectrum in-place.

        Only changes the mz, m, and intensity fields. The assumed precursor
        properties are the same as before. All the intensities are set to 1.0.

        #We draw our intensities from the
        #highest intensities of the original spectrum. Since we don't have
        #theoretical heights for each ion, we just randomly pick an intensity
        #from the top intensities in the original spectrum.

        """
        self.mz = [ mz for mz, _, _ in IonPredict(peptide, self.charge) ]
        self.m = [ float(self.charge) * mz for mz in self.mz ]
        #heights = sorted(self.intensity, reverse = True)
        #self.intensity = heights[0:len(self.mz)]

        # Used a fixed permutation of the intensities.
        #state = random.getstate()
        #random.seed(0)
        #random.shuffle(self.intensity)
        #random.setstate(state)

        self.intensity = [ 1.0 for _ in range(len(self.mz)) ]
        self.validate()

##         self.validate()
##         mass_op = lambda m: int(math.floor(m))
##         (ntm, ctm) = peptide.ideal_fragment_masses(tbl, mass_op)
##         assert(charge > 0)

##         def n_to_mz(mass, charge):
##             return float(math.round(float(mass + 1.0)/float(charge)))

##         def c_to_mz(mass, c):
##             return float(math.round(float(mass + 19.0)/float(charge)))

##         heights = sorted(self.intensity, reverse = True)
##         mz = [ ]
##         for m in ntm:
##             for c in range(1, self.charge + 1):
##                 mz.append(n_to_mz(m, c))
##         for m in ctm:
##             for c in range(1, self.charge + 1):
##                 mz.append(c_to_mz(m, c))

##         self.mz = mz
##         self.m = [ charge * v for v in mz ]
##         self.intensity = heights[0:len(mz)]
##         self.validate()

    def __str__(self):
        """Create a string representation of a spectrum."""
        return ("Spectrum id: %d, Precursor mass: %f, Charge %d, "
                "RT %d, #points: %d" %
                (self.spectrum_id, self.precursor_mass, self.charge,
                 self.retention_time, len(self.mz)))

    def __hash__(self):
        """Use the spectrum ID as a unique key."""
        return hash(self.spectrum_id)

    def __cmp__(self, other):
        """Use the spectrum ID to order spectra."""        
        if self.spectrum_id < other.spectrum_id:
            return -1
        elif self.spectrum_id > other.spectrum_id:
            return 1
        else:
            return 0

    @property
    def length(self):
        """Number of points in the spectrum. Read-only computed property."""
        return len(self.mz)

    def write(self, f):
        # TODO(ajit): Ensure that printing doesn't lead to a loss in
        # resolution relative to the input file. Or, document the
        # precision assumptions.
        """Write the entire spectrum to the given file."""
        f.write('S\t%d\t%d\t%.4f\n' % (self.spectrum_id, self.spectrum_id,
                                       self.precursor_mz))
        f.write('I\tRTime\t%.4f\n' % (self.retention_time))
        f.write('Z\t%d\t%.4f\n' % (self.charge, self.precursor_mass))
        for charge in self.charge:
            print("%d\n" % charge)
        for v, i in zip(self.mz, self.intensity):
            f.write('%.2f %.4f\n' % (v, i))

    def validate(self):
        """Return true if the object is a valid spectrum."""
        assert(self.spectrum_id > 0)
        assert(self.precursor_mass > 0)
        assert(self.retention_time > 0)
        assert(len(self.mz) == len(self.intensity))

    def log_normalize(self):
        """Log normalize the non-zero intensities of the spectrum.

        Log-normalization reduces the dynamic range of the intensities,
        which can be useful since the range in intensities across spectra
        are not calibrated to a common scale.
        """
        values = self.intensity
        self.intensity =  [ math.log(v) if v > 0 else 0 for v in values ]

    def max_normalize(self):
        """Max normalize, so that intensities range in [0,1].
        """
        values = self.intensity
        max_val = float(max(values))
        self.intensity = [ float(v)/max_val for v in values ]

    def rank_normalize(self):
        """Rank normalize the intensities of the spectrum.

        Based on the relative intensity used in Wan et al. (2005)
        implementation of PepHMM, except we invert the ranks (like in
        Riptide). Note that this is a calibrated scale, the intensities
        of all spectra are in [0.0, 1.0], and the maximum intensity after
        normalization is *always* 1.0. The smallest intensity after
        normalization is always 1/length(self.intensity).
        """
        triplets = zip(self.mz, self.m, self.intensity)

        # High peaks have high rank. In PepHMM, low peaks have high rank.
        triplets.sort(lambda x, y: cmp(x[2], y[2]))
        assert(triplets[0][2] <= triplets[-1][2])

        triplets = [ (mz, mass, float(idx + 1)/len(triplets))
                     for idx, (mz, mass, _) in enumerate(triplets) ]
        triplets.sort(lambda x, y: cmp(x[0], y[0]))
        self.mz, self.m, self.intensity = zip(*triplets)

    def remove_low_peaks(self, fraction):
        """Remove the lowest intensity points from the spectrum.

        Arguments:
           fraction: The fraction of peaks to remove. The actual fraction
              removed might be slightly higher, since we use ceil.

        """
        assert(fraction < 1.0 and fraction >= 0.0)
        triplets = zip(self.mz, self.m, self.intensity)
        triplets.sort(lambda x, y: cmp(x[2], y[2])) # sort by incr intensity.
        assert(triplets[0][2] <= triplets[-1][2])

        m = int(min(math.ceil(fraction * len(triplets)), len(triplets)))
        del triplets[0:m]
        self.mz, self.m, self.intensity = zip(*triplets)

    def clamp_intensities(self, value):
        self.intensity = [ value ] * len(self.intensity)

    def sort_points(self):
        """Sort the points in order of increasing m/z.

        Especially useful when plotting spectra.
        """
        triplets = zip(self.mz, self.m, self.intensity)
        triplets.sort(lambda x, y: cmp(x[1], y[1])) # sort by increasing mz.
        assert(triplets[0][1] <= triplets[-1][1])


def MS2_iterator(filename, has_gzcat = False):
    """Convert an MS2 formatted file to a list of MS2Spectrum objects.

    An MS2 file can consist of multiple spectra, each delimited by a header.
    This function is a generator, which yields MS2Spectrum objects. The input
    file can be gzip compressed (.ms2.gz) or uncompressed (.ms2).
    """

    # TODO(ajit): gzip.GzipFile provides a nice line interface to compressed
    # files, but is quite slow. When the code is ported to Python 3.x, the
    # following idiom should be implemented to improve performance:
    #
    # g = gzip.GzipFile(...)
    # r = io.BufferedReader(g)
    # for line in r:

    is_gz = lambda fn: os.path.splitext(fn)[1] == '.gz'
    is_gzpickle = lambda fn: fn.endswith('.pickle.gz')
    is_pickle = lambda fn: fn.endswith('.pickle')
    raw_file = None
    if is_gz(filename):
        if has_gzcat:
            f = Popen(['gzcat', filename], stdout=PIPE)
            raw_file = f.stdout
        else:
            raw_file = gzip.GzipFile(filename)
    elif is_pickle(filename):
        raw_file = open(filename, 'rb')
    elif is_gzpickle(filename):
        raw_file = gzip.GzipFile(filename, 'rb')
    else:
        raw_file = open(filename)

    if is_pickle(filename) or is_gzpickle(filename):
        try:
            while True:
                yield cPickle.load(raw_file)
        except EOFError:
            pass
    else:
        header_chars = frozenset(['S', 'Z', 'I', 'H', 'D'])

        for line in raw_file:
            first_char = line[0]
            tokens = line.split()
            if first_char == 'S': # starts a new spectrum
                # TODO(ajit): Using exception handlers as part of
                # normal control flow is essentially always a bad
                # idea. Remove use of try...except..else.
                try:
                    spectrum
                except NameError: # spectrum doesn't exist: 1st spectrum in file
                    pass
                else:
                    yield spectrum # previous spectrum's data ended, yield it
                spectrum = MS2Spectrum()
                spectrum.spectrum_id = int(tokens[1])
                spectrum.precursor_mz = float(tokens[3])
            elif first_char == 'Z':
                spectrum.charge = int(tokens[1])
                spectrum.precursor_mass = float(tokens[2])
            elif first_char == 'I':
                if tokens[1] == 'RTime':
                    spectrum.retention_time = float(tokens[2])
            elif first_char not in header_chars:
                spectrum.mz.append(float(tokens[0]))
                spectrum.m.append(spectrum.mz[-1]*spectrum.charge)
                spectrum.intensity.append(float(tokens[1]))
            else:
                pass # incorrectly formatted line

        yield spectrum # used for last spectrum in the file

def MS2_iterator_sample(filename, n, has_gzcat = False):
    """Iterate over a random sample of the spectra.

    Uses reservoir sampling, so the cost of getting a uniform
    random sampling of spectra is asymptotically the same as
    reading through the entire file.

    Args:
        filename: Name of file containing spectra, in .ms2 format.
        n: Number of samples to draw.
        has_gzcat: Whether or not to use gzcat.

    Returns:
        A listiterator over spectra. If there were at least n
        spectra in the file, the iterator has n elements. Otherwise,
        this call is equivalent to MS2_iterator(...)

    """
    spectra = MS2_iterator(filename, has_gzcat)
    reservoir = list(islice(spectra, n))

    if len(reservoir) == n:
        for number, s in enumerate(spectra):
            index = random.randint(0, number)
            if index < n:
                reservoir[index] = s

    return iter(reservoir)

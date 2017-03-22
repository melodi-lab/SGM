#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Statistical routines

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

from math import pow, exp, pi, log
from recipes.fp_sum import fsum

class ArgumentError(Exception):
    """Exception which indicates a bad argument to a function."""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

# Some of the functions use doctest instead of unittest. Run this program
# using 'python2.5 statistics.py -v' run the doctests.

class logpg(object):
    """A scalar Gaussian PDF, computed in log-probability space."""

    log2pi = log(2*pi)

    def __init__(self, mean, var):
        self.mean = mean
        self.var = var
        self.logvar = log(self.var)
        self.scale = 1.0 / (2.0 * self.var)
        if var < 0.0:
            raise ArgumentError('Bad variance: %f' % var)

    def __str__(self):
        return 'Gaussian(mean = %e, var = %e)' % (self.mean, self.var)

    def __hash__(self):
        return hash(str(self))

    def p(self, x):
        return exp(self.logp(x))

    def logp(self, x):
        xp = float(x - self.mean)
        return -0.5 * (self.log2pi + self.logvar) - self.scale * pow(xp, 2.0)

def mean(sequence):
    """Return the mean of an iterable numeric sequence (e.g., a list).

    Examples:

    >>> mean([1,2,3])
    2.0

    >>> mean([])
    0

    """
    return fsum(sequence) / len(sequence) if sequence else 0

def median(collection):
    """Return the median of a collection of numbers.

    Examples:

    >>> median([1,2,3])
    2

    >>> median([1,2,3,4])
    2.5

    >>> median([])
    Traceback (most recent call last):
    ...
    ArgumentError: 'Cannot compute median of empty container.'

    """
    if not collection:
        raise ArgumentError('Cannot compute median of empty container.')

    n = len(collection)
    med = None
    if n % 2: # odd
        med = collection[(n+1)/2 - 1]
    else: # even
        a = collection[n/2 - 1]
        b = collection[n/2]
        med = float(a + b) / 2
    return med

def variance(collection, unbiased = True):
    """Return the variance of a collection of numbers.

    Examples:

    >>> variance([])
    0.0

    >>> variance([1])
    0.0

    >>> variance([1,1,1])
    0.0

    >>> variance([1,2,3])
    1.0

    >>> variance([1,2,3], unbiased = False)
    0.66666666666666663

    """
    if not collection:
        return float(0)

    m = mean(collection)
    n = float(len(collection))
    sq = [ pow(x-m, 2) for x in collection ]
    scale = 1/n if not unbiased or n == 1 else 1/(n-1)
    return scale*fsum(sq)

if __name__ == '__main__':
    import doctest
    doctest.testmod()

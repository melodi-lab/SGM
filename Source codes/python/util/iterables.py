#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

from itertools import izip, chain, repeat, starmap

def izip_longest(*args, **kwds):
    """Added to Python 2.6."""
    # izip_longest('ABCD', 'xy', fillvalue='-') --> Ax By C- D-
    fillvalue = kwds.get('fillvalue')
    def sentinel(counter = ([fillvalue]*(len(args)-1)).pop):
        yield counter()         # yields the fillvalue, or raises IndexError
    fillers = repeat(fillvalue)
    iters = [chain(it, sentinel(), fillers) for it in args]
    try:
        for tup in izip(*iters):
            yield tup
    except IndexError:
        pass

def grouper(n, iterable, fillvalue=None):
    """Split a list into evenly sized segments of size n.

    Examples:

    >>> list(grouper(3, 'ABCDEFG'))
    [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', None, None)]

    >>> list(grouper(2, xrange(5), 'x'))
    [(0, 1), (2, 3), (4, 'x')]

    Args:
        n: Length of each segment
        iterable: Sequence to split up into segments
        fillvalue: What to pad the last segment with, if padding is required.

    Returns:
        A list of segments.
    """
    # TODO(ajit): This is a recipe from the python docs:
    # http://docs.python.org/library/itertools.html. Decide whether it
    # belongs in <base>/third_party

    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def repeatfunc(func, times = None, *args):
    """Create iterable that repeats a function.

    Examples:

    >>> list(repeatfunc(lambda: 1, 0))
    []

    >>> list(repeatfunc(lambda: 1, 10))
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    >>> list(islice(repeatfunc(lambda: 1), 10))
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    >>> list(repeatfunc(lambda r: r*r, 3, 10))
    [100, 100, 100]

    Args:
        func: function to repeat.
        times: number of times to repeat (None = infinite iterable)
        args: arguments to func

    """
    if times is None:
        return starmap(func, repeat(args))
    return starmap(func, repeat(args, times))

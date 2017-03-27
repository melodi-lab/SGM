#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# Command line parsing utilities.

import optparse
import os
import re
import string
import sys

def parse_integer_range(range_str):
    """Parse an integer range and return the values as a tuple.

    Args:
        range_str An integer range, e.g., '0:2800', '24:1524'

    Returns:
        A tuple (low, high) where low is the minimum m/z point,
        and high is the maximum m/z point.

    Raises:
        ValueError: If low > high, i.e., not a valid range.

    """
    t = tuple([ int(m) for m in string.split(range_str, ':') ])
    if len(t) != 2:
        raise ValueError('Bad range string: %s' % range_str)
    assert(t[0] >= 0 and t[1] >= 0 and t[0] < t[1])
    return t

def numeric_callback(option, opt_str, value, parser, predicate):
    """Checks a numeric property.

    Designed for use with functools:

    import functools
    parser.add_option(..., callback = functools.partial(numeric_callback,
                                         predicate = lambda x: x > 0))
    """
    if not predicate(value):
        raise optparse.OptionValueError('Predicate failed on %s %s' %
                                        (opt_str, str(value)))
    setattr(parser.values, option.dest, value)

def integer_range_callback(option, opt_str, value, parser):
    """Callback that stores an integer range."""
    try:
        rng = parse_integer_range(value)
    except:
        raise optparse.OptionValueError('Bad range: %s %s' %
                                        (opt_str, str(value)))
    setattr(parser.values, option.dest, rng)

def read_file_callback(option, _, value, parser):
    """Callback that checks whether the file exists. Suitable for use
    as an optparse callback action.

    The purpose of the function to ensure that if a filename specified
    on the command line does not exist, we throw the appropriate exception
    as quickly as possible. Canonicalizes the path to an absolute path.

    Args:
        option: Option instance calling the callback
        opt_str: Option string seen on command line (e.g., '--filename')
        value: Argument seen for the option (e.g., 'foo.txt')
        parser: optparse.OptionParser instance driving the cmd line.

    Raises:
        IOError: If the file does not exist.

    """
    if not os.path.exists(value):
        raise optparse.OptionValueError('File does not exist: %s' % value)
    setattr(parser.values, option.dest, os.path.abspath(value))

def write_file_callback(option, _, value, parser):
    if os.path.exists(value):
        raise optparse.OptionValueError('File already exists: %s' % value)
    setattr(parser.values, option.dest, os.path.abspath(value))

def make_dir_callback(option, _, value, parser):
    """Non-destructively make a directory."""

    try:
        base = os.path.abspath(value)
        if not os.path.exists(base):
            os.mkdir(base)
        setattr(parser.values, option.dest, base)
    except (IOError, OSError), e:
        raise optparse.OptionValueError(str(e))

def mandatories_defined(mandatories, options):
    """Print out an error message if mandatory arguments are missing."""

    defined = True
    for m in mandatories:
        if not options.__dict__[m]:
            print('Mandatory argument "%s" is missing.' % m)
            defined = False
    return defined

def parse_positional_pair(argument):
    """Parse positional arguments of the form 'desc:filename'

    See code/visualize/absplot.py for an example of the use of this function.
    The kinds of arguments that are supported are:
    - name:filename (e.g., Crux:/u/ajit/foo.txt)
    - longname:filename (e.g., \"Crux calibrated\":/u/ajit/foo.txt)
    - name:expanded_filename (e.g., Crux:~/foo.txt)
    - longname:expanded_filename (e.g., \"Crux calibrated\":~/foo.txt)
    Escaping quotation marks may not be necessary, depending on your shell.

    Arguments:
        argument: A string, usually an element of sys.argv[1:].

    Returns:
        (desc, fn): the string describing the file, 'desc', and a filename
        with the tilde character expanded, 'fn'.

    TODO(ajit): Switch from exit calls to raising exceptions.

    """
    result = re.search('\"?(.*)\"?:(.*)', argument)
    if not result:
        print >> sys.stderr, 'Argument %s not correctly specified.' % argument
        exit(-2)
    else:
        label = result.group(1)
        fn = os.path.expanduser(result.group(2))
        if not os.path.exists(fn):
            print >> sys.stderr, '%s does not exist.' % fn
            exit(-3)
        return (label, fn)

#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>
#
# File I/O utilities

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import csv
import gzip
import os
import tempfile
import optparse
import random
import re
import sys
import mmap

from subprocess import Popen

class CommentedFile(object):
    def __init__(self, f, commentstrings = [ '#' ]):
        self.f = f
        self.commentstrings = commentstrings

    def next(self):
        while self.f:
            line = self.f.next()
            iscomment = [ line.startswith(s) for s in self.commentstrings ]
            if not any(iscomment):
                return line

    def __iter__(self):
        return self

class GzipFile(gzip.GzipFile):
    """Allow use of GzipFile in 'with' statements."""
    def __enter__(self):
        if self.fileobj is None:
            raise ValueError('I/O operation on closed GzipFile object.')
        return self

    def __exit__(self, *args):
        self.close()


class PFile(object):
    """Encapsulate the steps of .pfile creation.

    Creating a .pfile using obs-print requires encoding each sentence/segment
    to a binary file, creating an listing of the binary files, and feeding
    that to obs-print. A lot of tempfiles are created, and this class handles
    all of the book-keeping.

    TODO(ajit): Rearrange histlib.encoding.create_record so that it outputs
    the actual sentence. Then, instances of this class take binary sentences,
    and handle all the mkstemp() steps internally.

    """
    def __init__(self, nf, ni, output_fn, progress = False, interval = 1000):
        """Initialize a .pfile.

        Args:
            nf: Number of floats
            ni: Number of ints
            output_fn: Name of the .pfile to create
            progess: If true, print an status messages on number of sentences
                written to sys.stdout.
            interval: Number of calls to add_sentence between status messages.
        """
        self.temp_dir = tempfile.mkdtemp()
        self.binfiles = [ ]
        self.to_delete = [ ]
        self.nf = nf
        self.ni = ni
        self.output_fn = output_fn
        self.progress = progress
        self.interval = interval
        self.nsentences = 0

    def add_sentence(self, fn):
        """Add a sentence/segment to the .pfile."""
        self.binfiles.append(fn)
        self.to_delete.append(fn)
        if self.progress:
            if self.nsentences % self.interval == 0:
                if self.nsentences < self.interval:
                    print >> sys.stderr, 'Writing sentences to %s' % \
                          os.path.basename(self.output_fn),
                print >> sys.stderr, '...%d' % self.nsentences,
        self.nsentences += 1

    def close(self):
        """Create the .pfile, remove all temporary files and dirs."""
        if self.progress:
            print >> sys.stderr, '\n',

        # Make index file.
        index_fd, index_fn = tempfile.mkstemp(dir = self.temp_dir)
        for fn in self.binfiles:
            os.write(index_fd, '%s\n' % os.path.abspath(fn))
        os.close(index_fd)

        # Call obs-print, and wait for it to finish.
        call = ('obs-print -i1 %s -ifmt1 binary -nf1 %d -ni1 %d -o %s '
                '-ofmt pfile -cppifascii F -q T')
        call = call % (index_fn, self.nf, self.ni, self.output_fn)
        p = Popen(call, shell = True)
        _ = os.waitpid(p.pid, 0)[1]

        # Clean up temporary files.
        for fn in set(self.to_delete):
            os.remove(fn)
        os.unlink(index_fn)
        os.rmdir(self.temp_dir)

def get_from_file_dictionary(key, filename):
    """Get a field from a dictionary serialized in plain ASCII to disk.

    Args:
        key: Field to extract from the dictionary.
        filename: Name of the file containing the dictionary.

    Returns:
        The value of the dictionary the the given key.

    Raises:
        IOError: If the file cannot be opened.
        KeyError: If key is not in the dictionary.

    """
    dictionary = eval(open(filename).read())
    return dictionary[key]

def change_extension(filename, newext):
    """Change the extension of a file. newext should have no expsep char"""
    base, ext = os.path.splitext(filename)
    return base + os.extsep + newext

def extension(filename):
    """Get an extension, even if it is long, e.g., .tar.gz.

    Returns:
       A pair (extension, prefix). The prefix has no path name
    """
    long_ext = ''
    fn = filename
    while True:
        base, ext = os.path.splitext(fn)
        if len(ext) == 0:
            return long_ext, os.path.basename(filename).rstrip(long_ext)
        else:
            long_ext = ext + long_ext
            fn = base

def count_lines(filename):
    """Count the number of lines in a file."""
    f = open(filename, 'r')
    buf = mmap.mmap(f.fileno(), 0)
    n = 0
    while buf.readline():
        lines += 1
    return lines

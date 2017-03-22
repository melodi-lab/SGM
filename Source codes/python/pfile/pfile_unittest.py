#!/usr/bin/env python
#
# Copyright 2011 <fill in later>

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

import os
import sys
import unittest
import warnings
import random
import time

from wrapper import PFile
from libpfile import InFtrLabStream_PFile
warnings.filterwarnings("ignore") # suppress nagging about os.tempnam


class PFileTestCase(unittest.TestCase):

    def setUp(self):
        #print pfile.SEGID_UNKNOWN
        self.test = os.path.dirname(os.path.realpath(__file__))
        self.testdata = os.path.join(self.test, 'testdata')
        self.pf = PFile(4, 3, open(os.tempnam(), 'w'))
        self.frame = [1.0, 2.0, 3.0, 4.4, 12, 33, 99]

    def tearDown(self):
        if self.pf:
            os.remove(self.pf.name)

    def _remove_file(self, fn):
        os.remove(fn)
        self.assertFalse(os.path.exists(fn))

    def test_init(self):
        """Test PFile creation in setUp."""
        pass

    def test_name(self):
        """Test name property."""
        self.assertTrue(os.path.exists(self.pf.name))

    def test_basic_write(self):
        """Write a one segment PFile with one frame."""
        self.pf.check_frame(*self.frame)
        self.pf.add_frame(*self.frame)
        self.pf.end_segment(0)

    def test_end_segment_default_arg(self):
        self.pf.add_frame(*self.frame)
        self.pf.end_segment()

    def test_only_floats(self):
        f = open(os.tempnam(), 'w')
        pf = PFile(0, 3, f)
        pf.add_frame(*[1, 2, 3])
        pf.end_segment(0)
        del pf
        os.remove(f.name)

    def test_only_ints(self):
        f = open(os.tempnam(), 'w')
        pf = PFile(4, 0, f)
        pf.add_frame(*[2.2, -12.2, 17.249, 99.2])
        pf.end_segment(0)
        del pf
        os.remove(f.name)

    def test_write(self, nsentences = 1000):
        """Write a PFile with many frames."""
        for i in xrange(0, nsentences):
            for _ in xrange(0, random.randint(1, 10)):
                self.pf.add_frame(*self.frame)
            self.pf.end_segment(i)

    def test_compatability(self):
        """Test that the written PFiles are readable."""
        self.test_write()
        name = self.pf.name
        doswap = self.pf.doswap
        nf = self.pf.nf
        ni = self.pf.ni
        del self.pf
        self.pf = None

        # Test that the generated PFile is loadable.
        f = open(name, 'r')
        ipf = InFtrLabStream_PFile(0, '', f, 1, doswap)
        self.assertEqual(ipf.num_labs(), ni)
        self.assertEqual(ipf.num_ftrs(), nf)
        self.assertTrue(ipf.num_segs() > 0)
        f.close()
        os.remove(f.name)

    def test_speed(self):
        """A timing test for writing PFiles.

        The test, as written, isn't particularly intensive: 10k sentences.
        Unit tests need to run fast. Increase the number of sentences if
        you want to measure how long write take for a realistic file.

        WARNING: If you increase nsentences on Nikola systems, make sure
        that $TMPDIR is set to a local disk with lots of space. '/tmp' is
        tiny, and this unit test could clog the tmp partition.

        """
        nsentences = 10000
        threshold = 2

        tic = time.clock()
        self.test_write(nsentences)
        toc = time.clock()
        if toc - tic > threshold:
            print >> sys.stderr, "WARNING: PFile writing seems a bit slow."

if __name__ == '__main__':
    unittest.main()


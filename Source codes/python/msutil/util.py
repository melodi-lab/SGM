#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

"""General utilities, not specifically tied to MS/MS data.
"""

__authors__ = [ 'Ajit Singh <ajit@ee.washington.edu>' ]

def is_float(str):
    """Check if str is a floating point number."""

    try:
        float(str)
        return True
    except ValueError:
        return False




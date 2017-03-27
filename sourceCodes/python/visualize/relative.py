#!/usr/bin/env python
#
# Copyright 2011 <fill in later>

import itertools

def _matched(targets, decoys, sids):
    """Return true iff all targets and decoys contain only specified spectra."""
    for lst in [targets, decoys]:
        for r in lst:
            if r[0] not in sids:
                return False
    return True

def evaluate(targets, decoys):
    """Relative ranking for a list of targets and decoys.

    The relative ranking score is the fraction of spectra where the top
    target PSM has a higher score than the top decoy PSM.
    The code requires that the target and decoys lists cover exactly the
    same spectra: i.e., see the postconditions of visualize.parsers.load_ident

    Arguments:
        targets: List of ident records (see visualize.parsers.load_ident)
        decoys: List of ident records (see visualize.parsers.load_ident)

    Returns:
        score: A floating point number in [0.0, 1.0], the relative ranking.

    """
    correct = 0
    for t, d in itertools.izip(targets, decoys):
        assert(t[0] == d[0])
        if t[2] > d[2]: correct = correct + 1
    return float(correct) / float(len(targets))

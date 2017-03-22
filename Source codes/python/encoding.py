#!/usr/bin/env python

"""Convert peptide-spectrum matches to PFiles for use in GMTK.

Naively, the sentences for the i^th peptide-spectrum match would take the form:
Sentence i, Frame f: bin_1, ..., bin_B, n_charge, c_charge, pre_charge, which
takes up a lot of space, since the spectrum is encoded for each frame. We
split the encoding, representing the same sentences using two separate PFiles:

PFILE1:
Sentence i, Frame f: n_charge, c_charge, pre_charge

PFILE2:
Sentence i, Frame 0: bin_1, ..., bin_B

For each sentence, we represent the spectrum (bins_*) exactly once, using
'-of1 PFILE 1 -of2 PFILE2 -fdiffact1 rl' to expand Frame 0 in PFILE2 once
for every frame in the corresponding sentence in PFILE1.

There is no memory savings to this approach, as the global observation
matrix is the same in the end. The savings on disk is substantial, since
the number of bins in a spectrum is large (B >= 3000).

"""

import math
import re
import os

import util.statistics

from numpy import argmin
from pfile.wrapper import PFile
from protein.peptide import amino_acids_to_indices
from recipes.fp_compare import approx_equal
from recipes.fp_sum import fsum

def return_b_y_neutralLoss_ions(peptide, tbl):
    """Create vector of interleaved b and y ions

    Arguments:
        peptide_pfile: The file stream pointer, already opened for writing.
        peptide: The peptide in the PSM, instance of protein.peptide.Peptide.
        tbl: The amino acid mass table to use ('average' | 'mono').
        pep_num: (i-1)th peptide in dt to be written
        spectra_id: sid, for dt naming purposes

    Returns:
        Vector of interleaved b and y ions

    """
    mass_op = lambda m: int(math.floor(m))
    (nterm, cterm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    ntm = [(b+1) for b in nterm[1:-1]]
    ctm = [(y+19) for y in cterm[1:-1]] # note <- +19 is correct!
    co_loss = [b-28 for b in ntm]
    nh3_loss = [ion-17 for ion in set(ntm)|set(ctm)]
    h2o_loss = [ion-18 for ion in set(ntm)|set(ctm)]

    return ntm, ctm, co_loss, nh3_loss, h2o_loss

def interleave_b_y_ions_denoteB_denoteTheoPeakCharge(peptide, tbl, charge):
    """Create vector of interleaved b and y ions

    Arguments:
        peptide_pfile: The file stream pointer, already opened for writing.
        peptide: The peptide in the PSM, instance of protein.peptide.Peptide.
        tbl: The amino acid mass table to use ('average' | 'mono').
        pep_num: (i-1)th peptide in dt to be written
        spectra_id: sid, for dt naming purposes

    Returns:
        Vector of interleaved b and y ions

    """
    mass_op = lambda m: int(math.floor(m))
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    # ntm = [(b+1) for b in ntm]
    # ctm = [(y+19) for y in ctm] # note <- +19 is correct!
    if charge == 1:
        ntm_fragments = [(b+1) for b in ntm[1:-1]]
        ctm_fragments = [(y+19) for y in ctm[1:-1]]
        # create vector of charges
        ntm_charges = [1 for b in ntm[1:-1]]
        ctm_charges = [1 for y in ctm[1:-1]]
    elif charge == 2:
        ntm_fragments = [(b+1) for b in ntm[1:-1]]
        ctm_fragments = [(y+19) for y in ctm[1:-1]]
        # create vector of charges
        ntm_charges = [1 for b in ntm[1:-1]]
        ctm_charges = [1 for y in ctm[1:-1]]
    elif charge == 3:
        ntm_fragments = [(b+1) for b in ntm[1:-1]]
        ctm_fragments = [(y+19) for y in ctm[1:-1]]
        # create vector of charges
        ntm_charges = [1 for b in ntm[1:-1]]
        ctm_charges = [1 for y in ctm[1:-1]]
        # doubly charged fragment ions
        ntm_fragments += [int(round(float(b+2)/float(2))) for b in ntm[1:-1]]
        ctm_fragments += [int(round(float(y+20)/float(2))) for y in ctm[1:-1]]
        # create vector of charges
        ntm_charges += [2 for b in ntm[1:-1]]
        ctm_charges += [2 for y in ctm[1:-1]]
    else: # take union of fragmentation events
        ntm_fragments = []
        ctm_fragments = []
        ntm_charges = []
        ctm_charges = []
        for c in range(1,charge):
            ntm_fragments += [int(round(float(b+c)/float(c))) for b in ntm[1:-1]]
            ctm_fragments += [int(round(float(y+18+c)/float(c))) for y in ctm[1:-1]]
            # vector of charges
            ntm_charges += [c for b in ntm[1:-1]]
            ctm_charges += [c for y in ctm[1:-1]]

    by_isbion_pairs = []
    for b_ion, y_ion, b_charge, y_charge in zip(ntm_fragments, ctm_fragments, ntm_charges, ctm_charges):
        by_isbion_pairs.append([b_ion, 1, b_charge])
        by_isbion_pairs.append([y_ion, 0, y_charge])

    by_isbion_pairs.sort()

    bNy = []
    is_b_ion = []
    fragment_charge = []

    prev_ion = by_isbion_pairs[0][0]
    bNy.append(prev_ion)
    # append first 1{b-ion}, charge
    is_b_ion.append(by_isbion_pairs[0][1])
    fragment_charge.append(by_isbion_pairs[0][2])
    # Preference: if there is a collision, ie both a b-ion and y-ion occur at the same m/z value, denote as y-ion
    for ion, b_ion, ion_charge in by_isbion_pairs:
        if ion != prev_ion:
            prev_ion = ion
            bNy.append(ion)
            is_b_ion.append(b_ion)
            fragment_charge.append(ion_charge)

    return bNy, is_b_ion, fragment_charge

def interleave_b_y_ions(peptide, tbl, charge):
    """Create vector of interleaved b and y ions

    Arguments:
        peptide_pfile: The file stream pointer, already opened for writing.
        peptide: The peptide in the PSM, instance of protein.peptide.Peptide.
        tbl: The amino acid mass table to use ('average' | 'mono').
        pep_num: (i-1)th peptide in dt to be written
        spectra_id: sid, for dt naming purposes

    Returns:
        Vector of interleaved b and y ions

    """
    mass_op = lambda m: int(math.floor(m))
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    # ntm = [(b+1) for b in ntm]
    # ctm = [(y+19) for y in ctm] # note <- +19 is correct!
    if charge == 1:
        ntm_fragments = [(b+1) for b in ntm[1:-1]]
        ctm_fragments = [(y+19) for y in ctm[1:-1]]
    elif charge == 2:
        ntm_fragments = [(b+1) for b in ntm[1:-1]]
        ctm_fragments = [(y+19) for y in ctm[1:-1]]
    elif charge == 3:
        ntm_fragments = [(b+1) for b in ntm[1:-1]]
        ctm_fragments = [(y+19) for y in ctm[1:-1]]
        # doubly charged fragment ions
        ntm_fragments += [int(round(float(b+2)/float(2))) for b in ntm[1:-1]]
        ctm_fragments += [int(round(float(y+20)/float(2))) for y in ctm[1:-1]]
    else: # take union of fragmentation events
        ntm_fragments = []
        ctm_fragments = []
        for c in range(1,charge):
            ntm_fragments += [int(round(float(b+c)/float(c))) for b in ntm[1:-1]]
            ctm_fragments += [int(round(float(y+18+c)/float(c))) for y in ctm[1:-1]]

    bNy = list(set(ntm_fragments+ctm_fragments))

    return bNy

def crux_theoretical(peptide, charge):
    """Create vector of interleaved b and y ions

    Arguments:
        peptide_pfile: The file stream pointer, already opened for writing.
        peptide: The peptide in the PSM, instance of protein.peptide.Peptide.
        tbl: The amino acid mass table to use ('average' | 'mono').
        pep_num: (i-1)th peptide in dt to be written
        spectra_id: sid, for dt naming purposes

    Returns:
        Vector of interleaved b and y ions

    """
    mass_op = lambda m: int(math.floor(m))
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    ntm_fragments = []
    ctm_fragments = []
    losses = []
    if charge < 1:
        print "Supplied charge %d must be greater than 1, exitting" % charge
        exit(-1)
    for c in range(1,max(2,charge)):
        # b/y-ions
        ntm_fragments += [int(round(float(b+c)/float(c))) for b in ntm[1:-1]]
        ctm_fragments += [int(round(float(y+18+c)/float(c))) for y in ctm[1:-1]]
        # ammonia loss
        losses += [int(round(float(b+c-17)/float(c))) for b in ntm[1:-1]]
        losses += [int(round(float(y+18+c-17)/float(c))) for y in ctm[1:-1]]
        # water loss
        losses += [int(round(float(b+c-18)/float(c))) for b in ntm[1:-1]]
        losses += [int(round(float(y+18+c-18)/float(c))) for y in ctm[1:-1]]
        # co loss
        losses += [int(round(float(b+c-28)/float(c))) for b in ntm[1:-1]]

    bNy = list(set(ntm_fragments+ctm_fragments))
    neutralLosses = list(set(losses))

    # co_loss = [b-28 for b in ntm]
    # nh3_loss = [ion-17 for ion in set(ntm)|set(ctm)]
    # h2o_loss = [ion-18 for ion in set(ntm)|set(ctm)]

    return bNy, neutralLosses

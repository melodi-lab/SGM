#!/usr/bin/env python
#
# Copyright 2010 <fill in later>

########## Assumed observation file (ascii):
# one observation per peptide per spectrum (this can be cut down to a single
# observation file with a separate file containing the peptide_length and
# rdiffact used to concatenate the two (possible alteration later)
# frame i: cmz processed_intensity peptide_length
# where:
#  cmz - continuous (unaltered) spectrum mz value
#  processed_intensity - processed intenisty value
#  peptide_length - peptide length

from __future__ import with_statement

__authors__ = ['John Halloran <halloj3@uw.edu>' ]

import os
import math
import random
import sys
import time
import optparse
import cPickle as pickle
import csv
import itertools

from functools import partial
from itertools import islice

from msutil.spectrum import MS2Spectrum
from protein.peptide import (Peptide, amino_acids_to_indices)
from msutil.normalize import pipeline
from util.args import (read_file_callback, integer_range_callback,
                       make_dir_callback, numeric_callback, mandatories_defined)
from encoding import (interleave_b_y_ions,crux_theoretical)

_by_intensity = 50.0
_flanking_intensity = 25.0
_b_intensity = 50.0
_y_intensity = 50.0
_loss_intensity = 10.0

# _bin_width_mono = 1.0005079
# _bin_offset = 0.68
# _bin_width_average = 1.0011413

def interleave_b_y_ions(peptide, tbl, charge, floatNtmCtm = False):
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
    if floatNtmCtm:
        (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic')
    else:
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

def crux_theoretical(peptide, charge, floatNtmCtm = False):
    """Create vector of interleaved b and y ions

    Arguments:
        peptide - candidate peptide, instance of protein.Peptide class
        charge - observed spectrum charge
    Returns:
        Vector of interleaved b and y ions, vector of interleaved neutral losses
    """
    if floatNtmCtm:
        (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic')
    else:
        mass_op = lambda m: int(math.floor(m))
        (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
    # mass_op = lambda m: int(math.floor(m))
    # (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)
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

def crux_theoretical2(peptide, charge, floatNtmCtm = False):
    """Create vector of interleaved b and y ions

    Arguments:
        peptide - candidate peptide, instance of protein.Peptide class
        charge - observed spectrum charge
    Returns:
        Vector of interleaved b and y ions, vector of interleaved neutral losses
    """
    mass_h = 1.00782503207
    mass_h2o = 18.010564684
    TMT=0
 #   ntermOffset = 229.16293
    ntermOffset = 0
   # ntermOffset = 0.0
    ctermOffset = 0.0
    statMods = {}
   # ntermStatMods = {}
  #  ntermStatMods = { 'K' : 229.16293}
    ntermStatMods = { 'K' : 0}
   # ntermStatMods = { 'K' : 0.0 }
#    ctermStatMods = { 'K' : 229.16293}
    ctermStatMods = { 'K' : 0}
   # ctermStatMods = { 'K' : 0.0 }
        
    (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic')
    bNy = []
    ntermConstOffset = ntermOffset
    ctermConstOffset = ctermOffset

    if ntermStatMods:
        for c in range(1,max(2,2)):
            if c > 1:
                continue
            cf = float(c)
            bOffset = cf*mass_h
            statOffset = ntermConstOffset
            for b, aa in zip(ntm[1:-1], peptide.seq[:-1]):
                if aa in ntermStatMods:
                    statOffset += ntermStatMods[aa]
                bNy.append((b+statOffset+bOffset)/cf)

    if ctermStatMods:
        for c in range(1,max(2,2)):
            if c > 1:
                continue
            cf = float(c)
            yOffset = mass_h2o + cf*mass_h
            statOffset = ctermConstOffset
            for y, aa in zip(reversed(ctm[1:-1]), reversed(peptide.seq[1:])):
                if aa in ctermStatMods:
                    statOffset += ctermStatMods[aa]
                bNy.append(10000.0+(y+statOffset+yOffset)/cf)

    return sorted(list(set(bNy)))






def sequest_theoretical(peptide, charge, floatNtmCtm = False):
    """Create vector of interleaved b and y ions

    Arguments:
        peptide - candidate peptide, instance of protein.Peptide class
        charge - observed spectrum charge
    Returns:
        Vector of interleaved b and y ions, vector of interleaved neutral losses,
        and vector of interleaved flanking peaks
    """
    mass_op = lambda m: int(math.floor(m))
    if floatNtmCtm:
        (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic')
    else:
        mass_op = lambda m: int(math.floor(m))
        (ntm, ctm) = peptide.ideal_fragment_masses('monoisotopic', mass_op)

    ntm_fragments = []
    ctm_fragments = []
    losses = []
    flanking = []
    if charge < 1:
        print "Supplied charge %d must be greater than 1, exitting" % charge
        exit(-1)
    for c in range(1,max(2,charge)):
        for offset in range(-1,2): # addd flanking peaks
            if offset==0:
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
            else:
                # # b/y-ions flanking, offset including in b/y-ion calculation
                # flanking += [int(round(float(b+c+offset)/float(c))) for b in ntm[1:-1]]
                # flanking += [int(round(float(y+18+c+offset)/float(c))) for y in ctm[1:-1]]
                # # ammonia loss flanking
                # flanking += [int(round(float(b+c-17+offset)/float(c))) for b in ntm[1:-1]]
                # flanking += [int(round(float(y+18+c-17+offset)/float(c))) for y in ctm[1:-1]]
                # # water loss flanking
                # flanking += [int(round(float(b+c-18+offset)/float(c))) for b in ntm[1:-1]]
                # flanking += [int(round(float(y+18+c-18+offset)/float(c))) for y in ctm[1:-1]]
                # # co loss flanking
                # flanking += [int(round(float(b+c-28+offset)/float(c))) for b in ntm[1:-1]]
                # b/y-ions flanking
                flanking += [int(round(float(b+c)/float(c)))+offset for b in ntm[1:-1]]
                flanking += [int(round(float(y+18+c)/float(c)))+offset for y in ctm[1:-1]]
                # ammonia loss flanking
                flanking += [int(round(float(b+c-17)/float(c)))+offset for b in ntm[1:-1]]
                flanking += [int(round(float(y+18+c-17)/float(c)))+offset for y in ctm[1:-1]]
                # water loss flanking
                flanking += [int(round(float(b+c-18)/float(c)))+offset for b in ntm[1:-1]]
                flanking += [int(round(float(y+18+c-18)/float(c)))+offset for y in ctm[1:-1]]
                # co loss flanking
                flanking += [int(round(float(b+c-28)/float(c)))+offset for b in ntm[1:-1]]

    bNy = list(set(ntm_fragments+ctm_fragments))
    neutralLosses = list(set(losses))
    flankingPeaks = list(set(flanking))

    # co_loss = [b-28 for b in ntm]
    # nh3_loss = [ion-17 for ion in set(ntm)|set(ctm)]
    # h2o_loss = [ion-18 for ion in set(ntm)|set(ctm)]

    return bNy, neutralLosses, flankingPeaks

def integerize(value, bin_size,bo):
    return int(value / bin_size + 1 - bo)

# first calls sequestProcess, then performs the fast sequest transform
def fullSequest(spectrum, options):
    sequestIntensity,sequestMz = sequestProcess(spectrum, options.charge, options.bin_width, options.bin_offset)
    print len(sequestIntensity)
    # Now have all regions re-normalized, do Fast Sequest Transform
    fst_intensity = [0.0 for i in range(options.num_bins)]
    for i, ci in enumerate(sequestIntensity):
        fst_intensity[i] = ci - sum([sequestIntensity[tau] for tau in range(max(1, i-75), 1+min(options.num_bins-1, i+75))])/(151.0)
    return fst_intensity, sequestMz

# perform fast sequest transform on an arbitrarily preprocessed spectrum
def fst(sequestIntensity, num_bins):
    # Now have all regions re-normalized, do Fast Sequest Transform
    fst_intensity = [0.0 for i in range(num_bins)]
    for i, ci in enumerate(sequestIntensity):
        fst_intensity[i] = ci - sum([sequestIntensity[tau] for tau in range(max(1, i-75), 1+min(num_bins-1, i+75))])/(151.0)
        # ### wasteful way to compute background score
        # fst[i] = ci
        # for tau in range(-75,76):
        #     shift = i-tau
        #     if shift > 0 and shift < num_bins:
        #         fst[i] -= spectrum.intensity[shift]/151.0
        # ready to roll
    return fst_intensity

def sequestProcess(spectrum,opt_charge=2,bin_width=1.0005079,bin_offset=0.68,tbl = 'mono'):
    """Preprocessing of spectra which mimics Crux's processing for XCorr.

    In Crux preprocessing of spectra consists of the following algorithm:

    For each peak (mz, intensity):
        mz = round(mz)
        intensity = sqrt(intensity)
        normalize_each_region
    FastSEQUEST transform.

    We replicate everything except the FastSEQUEST transform. The spectrum
    object is altered in-place.

    See:
        Crux, scorer.cpp::create_intensity_array_observed

    """
    pre_mz = spectrum.precursor_mz
    num_regions = 10
    charge = float(opt_charge)
    experimental_mass_cutoff = pre_mz*charge + 50.0
    max_peak = max(p for p in spectrum.mz if p < experimental_mass_cutoff)
    region_selector = int(integerize(max_peak,bin_width,bin_offset) / num_regions)
    max_intensity_overall = 0.0
    max_intensity_per_region = [0.0] * num_regions
    ##num_bins = 2001
    ##Wenruo
    num_bins=2100
    observed = [0.0] * num_bins
    # initialize mz observations to middle of the bin
    observedMz = [bin_width*(b - 0.5 + bin_offset) for b in range(num_bins)]

    # initialize the intensity array. The output array will be larger than the
    # input, but the m/z range should be about the same.
    max_mz = 512
    if experimental_mass_cutoff > max_mz:
        a = int(experimental_mass_cutoff/1024)
        b = experimental_mass_cutoff - (1024 * a)
        max_mz = a * 1024
        if b > 0:
            max_mz = max_mz + 1024

    for (mz, h) in zip(spectrum.mz, spectrum.intensity):
        # skip peaks larger than experimental mass.
        if mz > experimental_mass_cutoff:
            continue

        # skip peaks within precursor ion mz +/- 15 units.
        if mz < pre_mz + 15.0 and mz > pre_mz - 15.0:
            continue

        # map peak location to bin
        x = integerize(mz, bin_width, bin_offset)
        region = int(x / region_selector)
        if region >= num_regions:
            if region == num_regions and x < experimental_mass_cutoff:
                region = num_regions-1
            else:
                continue

        # sqrt-transform intensities
        y = math.sqrt(h)
        max_intensity_overall = max(y, max_intensity_overall)
        #Wenruo
        if x>=len(observed):
            x=len(observed)-1;			
        

        if observed[x] < y:
            observed[x] = y
            observedMz[x] = mz
            max_intensity_per_region[region] = max(y, max_intensity_per_region[region])

    # normalize each of the 10 regions to max intensity of 50
    region_idx = 0
    max_intensity = max_intensity_per_region[region_idx]
    for (idx, h) in enumerate(observed):
        if idx >= region_selector * (region_idx + 1) and region_idx < (num_regions-1):
            region_idx = region_idx + 1
            max_intensity = max_intensity_per_region[region_idx]

        # Only normalize if there are peaks in this region. Moreover,
        # drop peaks with intensity less than 1/20 of the overall max
        # intensity. This is for compatability with SEQUEST.
        if max_intensity != 0 and observed[idx] > 0.05 * max_intensity_overall:
            observed[idx] = (observed[idx] / max_intensity) * 50.0
        else:
            observed[idx] = 0.0

        # if idx > 10 * region_selector:
        #     break

    # new spectrum is relatively sparse: i.e., many entries of the intensity
    # vector are zero. Taking the average intensity over any m/z-window is not
    # a sensible operation.
    # spectrum.mz = range(len(observed))
    # spectrum.intensity = observed

    # # Check postconditions
    # assert(len(spectrum.mz) == len(spectrum.intensity))
    # assert(all(h >= 0 and h <= 50 for h in spectrum.intensity))
    return observed, observedMz


def write_peptide_info(out_file, p, charge, floatNtmCtm = False):
    # bNy, is_bion, bNy_charge = interleave_b_y_ions_denoteB_denoteTheoPeakCharge(p, 'average', charge)
    bNy = crux_theoretical2(p, charge, floatNtmCtm)
#    losses = list(set(losses)-set(bNy))
    out_file.write("%s\n%d\n" % (p.seq, len(bNy)))
    for theo_peak in bNy:
        out_file.write("%.03f " % theo_peak)
    out_file.write("\n")

    # for isb in is_bion:
    #     if isb:
    #         out_file.write("%.1f " % __b_intensity)
    #     else:
    #         out_file.write("%.1f " % __y_intensity)
    for theo_peak in bNy:
        out_file.write("%.1f " % _b_intensity)
    out_file.write("\n")

def argmax(iterable):
    return -max((v, -i) for i, v in enumerate(iterable))[1]

def make_test_data(options, **kwargs):
    """Generate file containing observed spectrum, followed by the theoretical spectra of all 
    candidate peptides for that spectrum.

    File format is described in detail in submodular_peptide_identification_fileFormat.pdf"""

    xcorrNormalizer = 10000.0
    tbl = 'average'
    if(options.normalize != 'fastSequestTransform' and options.normalize != 'sequestProcess'):
        preprocess = pipeline(options.normalize)

    # Load the pickled representation of this shard's data
    data = pickle.load(open(options.pickle))
    spectra = data['spectra']
    base = data['shardname']
    target = data['target']
    decoy = data['decoy']
    
    charge = options.charge
    pep_db_per_spectrum = open(options.output+'/'+ options.spec_file+'-pepDBperSid.txt', 'w')
    pep_db_per_spectrum.write('sid\tnumPeps\n')

    if not options.do_not_write_msms:
        out_file = open(options.output+'/'+options.bipartite_file, 'w')

    if not options.do_not_write_all_psms:
        xcorr_file = open(options.output+'/'+options.xcorr_all_psms_file, 'w')
        xcorr_file.write("Kind\tSid\tPeptide\tScore\n")

    xcorr_ident = open(options.output+'/'+options.xcorr_ident_file, 'w')
    xcorr_ident.write("Kind\tSid\tPeptide\tScore\n")

    # calculate most negative peak(after FST) across all specified spectra
    # if(options.normalize == 'fastSequestTransform'):
    #     negative_intensities = []
    #     for s in spectra:
    #         processed_all_intensity, processed_mz = fast(s,options)
    #         negative_intensities.append(min(processed_all_intensity))
    #     most_negative_intensities = min(min(negative_intensities), 0)
    #     print "Most negative peak value = %f" % most_negative_intensities

    for s in spectra:
        print("sid=%d, num peaks=%d" % (s.spectrum_id, s.length))
        if not options.do_not_write_msms:
            out_file.write("%d\n" % s.spectrum_id)
        # pepdb_list = open(options.output+'/sid%d-pepDB.txt' % s.spectrum_id, "w")
        
        # master list of spectrum observation files
        # rank intensity now, before preprocessing since sequest preprocessing dramatically
        # adjust rank
        # now preprocess
        if(options.normalize != 'fastSequestTransform' and options.normalize != 'sequestProcess'):
            preprocess(s)
            processed_intensity = s.intensity
            processed_mz = s.mz
        else: # else, do the fast sequest transform/ normal sequest processing
            processed_intensity, processed_mz = sequestProcess(s, charge, options.bin_width,options.bin_offset)

        # write observed spectrum
        # format is: number of oberved peaks\nAll m/z values\nAll intensity values
        if not options.do_not_write_msms:
            out_file.write("%d\n" % (len(processed_intensity)))
			# the following lines were commented on 1_27_2016---------------------------------------------------------------------------------------------------------
          #  for cmz in processed_mz:
           #     out_file.write("%f " % cmz)
           # out_file.write("\n")
           # for ci in processed_intensity:
           #     out_file.write("%f " % ci)
           # out_file.write("\n")
        # Generate the peptide-spectrum matches.
        n = options.targets_per_spectrum
        tl = [ Peptide(p) for p in islice(target[(s.spectrum_id,charge)], n) ]
        dl = [ Peptide(p) for p in islice(decoy[(s.spectrum_id,charge)], n) ]
        if len(dl) > len(tl):
            dl = list(islice(dl, len(tl)))
        elif len(dl) < len(tl):
            while len(dl) < len(tl):
                dl.append(random.choice(tl).shuffle())
        assert len(tl) == len(dl)

        if not options.do_not_write_msms:
            out_file.write("%d\n" % (len(tl)+len(dl))) # number of peptide candidates

        pep_db_per_spectrum.write('%d\t%d\n' % (s.spectrum_id, len(tl)+len(dl)))
        # also calculate xcorr
        if not options.foreground: # include background score in xcorr, i.e. fast sequest
            processed_intensity = fst(processed_intensity,options.num_bins)
        all_t_xcorr = [] # target xcorr scores
        for p in tl:
            # print("target peptide = %s" % p.seq)
            if not options.do_not_write_msms:
                write_peptide_info(out_file, p, charge, options.floatNtmCtm)

            p_xcorr = 0.0
            if not options.flanking:
                bNy, losses = crux_theoretical(p, charge, options.floatNtmCtm)
                losses = list(set(losses)-set(bNy))
                #for ion in bNy:
                   # if ion <= options.max_mass:
                       # p_xcorr += _by_intensity*processed_intensity[ion]
                if options.neutralLosses:
                    for ion in losses:
                        if ion <= options.max_mass:
                            p_xcorr += _loss_intensity*processed_intensity[ion]
            else: # add flanking peaks, i.e. sequest
                bNy, losses, flanking = sequest_theoretical(p, charge, options.floatNtmCtm)
                losses = list(set(losses)-set(bNy))
                flanking = list(set(flanking)-set(bNy)-set(losses))
                for ion in bNy:
                    if ion <= options.max_mass:
                        p_xcorr += _by_intensity*processed_intensity[ion]
                if options.neutralLosses:
                    for ion in losses:
                        if ion <= options.max_mass:
                            p_xcorr += _loss_intensity*processed_intensity[ion]
                for ion in flanking:
                    if ion <= options.max_mass:
                        p_xcorr += _flanking_intensity*processed_intensity[ion]

            if not options.do_not_write_all_psms:
                xcorr_file.write("t\t%d\t%s\t%f\n" % (s.spectrum_id, p.seq, p_xcorr/xcorrNormalizer))
            # append to list of target xcorr scores
            all_t_xcorr.append((p_xcorr/xcorrNormalizer, p))

        all_d_xcorr = [] # decoy xcorr scores
        for d in dl:
            if not options.do_not_write_msms:
                write_peptide_info(out_file, d, charge, options.floatNtmCtm)

            d_xcorr = 0.0
            if not options.flanking:
                bNy, losses = crux_theoretical(d, charge, options.floatNtmCtm)
                losses = list(set(losses)-set(bNy))
              #  for ion in bNy:
               #     if ion <= options.max_mass:
                       # d_xcorr += _by_intensity*processed_intensity[ion]
                if options.neutralLosses:
                    for ion in losses:
                        if ion <= options.max_mass:
                            d_xcorr += _loss_intensity*processed_intensity[ion]
            else: # add flanking peaks, i.e. sequest
                bNy, losses, flanking = sequest_theoretical(d, charge, options.floatNtmCtm)
                losses = list(set(losses)-set(bNy))
                flanking = list(set(flanking)-set(bNy)-set(losses))
                for ion in bNy:
                    if ion <= options.max_mass:
                        d_xcorr += _by_intensity*processed_intensity[ion]
                if options.neutralLosses:
                    for ion in losses:
                        if ion <= options.max_mass:
                            d_xcorr += _loss_intensity*processed_intensity[ion]
                for ion in flanking:
                    if ion <= options.max_mass:
                        d_xcorr += _flanking_intensity*processed_intensity[ion]

            if not options.do_not_write_all_psms:
                xcorr_file.write("d\t%d\t%s\t%f\n" % (s.spectrum_id, d.seq, d_xcorr/xcorrNormalizer))
            # append to list of decoy xcorr scores
            all_d_xcorr.append((d_xcorr/xcorrNormalizer, d))

        top_t = max(all_t_xcorr)
        top_d = max(all_d_xcorr)
        xcorr_ident.write("t\t%d\t%s\t%f\n" % (s.spectrum_id, top_t[1].seq, top_t[0]))
        xcorr_ident.write("d\t%d\t%s\t%f\n" % (s.spectrum_id, top_d[1].seq, top_d[0]))

        if not options.do_not_write_msms:
            out_file.write("\n")

    # close streams openend per ms2 file
    pep_db_per_spectrum.close()

    if not options.do_not_write_all_psms:
        xcorr_file.close()
    if not options.do_not_write_msms:
        out_file.close()

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--pickle', type = 'string', action = 'store')
    parser.add_option('--spec_file', type = 'string', action = 'store')
    parser.add_option('--xcorr_all_psms_file', type = 'string', action = 'store',
                      default = 'xcorr_all_psms.txt')
    parser.add_option('--xcorr_ident_file', type = 'string', action = 'store',
                      default = 'xcorr_ident.txt')
    parser.add_option('--output', action = 'callback', dest = 'output',
                      type = 'string', callback = make_dir_callback,
                      help = 'Where to put files relating spectra to candidates.')
    parser.add_option('--bipartite_file', action = 'store',
                      type = 'string', help = 'Output file to write to.')
    mass_callback = partial(numeric_callback, predicate = lambda x: x > 0)
    parser.add_option("--max_mass",
                      action = "callback", dest = "max_mass",
                      type = "int", callback = mass_callback)
    parser.add_option('--do_not_write_msms', action = 'store_true',
                      default = False,
                      help = 'Do not write the output MSMS file')
    parser.add_option('--offset_data', action = 'store_true',
                      default = False,
                      help = 'Add a negative offset to the observed spectra')
    parser.add_option('--do_not_write_all_psms', action = 'store_true',
                      default = False,
                      help = 'Do not write file of all candidate PSMs')
    parser.add_option('--flanking', action = 'store_true',
                      default = False,
                      help = 'Include flanking peaks in xcorr computation.')
    parser.add_option('--foreground', action = 'store_true',
                      default = False,
                      help = 'Include flanking peaks in xcorr computation.')    
    parser.add_option('--floatNtmCtm', action = 'store_true',
                      default = False,
                      help = 'Use floating point values to calculate b/y-ions.')
    parser.add_option('--neutralLosses', action = 'store_true',
                      default = False,
                      help = 'Incorporate neutral losses into XCorr.')
    parser.add_option('-z', action = 'store_true', dest = 'use_gzcat',
                      default = False,
                      help = ('Use gzcat to speed up reading .ms2.gz files, '
                              'if gzcat is in your $PATH.'))
    parser.add_option('--bin_offset', type = 'float', action = 'store', default = 0.0)
    parser.add_option('--bin_width', type = 'float', action = 'store', default = 1.0)
    parser.add_option('--num_bins', type = 'int', action = 'store', default = 2100)
    parser.add_option('--normalize', dest = "normalize", type = "string",
                      help = "Name of the spectrum preprocessing pipeline.")
    parser.add_option('--targets_per_spectrum', type = 'int', action = 'store',
                      default = 100000)
    parser.add_option('--charge', type = 'int', action = 'store', default = 2)
    (options, args) = parser.parse_args()

    req = [ 'output', 'pickle', 'max_mass', 'normalize', 'spec_file', 'bipartite_file']
    if not mandatories_defined(req, options):
        parser.print_help()
        exit(-1)

    make_test_data(options)

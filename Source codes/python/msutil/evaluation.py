#!/usr/bin/env python2.5
#
# Copyright 2010 <fill in later>

__authors__ = [ 'Adam Gustafson <amg81@ee.washington.edu>' ]

from operator import itemgetter
from protein.peptide import Peptide
from msutil.psm import PeptideSpectrumMatches

def rel_rank(pep_list, pep_scores, d_pep_list, d_pep_scores, true_pep, K, norm 
= True):
    """
    An evaluation script for relative ranking for Tandem Mass Spectrometry.

    Inputs:
        pep_list:       A list of Peptide objects (see protein.peptide.Peptide
                        class) from the true peptide database.
        pep_scores:     The scores (from gmtkJT, etc.) corresponding to
                        the peptides in pep_list.  Assumes scores are NOT 
                        normalized
                        for the length of the peptide
        d_pep_list:     A list of Peptide objects from the decoy peptide db.
        d_pep_scores:   The scores corresponding to the peptides in d_pep_list.
        true_pep:       A peptide object corresponding to the true peptide used 
                        as evidence in gmtk (the peptide which generated the
                        observed spectrum).
        K:              An integer in {1,2, ..., |{true database}|} used for
                        returning the performance in relative ranking (see 
                        outputs).
        norm:           A boolean variable indicating whether or not to 
                        normalize the scores by the number of frames in gmtk 
                        (1 + peptide length).
    Outputs:
        performance:    Equal to the number of peptides from the true database
                        within the top K scores, divided by K.  I.e., if the top
                        K = 10 scores are comprised of 9 peptides from the
                        true database and 1 decoy peptide, 0.9 is returned.
        true_pep_pos:   The rank (from 1 to |{true database}| + |{decoy 
                        database}|) of the true peptide used as evidence.
                        If not found, -1 is returned.
    """
    
    if not isinstance(pep_list[0], Peptide):
        raise(TypeError, 'pep_list must be a list of peptides')
    if not isinstance(pep_scores[0], (int, float)):
        raise(TypeError, 'pep_scores must be a list of numbers')        
    if not isinstance(d_pep_list[0], Peptide):
        raise(TypeError, 'd_pep_list must be a list of peptides')    
    if not isinstance(d_pep_scores[0], (int, float)):
        raise(TypeError, 'd_pep_scores must be a list of numbers')          
    if not isinstance(true_pep, Peptide):
        raise(TypeError, 'true_pep must be a peptide object')          
    if not isinstance(K, int):
        raise(TypeError, 'K must be an integer')
    if not isinstance(norm, bool):
        raise(TypeError, 'norm must be a boolean')
        
    # make list of tuples of peptide objects, scores, and flag for true peptide
    # note: fix scores by diding 1 + length of the peptide
    ranked_list = []
    # true peptides
    for i in range(len(pep_list)):
        if norm:
            pep_scores[i] = float(pep_scores[i])/(pep_list[i].length + 1)
        else:
            pep_scores[i] = float(pep_scores[i])
        ranked_list.append((pep_list[i], pep_scores[i], True))
    # decoy peptides
    for i in range(len(d_pep_list)):
        if norm:
            d_pep_scores[i] = float(d_pep_scores[i])/(d_pep_list[i].length + 1)
        else:
            pep_scores[i] = float(pep_scores[i])
        ranked_list.append((d_pep_list[i], d_pep_scores[i], False))

    # rank the list according to decreasing score order
    ranked_list = sorted(ranked_list, key = itemgetter(1), reverse = True)  

    # compute performance
    if K < 1 or K > len(pep_scores):
        raise(ValueError, 'incorrect value of K.')
    num_true = 0
    for i in range(0,K):
        current_record = ranked_list[i]
        if current_record[2] == True:
            num_true = num_true + 1
    performance = float(num_true)/K
    
    # find the position in the rankings of the peptide used as evidence.
    true_pep_pos = -1
    for i, record in enumerate(ranked_list):
        if record[0].seq == true_pep.seq and record[2] == True:
            true_pep_pos = i + 1

    return performance, true_pep_pos
    
    
def abs_rank(spectrum_IDs, guessed_peps, psm):
    """
    An evaluation script for absolute ranking for Tandem Mass Spectrometry.

    Inputs: 
        spectrum_IDs:   A list of IDs corresponding to the spectra used for
                        absolute ranking.
        guessed_peps:   A list of Peptide objects corresponding to the 
                        inferred peptide for the given spectrum_IDs.
        psm:            A PeptideSpectrumMatches object, which yields the 
                        bijection between spectrum_IDs and true peptides; i.e., 
                        the true matches.
    Outputs:
        performance:    The percentage of peptides that were correctly inferred.
        correct_IDs:    A list of spectrum_IDs corresponding to the peptides 
                        which were correctly inferred.
    """
    
    if not isinstance(spectrum_IDs[0], int):
        raise(TypeError, 'spectrum_IDs must be a list of intger IDs')
    if not isinstance(guessed_peps[0], Peptide):
        raise(TypeError, 'guessed_peps must be a list of peptides')    
    if not isinstance(psm, PeptideSpectrumMatches):
        raise(TypeError, 'psm must be a PeptideSpectrumMatches object')
    if len(spectrum_IDs) != len(guessed_peps):
        raise(RuntimeError, 'must have len(spectrum_IDs) == len(guessed_peps)')
    
    # find total number of spectra from which peptides were inferred
    N = len(spectrum_IDs)
    
    # create dictionary of spectrum_IDs and guessed_Peptides
    guessed_psm = {}
    for i in range(0,N):
        guessed_psm[spectrum_IDs[i]] = guessed_peps[i]
        
    # iterate over keys, checking to see if the peptide was guessed correctly
    correct_IDs = []
    for key in spectrum_IDs:
        if psm[key].seq == guessed_psm[key].seq:
            correct_IDs.append(key)
    performance = float(len(correct_IDs))/N
    
    return performance, correct_IDs
    

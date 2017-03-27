#!/usr/bin/env python

import itertools

def refine(scorelists):
    """Create a list of method scores which contain only shared spectra.

    Arguments:
       scorelists: List of [(targets,decoys)] pairs, where each of targets
           and decoys is itself an iterable of (sid, peptide, score) records.
           Each method is represented by a (targets,decoys) pair.

    Returns:
       newscorelists: List of [(targets,decoys)] pairs, where each of targets
           and decoys is a list of scores for the spectra that are scored in
           all of the methods.

    """
    # Find the sids common to all the methods.
    sids = [ ]
    for targets, decoys in scorelists:
        sids.append(set(r[0] for r in targets) & set(r[0] for r in decoys))
    final = sids[0]
    for ids in sids[1:]:
        final = final & ids

    # Filter the (sid, peptide, score) records to include those in sids.
    newscorelists = [ ]
    pred = lambda r: r[0] in final
    for targets, decoys in scorelists:
        newtargets = list(itertools.ifilter(pred, targets))
        newdecoys = list(itertools.ifilter(pred, decoys))
        #newscorelists.append( (newtargets, newdecoys) )
        newscorelists.append( (list(r[2] for r in newtargets),
                               list(r[2] for r in newdecoys)) )
    return newscorelists

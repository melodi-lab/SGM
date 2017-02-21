#!/usr/bin/env python2.5
#
# Adapted from Bill Noble's make-pq-plot.py.

def absolute_ranking(targets, decoys):
    """Get the points on a #positive matches vs. q-value plot.

    Args:
        targets: Each elt is the score of the best peptide-spectrum match
            from the target database, for a particular spectrum. We don't
            care which elt corresponds to which spectrum.
        decoys: Each elt is the score of the best peptide-spectrum match
            from the decoy database, for a particular spectrum. We don't
            care which elt corresponds to which spectrum. We do not require
            that elt i in decoys corresponds to the same spectrum as elt i
            in targets.

    Returns:
        qvalues: A vector whose length equals the number of spectra, where
            elt i is the q-value we get by returning the top i matches from
            a ranked list of the matches in union(target, decoys). An absolute
            ranking plot is produced by plotting range(len(qvalues)) on the
            x-axis against qvalues on the y-axis.
    """
    # Sort the scores in decreasing order.
    targets.sort(reverse = True)
    decoys.sort(reverse = True)

    # Factor to account for relative size of targets and decoys.
    ratio = float(len(targets)) / float(len(decoys))

    # Compute FDRs.
    numDecoys = [ ]
    fdrs = [ ]
    decoyIndex = 0
    decoy = decoys[decoyIndex]
    for targetIndex in xrange(0, len(targets)):
        target = targets[targetIndex]

        # Find the first decoy < this target score.
        while (decoy >= target):
            decoyIndex = decoyIndex + 1
            if (decoyIndex >= len(decoys)):
                break
            decoy = decoys[decoyIndex]

        # Estimate the FDR, with a ceiling of 1.0.
        fdr = ratio * (float(decoyIndex) / float(targetIndex + 1))
        if (fdr > 1.0):
            fdr = 1.0

        numDecoys.append(decoyIndex)
        fdrs.append(fdr)

    # Compute q-values.
    qvalues = []
    maxFDR = 0.0
    for fdrIndex in xrange(0, len(fdrs)):
        if (fdrs[fdrIndex] > maxFDR):
            maxFDR = fdrs[fdrIndex]
        qvalues.append(maxFDR)
    return qvalues

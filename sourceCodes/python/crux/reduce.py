#!/usr/bin/env python
#
# Copyright 2011 <fill in later>

"""Code for educing crux search-for-matches logs to top-1 identifications."""

import cPickle as pickle
import os
import csv
import itertools
import progress

from protein.peptide import Peptide

def load_crux(directory, cachefile = None, sids = None, progressbar = False,
              scorelbl = "xcorr score", maxfcn = max, convfcn = lambda x: x,
              all = False, norm = False):
    """Collect the best target and decoy matches from crux-output.

    Arguments:
        directory: The crux-output directory, which must contain two files:
            search.target.txt and search.decoy.txt.
        sids: If None, load all the data. If any iterable of sids, only return
            records for those spectra.
        cachefile: Name of a cachefile. If the file exists, we load the cached
            results from it. Otherwise, create a cache once the crux data has
            been loaded.
        progressbar: If true, print out a progress bar while parsing crux logs.
        scorelbl: See _load_crux
        maxfcn: See _load_crux
        convfcn: See _load_crux

    """
    if cachefile and os.path.exists(cachefile):
        return _load_crux_cachefile(cachefile, sids)

    files = [ 'search.target.txt', 'search.decoy.txt' ]
    output = ( [ ], [ ] )
    for fn, lst in zip(files, output):
        for r in _load_crux(os.path.join(directory, fn),
                            sids = None,
                            progressbar = progressbar, scorelbl = scorelbl,
                            maxfcn = maxfcn, convfcn = convfcn, all = all,
                            norm = norm):
            lst.append(r)

    if cachefile:
        _save_crux_cachefile(cachefile, output)
    if sids:
        lst = list(itertools.ifilter(lambda r: r[0] in sids, lst))

    return output

def _load_crux_cachefile(cachefile, sids = None):
    """Load crux results from a cachefile. See _save_crux_cachefile.

    Arguments:
        sids: If defined, only return records for these spectrum ids.

    """
    f = open(cachefile)
    output = pickle.load(f)
    assert(output[0])
    assert(output[1])
    assert(len(output) == 2)

    if sids:
        targets = list(itertools.ifilter(lambda r: r[0] in sids, output[0]))
        decoys = list(itertools.ifilter(lambda r: r[0] in sids, output[1]))
        assert(len(targets) == len(sids))
        assert(len(decoys) == len(sids))
        output = (targets, decoys)

    return output

def _save_crux_cachefile(cachefile, payload):
    """Save the payload to a cache file."""

    f = open(cachefile, 'wb')
    pickle.dump(payload, f, -1)

def _load_crux(filename, scorelbl = "xcorr score", maxfcn = max,
               convfcn = lambda x: x,
               sids = None, progressbar = False, all = False,
               norm = False):
    """Load a crux search-for-matches txt file. See load_crux.

    Arguments:
       filename: Crux TXT file to open.
       scorelbl: Column in 'filename' that represents the PSM score.
       maxfcn: Defines the best (maximum) score. Defaults to max(), but
         the user is free to define any function that select one floating
         point number out of an iterable of fp numbers. When the score is
         a p-value, set maxfcn = min: smaller p-values mean more significant.
       convfcn: Convert the score in the yielded scored PSM.
       sids: Select only the spectra in the given set (if None, no filtering)
       progressbar: Display an ASCII progressbar.

    """
    pb = None
    if progressbar:
        widgets = [ progress.Percentage(), progress.Bar(), progress.ETA() ]
        size = os.path.getsize(filename)
        pb = progress.ProgressBar(widgets = widgets, maxval = size)
        pb.start()

    f = open(filename)
    reader = csv.DictReader(f, delimiter = '\t')
    if scorelbl=='b/y ions matched':
        for sid, rows in itertools.groupby(reader, lambda r: int(r["scan"])):
            if not sids or sid in sids:
                if pb: pb.update(f.tell())
                if all:
                    for r in rows:
                        yield (sid, r["sequence"], convfcn(float(r[scorelbl])/(float(2*len(r["sequence"])-2))))
                else:
                    if(norm == False):
                        r = maxfcn(rows, key = lambda r: float(r[scorelbl]))
                        yield (sid, Peptide(r["sequence"]), convfcn(float(r[scorelbl])))
                    else:
                        r = maxfcn(rows, key = lambda r: float(r[scorelbl])/(float(2*len(r["sequence"])-2)))
                        yield (sid, Peptide(r["sequence"]), convfcn(float(r[scorelbl])/(float(2*len(r["sequence"])-2))))

        if progressbar:
            print '\n'
    else:
        for sid, rows in itertools.groupby(reader, lambda r: int(r["scan"])):
            if not sids or sid in sids:
                if pb: pb.update(f.tell())
                if all:
                    for r in rows:
                        yield (sid, r["sequence"], convfcn(float(r[scorelbl])))
                else:
                    r = maxfcn(rows, key = lambda r: float(r[scorelbl]))
                    yield (sid, Peptide(r["sequence"]), convfcn(float(r[scorelbl])))
                
        if progressbar:
            print '\n'


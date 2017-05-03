import csv
from math import *
import os.path
from os import listdir
from os.path import isfile, join



def readindent(f):
    #print "-------%s"%f
    input_file = csv.DictReader(open(f), delimiter='\t')
    spec2 = []
    for row in input_file:
        spec2.append((row["Kind"], int(row["Sid"]), row[
                     "Peptide"], float(row["Score"])))
    spec2sorted = sorted(spec2, key=lambda x: x[3], reverse=True)
    if len(spec2sorted) == 0:
        return [], 0, 1

    return spec2sorted, 1, 1 

    if len(spec2sorted) <= 10:
        maxs = spec2sorted[0][3]
        maxs = 1
        return [(i[0], i[1], i[2], i[3] / maxs) for i in spec2sorted], 1, maxs
    n_decoy = 0
    n_all = 0
    for row in spec2sorted:
        if row[0] == 't':
            n_all += 1
        else:
            n_all += 1
            n_decoy += 1
        if float(n_decoy) / float(n_all) > 0.01:
            break
    
    score2 = (spec2sorted[n_all][3])

    #print n_all, score2
    #if score2 < 1e-3:
    #    return spec2sorted, 0, score2
    score2 = 1

    return [(i[0], i[1], i[2], i[3] / score2) for i in spec2sorted], 1, score2


def merge_within(filenames):
    Spec1 = {}
    for f in filenames:
        if not os.path.isfile(f):
            continue
        spec, flag, rate = readindent(f)
        for i in spec:
            if not i[1] in Spec1:
                Spec1[i[1]] = i
            else:
                s1 = Spec1[i[1]][3]
                s2 = i[3]
                if s1 < s2:
                    Spec1[i[1]] = i
    return Spec1.values()


def merge_out(Spec):
    Spec1 = {}
    for spec in Spec:
        for i in spec:
            sid = i[1]
            while sid in Spec1:
                sid += 10000
            Spec1[sid] = (i[0], sid, i[2], i[3])
    return Spec1.values()


def printIdent(spec):
    spec2sorted = sorted(spec, key=lambda x: x[3], reverse=True)
    n_decoy = 0
    n_all = 0
    for row in spec2sorted:
        if row[0] == 't':
            n_all = n_all + 1
        else:
            n_all = n_all + 1
            n_decoy = n_decoy + 1
        if float(n_decoy) / float(n_all) > 0.01:
            break
    print "q=0.01: %d" % n_all


def writeSpec(spec, outputfile):
    f = open(outputfile, 'w')
    f.write('Kind\tSid\tPeptide\tScore\n')
    for row in spec:
        f.write('%s\t%s\t%s\t%.15f\n' % (row[0], row[1], row[2], row[3]))


#
import sys
Dir = sys.argv[1]
Org = sys.argv[2]
outputfile = Dir + Org + '.txt'
#Dir = "/home/ubuntu/SGM/sourceCodes/MSGFPlus/" 
#Dir = Dir + Org + "-no-isotopy/"
Spec = []
for idfolder in os.listdir(Dir):
    idfoldername = Dir + idfolder
    if not os.path.isdir(idfoldername):
        continue
    filenames = []    
    for chargefolder in os.listdir(idfoldername):
        chargefoldername = idfoldername + '/' + chargefolder
        if not os.path.isdir(chargefoldername):
            continue
        charge = int(''.join(chargefolder.split(Org + '-ch')))        
        if charge >= 6:
            continue
        target_file = chargefoldername + '/' + Org + "-ch%d-targets.tsv"%(charge)
        decoy_file = chargefoldername + '/' + Org + "-ch%d-decoys.tsv"%(charge)
        output = chargefoldername + '/' + Org + ".txt"
        filenames.append(output)

    spec = merge_within(filenames)
    Spec.append(spec)
S = merge_out(Spec)
printIdent(S)
writeSpec(S, outputfile)

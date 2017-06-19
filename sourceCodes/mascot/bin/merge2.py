import csv
from math import *
import os.path
from os import listdir
from os.path import isfile, join

def read_indent(FILE_NAME):
    input_file=csv.DictReader(open(FILE_NAME),delimiter='\t')
    
    target_decoy=[]
    for row in input_file:
        target_decoy.append((row["Kind"],int(row["Sid"]),row["Peptide"],float(row["Score"])))
    sorted_target_decoy=sorted(target_decoy,key=lambda x:-x[3])
#   print sorted_target_decoy[0:3]
    n_t=0;
    n_d=0;
    score=0;
    n=-1;
    for row in sorted_target_decoy:
        if row[0]=='t':
            n_t=n_t+1
        elif row[0]=='d':
            n_d=n_d+1
        if float(n_d)/float((n_t+n_d))>0.01:
            score=row[3]
            n=n_t+n_d
            break
    #print(score)
    #print(n)
    #print(n_t)
    #print(n_d)
    n_t=0;
    n_d=0;
    score=0;
    n=-1;
    for row in sorted_target_decoy:
        if row[0]=='t':
            n_t=n_t+1
        elif row[0]=='d':
            n_d=n_d+1
        if float(n_d)/float((n_t+n_d))>0.04:
            score=row[3]
            n=n_t+n_d
            break
    #print(score)
    #print(n)
    #print(n_t)
    #print(n_d)
    
    #print(sid)
    q_range=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
    q_range=range(101)
    q_range=[float(i)/1000 for i in q_range]
    n_range=[]
    score_range=[]
    for q in q_range:
        n_t=0;
        n_d=0;
        score=0;
        n=-1;
        for row in sorted_target_decoy:
            if row[0]=='t':
                n_t=n_t+1
            elif row[0]=='d':
                n_d=n_d+1
            if float(n_d)/float((n_t+n_d))>q:
                score=row[3]
                n=n_t+n_d
                break
        n_range.append(n_t)
        score_range.append(score)
#   print "q:\t",
#   for q in q_range:
#       print"%f\t"%(q),
    print "\nn:\n",
    print ','.join([str(n) for n in n_range])
    #for n in n_range:
#       print"%f\t"%(n),
#   print "\nscore:\t",
#   for s in score_range:
#       print"%f\t"%(s),
 
    return n_range
    #print(sid)

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
filenames = []    
for idfolder in os.listdir(Dir):
    idfoldername = Dir + idfolder
    if not os.path.isdir(idfoldername):
        continue
     
    for chargefolder in os.listdir(idfoldername):
        chargefoldername = idfoldername + '/' + chargefolder
        if not os.path.isdir(chargefoldername):
            continue
        charge = int(''.join(chargefolder.split(Org + '-ch')))        
        if charge >= 6:
            continue

        output = chargefoldername + '/' + Org + ".txt"
        filenames.append(output)


import numpy as np
results = np.array([0]*101)
for f in filenames:

    a = np.array(read_indent(f))
    results = results + a
print ','.join([str(i) for i in results])    

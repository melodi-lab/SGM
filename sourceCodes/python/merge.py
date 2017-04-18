import csv
from math import *
import os.path
from os import listdir
from os.path import isfile, join

# ch2file="test_scores_output_fast_yeast_01_ch2.txt.sumScoreIdent"
# ch3file="test_scores_output_fast_yeast_01_ch3.txt.sumScoreIdent"
# ch23file="test_scores_output_fast_yeast_01_ch23.txt.sumScoreIdent"
# ch2file="../output/test_scores_output_fast_bullseyeworm-ch4.txt.sumScoreIdent"
# ch3file="../output/test_scores_output_fast_bullseyeworm-ch3.txt.sumScoreIdent"
# ch3file="test_scores_output_fast_bullseyeworm-ch2.txt.sumScoreIdent2"
# ch23file="test_scores_output_fast_bullseyeworm-ch2345.txt.sumScoreIdent2"
# ch2file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent"
# ch3file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch4.txt.sumScoreIdent"
# ch23file="/s1/wrbai/codes/output/test_scores_output_fast_malaria-ch2345.txt.sumScoreIdent"
# ch2file="/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch2/malaria-ch2-ident.txt"
# ch2file="malaria-ch2345-ident.txt"
# ch3file="/s1/wrbai/codes/john_idents/7_23/malaria/malaria-ch3/malaria-ch3-ident.txt"
# ch23file="malaria-ch2345-ident.txt"


def readindent(f):
    print "-------%s"%f
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

    print n_all, score2
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
                if s1 > s2:
                    continue
                else:
                    Spec1[i[1]] = i
    return Spec1.values()


def merge_out(Spec):
    Spec1 = {}
    for spec in Spec:
        for i in spec:
            if not i[1] in Spec1:
                Spec1[i[1]] = i
            else:
                sid = i[1]
                while sid in Spec1:
                    sid += 10000
                Spec1[i[1]] = (i[0], sid, i[2], i[3])
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

# def merge_indent(CHA, CHB, CHAB, fix_r=0):
#     print("-------------------------")
#     print("CHA   :" + CHA)
#     print("CHB   :" + CHB)
#     print("CHAB  :" + CHAB)
#     input_file = csv.DictReader(open(CHA), delimiter='\t')
#     spec2 = []
#     for row in input_file:
#         spec2.append((row["Kind"], int(row["Sid"]), row[
#                      "Peptide"], float(row["Score"])))
#     spec2sorted = sorted(spec2, key=lambda x: x[3], reverse=True)
#     n_decoy = 0
#     n_all = 0
#     for row in spec2sorted:
#         if row[0] == 't':
#             n_all = n_all + 1
#         else:
#             n_all = n_all + 1
#             n_decoy = n_decoy + 1
#         if float(n_decoy) / float(n_all) > 0.01:
#             break
#     n_all_A = n_all
#     score2 = (spec2sorted[n_all][3])

#     input_file = csv.DictReader(open(CHB), delimiter='\t')
#     spec3 = []
#     for row in input_file:
#         spec3.append((row["Kind"], int(row["Sid"]), row[
#                      "Peptide"], float(row["Score"])))
#     spec3sorted = sorted(spec3, key=lambda x: x[3], reverse=True)
#     n_decoy = 0
#     n_all = 0
#     for row in spec3sorted:
#         if row[0] == 't':
#             n_all = n_all + 1
#         else:
#             n_all = n_all + 1
#             n_decoy = n_decoy + 1
#         if float(n_decoy) / float(n_all) > 0.01:
#             break
#     n_all_B = n_all
#     score3 = (spec3sorted[n_all][3])
#     # print(score2/score3)
#     # print(len(spec2sorted))
#     # print(len(spec3sorted))

#     sid_all = []
#     sid_2t = []
#     sid_2d = []

#     sid_3t = []
#     sid_3d = []
#     for row in spec2:
#         sid_all.append(row[1])
#         if row[0] == 't':
#             sid_2t.append(row[1])
#         else:
#             sid_2d.append(row[1])
#     for row in spec3:
#         sid_all.append(row[1])
#         if row[0] == 't':
#             sid_3t.append(row[1])
#         else:
#             sid_3d.append(row[1])

# # sid_not_2=list(set(sid_all)-set(sid_2))
# # sid_not_3=list(set(sid_all)-set(sid_3))

#     for row in list(set(sid_all) - set(sid_2t)):
#         spec2.append(('t', row, 'AAAA', -1000000000))
#     for row in list(set(sid_all) - set(sid_2d)):
#         spec2.append(('d', row, 'AAAA', -1000000000))
#     for row in list(set(sid_all) - set(sid_3t)):
#         spec3.append(('t', row, 'AAAA', -1000000000))
#     for row in list(set(sid_all) - set(sid_3d)):
#         spec3.append(('d', row, 'AAAA', -1000000000))

#     spec2 = sorted(spec2, key=lambda x: x[0], reverse=True)
#     spec2 = sorted(spec2, key=lambda x: x[1])
#     spec3 = sorted(spec3, key=lambda x: x[0], reverse=True)
#     spec3 = sorted(spec3, key=lambda x: x[1])
# # print(spec2[0])
# # print(spec2[1])
# # print(spec2[2])
# # print(spec2[3])
# # print(spec2[4])
# # print(spec2[5])
# # print(spec3[0])
# # print(spec3[1])
# # print(spec3[2])
# # print(spec3[3])
# # print(spec3[4])
# # print(spec3[5])
#     print(len(spec2))
#     print(len(spec3))
#     r1 = 1
#     r2 = 1
#     if score2 * score3 > 1e-10:
#         r1 = score2 / 1000
#         r2 = score3 / 1000
#     else:
#         r1 = 1
#         r2 = 1
#     if fix_r:
#         r1 = 1
#         r2 = 1
#     if len(spec3) < 100:
#         r1 = 1
#         r2 = 1
#     if len(spec2) < 100:
#         r1 = 1
#         r2 = 1
#     if r1 < 0.00000001:
#         r1 = 1
#         r2 = 1

#     if r2 < 0.00000001:
#         r1 = 1
#         r2 = 1
#     if r1 / r2 > 100:
#         r1 = 1
#         r2 = 1
#     if r2 / r1 > 100:
#         r1 = 1
#         r2 = 1
#     spec23 = []
#     print r1, r2
#     for row2, row3 in zip(spec2, spec3):
#         # if row2[0]!=row3[0]:
#         #     print "asdfasdf"
#         if row2[3] / r1 >= row3[3] / r2:
#             s2 = row2[3]
#             if s2 > -1000000:
#                 s2 = s2 / r1
#             row2r = (row2[0], row2[1], row2[2], s2)
#             spec23.append(row2r)

#         # spec23.append(row2);
#         if row2[3] / r1 < row3[3] / r2:
#             s3 = row3[3]
#             if s3 > -1000000:
#                 s3 = s3 / r2

#             row3r = (row3[0], row3[1], row3[2], s3)
#             spec23.append(row3r)

#     spec23sorted = sorted(spec23, key=lambda x: x[3], reverse=True)
#     n_decoy = 0
#     n_all = 0
#     for row in spec23sorted:
#         if row[0] == 't':
#             n_all = n_all + 1
#         else:
#             n_all = n_all + 1
#             n_decoy = n_decoy + 1
#         if float(n_decoy) / float(n_all) > 0.01:
#             break
# # print ('n_A:%d n_b:%d n_ab:%d'%(n_all_A,n_all_B,n_all))

#     f = open(CHAB, 'w')
#     f.write('Kind\tSid\tPeptide\tScore\n')
#     for row in spec23:
#         f.write('%s\t%s\t%s\t%.15f\n' % (row[0], row[1], row[2], row[3]))


# def merge_mutilple_file(allfilenames, outputfile, fix_r=0):
#     filenames = []
#     for row in allfilenames:
#         if os.path.isfile(row) and os.path.getsize(row) > 0:
#             filenames.append(row)
#         else:
#             print ("%s does not exist" % row)

#     if len(filenames) <= 1:
#         print "less than two files-----------------------"
#         return

#     merge_indent(filenames[0], filenames[1], outputfile, fix_r)
#     if len(filenames) == 2:
#         return
#     for i in filenames[2:]:
#         merge_indent(i, outputfile, outputfile, fix_r)


#
outputfile = "result.txt"
Spec = []
for i in range(20):
    id = i + 1
    filenames = []
    dir = "../MSGFPlus"
    for j in range(4):
        c = j + 2
        filenames.append("%s/plasm-%d/plasm-ch%d/plasm.txt" % (dir, id, c))
    spec = merge_within(filenames)
    Spec.append(spec)
S = merge_out(Spec)
printIdent(S)
writeSpec(S, outputfile)

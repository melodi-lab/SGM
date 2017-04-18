import csv
# target_file="Linfeng_112710_output/tide-search.target.txt"
import os.path
import sys


Dir = sys.argv[1]
Org = sys.argv[2]
print Dir
#Dir = "/home/ubuntu/SGM/sourceCodes/MSGFPlus/" 
#Dir = Dir + Org + "-no-isotopy/"

for idfolder in os.listdir(Dir):
    idfoldername = Dir + idfolder
    if not os.path.isdir(idfoldername):
        continue
    for chargefolder in os.listdir(idfoldername):

        chargefoldername = idfoldername + '/' + chargefolder
        if not os.path.isdir(chargefoldername):

            continue
        charge = int(''.join(chargefolder.split(Org + '-ch')))
        
        target_file = chargefoldername + '/' + Org + "-ch%d-targets.tsv"%(charge)
        decoy_file = chargefoldername + '/' + Org + "-ch%d-decoys.tsv"%(charge)

        output = open(chargefoldername + '/' + Org + ".txt", 'w')
#		output=open("kim-%d/kim-ch%d/kim.txt"%(id,charge),'w')

        if not os.path.isfile(target_file):
            
            continue
        if not os.path.isfile(decoy_file):
            
            continue
        # decoy_file="Linfeng_010211_output_2/tide-search.decoy.txt"
        target = []
        decoy = []
        name = "EValue"
    #	name="MSGFScore"

        with open(target_file) as csvfile:
            spamreader = csv.DictReader(csvfile, delimiter='\t')
            for row in spamreader:
                s = row["Peptide"]
                a = s.split("+229.163")
                b = ''.join(a)
                c = b.split("+57.021")
                d = ''.join(c)
                a1 = row["ScanNum"]
                a2 = row["Title"]
                a3 = a2.split(".")
                a4 = a3[-1]
                sid = int(a4)

                target.append(('t', sid, d, -float(row[name])))

        with open(decoy_file) as csvfile:
            spamreader = csv.DictReader(csvfile, delimiter='\t')
            for row in spamreader:
                s = row["Peptide"]
                a = s.split("+229.163")
                b = ''.join(a)
                c = b.split("+57.021")
                d = ''.join(c)
                a1 = row["ScanNum"]
                a2 = row["Title"]
                a3 = a2.split(".")
                a4 = a3[-1]
                sid = int(a4)
                target.append(('d', sid, d, -float(row[name])))

        Spec = {}
        for t in target:
            if not t[1] in Spec:
                Spec[t[1]] = t
            else:
                s1 = Spec[t[1]][3]
                s2 = t[3]
                if s2 > s1:
                    Spec[t[1]] = t   
        output.write('Kind\tSid\tPeptide\tScore\n')
        for row1 in Spec.values():
            output.write("%s\t%d\t%s\t%.15f\n" % (row1[0], row1[1], row1[2], row1[3]))


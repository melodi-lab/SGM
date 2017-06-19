import csv
import os
import subprocess
import sys
from multiprocessing import Pool
from multiprocessing import cpu_count

def runCommand(command):
    sys.stderr.write("%s\n"%command)
    subprocess.call(command, shell=True)

def parseMGFfilename(mgf):
    a = mgf.strip('.mgf')
    a = a.split('-')

    ORG = a[0]
    Charge = int(a[2][-1])
    ID1 = int(a[1])
    ID2 = 0
    if len(a) >= 4:
        ID2 = int(a[3])
    ID = 100 * ID1 + ID2
    return (ORG, Charge, ID)


def read_file(FILENAME):
    targetRead = open(FILENAME)
    status = 0
    peakDataPathName = "Peak list data path"
    proteinHitName = "Protein hits"
    databaseName = "Database"

    target = []
    for line in targetRead.readlines():
        if status == 0:
            if line[:len(peakDataPathName)] == peakDataPathName:
                splitLine = line.split('"')
                splitLine = splitLine[1].split('\\')
                targetMgfFilename = splitLine[-1]
            elif line[:len("Database")] == "Database":
                splitLine = line.split('"')
                splitLine = splitLine[1].split('\\')
                Database = splitLine[-1]
                if "target" in Database:
                    kind = 't'
                elif "decoy" in Database:
                    kind = 'd'
                else:
                    kind = 'unknown'
    return target_sid.values(), kind, parseMGFfilename(targetMgfFilename)                

def read_file_dat(FILENAME):
    kind = ''
    targetRead = open(FILENAME)
    line_faste = "COM="
    line_mgf = "FILE="
    status = 0
    status2 = 0
    target_sid = {}
  
    for line in targetRead.readlines():        
        if "fasta" in line:
            if "target" in line:
                kind = 't'
            elif "decoy" in line:
                kind = 'd'
            status += 1
        
        if line[:len(line_mgf)] == line_mgf:
            line = line.strip('\n')
            splitLine = line.split('\\')
            targetMgfFilename = splitLine[-1]
            status2 += 1
        if '_' in line:
            cc = line.index('_')
            if line[0] == 'q' and line[1:cc].isdigit():
                
                ss = line.split(",")
                if len(ss) <= 8:
                    continue
                if 'term' in line:
                    continue
                t = {}
                t['sid'] = int(line[1:cc])
                #print t['sid']
                t['pep_seq'] = ss[4]
                t['pep_score'] = float(ss[7])
                t['kind'] = kind
                if t["sid"] in target_sid:
                    if t['pep_score'] > target_sid[t["sid"]]['pep_score']:
                        target_sid[t["sid"]] = t
                else:
                    target_sid[t["sid"]] = t
   
    if 'merge'  not in targetMgfFilename:

        return kind, parseMGFfilename(targetMgfFilename),target_sid.values(), 0            
    return kind, [], [],1 

            


def file_len(fname):

    return sum(1 for line in open(fname))

  

def read_file_csv(FILENAME, kind):  
    if file_len(FILENAME) <= 1:
        return []
    targetRead = open(FILENAME)
    target_sid = {}
    for line in targetRead.readlines():
        splitLine = line.split(';')
        t = {}
        t['sid'] = int(splitLine[2].split('.')[-1])
        #print t['sid']
        t['pep_seq'] = splitLine[4]
        t['pep_score'] = float(splitLine[-3])
        t['kind'] = kind
        if t["sid"] in target_sid:
            if t['pep_score'] > target_sid[t["sid"]]['pep_score']:
                target_sid[t["sid"]] = t
        else:
            target_sid[t["sid"]] = t

    return target_sid.values()        




    


datfolder = "../result/6-12-2017-dat/1/"
csvfolder = "../result/6-12-2017-dat/1/csv/"
Outfolder = "../result/6-12-2017-dat-yeast/"

datfolder = "../result/5-24-2017-dat/2/"
csvfolder = "../result/5-24-2017-dat/2/csv/"
Outfolder = "../result/5-24-2017-dat-human/"

datfolder = "../result/5-24-2017-dat/106/"
csvfolder = "../result/5-24-2017-dat/106/csv/"
Outfolder = "../result/5-24-2017-dat-plams4/"


allfile= []
for file in os.listdir(datfolder):
    if file[-4:] == '.dat':

        f = file.strip('.dat')
        allfile.append(f)

if False:  

    os.makedirs(os.path.dirname(csvfolder))
    commands = []
    for f in allfile:
        command = "java -jar ../datParser/mascotdatfile-3.6.0/mascotdatfile-3.6.0.jar 1000 %s%s.csv %s%s.dat"%(csvfolder,f,datfolder,f)
        print command
        commands.append(command)
    pool = Pool()  
    pool.map(runCommand, commands)      


TARGET = {}
DECOY = {}
Names = set()
for filename in sorted(allfile):
    kind, name, target, merged = read_file_dat(datfolder+filename+'.dat')
    print kind, name
    print csvfolder+filename+'.csv'
    if merged:
        continue
    #target = read_file_csv(csvfolder+filename+'.csv', kind)
    print kind
    if kind == 't':
        TARGET[name] = target
    elif kind == 'd':
        DECOY[name] = target
    Names.add(name)
#print DECOY
for name in Names:
    # print name
    target = []
    decoy = []
    if name in TARGET:
        target = TARGET[name]
    else:
        target = []

    if name in DECOY:
        decoy = DECOY[name]
    else:
        decoy = []
    #print decoy

    #print target
    #print decoy
    #print "==========="
    # if not name in DECOY:
    
    #    print target[0]['mgf']
    #    print "======================================================================================================================="
    #    continue
    # decoy = DECOY[name]

    t_and_d = {}

    for t in target:
        if t["sid"] in t_and_d:
            if t['pep_score'] > t_and_d[t["sid"]]['pep_score']:
                t_and_d[t["sid"]] = t
        else:
            t_and_d[t["sid"]] = t

    for t in decoy:
        if t["sid"] in t_and_d:
            if t['pep_score'] > t_and_d[t["sid"]]['pep_score']:
                t_and_d[t["sid"]] = t
        else:
            t_and_d[t["sid"]] = t    

    t_and_d = [(t["kind"], t["sid"], t["pep_seq"], t["pep_score"]) for t in t_and_d.values()]    
    t_and_d = sorted(t_and_d, key = lambda x : -x[3])    
    outputname = "%s/%s-%d/%s-ch%d/%s.txt"%(Outfolder,name[0],name[2],name[0],name[1],name[0]) 
    if not os.path.exists(os.path.dirname(outputname)):
        try:
            os.makedirs(os.path.dirname(outputname))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    output = open(outputname, 'w')
    output.write('Kind\tSid\tPeptide\tScore\n')
    for row1 in t_and_d:
        output.write("%s\t%d\t%s\t%.15f\n" % (row1[0], row1[1], row1[2], row1[3]))  
    #print "================="
    #print "processing target and decoy"
    #print "ORG: %s, id: %d-%d, charge: %d"%( name[0], name[2]/100, name[2]%100, name[1])    
    #print "Writing:", outputname
# print [t["pep_seq"]for t in target]

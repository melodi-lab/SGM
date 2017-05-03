import csv
import os
import subprocess
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
    return

def read_file_csv(FILENAME):  
    return  


    
datfolder = "test/"
csvfolder = "test_csv/"
Outfolder = "plasm11"
allfile= []
for file in os.listdir(datfolder):
    f = file.strip('.dat')
    allfile.append(f)
import shutil
if os.path.exists(os.path.dirname(csvfolder)):   
    shutil.rmtree(csvfolder)
if not os.path.exists(os.path.dirname(csvfolder)):   
    os.makedirs(os.path.dirname(csvfolder))
    commands = []
    for f in allfile:
        command = "java -jar mascotdatfile-3.6.0/mascotdatfile-3.6.0.jar 10000 %s%s.csv %s%s.dat"%(csvfolder,f,datfolder,f)
        print command
        commands.append(command)
    pool = Pool()  
    pool.map(runCommand, commands)      



exit()    

TARGET = {}
DECOY = {}
Names = set([])
for filename in allfile:
    target, kind, name = read_file(filename)
    if kind == 't':
        TARGET[name] = target
    elif kind == 'd':
        DECOY[name] = target
    Names.add(name)

for name in Names:
    # print name
    if name in TARGET:
        target = TARGET[name]
    else:
        target = []
    if name in DECOY:
        decoy = DECOY[name]
    else:
        decoy = []
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
    print "================="
    print "processing target and decoy"
    print "ORG: %s, id: %d-%d, charge: %d"%( name[0], name[2]/100, name[2]%100, name[1])    
    print "Writing:", outputname
# print [t["pep_seq"]for t in target]

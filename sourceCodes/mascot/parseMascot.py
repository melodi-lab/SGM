import csv
import os
def parseMGFfilename(mgf):
    a = mgf.strip('.mgf')
    a = a.split('-')
    
    ORG = a[0]
    Charge = int(a[2][-1])
    ID1 = int(a[1])
    ID2 = 0
    if len(a) >= 4:
        ID2 = int(a[3])
    ID = 100*ID1 + ID2
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
 
        
            elif line[1:len(proteinHitName) + 1] == proteinHitName:
                status += 1
        elif status == 1:
            status += 1
        elif status == 2:
            line = line.strip()
            csvEntry = line.split(',')
            status += 1
        elif status == 3:
            line = line.strip()
            start = 0
            end = 0
            quote = 0
            entry0 = []
            while end < len(line):
                if line[end] == ',' and quote == 0:
                    entry0.append(line[start:end])
                    start = end + 1
                    end += 1

                elif line[end] == '\"':
                    quote = 1 - quote
                    end += 1
                else:
                    end += 1
            entry0.append(line[start:end])

            # print ",".join(entry0)== line

            #entry = line.split(',')
            t = {}
            for key, e in zip(csvEntry, entry0):
                t[key] = e
            t['all'] = line
            t['database'] = Database
            t['mgf'] = targetMgfFilename
            t['kind'] = kind
            t['pep_score'] = float(t['pep_score'])
            t['sid'] = int(t['pep_scan_title'].split('.')[-1][:-1])
            target.append(t)
    target_sid = {}
    for t in target:
        if t["sid"] in target_sid:
            if t['pep_score'] > target_sid[t["sid"]]['pep_score']:
                target_sid[t["sid"]] = t
        else:
            target_sid[t["sid"]] = t
    #print csvEntry
    print "----------------------------"
    print "proccesing: ", FILENAME   

    print "mgf_file: ", targetMgfFilename 
    print "parse filename: ", parseMGFfilename(targetMgfFilename)
    print "database: ", Database
    if kind == 't':        
        print "kind: target"
    elif kind == 'd':
        print "kind: decoy"
    print "----------------------------"

    return target_sid.values(), kind, parseMGFfilename(targetMgfFilename)
allfilefolder = "/home/ubuntu/wenruo_aws/result_backup/5_2_2017_plasm_mascot/"
Outfolder = "plasm5_2"
allfile= []
for file in os.listdir(allfilefolder):
    allfile.append(allfilefolder + file)

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
    #print name
    if name in TARGET:
        target = TARGET[name]
    else:
        target = []
    if name in DECOY:
        decoy = DECOY[name]
    else:
        decoy = []
    #if not name in DECOY:
    
    #    print target[0]['mgf']
    #    print "======================================================================================================================="
    #    continue
    #decoy = DECOY[name]

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
#print [t["pep_seq"]for t in target]

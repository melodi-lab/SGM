import csv
targetfile = "test/F004920.csv"
targetRead = open(targetfile)
status = 0
peakDataPathName = "Peak list data path"
proteinHitName = "Protein hits"

target = []
for line in targetRead.readlines():
    if status == 0:
        if line[:len(peakDataPathName)] == peakDataPathName:

            splitLine = line.split('"')
            splitLine = splitLine[1].split('\\')
            targetMgfFilename = splitLine[-1]
            status += 1
    elif status == 1:
        if line[1:len(proteinHitName)+1] == proteinHitName:
            status += 1
    elif status == 2:
        status += 1
    elif status == 3:
    	line = line.strip()    	
        csvEntry = line.split(',')
        status += 1
    elif status == 4:
    	line = line.strip()
        entry = line.split(',')
        t = {}
        for key, e in zip(csvEntry, entry):
        	t[key] = e

        t['pep_score'] = float(t['pep_score'])
        t['sid'] = t['pep_scan_title'].split('.')[-1][:-1]
        target.append(t)


print [(t["sid"], t["pep_scan_title"]) for t in target  ]      
#print target   
print csvEntry

import csv
targetfile = "test/F004920.csv"
targetRead = open(targetfile)
status = 0
peakDataPathName = "Peak list data path"
proteinHitName = "Protein hits"

csvEntryId = {}
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
        csvEntry = line.split(',')
        id = 0
        for i in csvEntry:
            csvEntryId[i] = id
            id += 1
        status += 1
    elif status == 4:
        entry = line.split(',')

        target.append((entry[csvEntryId['prot_hit_num']]))
print target    
#print csvEntryId.keys()

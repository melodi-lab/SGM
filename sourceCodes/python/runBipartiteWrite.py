import os
import sys
import subprocess
from multiprocessing import Pool
from multiprocessing import cpu_count

inputFolder = "../../../ms2data/malaria"
outputFolder = "../data/ms2file"

inputfiles = []
for file in os.listdir(inputFolder):
    if file.endswith(".ms2"):
        inputfiles.append("%s/%s" % (inputFolder, file))
inputfiles.sort()


if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
commands = []
index = 1    
for f in inputfiles:
    for charge in [2,3,4,5,6,7]:
        outputfile = "%s/%d-ch%d.ms2" % (outputFolder, index, charge)
        if os.path.exists(outputfile):
            os.remove(outputfile)
        command = "python filter_spectra.py -z -i %s -c %d > %s" % (f, charge, outputfile)
        commands.append(command)
        #sys.stderr.write("%s\n"%command)
        #subprocess.call(command, shell=True)
    index += 1
def runCommand(command):
    subprocess.call(command, shell=True)

pool = Pool()  
pool.map(runCommand, commands)    



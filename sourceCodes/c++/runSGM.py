import os
import sys
import subprocess
from multiprocessing import Pool
from multiprocessing import cpu_count
import shutil

inputFolder = "../data/ms2file"

outputFolder = "../data/result"

if os.path.exists(outputFolder):
    shutil.rmtree(outputFolder)
os.makedirs(outputFolder)


inputfiles = []
inputfilenames = []
for file in os.listdir(inputFolder):
    if file.endswith(".ms2"):
        #num_lines = sum(1 for line in open('myfile.txt'))
        inputfiles.append("%s/%s" % (inputFolder, file))
        inputfilenames.append(file)
#inputfiles.sort()
#inputfilenames.sort()


commands = []
def getCharge(filename):
    a = '\t'.join(filename.split('.ms2'))
    a = '\t'.join(a.split('-ch'))
    a = a.strip()
    b = a.split('\t')

    return b[0], b[1]


for f, g in zip(inputfiles, inputfilenames):
    #
    index, charge = getCharge(g)
    command = "./runSGM.sh %s %s %s"%(index, charge, outputFolder)
    commands.append(command)
    #sys.stderr.write("%s\n"%command)
    #subprocess.call(command, shell=True)

def runCommand(command):
    sys.stderr.write("%s\n"%command)
    subprocess.call(command, shell=True)

pool = Pool()  
pool.map(runCommand, commands)    



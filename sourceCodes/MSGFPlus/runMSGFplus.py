import os
import sys
import subprocess
from multiprocessing import Pool
from multiprocessing import cpu_count

commands = []
for id in range(1,21):
    for charge in [2,3,4,5]:
        command = "./runPlasm.sh %d %d"%(id,charge)
        commands.append(command)

        command = "./runPlasm_isotopy.sh %d %d"%(id,charge)
        commands.append(command)


def runCommand(command):
    subprocess.call(command, shell=True)

#pool = Pool()  
#pool.map(runCommand, commands)    
def postPressing(Dir, Org):
    subprocess.call("python bin/convert_into_ident.py %s %s > log.txt"%(Dir, Org), shell=True)
    subprocess.call("python bin/merge.py %s %s > log.txt"%(Dir, Org), shell=True)
    subprocess.call("python bin/readident.py %s %s > log.txt"%(Dir, Org), shell=True)

Dir = "/home/ubuntu/SGM/sourceCodes/MSGFPlus/plasm-no-isotopy/"
Org = "plasm"
postPressing(Dir, Org)

Dir = "/home/ubuntu/SGM/sourceCodes/MSGFPlus/plasm-with-isotopy/"
Org = "plasm"
postPressing(Dir, Org)




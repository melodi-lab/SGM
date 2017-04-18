import os
import sys
import subprocess
from multiprocessing import Pool
from multiprocessing import cpu_count

commands = []
for id in range(1,21):
    for charge in [2,3,4,5,6,7]:
        command = "./run.sh %d %d"%(id,charge)
        commands.append(command)
def runCommand(command):
    subprocess.call(command, shell=True)

pool = Pool()  
pool.map(runCommand, commands)    



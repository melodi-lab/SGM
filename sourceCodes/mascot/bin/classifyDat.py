import sys
import os
import shutil

inputfolder = "../result/5-24-2017-dat/"
outputfolder = "../result/5-24-2017-dat-yeast/"

DB = {""}
start = 7028
end = 7137
for i in range(start,end+1):
    fn = inputfolder + "F00%d.dat"%i
    shutil.copyfile(fn, outputfolder + "F00%d.dat"%i )


exit()


for fn in sorted(os.listdir(inputfolder)):
    if os.path.isfile(inputfolder+fn):
        database = ""
        search_title = ""
        count = 0
        for line in open(inputfolder+fn):
            count += 1
            if count > 100:
                continue

            if line[:3] == "DB=":
                database = (line[3:]).strip("\n")
            if line[:4] == "COM=":
                s = (line[4:]).strip("\n")
                s1 = s.split()
                search_title = s1[-1]
        print fn, search_title        


        
          


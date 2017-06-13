import sys
import os
import shutil

def remove_td(s):
    ss = s.split("target")
    s1 = '[t or d]'.join(ss)
    ss = s1.split("decoy")
    s1 = '[t or d]'.join(ss)
    return s1

inputfolder = "../result/5-24-2017-dat/"

searches = {}
folder_id = 1
cc = 0
for fn in sorted(os.listdir(inputfolder)):
    sys.stderr.write("%d"%cc)
    cc += 1
    if os.path.isfile(inputfolder+fn):
        database = ""
        search_title = ""
        fid = 0
        count = 0
        for line in open(inputfolder+fn):
            count += 1
            if count > 100:
                break

            if line[:3] == "DB=":
                #database = (line[3:]).strip("\n")
            #if line[:4] == "COM=":
                if not remove_td(line) in searches:
                    searches[remove_td(line)] = folder_id
                    folder_id += 1
                fid = searches[remove_td(line)]
        folder = "%s%d/"%(inputfolder,fid)
        if not os.path.exists(folder):
            os.makedirs(folder)
        shutil.copy(inputfolder+fn,folder)
for key in searches:
    print searches[key],key                   




        
          


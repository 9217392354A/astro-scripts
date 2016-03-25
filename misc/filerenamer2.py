#Chris Fuller, Changes the name of files

import numpy as np 
import os
from os.path import join as pj
import shutil as s

dst = '/Users/chrisfuller/Documents/Personal/The-Old-Republic-Annihilation/'

files = os.listdir(dst)
files = [x for x in files if x[0] != "."]

for f in files:
    print f.split("-")[2]



        
        
        #outname = str(count) + "-StarWars-Annihilation" + "disk-" + disk + "-track-" + track +".mp3"
        #os.rename(src, pj(sourcedir,outname))        
        #command = "mv " + src + " " + pj(dst,outname)
        #print command
        #os.system(command)
        
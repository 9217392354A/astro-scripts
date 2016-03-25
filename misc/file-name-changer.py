#Chris Fuller, Changes the name of files

import numpy as np 
import os
from os.path import join as pj
import shutil as s


rootdir = "/Users/chrisfuller/Documents/Personal/Star-Wars-The-Old-Republic-Annihilation/"
dst = '/Users/chrisfuller/Documents/Personal/The-Old-Republic-Annihilation/'

topDirs = os.listdir(rootdir)
topDirs = [x for x in topDirs if x[0] != "."]
count = 0
for folder in topDirs:
    sourcedir = pj(rootdir,folder)
    files = os.listdir(sourcedir)
    #files = [x for x in topDirs if x[-3:] != "mp3"]
    for f in files: 
        count += 1
        src = pj(sourcedir,f)
        s.copy(src,dst)
        
        
        
        #outname = str(count) + "-StarWars-Annihilation" + "disk-" + disk + "-track-" + track +".mp3"
        #os.rename(src, pj(sourcedir,outname))        
        #command = "mv " + src + " " + pj(dst,outname)
        #print command
        #os.system(command)
        
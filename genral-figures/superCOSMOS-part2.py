#program to take folders full of supercosmos images and then re-name them back to thier object names
#Chris Fuller, Feb 2013
#part 2 of supercosmos program

import os
import numpy as np
from os.path import join as pj
import shutil as s

###################### inputs ######################
rootfolder = "/Users/chrisfuller/Desktop/superCOSMOS/"
outfolder = "/Users/chrisfuller/Desktop/superCOSMOS/fits/"
try: os.mkdir(outfolder)
except: print "outfolder aready exists..."
##################### functions #####################

##################### main program ##################

#make a list of root folders
folders = [x for x in os.listdir(rootfolder) if x[:6] == "super-"]

for folder in folders:
    print "starting folder ",folder
    #open list of each object name in each folder
    sec = (folder.split("-")[1])
    objs = np.loadtxt(pj(rootfolder,folder,"record-"+sec+".txt"), dtype=str,usecols=[0],  skiprows=0, delimiter=" ", unpack=False)
    
    #get list of fits files
    fits = os.listdir(pj(rootfolder,folder,"fits"))
    
    #loop through all the fits file and rename and copy the files to unified dir
    for i in range(0,len(fits)):
        outname = objs[i]+".fits"
        src = pj(rootfolder,folder,"fits",fits[i])
        dst = pj(rootfolder,outfolder,outname)
        s.copy(src,dst)
        print "copying...",src," to ",dst
        
print "program finished, onwards and upwards..."
        
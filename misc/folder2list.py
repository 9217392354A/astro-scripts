#makes a csv file from a folder
#chris fuller, march 13

#detection rate caculator, Chris Fuller March 2013

import os
import numpy as np  
from os.path import join as pj
import glob as g
import shutil as sh
###################### inputs ######################
folder = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/output/bad"
srcFolder = "/Users/chrisfuller/Desktop/current/PSW/"
destFolder = "/Users/chrisfuller/Desktop/current/check-coords/"

def folder2list(folder):
    galaxies = [s for s in os.listdir(folder) if "postscript" in s]
    l =[]
    for gal in galaxies:
        l.append(int(gal.split("-")[0]))
    l = np.array(l)
    l.sort()
    return l
 # for saving an array to a text file .csv
def savetext(folder,filename,array):
    outfile = open(pj(folder, filename), 'w')
    for element in array:
        x = str(element)
        outfile.write(x)
        outfile.write("\n")
    outfile.close()
    
def movefiles(srcList,dst):
    for src in srcList:
        print src.split("/")[-1]
        sh.copy(src,pj(dst,src.split("/")[-1]))
        
    
a = folder2list(folder)
savetext("/Users/chrisfuller/Desktop/","bad-gals.txt",a)

#for gal in a:
    #src = g.glob(pj(srcFolder,str(gal)+"*"))
    #movefiles(src,destFolder)
# Program to create region file for apature photomerty of galaxys with funtools in ds9
#Chris Fuller
#dec 2011
#edited may2012


# import stuff
from numpy import *
import numpy
import scipy
import math
import sys
import os
from scipy.optimize import fmin
from os.path import join as pj

aloc={}
beam = 13.64
rap = beam

outer = sqrt(2*(beam**2))

#file stuff
firstFile = "bigcoma.csv"
folder = "/Users/chrisfuller/Dropbox/coma/Catalogues"
outfolder = "/Users/chrisfuller/Dropbox/coma/regions"
# outfiles
outName = "App160.reg"

outBack = "App160back.reg"
# write out overall list at end as a .reg file
outFile = open(pj(outfolder,outName), 'w')
outFile.write("# Region file format: DS9 version 4.1\n")
outFile.write("# Filename:  PLWmap-v51.fits[image]\n")
outFile.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
outFile.write("fk5\n")
outFile.write("\n")
out = open(pj(outfolder,outBack), 'w')
out.write("# Region file format: DS9 version 4.1\n")
out.write("# Filename:  PLWmap-v51.fits[image]\n")
out.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
out.write("fk5\n")
out.write("\n")
firstFile = open(pj(folder,firstFile), 'r')

# start count of sorces at zero
count = 0
# read in all lines
for line in firstFile.readlines():
    # skip lines that are not useful
    if line[0] == "o":
        continue
    # sepearte with deliminator
    info = line.split(",")    
    a = float(info[9])*60.0
    #rap = sqrt(a**2 + beam**2)
    count += 1    
    aloc[count] = {"ra":(info[1]), "dec":(info[2])}        
    # make new line
    line = "circle(" + str(aloc[count]["ra"]) + "," + str(aloc[count]["dec"])+ "," + str(rap) +"\")" +"\n"
    background = "annulus(" + str(aloc[count]["ra"]) + "," + str(aloc[count]["dec"])+ "," + str(rap) + "\"," + str(rap*2) +"\")" +"\n"
    # write new line
    outFile.write(line)
    out.write(background)        
# close output
outFile.close()
out.close()
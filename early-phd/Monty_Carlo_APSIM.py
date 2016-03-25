#Program by Chris 'Hacker' Fuller Feb 2012
#Program to randomly lay down apatures over the ngp atlas field

# import stuff
from numpy import *
import numpy
import scipy
import math
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt
import pyfits
import pywcs
import random

#File Stuff
folder = "/Users/chrisfuller/Dropbox/coma/regions/"
outName = "coma-cluster.csv"
regName = "monty_carlo_250_a.reg"

#Inputs
beam = 18.2/3600.0
numApp = 752


#####################################################################################################




# write out overall list at end
regFile = open(pj(folder,regName), 'w')
regFile.write("# Region file format: DS9 version 4.1\n")
regFile.write("# Filename:  PLWmap-v51.fits[image]\n")
regFile.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
regFile.write("physical\n")
regFile.write("\n")


#Read in fits file
hdulist = pyfits.open("/Users/chrisfuller/Documents/coma/PSW.fits")

header = hdulist[0].header
wcs =  pywcs.WCS(hdulist[0].header)
wcs.wcs.print_contents()

#Working out beam size
pixScale = header['CDELT2']
appRad = beam/pixScale


print ' reading in fits data'

ngpMap = hdulist[0].data

print 'Making Random Apature Location'

count = 0
while count < numApp:
       
    X = random.randint(0,ngpMap.shape[1])
    Y = random.randint(0,ngpMap.shape[0])
    
    if ngpMap[Y,X] == 0.0:
        continue
    else:
        # write new line
        region = "circle(" + str(X) + "," + str(Y) + "," + str(appRad) +")" +"\n"
        regFile.write(region)
    
        count += 1
        
regFile.close()
print 'finished'
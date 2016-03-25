
#Program written by Chris Fuller to match and find fluxes from a cataloge of 
#optically selected galaxies from the NGP atlas maps 
#Feb 2012


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

#File Stuff
folder = "/export/home/spx6cff/hefocs/coma/"
catalogeFile = "ngpfluxcat.csv"
sourceList = "coma-goldmine.csv"
outName = "fluxcat.csv"
#Varibles

threshold = float(5.0/3600.0)

##########################################################################################################################################

# write out overall list 
outFile = open(pj(folder,outName), 'w')

outFile.write("ra,dec,F250,F350,F500\n")


#into an array

fluxCat = loadtxt((folder+catalogeFile), dtype=float, delimiter=",", skiprows=1, usecols=(4,5,7,8,9), unpack=False)

sourceCat = loadtxt((folder+sourceList), dtype=float, delimiter=",", skiprows=1,usecols=(1,2), unpack=False)

badgal = 0
count = 0
nogal = 0


finalfluxcat = numpy.array(([0,0,0,0,0], [0,0,0,0,0]),  dtype = float)

for row in fluxCat:
    
    selection = numpy.where(((abs(row[0] - sourceCat[ : ,0])*math.cos(row[1]))**2.0 + abs(row[1] - sourceCat[ : ,1])**2.0 < (threshold)**2.0))
    count += 1


    ### see if we find match or not
        # see if have more than one
    if len(selection[0]) > 1:
        badgal += 1             
        print "Multiple object found for ", row[0:1], " trying halving threshold"
        selection = numpy.where(((abs(row[0] - sourceCat[ : ,0])*math.cos(row[1]))**2.0 + abs(row[1] - sourceCat[ : ,1])**2.0 < (threshold)**2.0))
        
        if len(selection[0]) != 1:
            print "Still getting multiple objects for ", row[0:1]
            continue
    
    # if only one match continue
    if len(selection[0]) == 1:        
        finalfluxcat = concatenate((finalfluxcat,[row]), axis=0)
    # if no matches add source to list
    elif len(selection[0]) == 0:
        nogal += 1
        continue

finalfluxcat = numpy.delete(finalfluxcat, (0), axis=0)

savetxt("/export/home/spx6cff/hefocs/coma/finalfluxcatsavetext.csv", finalfluxcat, delimiter = ",")

for line in finalfluxcat:

    # make new line
    line = str(line[0]) + " , " + str(line[1]) + " , " + str(line[2]) + " , " + str(line[3]) + " , " + str(line[4]) + "\n"
    
    # write new line
    outFile.write(line)

# close output
outFile.close()

print "program complete"    
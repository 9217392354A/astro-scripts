
# Program to compare points source catalogs from different scans
# written by Chris Fuller & Matt Smith Nov 2011

# import stuff
from numpy import *
import numpy
import scipy
import math
import sys
import os
from os.path import join as pj



# define threshold
threshold = 15.0 / 3600.0

# outfile

folder = "/Users/chrisfuller/Dropbox/phd/herchel/virgo/offsets/V1"
field = folder.split("/")
outName = "master-catlog-" + field[-1] + ".csv"
##################################################################################

# search folders
files = os.listdir(folder)

# Check all maps are fits files
files = [x for x in files if x[-10:] == "-input.csv"]

# open first file
firstFile = open(pj(folder,files[0]), 'r')

# make an empty dictionary to store data
master = {}

# start count of sorces at zero
count = 0
badgal = 0
# read in all lines
for line in firstFile.readlines():
    # skip lines that are not useful
    if line[0] == "r" or line[0] == "D" or line[0] == "," or line[0] == 'i' or line[0] == '#' or line[0] == 'd':
        continue
    
    # sepearte with deliminator
    info = line.split(",")
    
    # save info
    count += 1
    master[count] = {"ra":double(info[0]), "dec":double(info[1])}

# close file
firstFile.close()

# loop through and make an array of ra and dec
ra = []
dec = []
for key in master.keys():
    ra.append(master[key]["ra"])
    dec.append(master[key]["dec"])
# convert arrays to numpy array
ra = numpy.array(ra)
dec = numpy.array(dec)

# loop over rest of files
for i in range(1,len(files)):
    # open file
    newFile = open(pj(folder,files[i]), 'r')
    
    # read in all lines
    for line in newFile.readlines():
        # skip lines that are not useful
        if line[0] == "r" or line[0] == "D" or line[0] == "," or line[0] == 'i' or line[0] == '#' or line[0] == 'd':
            continue
    
        # sepearte with deliminator
        info = line.split(",")
        
        lineRA = double(info[0])
        lineDEC = double(info[1])
        
        # see if exist in current list
        selection = numpy.where(((abs(lineRA - ra)*math.cos(math.radians(lineDEC)))**2.0 + abs(lineDEC - dec)**2.0 < (threshold)**2.0))
        
        ### see if we find match or not
        # see if have more than one
        if len(selection[0]) > 1:
            badgal += 1             
            print "Multiple object found for ", lineRA, lineDEC, " trying halving threshold"
            selection = numpy.where(((abs(lineRA - ra)*math.cos(math.radians(lineDEC)))**2.0 + abs(lineDEC - dec)**2.0 < (threshold/2.0)**2.0))
            if len(selection[0]) != 1:
                #raise "Still getting multiple objects for ", lineRA, lineDEC
                print "Still getting multiple objects for ", lineRA, lineDEC
                count +=1
                continue
        # if only one match continue
        if len(selection[0]) == 1:
            continue
        # if no matches add source to list
        elif len(selection[0]) == 0:
            # save info
            count += 1
            master[count] = {"ra":double(info[0]), "dec":double(info[1])}

    # close file
    newFile.close()
    
    # make new overall ra and dec files
    ra = []
    dec = []
    for key in master.keys():
        ra.append(master[key]["ra"])
        dec.append(master[key]["dec"])
    # convert arrays to numpy array
    ra = numpy.array(ra)
    dec = numpy.array(dec)
            
# write out overall list at end
outFile = open(pj(folder,outName), 'w')

outFile.write("id, ra, dec \n")

# loop over each line in master dictionary
keys = master.keys()
keys.sort()
for key in keys:

    # make new line
    line = str(key) + " , " + str(master[key]["ra"]) + " , " + str(master[key]["dec"]) + " \n"
    
    # write new line
    outFile.write(line)

# close output
outFile.close()

print "Program Finished Successfully"
print "Number of galaxies found ", len(master)
print "Number of multiple objects found = ", badgal
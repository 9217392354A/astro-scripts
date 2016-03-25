# Program to compute offsets for catalog created by CatComp
# written by Chris Fuller & Matt Smith Nov 2011

# import stuff
from numpy import *
import numpy
import scipy
import math
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt

# seclect folder
#folder = "/Users/chrisfuller/Dropbox/phd/herchel/virgo/offsets/V4"
folder = "/Users/chrisfuller/Dropbox/phd/herchel/virgo/offsets/V4"
field = folder.split("/")

sdssfile = "SDSS-" + field[-1] + ".csv" 

#name_file = "e.csv" #if file name Brigade.csv
name_file = "s.csv" #if file name fits.csv



# define threshold
threshold = (15.0 / 3600.0)

# outfile
outName = "offsets" + field[-1] + ".txt"

# write out 
outFile = open(pj(folder,outName), 'w')
outFile.write("scan, mean ra, mean dec, median ra, median dec \n")



############################################################################

#Import and create dictionary of sdss locations 

#create empty dictionary to store sdss locatins 
sdsslocations = {}

sdssMaster = open(pj(folder,sdssfile), 'r')

count = 0 

# read in all lines
for line in sdssMaster.readlines():
    # skip lines that are not useful
    if line[0] == "R":
        continue
    
    # sepearte with deliminator
    info = line.split(",")
    
    count += 1
    
    sdsslocations[count] = {"ra":double(info[0]), "dec":double(info[1])}

# close file
sdssMaster.close()

print 'sdss locations read in; len = ', len(sdsslocations)

############################################################################
offset = {}
k=0
badgal = 0

offSetVal_RA =  double()
offSetVal_DEC = double()

# loop through and make an array of ra and dec
ra = []
dec = []
for k in sdsslocations:
    ra.append(sdsslocations[k]["ra"])
    dec.append(sdsslocations[k]["dec"])
# convert arrays to numpy array
raSdss = numpy.array(ra)
decSdss = numpy.array(dec)


#read in the scan files
# search folder
files = os.listdir(folder)
fig = plt.figure(figsize = (11,8))
# Check all maps are fits files
files = [x for x in files if x[-5:] == name_file]

# loop over rest of files
c = 1
for i in range(0,len(files)):
    # open file
    newFile = open(pj(folder,files[i]), 'r')
    scanname = files[i].split("-")[2]
    offSet_RA = []
    offSet_DEC = []
    count2 = 0
    # read in all lines
    
    for line in newFile.readlines():
        # skip lines that are not useful
        if line[0] == "r" or line[0] == "D" or line[0] == "," or line[0] == 'i':
            continue
    
        # sepearte with deliminator
        info = line.split(",")
        
        lineRA = double(info[0])
        lineDEC = double(info[2])
        
        # see if exist in current list
        selection = numpy.where(((abs(lineRA - raSdss)*math.cos(math.radians(lineDEC)))**2.0 + abs(lineDEC - decSdss)**2.0 < (threshold)**2.0))
        
        ### see if we find match or not
        # see if have more than one
        if len(selection[0]) > 1:
            print "Multiple object found for ", lineRA, lineDEC, " trying halving threshold"
            selection = numpy.where(((abs(lineRA - ra)*math.cos(math.radians(lineDEC)))**2.0 + abs(lineDEC - dec)**2.0 < (threshold/2.0)**2.0))
            if len(selection[0]) != 1:
                print "Still getting multiple objects for ", lineRA, lineDEC
                continue
        # if only one match continue
        if len(selection[0]) == 1:
            offSetVal_RA =  ((raSdss[selection][0] - lineRA)*3600.0)
            offSetVal_DEC = ((decSdss[selection][0] - lineDEC)*3600.0)
            
            if abs(offSetVal_RA) < (threshold*3600) and abs(offSetVal_DEC) < (threshold*3600) :
                # save info
                count2 += 1
                offSet_RA.append((raSdss[selection][0] - lineRA)*3600.0)
                offSet_DEC.append((decSdss[selection][0] - lineDEC)*3600.0)
                continue
            elif abs(offSetVal_RA) < (threshold*3600):
                print "error offset greater than threshold in RA offset size = ", offSetVal_RA, "at ", selection
                continue
        # if no matches add source to list
        elif len(selection[0]) == 0:
            # save info
            #print "source removed in pre-processing", lineRA, lineDEC
            badgal = badgal + 1
   
    # first create a single histogram
    meanRA = mean(offSet_RA)
    meanDEC = mean(offSet_DEC)
    medianRA = median(offSet_RA)
    medianDEC = median(offSet_DEC)
    
    #Plotting Graph's   
    #f1 = plt.axes([0.1,0.1,0.85,0.85]) 
    #f2 = plt.axes([0.1,0.1,0.85,0.85]) 
    f1 = plt.subplot(4,4,c)
    f2 = plt.subplot(4,4,(c+1))
    c += 2
    # overplot mean and median
    n, bins, patches = f1.hist(offSet_RA, 40, normed=False)
    n, bins, patches = f2.hist(offSet_DEC, 40, color = "green", normed=False)
    #plt.xlabel('Offset RA, [arcsec] ')
    #plt.title("Scan = V1-" + scanname)
    f1.plot([meanRA, meanRA], [0,100], 'r--', linewidth=3.0)
    f1.plot([medianRA, medianRA], [0,100 ], '--', color='black', linewidth=3.0) 
    f2.plot([meanDEC, meanDEC], [0,100], 'r--', linewidth=3.0)
    f2.plot([medianDEC, medianDEC], [0,100 ], '--', color='black', linewidth=3.0)     
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    f1.set_ylim(0,35)
    f2.set_ylim(0,35)
    f1.set_xlim(-5,5)
    f2.set_xlim(-5,5)  
    f2.tick_params(axis='y', labelleft='off', labelright='off')
    if c != 3 and c != 7 and c != 11 and c != 15:
        f1.tick_params(axis='y', labelleft='off', labelright='off')
    f1.text(2.3,19.0,scanname, fontsize=15)
    # make new line
    line = (scanname + ", " + str(meanRA) + ", " + str(meanDEC) + ", " + str(medianRA) + ", " + str(medianDEC) + "\n")
    # write new line
    outFile.write(line)
    
    # close file
    newFile.close()
print "number of sources that where removed in pre-proccessing = ", badgal    
# close output
outFile.close()
fig.text(0.43,0.04,'Offset (arcsec)',fontsize=15,color='black')
fig.text(0.04,0.55,'Number in Bin',fontsize=15,color='black', rotation='vertical')
fig.savefig(pj(folder,"offset-" + field[-1] + ".png"))
plt.show()


print "Program finished!"
#print sdsslocations
#print sdssMaster

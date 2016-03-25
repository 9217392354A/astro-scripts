#Program written by Chris Fuller to find if objects are in a fits field.
#Jan 2012


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
#imagefolder = "/home/herdata/spxrra/myoutput/atlas/ngp/maps/ngp-mosaic-240311-PSWmap-v51.fits
#fitsfile = "ngp-mosaic-240311-PSWmap-v51.fits"
catalogefile = "comaclustersuper.csv" 

outName = "coma-cluster.csv"
regName = "comaregionsPLW.reg"
fulcat = "coma-cluster-gmformat.csv"



beam = 36.3/60.0
frac = 0.3 #fraction of semi-maiour axis

############################################################################
cat = loadtxt(pj(folder,catalogefile), dtype=str, delimiter=",", unpack=False)
#catHeader =loadtxt(pj(folder,catalogefile), dtype=str, delimiter=",", unpack=False)


finalcat = numpy.array(cat[0])
finalcat = numpy.reshape(finalcat, (1,24))


#Import and create dictionary of objects locations 





#create empty dictionary to store  locations 
locations = {}

master = open(pj(folder,catalogefile), 'r')

count = 0 

# read in all lines
for line in master.readlines():
    # skip lines that are not useful
    if line[0] == "O":
        continue
    
    # seperate with deliminator
    info = line.split(",")

    
    locations[count] = {"object":str(info[0]), "ra":float(info[6]), "dec":float(info[7]), "a":str(info[9]), "b":float(info[10]), "distance":float(info[11]),"type":int(info[13])}
    count += 1
# close file
master.close()

print 'locations read in; len = ', len(locations)
print 'idiot test; first location = ', locations[15]

# loop through and make an array of ra and dec
ra = []
dec = []
galname = []
a = []
b = []
distance = []
galtype = []


for key in locations.keys():
    ra.append(locations[key]["ra"])
    dec.append(locations[key]["dec"])
    galname.append(locations[key]["object"])
    a.append(locations[key]["a"])
    b.append(locations[key]["b"])
    distance.append(locations[key]["distance"])    
    galtype.append(locations[key]["type"])
# convert arrays to numpy array
ra = numpy.array(ra)
dec = numpy.array(dec)
galname = numpy.array(galname)
a = numpy.array(a)
b = numpy.array(b)
distance = numpy.array(distance)
galtype = numpy.array(galtype)

print 'reading in fits file....'

#Read in fits file

hdulist = pyfits.open("/home/herdata/spxrra/myoutput/atlas/ngp/maps/ngp-mosaic-240311-PSWmap-v51.fits")
#hdulist = pyfits.open("/export/home/spx6cff/hefocs/PLW-mask.fits")

header = hdulist[0].header
wcs =  pywcs.WCS(hdulist[0].header)
wcs.wcs.print_contents()

mask = hdulist[0].data

found = 0
notfound = 0
count = 0 
pointcount = 0
print "starting galaxy check......"

# write out overall list at end
outFile = open(pj(folder,outName), 'w')
fulcatFile = open(pj(folder,fulcat), 'w')

outFile.write("#object, ra, dec, a, b, distance, type \n") 

# write out overall list at end
regFile = open(pj(folder,regName), 'w')
#ngp-mosaic-240311-
regFile.write("# Region file format: DS9 version 4.1\n")
regFile.write("# Filename:  PLWmap-v51.fits[image]\n")
regFile.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
regFile.write("fk5\n")
regFile.write("\n")

for count in range(len(ra)):
    skyRA = ra[count]
    skyDEC = dec[count]
    galName = galname[count]
    agal = a[count]
    bgal = b[count]
    distancegal = distance[count]
    GALTYPE = galtype[count]
    
 
    skycrd = numpy.array([[skyRA,skyDEC]],numpy.float_)
          
    pixcrd = wcs.wcs_sky2pix(skycrd, 1)
    count += 1 
    
    X = pixcrd[0,0]
    Y = pixcrd[0,1]
    
    if Y > ((-175.0 / 1248.0)*X + 1394.3) and Y > ((9072.0 / 1114.0)*X -81754.4): 
        found += 1
        if 0.5*(float(agal)) > beam:
            point = 0
            
        else:
            point = 1
            pointcount +=1
        
        x = (float(agal) + (float(bgal)))/2.0
        line = galName + "," + str(skyRA) + "," + str(skyDEC) +"," + str(agal) + "," + str(bgal) + "," + str(distancegal) +"," + str(GALTYPE) +"\n"
        # write new line
        region = "circle(" + str(skyRA) + "," + str(skyDEC) + "," + str(beam) +"\') # text={" + str(found) + "}" +"\n"
        
        outFile.write(line)
        
        regFile.write(region)
        catout = numpy.reshape(cat[count], (1,24))

        
        finalcat = concatenate((finalcat,catout), axis=0)  
    else:
        notfound += 1
        
fraction = (float(found)/float(count))*100
# close output
outFile.close() 
regFile.close() 
#fulcatFile.close()

#savetxt("/export/home/spx6cff/hefocs/coma/finalcatsavetext.txt", finalcat)
#fulcatFile.write(catout)

for row in finalcat:
    
    for element in reshape(row, 24):
        x = element + ","
        fulcatFile.write(x)
    y = "\n"
    fulcatFile.write(y)  
        
fulcatFile.close() 

print "Program Finished Successfully"
print 'total = ', count, 'found = ', found, 'not found = ', notfound, 'number of point sources = ', pointcount     
print round(fraction) , '% in field'
hdulist.close()

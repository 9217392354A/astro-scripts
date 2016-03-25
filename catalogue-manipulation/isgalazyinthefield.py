#Program written by Chris Fuller to find if objects are in a fits field.
#Jan 2012
# import stuff
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
folder = "/export/home/spx6cff/hefocs/"
fitsfile = "PLW-mask.fits"
catalogefile = "a1367withvelcut.csv" 
############################################################################
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
    locations[count] = {"ra":float(info[5]), "dec":float(info[6])}
    count += 1
# close file
master.close()

# loop through and make an array of ra and dec
ra = []
dec = []
for key in locations.keys():
    ra.append(locations[key]["ra"])
    dec.append(locations[key]["dec"])
# convert arrays to numpy array
ra = numpy.array(ra)
dec = numpy.array(dec)

#Read in fits file
hdulist = pyfits.open("/export/home/spx6cff/hefocs/PLW-mask.fits")
header = hdulist[0].header
wcs =  pywcs.WCS(hdulist[0].header)
#wcs.wcs.print_contents()
mask = hdulist[0].data
count = 0 
print "starting galaxy check......"

for count in range(len(ra)):
    skyRA = ra[count]
    skyDEC = dec[count] 
    skycrd = numpy.array([[skyRA,skyDEC]],numpy.float_)          
    pixcrd = wcs.wcs_sky2pix(skycrd, 1)
    count += 1     
    pixval = mask[pixcrd[0,1],pixcrd[0,0]]    
    print pixval
print count    
hdulist.close()
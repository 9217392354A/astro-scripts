#Program written by Chris Fuller to find if objects are in a fits field.
#Nov 2013
# import stuff
print 'importing stuff....'
import numpy as np
from os.path import join as pj
import pyfits 
import pywcs 
from atpy import Table
import scipy.signal as s
import matplotlib.pyplot as plt
from math import isnan

#Inputs
folder = "/Users/chrisfuller/Desktop/"
map_name = 'NGP-entire-pacs160-mosaic-20120512.fits'
cat  = Table(pj(folder,"ngp_input_ready.fits" ))
outname = "ngp_input_pacsready.fits"

 

delta = 100

# extract ra
RA = cat.GRA2000
DEC = cat.GDEC2000
OBJ = cat.OBJECT

ra = np.array(RA)
dec = np.array(DEC)


#Read in fits file
print 'reading in data...'
hdulist = pyfits.open(pj(folder, map_name))
header = hdulist[0].header
wcs =  pywcs.WCS(hdulist[0].header)
#wcs.wcs.print_contents()
im = hdulist[0].data

#run through table and check if galaxies are inside the map.

print "starting galaxy check......"
counter  = 0 
#create list to store the values in
in_pacs = []
for i in range(0,len(ra)):
    #put sky cordinates into a single array
    skycrd = np.array([[ra[i],dec[i]]],np.float_)

    #get pixel cordinates from sky coordinates  
    pixcrd = wcs.wcs_sky2pix(skycrd, 1)
    x = pixcrd[0,0]
    y = pixcrd[0,1]

    #extract array
    pixval_array = im[ y-delta : y+delta , x-delta : x+delta ]

    #flatten array
    pixval_array_flat = pixval_array.flatten()

    #select all pixels that are nans
    w = np.where(np.isnan(pixval_array_flat))[0]

    val = 1    
    # print out full details of any galaxies that have nanage
    if len(w) != 0:
        counter += 1
        print OBJ[i], ' Extracted array size: ' , len(pixval_array_flat), pixval_array.shape , ' No. NAN: ', len(w), ' Sky: ' , skycrd , ' Pix: ' , pixcrd
        
        #change HEVICS_PLW to 0
        val = 0

    in_pacs.append(val)



cat.add_column('in_pacs', np.array(in_pacs, dtype=np.int16), dtype=np.int16)

print 'number of galaxies outside map: ', counter
cat.write(pj(folder, outname), overwrite=True)
hdulist.close()
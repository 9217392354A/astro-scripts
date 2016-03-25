# program to create mask from fits file
# Chris Fuller, August 2013

from astropy.io import fits
from os import chdir
import numpy as np

#change to sss dir
chdir('/Users/chrisfuller/Documents/maps/fornax/')

#read in fits files
hdulist = fits.open('HeFoCS-All-PLWmap-mosaic_MS-DR35.fits')

#extract image and hdr data
image_nan = hdulist[1].data
hdr = hdulist[1].header

# now add 1000
image = np.nan_to_num(image_nan)

#now find all values that are not equal to 0
w1 = np.where(image != 0.)[0]
print np.sqrt(len(w1))

print np.mean(image), np.mean(image[w1])
#fix all pixels that arn't 0 to 1
image[w1] = 1





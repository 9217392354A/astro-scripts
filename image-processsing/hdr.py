#header finder

import pyfits
import os
import numpy as np  
from os.path import join as pj

folder = '/home/scuba2/spx6cff/coma/maps/'

imps = [s for s in os.listdir(folder) if ".fits" in s]

for im in imps:

	print '-'*80.0
	#Read in fits file
	hdulist = pyfits.open(pj(folder,"HeFoCS-All-PLWmap-mosaic_MS-DR35.fits"))
	header = hdulist[0].header
	
	hdulist.close()

	print header

#program that reads in results from Matt's SED fitter and 
# then creates a unfied dust mass cloumn

import numpy as np  
from os.path import join as pj
from atpy import Table
import matplotlib.pyplot as plt
import pylab as pl

cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/sed-fits/sed-all.fits")

bands = ["F500","F350","F250","F160","F100"]
""" Create detected column """
cat.remove_columns(['DETECTED']) 

""" Loop through cat bands and add number to detected col """

for i in range(len(bands)):
	band = bands[i]

	cat.add_column('D' + band[-3:], [0]*len(cat))

	flux = np.nan_to_num(cat[band])

	w = np.where(flux != 0.0)[0]

	cat['D' + band[-3:]][w] = 1

total = cat.D500 + cat.D350 + cat.D250 + cat.D160 + cat.D100


five_cluster = cat.where((total == 5) & (cat.RADIUS_VIR <= 1.0))
five_filament = cat.where((total == 5) & (cat.RADIUS_VIR > 1.0))

""" create new columum dmass and dmass_type """

dmass, dmass_type = np.array([0.0]*len(cat)) , np.array([0]*len(cat))

w_five = np.where(total == 5)[0]
w_250 = np.where((total != 5) & (total != 0))[0]

dmass[w_five] = cat.DMASS_SED[w_five]
dmass_type[w_five] = 2
dmass[w_250]  = cat.DMASS_250[w_250]
dmass_type[w_250] = 1


cat.add_column('DMASS', dmass)
cat.add_column('DMASS_TYPE', dmass_type)

cat.write('/Users/chrisfuller/Dropbox/phd/herchel/coma/sed-fits/test.fits',  overwrite=True)

cluster = cat.






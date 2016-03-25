# program to work out the counts in a bin from SExtractor for background-contamination-SPIRE.py
# Chris Fuller, August 2013

import numpy as np
from os.path import join as pj
from os import listdir
from atpy import Table
import matplotlib.pyplot as plt
import pylab as pl

######################################### i/o ####################################################
folder = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/SPIRE-contamination/contamination/'

psw = Table(pj(folder,'250.fits'))
pmw = Table(pj(folder,'350.fits'))
plw = Table(pj(folder,'500.fits'))

cats = [psw, pmw, plw]
bin_s = [15.0,  30.0, 100.0, 1000.0]
###################################### functions ##################################################
def my_print(x):
	Np = 'Np = np.array(['
	for i in range(0,len(x)):
		if i != len(x)-1: Np = Np + str(x[i]) + ','
		else: Np = Np + str(x[i])
	print Np + '])'
###################################### control ####################################################

#conversion factors
#X_jy = [ (1.0**2)/423., (1.0**2)/751., (1.0**2)/1587.]

X_jy = [ (6.0**2)/423., (8.0**2)/751., (12.0**2)/1587.]
min_size = [12.0,]

#X_jy = [ 1., 1., 1.]
#loop through all 3 cats
for i in range(0, len(cats)):
	#print i, range(0, len(cats))
	cat = cats[i]

	#find out where values are above and area of 12
	w1 = np.where(cat.ISOAREA_IMAGE >= 12.)[0]
	#print "number greater than 12 pixels in area: ", len(w1)
	
	#select all fluxes that have and area of above 12 
	flux = cat.FLUX_AUTO[w1]*X_jy[i]*1000.0

	#bin it up according to binnarge
	hist, bin_edges = np.histogram(flux, bins= bin_s)

	#hist to show visually whats going on!
	#pl.hist(flux, bins= bin_s)#np.logspace(1,2))
	#pl.semilogx()
	#pl.show()
	
	# convert into per square degree
	hist = hist / 20.0


	#print bin_edges
	print my_print(hist)
	#print hist





	

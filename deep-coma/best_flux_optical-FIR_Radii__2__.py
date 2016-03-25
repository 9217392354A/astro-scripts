# program for removing sources bassed on the beam size, to be run after best_flux_point-vs-extended.py
# program written by Chris Fuller, Jan 2014

#import modules
import numpy as np
from os.path import join as pj
import atpy as at
import matplotlib.pyplot as plt
from copy import copy, deepcopy

#folder stuff
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/catalogue_creation_phase/'
catname = 'ngp+20140331__best-flux__1.fits'
outname = 'ngp+20140331__best-flux-3.5arcsec__1-2__.fits'
#load catalogues one where all are forced as point to serve as the base catalogue
cat = at.Table(pj(folder,catname),type='fits')

#bands
#bands = ['500','350','250','160','100']
bands = ['250']
beams =  [36.0, 24.5, 18.2, 13.4, 9.4] #beams in arcsec
colour = ['b', 'r', 'k', 'g', 'c']

# # # # # # # # # # # # # # # # # # # # # # # # Functions # # # # # # # # # # # # # # # # # # # # # # # #
def remove_single_from_cat(index,t,band):
	cols = ['F', 'SN', 'R']
	#now times the error by 3 for all galaxies removed
	if t['F'+band][index] != 0.0:
		#print 'OBJECT: ', t['OBJECT'][index], 'BAND: ', band, 'FLUX: ', t['F'+band][index], 'ERROR: ', t['E'+band][index]
		t['E'+band][index] = t['E'+band][index]*3.0

	#loop through each col
	for col in cols:
		#set col to zero
		t[col+band][index] = 0.0 

# # # # # # # # # # # # # # # # # # # # # # # # Functions # # # # # # # # # # # # # # # # # # # # # # # #
new = deepcopy(cat)
counts = []

#repeat analsis for each band
for j in range(0, len(bands)):
	band = bands[j]
	beam = beams[j]
	print band
	print ''
	count = 0
	#work rowwise down each band
	for i in range(0,len(cat)):
		#fluxes
		F_cat = cat['F'+band][i]
		
		#extendedness
		sep = cat['Separation'][i]

		if sep == '': 
			sep = 999.9
			new['Separation'][i] = 999.9
		
		#if detected and if source is greater than beam/2.0 from source remove from cat
		if (F_cat != 0.0) and (sep >= 3.5):
			remove_single_from_cat(i,new,band)
			count +=1
	counts.append(count)

print bands
print counts
new.write(pj(folder, outname), overwrite=True)

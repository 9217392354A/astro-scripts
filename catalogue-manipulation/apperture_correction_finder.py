#program to find the median aperture correction

import os
import numpy as np  
from os.path import join as pj
import atpy as at

folder = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/'
cat = at.Table(pj(folder, 'stellar-mass-fornax_final.fits'))

bands = [500, 350, 250, 160, 100]

for i in range(0,len(bands)):
	band = bands[i]
	temp = []
	for j in range(0,len(cat)):
		if cat['EXTENDEDNESS'][j][i] == 'P' or cat['R'+str(band)][j] == 0.0: continue
		temp.append(cat['R'+str(band)][j])

	print band, ' mean: ', np.mean(np.array(temp, dtype = np.float))*60.0, ' median: ', np.median(np.array(temp, dtype = np.float))*60.0, ' number of galaxies extended: ', len(np.array(temp, dtype = np.float))








# program to take my cat and make it into a format that Matt's sed fitter likes
# Chris Fuller',' January 2014

import numpy as np
from os.path import join as pj
import atpy as at

num_bands = 5.0 #number of bands that we require a detection in

#headers for new cat
headers = np.array(['OBJECT','z','z_error','distance','distance_error','PACS100_w','PACS100','PACS100_e',\
'PACS160_w','PACS160','PACS160_e','SPIRE250_w','SPIRE250','SPIRE250_e','SPIRE350_w','SPIRE350','SPIRE350_e',\
'SPIRE500_w','SPIRE500','SPIRE500_e'], dtype = np.str).reshape(1, 20)

#open old cat
cat = at.Table('/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/cluster+filament-mybg-130114-3.1arcsec.fits', type='fits')

#bands
bands =['F500', 'F350', 'F250', 'F160', 'F100']


#loop through all the galaxies in the cat
for i in range(0, len(cat)):
	dtest = 0.
	#find if galaxy is detected in at least 3 bands
	for band in bands: 
		if cat[band][i] > 0.: dtest += 1.

	if dtest >= num_bands: 
		#add
		print 'add', dtest

	else:
		continue 


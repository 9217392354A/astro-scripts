#program to move a catalogue to the distance of coma
# Chris Fuller March 2014

#import modules
import numpy as np
from os.path import join as pj
import atpy as at
from copy import copy, deepcopy
import matplotlib.pyplot as plt
from pylab import bar
import pdb

#pdb.set_trace()

#i/o
print 'reading in cat . . .'
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data'
cat = at.Table(pj(folder, 'fornax_input.fits'),type='fits')
output = 'fornax_at_100mpc-030314.fits'

#key parameters
#coor_names = ['RA (2000)', 'DEC (2000)'] # these are the colum names that containe ra and dec ### Virgo  ###
coor_names = ['GRA2000', 'GDEC2000'] # these are the colum names that containe ra and dec ####### Fornax ###
optical_col = 'BTmag_1'
flux_cols = ['F100', 'F160', 'F250', 'F350', 'F500' ]

optical_lim = 14.89 # faintest magnitude that is possible to select at the distance of the coma cluster
x = 0.30 #conversion between deg and mpc
dist_x = 0.0289#scale fluxes



# conversion between degrees to mpc 
#coma x = 1.77
#virgo x= 0.25
#fornax x=0.30

#flux scales
#coma   = 1.0
#virgo  = 0.0196
#fornax = 0.0289

# # # # # # # # # # # # # # # Function # # # # # # # # # # # # # # # # # # # # # # # # # #
#function to produce new cat with column added for the nth nearest neigboure
def nth_nearest_neighbour(t, coor_names):
	print 'nth_nearest_neighbour....'

	#add columnd for D1,D5, and D10
	t.add_empty_column('D1', dtype = np.float)
	t.add_empty_column('D5', dtype = np.float)
	t.add_empty_column('D10', dtype = np.float)
	t.add_empty_column('SIGMA1', dtype = np.float)
	t.add_empty_column('SIGMA5', dtype = np.float)
	t.add_empty_column('SIGMA10', dtype = np.float)
	
	###### part 2 #######
	# find nearest neighbours

	#ra1 and dec1
	ra_1 = t[coor_names[0]]
	dec_1 = t[coor_names[1]]

	#loop through all members of catA
	for i in range(0, len(t)):

		ra = t[coor_names[0]][i]
		dec = t[coor_names[1]][i]

		#ra2 and dec2
		ra_2 = np.array([ra]*len(ra_1), dtype=np.float)
		dec_2 = np.array([dec]*len(ra_1), dtype=np.float)
		
		#caculate distance to all sources from ra1 and dec1
		radius = np.sort(distance(ra_1, dec_1, ra_2, dec_2 ))

		#print radius[1]*1.77*1000.0, np.min(radius)

		#add values to table
		t['D1'][i] = radius[1] * x
		t['D5'][i] = radius[5] * x
		t['D10'][i] = radius[10]* x 
		t['SIGMA1'][i] = np.log10(1.0 / (np.pi*(radius[1]*x)**2.0) )
		t['SIGMA5'][i] = np.log10(5.0 / (np.pi*(radius[5]*x)**2.0) )
		t['SIGMA10'][i] = np.log10(10.0 / (np.pi*(radius[10]*x)**2.0)) 


	return t

#distance equation designed to do arraywise caculations
def distance(ra1, dec1, ra2, dec2):
	delta_ra = (ra1 - ra2) * np.cos(np.radians((dec1+dec2)/2.0))
	delta_dec = (dec1 - dec2)

	return np.sqrt(delta_ra**2.0 + delta_dec**2.0)

 # # # # # # # # # # # # # # # Main Program # # # # # # # # # # # # # # # # # # # # # # # #

#scale fluxes so as to appear at the distance of coma
for i in range(len(flux_cols)):
	col = flux_cols[i]

	#scale to the distance of coma
	cat[col] = cat[col]*dist_x 

	#if less than 15mjy then set to 0
	w = np.where(cat[col] < 0.015)[0]
	cat[col][w] = 0.0


#make an optical selection for the cluster 
optical = cat.where((cat[optical_col] <= optical_lim) & (np.nan_to_num(cat[optical_col]) != 0.0))
new_cat  = nth_nearest_neighbour(optical, coor_names)






new_cat.write(pj(folder, output), overwrite=True)


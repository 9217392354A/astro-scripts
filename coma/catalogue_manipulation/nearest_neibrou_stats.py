# program written to find the angular contamination by using a simulated optical catalogue

# Chris Fuller, Jan 2013

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
cat = at.Table(pj(folder, 'large_coma_spx6cff.fits'),type='fits')


#key parameters
coor_names = ['ra', 'dec'] # these are the colum names that containe ra and dec
x = 1.77 #conversion between deg and mpc


# # # # # # # # # # # # # # # Function # # # # # # # # # # # # # # # # # # # # # # # # # #
#function to produce new cat with column added for the nth nearest neigboure
def nth_nearest_neighbour(t, coor_names):
	print 'nth_nearest_neighbour....'

	#add columnd for D1,D5, and D10
	t.add_empty_column('D1_1000', dtype = np.float)
	t.add_empty_column('D5_1000', dtype = np.float)
	t.add_empty_column('D10_1000', dtype = np.float)
	t.add_empty_column('SIGMA1_1000', dtype = np.float)
	t.add_empty_column('SIGMA5_1000', dtype = np.float)
	t.add_empty_column('SIGMA10_1000', dtype = np.float)
	
	###### part 2 #######
	# find nearest neighbours

	

	#loop through all members of catA
	for i in range(0, len(t)):

		#find velocity
		vel = t.velocity[i]

		#select only galaxies within 1000kms of galaxy
		tSub = t.where((t.velocity > vel-500.0) & (t.velocity < vel+500.0))

		ra = t[coor_names[0]][i]
		dec = t[coor_names[1]][i]

		#ra1 and dec1
		ra_1 = tSub[coor_names[0]]
		dec_1 = tSub[coor_names[1]]

		#ra2 and dec2
		ra_2 = np.array([ra]*len(ra_1), dtype=np.float)
		dec_2 = np.array([dec]*len(ra_1), dtype=np.float)
		
		#caculate distance to all sources from ra1 and dec1
		radius = np.sort(distance(ra_1, dec_1, ra_2, dec_2 ))

		#print radius[1]*1.77*1000.0, np.min(radius)

		#add values to table
		t['D1_1000'][i] = np.log10(radius[1] * x)
		t['D5_1000'][i] = np.log10(radius[5] * x)
		t['D10_1000'][i] = np.log10(radius[10]* x) 
		t['SIGMA1_1000'][i] = np.log10(1.0 / (np.pi*(radius[1]*x)**2.0))
		t['SIGMA5_1000'][i] = np.log10(5.0 / (np.pi*(radius[5]*x)**2.0)) 
		t['SIGMA10_1000'][i] = np.log10(10.0 / (np.pi*(radius[10]*x)**2.0)) 


	return t

#distance equation designed to do arraywise caculations
def distance(ra1, dec1, ra2, dec2):
	delta_ra = (ra1 - ra2) * np.cos(np.radians((dec1+dec2)/2.0))
	delta_dec = (dec1 - dec2)

	return np.sqrt(delta_ra**2.0 + delta_dec**2.0)

def createArray(table):
	array = np.array([table.ra,table,dec])

 # # # # # # # # # # # # # # # Main Program # # # # # # # # # # # # # # # # # # # # # # # #
new_cat  = nth_nearest_neighbour(cat, coor_names)

new_cat.write(pj(folder, 'ngp_optical_nearest-stats-v2.fits'), overwrite=True)


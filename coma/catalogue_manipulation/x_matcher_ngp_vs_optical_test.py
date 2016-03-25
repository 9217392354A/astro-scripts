# Program written to caculate what the flux would be at 250um from each galaxy if they
# were found using the ngp cat
# Chris Fuller, Mar 2013

#import modules
import numpy as np
from os.path import join as pj
import atpy as at
from copy import copy, deepcopy
import matplotlib.pyplot as plt
from random import random
from pylab import bar
from time import time
import datetime


#i/o
print 'reading in cats . . .'
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/positions/'
cat_mine = at.Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12.fits",type='fits')
catA = at.Table(pj(folder, 'ngp_psf.fits'),type='fits')

#key parameters
r_search = 300.0 / 3600.0 # search radius in degrees
coor_names_A = ['GRA2000', 'GDEC2000'] # these are the colum names that containe ra and dec


# # # # # # # # # # # # # # # Function # # # # # # # # # # # # # # # # # # # # # # # # # #
#function to produce new cat with matched columns
def x_matcher(catA, coor_names_A, catB, coor_names_B, r_search):
	print 'starting x matcher . . .'
	###### part 1 ######
	# first create new cat structure to hold old cat and new matches
	t = at.Table() #empty cat object

	colsA =  catA.columns #get columns of the first table
	colsB =  catB.columns #get columns of the second table

	cols_names_b = colsB.keys #col names for cat b

	#loop through each column adding it to the table 
	for j in range(0, len(colsA)): t.add_column(colsA.keys[j],catA[colsA.keys[j]],dtype = colsA[colsA.keys[j]])

	

	#add column for seporation and ngp250 flux
	t.add_empty_column('radial_sep_test', dtype = np.float)
	t.add_empty_column('F250_ngp_test', dtype = np.float)

	#if verbous: print t.columns #print list of new cols

	###### part 2 #######
	# find nearest neighbours

	#create arrays of ra and dec for cat B, save making new ones each loop
	ra_2_long = catB[coor_names_B[0]]
	dec_2_long = catB[coor_names_B[1]]

	#loop through all members of catA
	for i in range(0, len(t)):

		####  wiget #####
		if i%100 == 0: print np.round(i*100.0/len(t), decimals=2), '%'

		ra = t[coor_names_A[0]][i]
		dec = t[coor_names_A[1]][i]
		
		#create a search box around the source
		temp_cat = search_box(ra, dec, ra_2_long , dec_2_long, r_search, catB)

		#create new list of ra and dec
		ra_2 = temp_cat[coor_names_B[0]]
		dec_2 = temp_cat[coor_names_B[1]]

		#ra1 and dec1
		ra_1 = np.array([t[coor_names_A[0]][i]]*len(ra_2), dtype=np.float)
		dec_1 = np.array([t[coor_names_A[1]][i]]*len(ra_2), dtype=np.float)

		#if verbous: print 'caculating distance for all soruces'
		#caculate distance to all sources from ra1 and dec1
		radius = distance(ra_1, dec_1, ra_2, dec_2)

		#find where radius is closer than search radius
		w = np.where(radius < r_search)[0]	
	
		#now we need to find the nearest neigbourgh
		if len(w) == 0: #no source found that is whithin search radius
			#if verbous: print 'No source found within ', np.around(r_search*3600.0), ' moving on to next source'
			t['radial_sep'][i] = -1.0

		elif len(w) >= 1: # a single or multipul source found
			#add all the seps to a temp cat
			distances = np.array(radius[w], dtype=np.float)
			keys = np.array(w, dtype=np.int)

			d = np.sort(distances)[0]
			 
			#reshape into an array
			temp = np.append(distances.reshape(len(keys),1) , keys.reshape(len(keys),1),1)

			#sort to find the nearest source
			temp = deepcopy(temp[np.argsort(temp[:,0])])

			#now the top row second element will hold the key for the galaxy
			# from catB
			key = temp[0,1]
			rad = temp[0,0]
		
			#now add source to new cat
			t['radial_sep_test'][i] = d*3600.0
			t['F250_ngp_test'][i] = temp_cat['f250'][key]
			#loop through each column adding to new cat
			if False: 
				for j in range(0, len(cols_names_b)): t[cols_names_b[j]][i] = temp_cat[cols_names_b[j]][key]
	return deepcopy(t.where(t['radial_sep_test']>0.0)) #.write(pj(folder, 'test.fits'), overwrite=True)

#distance equation designed to do arraywise caculations
def distance(ra1, dec1, ra2, dec2):
	delta_ra = (ra1 - ra2) * np.cos(np.radians((dec1+dec2)/2.0))
	delta_dec = (dec1 - dec2)

	return np.sqrt(delta_ra**2.0 + delta_dec**2.0)

#returns a cat with galaxies inside a search box
def search_box(box_ra, box_dec, box_ra_2 , box_dec_2, r, box_cat):
	#factor of box greater than search radius
	factor = 2.0
	
	max_ra = box_ra + r*factor
	min_ra = box_ra - r*factor

	max_dec = box_dec + r*factor
	min_dec = box_dec - r*factor

	new_cat = box_cat.where((box_ra_2 < max_ra) & (box_ra_2 > (min_ra)) & (box_dec_2 < max_dec) & (box_dec_2 > min_dec))

	return new_cat

#first only select atlas sources that are above sn = 3
print 'original ngp cat: ', len(catA)
ngp = catA.where((catA.f250 / catA.e250) >= 3.0) 
print 'sources greater than sn3: ', len(ngp)

#find nearest hatlas neibours
new_cat_mine = x_matcher(cat_mine, ['GRA2000', 'GDEC2000'], catA, ['ra', 'dec'], r_search )
#now set all f250 test fluxes to 0 if radialsep is greater than 3.1"

w = np.where(new_cat_mine['radial_sep_test'] > 3.0)[0] #make selection of sources within threshold

new_cat_mine['F250_ngp_test'][w] = 0.0

new_cat_mine.write(pj('/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/', 'test.fits'), overwrite=True)


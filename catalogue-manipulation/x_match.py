# program written to find the nearest neighbor 

# Instructions 
# 1 # Change 'folder' to wherever you have your catalogs stored 
# 2 # Change the names for cata and catb
# 3 # select the threshold to search in
# 4 # change coor_names A and B to the column names in cat A and cat B for ra and dec
# 5 # make sure ra and dec are in degrees
# 6 # run! 

# Chris Fuller, Jan 2013

#import modules
import numpy as np
from os.path import join as pj
import atpy as at
from copy import copy, deepcopy
import matplotlib.pyplot as plt


#i/o
print 'reading in cats . . .'
folderA = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/"
catA = at.Table(pj(folderA, 'stellar-mass-fornax_final.fits'),type='fits')

folderB = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/SPIRE-contamination/contamination/'
catB = at.Table(pj(folderB, '0.9_full.FIT'),type='fits')

outfolder = '/Users/chrisfuller/Desktop/'
outcatname = 'test.fits'

#key parameters
r_search = 100.0 / 3600.0 # search radius in degrees
verbous = False # gives increased terminal bumf
auxdata = True # gives aux data and saves it into new cat (not woriking)

coor_names_A = ['GRA2000', 'GDEC2000'] # these are the colum names that containe ra and dec
coor_names_B = ['ALPHA_J2000', 'DELTA_J2000']

# # # # # # # # # # # # # # # Function # # # # # # # # # # # # # # # # # # # # # # # # # #
#function to produce new cat with matched columns
def x_matcher(catA, coor_names_A, catB, coor_names_B, r_search):
	###### part 1 ######
	# first create new cat structure to hold old cat and new matches
	t = at.Table() #empty cat object

	colsA =  catA.columns #get columns of the first table
	colsB =  catB.columns #get columns of the second table

	cols_names_b = colsB.keys #col names for cat b

	#loop through each column adding it to the table 
	for j in range(0, len(colsA)): t.add_column(colsA.keys[j],catA[colsA.keys[j]],dtype = colsA[colsA.keys[j]])

	#loop through each column adding an empty one to the new empty table
	for k in range(0, len(colsB)): t.add_empty_column(colsB.keys[k],dtype = colsB[colsB.keys[k]])

	#add column for seporation
	t.add_empty_column('Separation', dtype = np.float)

	#if verbous: print t.columns #print list of new cols

	###### part 2 #######
	# find nearest neighbours

	#create arrays of ra and dec for cat B, save making new ones each loop
	ra_2_long = catB[coor_names_B[0]]
	dec_2_long = catB[coor_names_B[1]]

	#loop through all members of catA
	for i in range(0, len(t)):

		####  wiget #####
		if i%50 == 0: print np.round(i*100.0/len(t), decimals=2), '%'

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
			t['Separation'][i] = -1.0

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
			t['Separation'][i] = d*3600.0

			#loop through each column adding to new cat
			if auxdata == True: 
				for j in range(0, len(cols_names_b)): t[cols_names_b[j]][i] = temp_cat[cols_names_b[j]][key]
	return deepcopy(t) #.write(pj(folder, 'test.fits'), overwrite=True)

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

 # # # # # # # # # # # # # # # Main Program # # # # # # # # # # # # # # # # # # # # # # # #

print 'starting x matcher . . .'
#call x_matcher and give it both cats and the columnnames for ra and dec of each 
# as well as your search threshold
new_cat = x_matcher(catA, coor_names_A, catB, coor_names_B, r_search)

new_cat.write(pj(outfolder, outcatname), overwrite=True)
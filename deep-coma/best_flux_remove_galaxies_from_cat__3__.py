# program removes the following from our catalogue;
	# - Bad galaxies from the list bad gal
	# - Galaxies not detected at 250
	# - Galaxies not in both pacs and spire bands

# program written by Chris Fuller, November 2013

#import numpy as np
print 'importing modules...'
from numpy import where, logical_and
from os.path import join as pj
from atpy import Table
#inputs
print 'opening tables...'
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/catalogue_creation_phase/'
catname = 'ngp+20140331__best-flux-3.5arcsec__1-2__.fits'
#load catalogues one where all are forced as point to serve as the base catalogue
cat = Table(pj(folder,catname),type='fits')

badgal = Table('/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/all-detections/bad-gal-ngp.fits',type='fits')
outname = 'ngp+20140331__best-flux-3.5arcsec__1-2-3__.fits'



############################## functions ####################################
#function that removes a galaxy from the catalogue
# -set all bands to 0 for SN, R, and F
# -3x error
# -set EXTENDEDNESS to PPPPP

def remove_single_from_cat(index,t):
	bands = ['500', '350', '250', '160', '100']
	cols = ['F', 'SN', 'R']

	#loop through each band
	for band in bands:
		#now times the error by 3 for all galaxies removed
		if t['F'+band][index] != 0.0:
			#print 'OBJECT: ', t['OBJECT'][index], 'BAND: ', band, 'FLUX: ', t['F'+band][index], 'ERROR: ', t['E'+band][index]
			t['E'+band][index] = t['E'+band][index]*3.0

		#loop through each col
		for col in cols:
			#set col to zero
			t[col+band][index] = 0.0 

		
	#set extendedness to PPPPP for all removed
	t['EXTENDEDNESS'][index] = "PPPPP"
	return t

#test if any galaxies is detected at 250 but not at any other bands
def needs_remove_no250(jdex,t):
	bands = ['500', '350', '250', '160', '100']

	if t.F250[jdex] == 0.0:
		for band in bands:
			if t['F'+band][jdex] != 0.0: return True

	return False
########################## master control program ###########################
counter1, counter2, counter3 = 0,0,0


#1 set all band to zero for all galaxies that are not detected at 250
print 'seting all bands to zero for galaxies not detected at 250...'
#w1 = where(cat['F250'] == 0.0)[0] #find where galaxies are not detected
#cat1 = remove_from_cat(w1,cat)
for j in range(0,len(cat['OBJECT'])):
	if needs_remove_no250(j,cat): 
		cat = remove_single_from_cat(j,cat)
		counter1 += 1

#2 set all band to zero for all galaxies that are listed as badgal
print 'seting all bands to zero for all galaxies that are listed as badgal...'
#loop through all bad gal
if True:
	for i in range(0,len(cat['OBJECT'])): 
		if len(where(cat['OBJECT'][i] == badgal['OBJECT'])[0]) == 1: 
			cat = remove_single_from_cat(i,cat)
			counter2 += 1


print 'Number of galaxies removed due to no 250 = ',counter1
print 'Number of galaxies removed due to badgal = ',counter2
cat.write(pj(folder, outname), overwrite=True)



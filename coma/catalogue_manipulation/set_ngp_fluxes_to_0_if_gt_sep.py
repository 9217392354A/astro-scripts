#code the reads in the main table and then sets ngp fluxes to 0 if sep gt 3.5
#Chris Fuller, March 2014

#import modules and functions
print 'importing modules and functions...'
from atpy import Table
from numpy import where, float
from os.path import join as pj

#inputs 
print 'reading in data...'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,cat_name))

#bands 
bands = ['500', '350', '250'] #f250_1a

#find where sep is greater than 3.5"
w_sp = where(cat.Separation_mine_npg > 3.5)[0]

#loop through bands
for band in bands:
	print 'starting ' + band

	#extract fluxes
	new_fluxes = cat['f' + band + '_1a']
	w_sn = where((new_fluxes / cat['e' + band + '_1a']) < 3.0)[0]

	#set all fluxes greater than 3.5 to 0
	new_fluxes[w_sp] = 0.0

	#set all fluxes to 0 that are less than s/n 3
	new_fluxes[w_sn] = 0.0

	#add new col to table with fluxes
	cat.add_column('NGPFLUX' + band, new_fluxes, unit='Jy', dtype=float)

cat.write(pj(folder, 'test.fits'), overwrite=True)
print 'program complete'
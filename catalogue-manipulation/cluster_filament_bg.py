# program to take level 3 catalouge and then asign coma_cluster, filament, and background.

from os.path import join as pj
from atpy import Table
from numpy import where
from copy import copy, deepcopy

# i/o
catname = 'ngp+20140109__best-flux-3.1arcsec__1-2-3__.fits'
outname = 'ngp+20140109__best-flux-3.1arcsec__1-2-3__FINAL.fits'

folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/catalogue_creation_phase/'


print 'opening table...'
cat = Table(pj(folder,catname),type='fits')


# vary
vmin = 4268.8 #kms
vmax = 9700.2 #kms
rmax = 6012.0 # arcsec

# find key cols
vel = cat['VELOCITY']
rad = cat['RADIUS']

#loop through all gals
for i in range(0,len(cat)): 
	#cluster
	if (vel[i] > vmin) and (vel[i] < vmax) and (rad[i] < rmax):
		cat['TYPE'][i] = 1

	#filament
	elif (vel[i] > vmin) and (vel[i] < vmax) and (rad[i] > rmax):
		cat['TYPE'][i] = 2

	#bg
	elif (vel[i] < vmin) or (vel[i] > vmax):
		cat['TYPE'][i] = 3

	else: raise 'error'

cat.write(pj(folder, outname), overwrite=True)

##### functions #####
def catsplit(val, name):
	newcat = deepcopy(cat.where(cat['TYPE'] == val))
	newcat.write(pj(folder,name), overwrite=True)

def catnot(val, name):
	newcat = deepcopy(cat.where(cat['TYPE'] != val))
	newcat.write(pj(folder,name), overwrite=True)
# split catalouge into filmament and cluster
if True:
	catsplit(1, 'coma_cluster-mgbg-130114-3.1arcsec.fits')
	catsplit(2, 'filament-mgbg-130114-3.1arcsec.fits')
	catnot(3,'cluster+filament-mybg-130114-3.1arcsec.fits')


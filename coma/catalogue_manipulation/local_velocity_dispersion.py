# programs to find veloisty dispersion of the n nearest neigbours
# Chris Fuller, March 2014

#import moduals
from atpy import Table
import numpy as np
from os.path import join as pj
from scipy import spatial
from pylab import hist, show
from multiprocessing import Pool

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
ngp = Table(pj(folder, 'coma_supercluster_cal12.fits'),type='fits')

#create a kdtree
tree_npg = spatial.cKDTree(zip(ngp['GRA2000'], ngp['GDEC2000']))

def mod(x):
	return np.sqrt(x**2)

#function to caculate veldis for nth nearst neigbours
def vel_dis_n(n):
	#query kdtree
	distance, ind  = tree_npg.query(zip(ngp['GRA2000'], ngp['GDEC2000']), k=n)

	#create empty list of velocity dispersions
	vel_dis = []

	#loop through cat
	for i in range(len(ind)): vel_dis.append(np.mean(mod(ngp['velocity'][ind[i]]-7000.0)))

	#add empty col to cat
	ngp.add_column('dispersion_'+str(n), vel_dis, unit='km/s', null='', dtype=np.float)

#pool = Pool(4)
#pool.map(vel_dis_n, [5,10,20,30,40,50,100])

for n in [5,10,20,30,40,50,100]: vel_dis_n(n)

ngp.write(pj('/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/', 'test.fits'), overwrite=True)

# program written to find the angular contamination by using a simulated optical catalogue

# Chris Fuller, Jan 2013

#import modules
import numpy as np
from os.path import join as pj
import atpy as at
from astropy.io import fits 
from copy import copy, deepcopy
import matplotlib.pyplot as plt
from random import random
from pylab import bar
from time import time
import datetime
from lmfit import minimize, Parameters,report_fit
from scipy import spatial
import multiprocessing 

#i/o
print 'reading in cats . . .'

if False:
	folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/positions/'
	cat_mine = at.Table(pj(folder, 'ngp+20140109__best-flux__.fits'),type='fits')
	catA = at.Table(pj(folder, 'ngp_psf.fits'),type='fits')

if True:
	folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/positions/'
	cat_mine = fits.open(pj(folder, 'ngp+20140109__best-flux__.fits'))[1].data
	catA = fits.open(pj(folder, 'ngp_psf.fits'))[1].data	

#key parameters
r_search = 10.0 / 3600.0 # search radius in degrees
coor_names_A = ['GRA2000', 'GDEC2000'] # these are the colum names that containe ra and dec
n= 100000 #number of interations of monty carlo
n_bins = 35 #number of radial bins

# # # # # # # # # # # # # # # Function # # # # # # # # # # # # # # # # # # # # # # # # # #
#distance equation designed to do arraywise caculations
def distance(ra1, dec1, ra2, dec2):
	delta_ra = (ra1 - ra2) * np.cos(np.radians((dec1+dec2)/2.0))
	delta_dec = (dec1 - dec2)

	return np.sqrt(delta_ra**2.0 + delta_dec**2.0)


def random_cat(cen_x, cen_y, size, length):
	rand1, rand2 =[], [] #create empty lists

	#create list of random numbers
	for i in range(0,length):
		rand1.append(random())
		rand2.append(random())

	# turn them into numpy arrays
	rand1 = np.array(rand1, dtype= np.float)
	rand2 = np.array(rand2, dtype= np.float)

	#create lists of ra and dec
	ra_rand =  cen_x + (rand1 - 0.5)*size
	dec_rand = cen_y + (rand2 - 0.5)*size

	# create new cat structure
	t = at.Table() #empty cat object

	#add cols for ra dec
	t.add_column('ra_rand',ra_rand, dtype=np.float)
	t.add_column('dec_rand', dec_rand, dtype=np.float)

	return t


def residual(params, x, y_data, y_error):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    e = params['e'].value
   
    y_model = a*x**4 + b*x**3 + c*x**2 + d*x + e

    return (y_data-y_model)/y_error

def lmfitter(x , y, y_error):
    params = Parameters()
    params.add('a', value=-1., vary=True)
    params.add('b', value=1., vary=True)
    params.add('c', value=1., vary=True)
    params.add('d', value=1., vary=True)
    params.add('e', value=0.0, vary=False)


    # remove inf values in errors
    y_error[y_error ==  np.float('inf')] = 1.0  
    out = minimize(residual, params, args=(x, y, y_error))
    report_fit(params)
    return out

def model(params, x):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    e = params['e'].value
   
    return a*x**4 + b*x**3 + c*x**2 + d*x + e


def rmodel(params, x):
    x_long = np.round(np.arange(0.,8.0, 0.001), decimals = 3)
    y_long = np.round(model(params, x_long), decimals = 3)

    #round all to the same decima  place
    sol = x_long[y_long == np.round(x, decimals=3)][0]
    print x_long[y_long == np.round(x, decimals=3)]
    return sol

def monty_carlo(iii):

	#create random cat
	rand_cat = random_cat(200.0, 29.0, 11.0, len(cat_mine))

	#find nearest neigbours of random cat to kd tree
	seps_list, seps_indexs = tree_npg.query(zip(rand_cat.ra_rand, rand_cat.dec_rand))

	#caculate distance in degrees
	seps_list = distance(ngp.ra[seps_indexs], ngp.dec[seps_indexs], rand_cat.ra_rand, rand_cat.dec_rand)*3600.0

	#create hist
	hist, binedges = np.histogram(seps_list, bins=np.arange(0.0,r_search*3600.0, r_search*3600.0/n_bins) , range=None, normed=False, weights=None, density=None)

	#add to list of hists
	#rand_hist = rand_hist + hist


	return hist

# # # # # # # # # # # # # # # Main Program # # # # # # # # # # # # # # # # # # # # # # # #

####### part 1 ###### create plot of random radial distribution ####### part 1 ###### 
start = time()

#select only sources that are above a sn of 3
ngp = catA[np.where((catA['f250'] / catA['e250']) >= 3.0)[0]]

#create kd tree
tree_npg = spatial.cKDTree(zip(ngp['ra'], ngp['dec']))

#create array for rand_hist
rand_hist = np.zeros(n_bins-1)

"""
print 'simulation random distribution'
for i in range(0,n):
	hist, bin = monty_carlo()
	rand_hist = rand_hist + hist

"""

pool = multiprocessing.Pool(4)
rand_hist =  pool.map(monty_carlo, range(n))

master = np.zeros(n_bins-1)
#add rand_hist together
for count in range(len(rand_hist)):
	master += rand_hist[count]

rand_hist = master

print 'completion time: ', str(datetime.timedelta(seconds=(time() - start)))

	
	#now do mean of sample
	#print 'single loop ', str(datetime.timedelta(seconds=loop_time))
	#new = time()
	#loop_time_current = new - old
	#old = new
	#print 'estimated completion time remaining: ', str(datetime.timedelta(seconds=loop_time_current*(n-i)))
	#loop_time = time() - start
	

#find average random hist
rand_hist = rand_hist.astype(np.float)/np.float(n)




####### part 2 ###### perform x_match on my optical cat ####### part 2 ###### 
print 'caculating my distribution'
#new_cat_mine = x_matcher(cat_mine, ['GRA2000', 'GDEC2000'], catA.where((catA.f250 / catA.e250) >= 3.0), ['ra', 'dec'], r_search )

seps_list_mine, seps_indexs = tree_npg.query(zip(cat_mine['GRA2000'], cat_mine['GDEC2000']))
seps_list = distance(ngp['ra'][seps_indexs], ngp['dec'][seps_indexs], cat_mine['GRA2000'], cat_mine['GDEC2000'])
seps_list_mine = seps_list_mine*3600.0
#find sources with matches
mine_hist, binedges = np.histogram(seps_list_mine, bins=np.arange(0.0,r_search*3600.0, r_search*3600.0/n_bins) , range=None, normed=False, weights=None, density=None)






####### part 3 ###### caculate the % contamination per bin ####### part 3 ###### 
conta_hist_fraction = rand_hist.astype(np.float)*100./ mine_hist.astype(np.float)


#fit polynomial to contamination plot
if False: 
	fit = np.polyfit(binedges[:-1]+r_search*3600.0/n_bins,conta_hist_fraction,2.)
else: 
	xx = binedges[:-1]+(r_search*3600.0/n_bins)*0.5
	yy = conta_hist_fraction
	try:
		yy_error = yy *  (1.0 / np.sqrt(mine_hist.astype(np.float)))
	except:
		yy_error = np.array([1.0]*len(yy), dtype= np.float)
		print 'error'



#fit = lmfitter(xx , yy,  np.array([1.0]*len(yy), dtype= np.float))
fit = lmfitter(xx , yy,  yy_error)


X_con = xx


#save all relevent details
a = rand_hist.reshape(len(conta_hist_fraction),1)
b = mine_hist.reshape(len(conta_hist_fraction),1)
c = conta_hist_fraction.reshape(len(conta_hist_fraction),1)
d = binedges[:-1].reshape(len(conta_hist_fraction),1) + r_search*3600.0/n_bins

a = np.append(a, b, 1)
a = np.append(a, c, 1)
a = np.append(a, d, 1)

header = np.array(['random_histogram', 'my_optical_hist', 'contamination_fraction', 'angular_distance_arcsec'], dtype=np.str).reshape(1,4)

a = np.append(header,a,0)


####### part 4 ####### plot results ###### part 4 #########

fig = plt.figure(figsize = (4.5,4.5),facecolor='w',edgecolor='w')
sub1 = fig.add_subplot(2,1,1)
sub2 = fig.add_subplot(2,1,2)

sub1.bar(xx, rand_hist, align='center', width=r_search*3600.0/n_bins )#, fc ='None',hatch='/')
sub1.bar(xx, mine_hist, align='center', width=r_search*3600.0/n_bins, fc ='None',hatch='/', yerr=np.sqrt(mine_hist.astype(np.float)))
sub1.tick_params(axis='x', labelbottom='off')
sub1.set_ylabel('Number of galaxies')

sub2.bar(xx, yy, align='center', width=r_search*3600.0/n_bins, fc ='b', yerr=yy_error)
sub2.plot(xx, model(fit.params,xx), 'r-')
sub2.set_xlabel('Angular Separation, (arcsec)')
sub2.set_ylabel('Contamination (%)')
sub2.set_ylim(0, 115)

sub1.set_xlim(0,r_search*3600.0)
sub2.set_xlim(0,r_search*3600.0)

plt.subplots_adjust(left=0.15, bottom=0.13, right=0.97, top=0.97,hspace=0.0)

try:
	contam = rmodel(fit.params, 5.0)

except: 
	contam = 999.
fig.savefig('/Users/chrisfuller/Dropbox/phd/thesis/main/chapter-coma/angular_contamination-v2.pdf')


print '5% contamination reached at ', contam
plt.show()

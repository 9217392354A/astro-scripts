# Program to make figures for velocity dispersions and fit
# Chris Fuller, April 2014

import numpy as np
import atpy as at
import matplotlib.pyplot as plt
from os.path import join as pj
from lmfit import minimize, Parameters, report_fit

#Inputs
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'
c1 = 'big_ngp.fits' #input name

#read in cat
cat = at.Table(pj(folder2,c1))

#distance equation designed to do arraywise caculations
def distance(ra, dec, ra1, dec1):

	ra2 = np.array([ra]*len(ra1), dtype= np.float) 
	dec2 = np.array([dec]*len(ra1), dtype= np.float)


	delta_ra = (ra1 - ra2) * np.cos(np.radians((dec1+dec2)/2.0))
	delta_dec = (dec1 - dec2)

	return np.sqrt(delta_ra**2.0 + delta_dec**2.0)

def residual(params, x, y_data, y_error):
	y_model = model(params, x)
	return (y_data-y_model)/y_error

def lmfitter(x, y, y_error):
    params = Parameters()
    params.add('sigma', value=905.2, vary=False)
    params.add('mean', value=6984.5, vary=False)
    params.add('hmax', value=1000., vary=True)
    params.add('hmin', value=1., vary=True)


    # remove inf values in errors
    out = minimize(residual, params, args=(x, y, y_error))
    report_fit(out.params)
    return out

def model(params, x):
    sigma = params['sigma'].value
    mean = params['mean'].value
    hmax = params['hmax'].value
    hmin = params['hmin'].value
	
    return hmax*np.exp(-(x-mean)**2/(2*sigma**2)) +hmin
    

##### main program #####

#caculate radius
rad = distance(194.9531, 27.9807, cat.ra, cat.dec) 
cat.add_column('RADIUS', rad, unit='Mpc', dtype=np.float)

#select cat where radius is less than 3.1Mpc
#cat = cat.where(cat.RADIUS <= 1.67)
cat = cat.where((cat.velocity > 4268.8 ) & (cat.velocity < 9700.2) & (cat.RADIUS <= 1.67))

#extract velocity
vel = cat.velocity
hist, binedges = np.histogram(vel,bins=25)


#fit velocity dispersion with gaussian
fit = lmfitter(binedges[:-1] , hist , np.sqrt(hist+1.))


#caculate gaussian 
xx = np.linspace(0, 16000.0, 10000)
yy = model(fit.params, xx)


#####plot result######

#create figure and subfigure
fig = plt.figure(figsize = (8.0,8.0),facecolor='w',edgecolor='w')
sub = fig.add_subplot(1,1,1)

#plot historgram
wid = binedges[1]-binedges[0]
sub.bar((binedges[:-1]+wid*0.5)/1000.0, hist, width=wid/1000.0, fc ='None',hatch='/')
sub.plot((xx+wid*0.5)/1000.0, yy, '-k')
sub.axvline(x=4268.8/1000.0, ls='--', c='k')
sub.axvline(x=9700.2/1000.0, ls='--', c='k')
sub.set_xlim(0,14)


#label
sub.set_ylabel('N')
sub.set_xlabel('Velocity (10$^{3}$km s$^{-1}$)')
plt.subplots_adjust(left=0.08, bottom=0.08, right=0.97, top=0.98)



#save and show 
#fig.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/coma','velocity_dispersion.pdf'))
#fig.savefig(pj('/Users/chrisfuller/Dropbox/phd/thesis/main/chapter-coma/','velocity_dispersion.pdf'))

plt.show()




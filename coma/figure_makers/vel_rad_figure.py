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
  

##### main program #####

#caculate radius
rad = distance(194.9531, 27.9807, cat.ra, cat.dec) * (1.66 / 3.01)  
cat.add_column('RADIUS', rad, dtype=np.float)  

#####plot result######

#create figure and subfigure
fig = plt.figure(figsize = (8.0,4.5),facecolor='w',edgecolor='w')
sub = fig.add_subplot(1,1,1)

sigma = 905.2/1000.0

rect = plt.Rectangle((0,4268.8/1000.0),1.0, sigma*6, color="red",  fc='None', hatch = '/', alpha=0.8)
plt.gca().add_patch(rect)

rect = plt.Rectangle((1,4268.8/1000.0),11.5-1., sigma*6, color="blue",  fc='None', hatch = 'X', alpha=0.8)
plt.gca().add_patch(rect)

sub.scatter(cat.RADIUS, cat.velocity/1000.0, s=10, c='k', marker='o')


#sub.plot((xx)/1000.0, yy, '-k')
sub.axhline(y=4268.8/1000.0, ls='--', c='k')
sub.axhline(y=9700.2/1000.0, ls='--', c='k')
sub.axvline(x=1.0, ls='--', c='k')
sub.set_xlim(0,6.5)
sub.set_ylim(2.5,11.5)

#label
sub.set_xlabel('Projected Cluster Radius (R/R$_{Virial}$)')
sub.set_ylabel('Velocity (10$^{3}$km s$^{-1}$)')
plt.subplots_adjust(left=0.1, bottom=0.12, right=0.97, top=0.98)



#save and show 
#fig.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/coma','velocity_radius.pdf'))
fig.savefig(pj('/Users/chrisfuller/Dropbox/phd/thesis/main/chapter-coma/','velocity_radius.pdf'))
plt.show()




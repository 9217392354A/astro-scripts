# programs to plot sigma 1 sigma 5 sigma 10 plots
# Chris Fuller, March 2014

#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#Inputs
folder = "/Users/chrisfuller/Documents/phd/herchel/coma/final_outputs/" # input/output folder
folder2 = '/Users/chrisfuller/Documents/phd/herchel/coma/aux_data/'

c1 = 'coma_supercluster_cal12.fits' #input name
c1 = 'coma_supercluster_cal12_pacscorrected.fits'
c3 = 'fornax_at_100mpc.fits'
c2 = 'virgo_at_100mpc.fits'

N = 3 #number of sigma and radial density bins between logmin and logmax
colours = ['Green', 'Orange', 'Red']

cats = [Table(pj(folder,c1)), Table(pj(folder2,c2)), Table(pj(folder2,c3))] 
cat_names = ['Coma', 'Virgo', 'Fornax']
morphsys = ['zoo', 'goldmine', 'goldmine']

#quanties on x axcis
#cols  = ['SIGMA1', 'SIGMA5', 'SIGMA10'] #, 'RADIUS_VIR']
cols  = ['SIGMA1_1000', 'SIGMA5_1000', 'SIGMA10_1000']
labels  = ['$\Sigma$ 1', '$\Sigma$ 5', '$\Sigma$ 10']

#xlabs = ['$Log_{10}(\Sigma 10) ($Mpc$^{-2})$', 'Projected Cluster Radius (Mpc)']
#xlims = [[-1.2,3.3], [0.0, 7.1]]


#create figure and subplots
fig, subs = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=False, squeeze=True, figsize = (6.*(3./4),4.*(3./4)), facecolor='w',edgecolor='w')

#loop throgh cats
for i in range(1):
	cat = cats[i]

	#loop through cols
	for j in range(len(cols)):
		para = cols[j]
		#create density bins
		density_bins = np.linspace(min(cat[cols[2]]), max(cat[cols[2]]), N+1)

		print density_bins
		sub = subs[0,j]

		#loop through density bins
		for k in range(len(density_bins)-1):
			lower = density_bins[k]
			upper = density_bins[k+1]

			#select cat
			temp = cat.where((cat[para] >= lower) & (cat[para] <= upper))



			sub.scatter(temp.GRA2000, temp.GDEC2000, s=15, marker='o', c=colours[k], edgecolor='None', alpha=0.6)
		sub.text(0.6, 0.98, labels[j] , transform=sub.transAxes, fontsize=16, verticalalignment='top')
	
		sub.autoscale_view(tight=True, scalex=True, scaley=True)
		if j != 0: 
			
			sub.tick_params(axis='y', labelleft='off')
		sub.tick_params(axis='x', labelbottom='off')
		#sub.invert_xaxis()

		#sub.set_xlabel('RA (J2000)')
		if j == 0: sub.set_ylabel('Dec (J2000)')

	#loop through cols
	for j in range(len(cols)):
		para = cols[j]
		#create density bins
		#density_bins = np.linspace(min(cat[para]), max(cat[para]), N+1)

		print density_bins

		sub = subs[1,j]

		#loop through density bins
		for k in range(len(density_bins)-1):
			lower = density_bins[k]
			upper = density_bins[k+1]

			#select cat
			temp = cat.where((cat[para] >= lower) & (cat[para] <= upper))



			sub.scatter(temp.GRA2000, temp.VELOCITY_1/1000.0, s=20, marker='o', c=colours[k], edgecolor='None', alpha=0.6)
			
		sub.autoscale_view(tight=True, scalex=True, scaley=True)
		if j != 0: 
			sub.tick_params(axis='y', labelleft='off')
		
		sub.invert_xaxis()

		sub.set_xlabel('RA (J2000)')
		if j == 0: sub.set_ylabel('Velocity, 10$^{3}$ km s$^{-1}$')
##### plotting options
plt.subplots_adjust(left=0.14, bottom=0.16, right=0.96, top=0.99, wspace=0.0, hspace=0.0)


fig.savefig('/Users/chrisfuller/Desktop/ra_vs_dec_sigmas.pdf')
plt.show()





		

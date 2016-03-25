# radius vs sigma stats
# Chris Fuller, March 2014

#import modules
import numpy as np
from atpy import Table
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'

c1 = 'coma_supercluster_cal12.fits' #input name
c3 = 'fornax_at_100mpc.fits'
c2 = 'virgo_at_100mpc.fits'

cats = [Table(pj(folder,c1)), Table(pj(folder2,c2)), Table(pj(folder2,c3))] 

viral = [3.01, 1.68, 0.7,]
#sigma stats
cols = ['SIGMA1', 'SIGMA5', 'SIGMA10']
colours = ['k', 'b', 'r' ]
cluster = ['Coma', 'Virgo', 'Fornax']
names = ['\Sigma1','\Sigma5','\Sigma10'] 
#create figure and subplots
fig, subs = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False, squeeze=False, figsize = (6.5,10.5),facecolor='w',edgecolor='w')




for j in range(len(cats)):
	#read in catalogue
	cat = cats[j]

	# loop through cols
	for i in range(len(cols)):

		col = cols[i] #asign col

		#asign x and y
		x= cat['RADIUS']/viral[j]
		y = cat[col]

		#plot x vs y
		sub = subs[i,0]
		sub.scatter(x,y, marker = 'o', color = colours[j], s=30, alpha=1., label = cluster[j])

		sub.xaxis.set_major_locator(MaxNLocator(4))
		sub.yaxis.set_major_locator(MaxNLocator(4))

		sub.set_ylabel('$\log_{10}(' + names[i] + ')$ $($Mpc$^{-2}$$)$')


sub.set_xlim(0, 6.8)
#sub.set_ylim(ymin = 0.0, ymax=4.1)
subs[2,0].set_xlabel('$($R$/$R$_{Virial})$')
subs[0,0].legend()
plt.subplots_adjust(left=0.13, bottom=0.11, right=0.98, top=0.99, wspace=0.0, hspace=0.0)

fig.savefig('/Users/chrisfuller/Desktop/radius_vs_sigma_stats_sub.pdf')

plt.show()
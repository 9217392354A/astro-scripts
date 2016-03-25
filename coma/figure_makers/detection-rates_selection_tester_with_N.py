# Program to make a figure of detection rate vs any parameter of choice
# parameters will be given as a list
# Chris Fuller, March 2014

#import
print 'importing modules...'
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
#remove numpy runtime warings
#np.seterr(invalid='ignore')

#Inputs
print 'reading in cats'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,cat_name))
cluster = cat.where(cat['RADIUS_VIR'] <=  1.0)
filament = cat.where(cat['RADIUS_VIR'] >  1.0)

#User variables
params = ['SMASS','Mr'] #, 'g', 'r', 'i', 'z'] #parameters to check against detection rates
detCols = ['F250']#,'NGPFLUX250']#, 'RERUNF250'] #detection rate columns
N = 5 #number of bins
morphological_types = ['early', 'late', 'inter'] #morphological types
morphological_colours = ['r', 'b', 'g'] #colous for each line
xlabs = ['$log_{10}(M_{star}/$M$_{\odot})$','$(M_{r})$']
#texts = ['(This paper)', 'NGP H-ATLAS Catalogue', 'RERUN FLUX250'] 

### Functions ###
def plot_sub(sub, dcol, para, bins, subLower):

	lines = []
	labels = []

	#loop through morphological types
	for ii in range(len(morphological_types)):
		morph = morphological_types[ii]
		col = morphological_colours[ii]

		cluster = cat.where((cat['RADIUS_VIR'] <=  1.0) & (cat[morph] == 1))
		filament = cat.where((cat['RADIUS_VIR'] >  1.0) & (cat[morph] == 1))

		""" Caculate Cluster Fraction Detected """


		clusterDet = cluster[para][np.where(np.nan_to_num(cluster[dcol]) > 0.0)[0]]
		clusterTot = cluster[para]

		clusterDet_hist, _ = np.histogram(clusterDet, bins=bins)
		clusterTot_hist, _ = np.histogram(clusterTot, bins=bins)

		clusterFraDet = clusterDet_hist.astype(dtype=np.float) / clusterTot_hist



		""" Caculate Filiament Fracion Detected """

		filamentDet = filament[para][np.where(np.nan_to_num(filament[dcol]) > 0.0)[0]]
		filamentTot = filament[para]

		filamentDet_hist, _ = np.histogram(filamentDet, bins=bins)
		filamentTot_hist, _ = np.histogram(filamentTot, bins=bins)

		filamentFraDet = filamentDet_hist.astype(dtype=np.float) / filamentTot_hist



		""" Plot each """
		delta = bins[1] - bins[0]

		l1 = sub.errorbar(bins[:-1]+delta*0.5, np.nan_to_num(clusterFraDet), yerr=er(clusterDet_hist,clusterTot_hist),  color = col, label='Cluster ' + morph, ls = '-')
		l2 = sub.errorbar(bins[:-1]+delta*0.5, np.nan_to_num(filamentFraDet), yerr=er(filamentDet_hist,filamentTot_hist),  color = col, label='Filament ' + morph, ls='--')



		subLower.errorbar(bins[:-1]+delta*0.5, clusterTot_hist, yerr=np.sqrt(clusterTot_hist),  color = col, ls = '-')
		subLower.errorbar(bins[:-1]+delta*0.5, filamentTot_hist, yerr=np.sqrt(filamentTot_hist),  color = col, ls = '--')

		#sub.plot(bins[:-1]+delta*0.5, clusterFraDet)
		sub.xaxis.set_major_locator(	MaxNLocator(6)	)
		sub.yaxis.set_major_locator(	MaxNLocator(6)	)

		subLower.xaxis.set_major_locator(	MaxNLocator(6)	)
		subLower.yaxis.set_major_locator(	MaxNLocator(6)	)

		lines.append(l1)
		lines.append(l2) 
		labels.append('Cluster ' + morph) 
		labels.append('Filament ' + morph)
	return lines, labels



def er(a,b):
    return np.sqrt(a) * 1. / b


### Main Program ###

#create figure and subplots
fig, subs = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, squeeze=False, figsize = (4.,4.), facecolor='w',edgecolor='w')


#loop through parameters 
for j in range(len(params)):
	parameter = params[j] #select parameter of intrest
	_, bins = np.histogram(cat[parameter], bins=N)

	#loop through detection rate columns
	for i in range(len(detCols)):
		detCol = detCols[i] #select detection column

		#select subplot from subs
		sub = subs[i,j]
		subLower = subs[i+1,j]
		
		sub.set_ylim(0,0.95)
		subLower.set_ylim(0,380)

		#plot detection rates vs parameter
		lines, labels = plot_sub(sub, detCol, parameter, bins,subLower)

		if i == len(params) - 1: sub.set_xlabel(xlabs[j])

		

		sub.tick_params(axis='x', labelbottom='off')

		#if j == 0: sub.text(0.02, 0.9, texts[i] , transform=sub.transAxes, fontsize=14, verticalalignment='top')

		if params[j] == 'SMASS': 
			sub.axvline(x=9.1, ls = '--', c='k')
			sub.axvspan(5.0, 9.1, facecolor='k', hatch='X', alpha=0.1)
			sub.text(9.1, 0.90, '$10^{-3}$' , fontsize=12, verticalalignment='top', family='sans-serif')
			sub.set_xlim(min(cat.SMASS), max(cat.SMASS))

			subLower.axvline(x=9.1, ls = '--', c='k')
			subLower.axvspan(5.0, 9.1, facecolor='k', hatch='X', alpha=0.1)
			subLower.text(9.1, 0.90, '$10^{-3}$' , fontsize=12, verticalalignment='top', family='sans-serif')
			subLower.set_xlim(min(cat.SMASS), max(cat.SMASS))

		if params[j] == 'Mr': 
			
			sub.set_xlim(min(cat.Mr), max(cat.Mr))
			#sub.invert_xaxis()

			subLower.set_xlim(min(cat.Mr), max(cat.Mr))
			#subLower.invert_xaxis()



		#subLower.set_ylim(0,1)
		#sub.set_ylim(0,0.95)


#fig.legend( lines, labels, loc = 8, ncol=6, fontsize=8 )


subs[0,0].set_ylabel('Fraction detected')
subs[1,0].set_ylabel('N')

subs[0,1].tick_params(axis='both', labelleft='off', labelbottom='off')
subs[1,1].tick_params(axis='both', labelleft='off')


subs[1,0].set_xlabel(xlabs[0])
subs[1,1].set_xlabel(xlabs[1])
plt.subplots_adjust(left=0.15, bottom=0.11, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
#fig.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/coma/','detection_rates.pdf'))
plt.show()

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
detCols = ['F250','NGPFLUX250']#, 'RERUNF250'] #detection rate columns
N = 4 #number of bins
morphological_types = ['early', 'late', 'inter'] #morphological types
morphological_colours = ['r', 'b', 'g'] #colous for each line
xlabs = ['$log_{10}(M_{star}/$M$_{\odot})$','$(M_{r})$']
texts = ['(This paper)', 'NGP H-ATLAS Catalogue', 'RERUN FLUX250'] 

### Functions ###
def plot_sub(sub, dcol, para, bins):

	lines = []
	labels = []

	#loop through morphological types
	for ii in range(len(morphological_types)):
		morph = morphological_types[ii]
		col = morphological_colours[ii]

		cluster = cat.where((cat['RADIUS_VIR'] <=  1.0) & (cat[morph] == 1))
		filament = cat.where((cat['RADIUS_VIR'] >  1.0) & (cat[morph] == 1))

		""" Caculate Cluster Fraction Detected """
		print morph, col, dcol

		clusterDet = cluster[para][np.where(np.nan_to_num(cluster[dcol]) > 0.0)[0]]
		clusterTot = cluster[para]

		clusterDet_hist, _ = np.histogram(clusterDet, bins=bins)
		clusterTot_hist, _ = np.histogram(clusterTot, bins=bins)

		clusterFraDet = clusterDet_hist.astype(dtype=np.float) / clusterTot_hist

		print 'Cluster Det Hist', clusterDet_hist
		print 'Cluster Tot Hist', clusterTot_hist

		print 'Cluster per Hist', clusterFraDet

		""" Caculate Filiament Fracion Detected """

		filamentDet = filament[para][np.where(np.nan_to_num(filament[dcol]) > 0.0)[0]]
		filamentTot = filament[para]

		filamentDet_hist, _ = np.histogram(filamentDet, bins=bins)
		filamentTot_hist, _ = np.histogram(filamentTot, bins=bins)

		filamentFraDet = filamentDet_hist.astype(dtype=np.float) / filamentTot_hist

		print 'filament Det Hist', filamentDet_hist
		print 'filament Tot Hist', filamentTot_hist

		print 'filament per Hist', filamentFraDet

		""" Plot each """
		delta = bins[1] - bins[0]

		l1 = sub.errorbar(bins[:-1]+delta*0.5, np.nan_to_num(clusterFraDet*100.0), yerr=er(clusterDet_hist,clusterTot_hist),  color = col, label='Cluster ' + morph, ls = '-')
		l2 = sub.errorbar(bins[:-1]+delta*0.5, np.nan_to_num(filamentFraDet*100.0), yerr=er(filamentDet_hist,filamentTot_hist),  color = col, label='Filament ' + morph, ls='--')

		#sub.plot(bins[:-1]+delta*0.5, clusterFraDet)
		sub.xaxis.set_major_locator(	MaxNLocator(6)	)
		sub.yaxis.set_major_locator(	MaxNLocator(6)	)

		lines.append(l1)
		lines.append(l2) 
		labels.append('Cluster ' + morph) 
		labels.append('Filament ' + morph)
	return lines, labels



def er(a,b):
    return np.sqrt(a) *100.0 / b


### Main Program ###

#create figure and subplots
fig, subs = plt.subplots(nrows=len(params), ncols=len(detCols), sharex=False, sharey=True, squeeze=True, figsize = (4.,4.), facecolor='w',edgecolor='w')


#loop through parameters 
for j in range(len(params)):
	parameter = params[j] #select parameter of intrest
	_, bins = np.histogram(cat[parameter], bins=N)

	#loop through detection rate columns
	for i in range(len(detCols)):
		detCol = detCols[i] #select detection column

		#select subplot from subs
		try: sub = subs[i,j]
		except: sub = subs

		#plot detection rates vs parameter
		lines, labels = plot_sub(sub, detCol, parameter, bins)

		if i == len(params) - 1: sub.set_xlabel(xlabs[j])

		if j == 0: sub.set_ylabel('Detected Rate (%)')

		if i != len(detCols) - 1: sub.tick_params(axis='x', labelbottom='off')

		#if j == 0: sub.text(0.02, 0.9, texts[i] , transform=sub.transAxes, fontsize=14, verticalalignment='top')

		if params[j] == 'SMASS': 
			sub.axvline(x=9.1, ls = '--', c='k')
			sub.axvspan(5.0, 9.1, facecolor='k', hatch='X', alpha=0.1)
			sub.text(9.1, 90.0, '$10^{-3}$' , fontsize=12, verticalalignment='top', family='sans-serif')
			sub.set_xlim(min(cat.SMASS), max(cat.SMASS))

		if params[j] == 'Mr': 
			
			sub.set_xlim(min(cat.Mr), max(cat.Mr))
			sub.invert_xaxis()
		sub.set_ylim(0,120)

#fig.legend( lines, labels, loc = 8, ncol=6, fontsize=10 )
fig.savefig(pj('/Users/chrisfuller/Desktop/','detection_rates.pdf'))


plt.subplots_adjust(left=0.08, bottom=0.18, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
plt.show()

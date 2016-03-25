#program to plot FIR properties 
# Chris Fuller, April 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MaxNLocator


#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(pj(folder,fname))

selection = cat.where(cat.DMASS_TYPE == 2)
#selection = cat.where((cat.DMASS_TYPE == 2) & (cat.F160 > (43.8*10**-3)*3.0) & (cat.F100 > (24.0*10**-3)*3.0))


types = ['early', 'inter', 'late'] 
labels = ['early', 'uncertain', 'late'] 
colours = ['r', 'g', 'b']
ylims = [10, 39, 78]
#quanties on the yaxis
xlabs = ['$Log_{10} (M_{dust} / $M$_{\odot}$)', 'Dust Temp. (K)','$Log_{10} (M_{dust} / M_{stars}$)']

	
#create figure and subplots
fig, subs = plt.subplots(nrows=3, ncols=3, sharex=False, sharey=False, squeeze=False, figsize = (8., 4.), facecolor='w',edgecolor='w')



#dust mass histogram
for j in range(3):
	t = types[j]
	print ''
	print 'Type: ', t

	subL = subs[j,0]
	subM = subs[j,2]
	subR = subs[j,1]
	

	subCat = selection.where(selection[t] == 1)

	#left --- Dust Mass ---
	hist_mass, bin_edges_mass = np.histogram(subCat.DMASS, bins=5)
	massMean, massErr = np.mean(subCat.DMASS), np.std(subCat.DMASS)/np.sqrt(len(subCat.DMASS)*1.0)

	print 'Dust Mass Mean: ', np.round(massMean, decimals=2), '$\\pm$', np.round(massErr, decimals=2)
	print 'Range: ', np.round(np.min(subCat.DMASS), decimals=2), ' - ', np.round(np.max(subCat.DMASS), decimals=2)

	w_mass = bin_edges_mass[1] - bin_edges_mass[0]

	
	subL.axvspan(massMean - massErr, massMean + massErr, color = colours[j], alpha=0.1)
	subL.axvline(x=massMean, c='k', ls='--')
	subL.bar(bin_edges_mass[:-1], hist_mass, w_mass, edgecolor=colours[j], fc ='None',hatch='\\')


	

	#mid --- Dust to STARs ---
	hist_d2s, bin_edges_d2s = np.histogram(subCat.DUST_STARS, bins=5)
	d2sMean, d2sErr = np.mean(subCat.DUST_STARS), np.std(subCat.DUST_STARS)/len(subCat)

	w_d2s = bin_edges_d2s[1] - bin_edges_d2s[0]

	subR.axvspan(d2sMean - d2sErr, d2sMean + d2sErr, color = colours[j], alpha=0.1)
	subR.axvline(x=d2sMean, c='k', ls='--')
	subR.bar(bin_edges_d2s[:-1], hist_d2s, w_d2s, edgecolor=colours[j], fc ='None',hatch='\\')

	print 'Dust-to-stars Mean: ', np.round(d2sMean, decimals=2), '$\\pm$', np.round(d2sErr, decimals=2)
	print 'Range: ', np.round(np.min(subCat.DUST_STARS), decimals=2), ' - ', np.round(np.max(subCat.DUST_STARS), decimals=2)
	#right --- Dust Temp ---
	hist_temp, bin_edges_temp = np.histogram(subCat.DTEMP, bins=5)
	tempMean, tempErr = np.mean(subCat.DTEMP), np.std(subCat.DTEMP)/len(subCat)

	w_temp = bin_edges_temp[1] - bin_edges_temp[0]

	subM.axvspan(tempMean - tempErr, tempMean + tempErr, color = colours[j], alpha=0.1)
	subM.axvline(x=tempMean, c='k', ls='--')
	subM.bar(bin_edges_temp[:-1], hist_temp, w_temp, edgecolor=colours[j], fc ='None',hatch='\\')


	print 'Dust Temp Mean: ', np.round(tempMean, decimals=2), '$\\pm$', np.round(tempErr, decimals=2)
	print 'Range: ', np.round(np.min(subCat.DTEMP), decimals=2), ' - ', np.round(np.max(subCat.DTEMP), decimals=2)
	subR.tick_params(axis='both', labelleft='off')
	subM.tick_params(axis='both', labelleft='off')

	font = {'weight' : 'normal', 'size'   : 12}
	subL.text(0.03, 0.98, labels[j] , transform=subL.transAxes, fontdict=font, verticalalignment='top')

	subR.set_ylim(0, ylims[j])
	subM.set_ylim(0, ylims[j])
	subL.set_ylim(0, ylims[j])

	if j == 2: 
		subL.set_xlabel(xlabs[0])
		subM.set_xlabel(xlabs[1])
		subR.set_xlabel(xlabs[2])

	subL.set_ylabel('N')

	if j != 2: 
		subL.tick_params(axis='both', labelbottom='off')
		subM.tick_params(axis='both', labelbottom='off')
		subR.tick_params(axis='both', labelbottom='off')

	subL.xaxis.set_major_locator(MaxNLocator(4))
	subL.yaxis.set_major_locator(MaxNLocator(4))

	subR.xaxis.set_major_locator(MaxNLocator(4))
	subR.yaxis.set_major_locator(MaxNLocator(4))

	subM.xaxis.set_major_locator(MaxNLocator(4))
	subM.yaxis.set_major_locator(MaxNLocator(4))

	subL.yaxis.set_label_coords(-0.16, 0.5)
	subM.yaxis.set_label_coords(-0.16, 0.5)
	subR.yaxis.set_label_coords(-0.16, 0.5)

	subL.set_xlim(6.2, 9.6)
	subR.set_xlim(-5., -0.6)
	subM.set_xlim(8, 29)




plt.subplots_adjust(left=0.08, bottom=0.13, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/sed_results.pdf')
plt.show()




#program to create functions


#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MaxNLocator
import scipy

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))
cat = cat.where(cat.SMASS > 6.8)


#cats
cats = [cat.where(cat.RADIUS_VIR <= 1.0), cat.where(cat.RADIUS_VIR > 1.0)]
cats_names = ['Cluster&', 'Filament&']

N = 12


#columns to plot
plots = [	['SMASS', 'HI_ALL2', 'DMASS']] 
lines = [9.1, 9.1, 6.4]
limits = [[8,12], [8, 11], [6,9]]

plots = np.array(plots)

#create dict for lables
labels = {	'DMASS' 		: '$\log_{10} (M_{dust} / $M$_{\odot}$)', 
			'SMASS' 		: '$\log_{10} (M_{stars} / $M$_{\odot}$)',
			'DUST_STARS' 	: '$\log_{10} (M_{dust} / M_{star}$)',
			'SRF' 			: '$\log_{10}$(SFR) (M$_{\odot}$ yr$^{-1}$)',
			'sSFR'			: '$\log_{10}$(sSFR) (yr$^{-1}$)',
			'METAL'			: '$\log_{10}$(O/H) + 12',
			'colour'		: 'g-r',
			'Mr'			: 'M$_{r}$',
			'r'				: 'm$_{r}$',
			'HI_ALL2'		: '$\log_{10} (M_{gas} / $M$_{\odot}$)',
			'my_morph'		: '$\phi$'}

method = {	'DMASS' 		: 'FIR', 
			'SMASS' 		: 'non-FIR',
			'DTEMP'			: 'FIR',
			'DUST_STARS' 	: 'FIR',
			'SRF' 			: 'SF',
			'sSFR'			: 'SF',
			'METAL'			: 'SF',
			'colour'		: 'non-FIR',
			'my_morph'		: 'non-FIR',
			'Mr'			: 'non-FIR',
			'r'				: 'non-FIR',
			'HI_ALL2'		: 'gas'}

select = {	'DMASS' 		:  8.2,
			'DTEMP'			:  .0, 
			'SMASS' 		:  7.6,
			'DUST_STARS' 	: -7.0,
			'SRF' 			: -100.,
			'sSFR'			: -100.,
			'METAL'			: -100.,
			'colour'		: -100.,	
			'my_morph'		: -100.,
			'HI_ALL2'		: 9.2	}

#morphological types
types = ['late', 'inter', 'early']
colours = ['b', 'g', 'r']
types_lab = ['late', 'intermediate', 'early']


colours_sam = [ 'r', 'b']
#set i and j range
i_range = plots.shape[1]
j_range = plots.shape[0]

#create figure and subplots
fig, subs = plt.subplots(nrows=j_range, ncols=i_range, 
						sharex=False, sharey=True, squeeze=False, 
						figsize = (8.5, 3.5), 
						facecolor='w',edgecolor='w')


# # # # Functions # # # #


#function to plot function
def functionPlotter(x, bins, colour, subplot):
	#create histogram of data
	hist, bins = np.histogram(x, bins=bins)

	#caculate errors
	hist_errors = np.sqrt(hist) 

	#cacuate function
	y  = hist
	y_err = hist_errors

	#width of bin
	w = bins[1] - bins[0]

	#find x
	x = bins[:-1] + w*0.5



	y_log = np.log10(y)

	y_errorUp = abs(np.log10(y) - np.log10(y+y_err))
	y_errorDown = abs(np.log10(y) - np.log10(y-y_err))

	y_errors = np.append(y_errorDown.reshape(1,len(y_err)), y_errorUp.reshape(1,len(y_err)), axis=0)

	#plot results
	subplot.errorbar(x, y_log, yerr=y_errors,  c=colour, ls='-')
	#subplot.plot(x, y_log, ls=lineStyle, c=colour)
	subplot.scatter(x, y_log, marker='o', s=10, facecolor=colour, edgecolor = colour)
	#plot schechter line

#loop throught plots
for j in range(j_range):
	for i in range(i_range):
		sub = subs[j,i]
		para = plots[j,i]
		print para

		if method[para] == 'FIR': bin_cat = cat.where(cat.D250 ==1)
		elif method[para] == 'gas': bin_cat = cat.where(np.nan_to_num(cat.HI_ALL2) > 4.2)
		else: bin_cat = cat 

		#create bins
		bins = np.linspace(min(bin_cat[para]),max(bin_cat[para]), N+1)

		
		#loop through cats
		selection = cat
		if method[para] == 'FIR': selection = selection.where(selection.D250 ==1)
		elif para == 'HI_ALL2': selection = selection.where(selection[para] > 4.2 )
		elif para == 'SMASS': selection = selection.where(selection[para] > 4.0 )
		#elif para == 'SRF': selection = selection.where(selection.bptclass == 1)
		#for l in range(len(types)):
		#	subCat = selection.where(selection[types[l]]==1)
		#	functionPlotter(subCat[para], volumes[k], bins, line_syles[k], colours[l], sub, guesses[i])

		functionPlotter(selection[para], bins, 'k', sub)


		sub.yaxis.set_major_locator(MaxNLocator(5))
		sub.xaxis.set_major_locator(MaxNLocator(4))
		#ploting options
		sub.set_xlabel(labels[para])
		sub.axvline(x=lines[i], ls='--', c='k')
		sub.set_xlim(limits[i][0], limits[i][1])
		#sub.tick_params(axis = 'x', labelbottom = 'off')


subs[0,0].set_ylabel('$\log_{10}$(N)')

plt.subplots_adjust(left=0.1, bottom=0.15, right=0.97, top=0.97, wspace=0.15, hspace=0.3)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/mass_hists.pdf')
plt.show()



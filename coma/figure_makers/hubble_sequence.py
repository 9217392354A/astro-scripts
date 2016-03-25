# various parameters across the hubble sequence
# using the goldmine morphology
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
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

#convert radius vir to log10
cat.add_column('G2S', cat.GMASS - cat.SMASS)
cat.add_column('G2D', cat.GMASS - cat.DMASS)
cat.add_column('SFR2D', cat.SRF - cat.DMASS)

w1 = np.where(cat.goldmine == 13)[0]
cat.goldmine[w1] = 9
w2 = np.where(cat.goldmine == 18)[0]
cat.goldmine[w2] = 8

#select detected galaxies
#cat = cat.where(np.nan_to_num(cat.K) != 0.0)
cat = cat.where((np.nan_to_num(cat.goldmine+10) != 0.0) & (cat.goldmine > -1000.0) & (cat.D250 == 1))

cols = ['SMASS', 'DUST_STARS', 'DMASS', 'sSFR','SRF', 'SFR2D']
labs = ['$\log_{10} (M_{stars} / $M$_{\odot}$)', '$\log_{10} (M_{dust} / M_{star}$)', '$\log_{10} (M_{dust} / $M$_{\odot}$)','$\log_{10}$(sSFR) (Yr$^{-1}$)','$\log_{10}$(SFR) (Yr$^{-1}$)', '$\log_{10}$($SFR/M_{Dust}$)']



#create figure and subplots
fig, subs = plt.subplots(nrows=len(cols), ncols=1, sharex=True, sharey=False, squeeze=False, figsize = (8., 9.5), facecolor='w',edgecolor='w')

##################################################################################
#function that bins up something and then finds errors that are error in mean
def mean_binup(x, y, bins):
	#intialise list to hold y values
	y_mean, x_mean, y_error, y_n, y_std = [], [], [], [], []

	if True:
		#loop through bins and bin up
		for i in range(len(bins)-1):
			#limits
			lower, upper = bins[i], bins[i+1]
			width = upper-lower
			#select y values that fall in bins
			y_vals = y[(x > lower ) & (x < upper)]

			#mean, error in the mean, and midbin val
			y_n.append(np.float(len(y_vals)))
			y_mean.append(np.mean(y_vals))
			y_error.append(np.std(y_vals)/ np.sqrt(len(y_vals)))
			x_mean.append((lower+upper)/2.0)
			y_std.append(np.std(y_vals))


	return np.array(x_mean), np.array(y_mean), np.array(y_error), np.array(y_n), np.array(y_std)
##################################################################################



cats = [cat.where(cat.RADIUS_VIR <= 1.0), cat.where(cat.RADIUS_VIR > 1.0)]
fcs = ['None', 'k']
colours = ['r', 'b']
mkrs = ['o', '*']

for j in range(len(cols)):
	sub = subs[j,0]
	col = cols[j]


	for ll in range(2):
		selection = cats[ll]
		sub_cat = selection.where((np.nan_to_num(selection[cols[j]]) != 0.0) & (selection[cols] > -100.0))

		x_raw = sub_cat.goldmine
		y_raw = sub_cat[cols[j]]

		#find mean and errors for each point
		x, y, y_err, y_n, y_std = mean_binup(x_raw, y_raw, np.arange(-0.5,10.5,1.))
		sub.scatter(x_raw, y_raw, s=40, marker=mkrs[ll], facecolor='k', edgecolor='None')

		sub.plot(x, y+y_std, ls = '--', c=colours[ll]) #this is tos
		sub.plot(x, y-y_std, ls = '--', c=colours[ll]) #this is tos

		sub.errorbar(x, y, yerr=y_err, c=colours[ll], ls='-')
		sub.scatter(x, y, marker=mkrs[ll], s=60, c=colours[ll])

	sub.xaxis.set_major_locator(MaxNLocator(12))
	sub.yaxis.set_major_locator(MaxNLocator(4))
	sub.set_xlim(-1,10)
	sub.set_ylabel(labs[j])

	sub.yaxis.set_label_coords(-0.08, 0.5)

	sub.grid()

#sub.set_ylim()
sub.set_xticklabels( ('','E', 'S0', 'S0/Sa', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'S', 'Pec','') )
plt.subplots_adjust(left=0.12, bottom=0.03, right=0.99, top=0.99, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/hubble_sequence.pdf')
plt.show()



#morph converter
#This program takes surface density and then looks at the fraction detected for a given morphological type
# Chris Fuller, Feb 2014


#import moduals
from atpy import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#switches
no_morph = 3 #number of morphology bins between 0 and 1
no_density = 4 #number of sigma and radial density bins between logmin and logmax
sfactor = 0.0
colours = ['r', 'g', 'b']
types = ['early', 'inter', 'late']
from lmfit import minimize, Parameters,report_fit

detect_col = 'F250' #'F250' 

#Inputs 
cat = Table("/Users/chrisfuller/Documents/phd/herchel/coma/final_outputs/coma_supercluster_cal12_pacscorrected.fits")
#cat = cat.where(cat.SMASS > 9.1)

cat.RADIUS_VIR = np.log10(cat.RADIUS_VIR)

###################### functions ######################
def er(a,b):
    return np.sqrt(a)  / b

def super_histogram(a, bins, sfactor):
	#create temp array to hold histogram
	hist = []

	#loop through bins
	for ii in range(0, len(bins)-1):
		width = bins[ii+1] - bins[ii]
		lower = bins[ii] - width*sfactor
		upper = bins[ii+1] + width*sfactor

		#select bins
		w = np.where((a >= lower) & (a <= upper))[0]
		hist.append(len(a[w]))

	return np.array(hist)


def s(x):
	return str(abs(np.round(x, decimals=3)))

#function for finding errors in a straight line of two arras x and y
# assuming least squares and equal error on each point

def residual(params, x, y_data, y_error):
	y_model = model(params, x)

	return (y_data-y_model)/y_error

def lmfitter(x , y, y_error):
	params = Parameters()
	params.add('m', value=1., vary=True)
	params.add('c', value=0., vary=True)


	# remove inf values in errors
	out = minimize(residual, params, args=(x, y, y_error))
	#report_fit(params)
	return out

def model(params, x):
	m = params['m'].value
	c = params['c'].value

	return m*x + c
################### end functions #####################

#create figure and subplots
fig, subs = plt.subplots(nrows=2, ncols=4, sharex=False, sharey=False, squeeze=False, figsize = (8.,4.8),facecolor='w',edgecolor='w')

#density colums that we care about
cols = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']
#cols  = ['SIGMA1_1000', 'SIGMA5_1000', 'SIGMA10_1000', 'RADIUS_VIR']
xlabs = ['$Log_{10}(\Sigma 1)($Mpc$^{-2})$','$Log_{10}(\Sigma 5)($Mpc$^{-2})$','$Log_{10}(\Sigma 10) ($Mpc$^{-2})$', '$\log_{10}$(R/R$_{Virial}$)']



#loop through density estimators
for i in range(0,len(cols)):
	#define cols 
	col = cols[i]

	#define density bins ****** these are created in log space ********
	bins_dens  = np.linspace(np.min(cat[col]),np.max(cat[col]),no_density+1)

	#define morphological bins  ****** these are created in log space ********
	bins_morph = np.linspace(0,1,no_morph+1)


	#Add title to subplot
	#subs[0,i].text(0.05, 0.95, col, transform=subs[0,i].transAxes, fontsize=12, verticalalignment='top')

	#create empty array to hold fraction of morphological types
	temp_totals = []
	
	print col
	#loop through morphology
	for j in range(0, no_morph):
		
		#print morph_lower, morph_upper
		#selected galaxies that are in morphological bin
		selection_detec = cat.where((cat[types[j]] == 1) & (cat[types[j]] == 1) & (np.nan_to_num(cat[detect_col]) != 0.0)) #total
		selection_total = cat.where((cat[types[j]] == 1) & (cat[types[j]] == 1)) #total

		#print bins_dens, len(selection_detec), len(selection_total)
		#create temp array to hold fraction detected
		if False:
			detected, _ = np.histogram(selection_detec[col], bins = bins_dens)
			totals  , _ = np.histogram(selection_total[col], bins = bins_dens)

		if True:
			detected = super_histogram(selection_detec[col], bins_dens, sfactor)
			totals = super_histogram(selection_total[col], bins_dens, sfactor)

		#fraction detected of a given morphological type over a range of densityis
		per_det = detected / totals.astype(dtype=np.float)
		print per_det
		# # # # #plot fraction detected # # # #

		x1 = bins_dens[:-1]+(bins_dens[1]-bins_dens[0])/2.0
		y1 = np.nan_to_num(detected / totals.astype(dtype=np.float))
		yerr1 = np.nan_to_num(er(detected,totals.astype(dtype=np.float)))
		
		#y1 = y[y1==np]
		yerr1[yerr1==0.0] = 0.5

		#fit straight line
		fit = lmfitter(x1 , y1, yerr1)
		m, dm = fit.params['m'].value, fit.params['m'].stderr
		c, dc = fit.params['c'].value, fit.params['c'].stderr
		x_line = x1
		y_line = model(fit.params, x_line)

		if col == "RADIUS_VIR": 
			linestyle1 = ' '
			linestyle2 = '--'
			alpha1 = 0.5
			alpha2 = 1.0

		else:
			linestyle1 = ' '
			linestyle2 = '--'
			alpha1 = 0.5
			alpha2 = 1.0


		subs[0,i].errorbar(x1, y1, yerr=yerr1, color=colours[j], label=types[j], ls =linestyle1, alpha=alpha1)
		subs[0,i].scatter(x1,y1, color=colours[j], edgecolor=colours[j], s=60, alpha=0.5)
		subs[0,i].plot(x_line, y_line, color=colours[j], ls=linestyle2, alpha=alpha2)


		#options 
		#subs[0,i].yaxis.set_major_locator(MaxNLocator(5))
		subs[0,i].xaxis.set_major_locator(MaxNLocator(5))
		subs[0,i].set_ylim(0.0, 0.95)

		subs[0,i].tick_params(axis='x', labelbottom='off')

		if i != 0: 
			subs[0,i].tick_params(axis='y', labelleft='off')

		if i!=  len(cols)-1:
			subs[0,i].set_xlim(-1.4, 4.1)
			#subs[0,i].invert_xaxis()

		if i == len(cols)-1:
			subs[0,i].set_xlim(-2.7,1.5)

		#append temp_totals for next phase
		temp_totals.append(totals.astype(dtype=np.float))		

	#plot fraction of any given type dressler 1980 M-D relation
	#loop through each morphology and plot fraction of any morphological type on sub2
	for k in range(0,no_morph):
		focus = temp_totals[k]

		temp = np.zeros(no_density)
		#loop through temp_totals to create 
		for l in range(0,no_morph): temp += temp_totals[l]

		x2 = bins_dens[:-1]+(bins_dens[1]-bins_dens[0])/2.0
		y2 = focus / temp
		yerr2 = er(focus , temp)

		yerr2[yerr2==0.0] = 0.5

		#fit straight line
		fit = lmfitter(x2 , y2, yerr2)
		m, dm = fit.params['m'].value, fit.params['m'].stderr
		c, dc = fit.params['c'].value, fit.params['c'].stderr
		x_line_2 = x2
		y_line_2 = model(fit.params, x_line_2)

		if col == "RADIUS_VIR": 
			linestyle1 = ' '
			linestyle2 = '--'
			alpha1 = 0.5
			alpha2 = 1.0

		else:
			linestyle1 = ' '
			linestyle2 = '--'
			alpha1 = 0.5
			alpha2 = 1.0


		subs[1,i].errorbar(x2, y2, yerr=yerr2, color=colours[k], label=types[k], ls = linestyle1, alpha = alpha1)
		subs[1,i].scatter(x2,y2, color=colours[k], edgecolor=colours[k], s=60, alpha=0.5)
		subs[1,i].plot(x_line_2, y_line_2, color=colours[k], ls = linestyle2, alpha = alpha2)


		
		#subs[1,i].yaxis.set_major_locator(MaxNLocator(5))
		subs[1,i].xaxis.set_major_locator(MaxNLocator(5))
		subs[1,i].set_ylim(0.0, 0.78)

		if i ==0:
			subs[0,0].set_ylabel('Fraction Detected')
			subs[1,0].set_ylabel('Fraction of Population')

		subs[1,i].set_xlabel(xlabs[i])

		if i != 0: 
			subs[1,i].tick_params(axis='y', labelleft='off')
	
		if i!= len(cols)-1:
			subs[1,i].set_xlim(-1.4, 4.1)
			#subs[1,i].invert_xaxis()
		if i == len(cols)-1:
			subs[1,i].set_xlim(-2.7,1.5)

	#control what axis are hidden, and x and y lim

	#hide all x axies for sub[0,i]
	#subs.yaxis.set_major_locator(MaxNLocator(5))
	

#subs[1,i].legend()
plt.subplots_adjust(left=0.07, bottom=0.12, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
#fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/density_vs_radius_log.pdf')
plt.show()
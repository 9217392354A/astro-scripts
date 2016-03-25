#This program takes surface density and then looks at the fraction detected for a given morphological type 
# in addition to sigma_vs_radius_fdetected.py it adds fornax and virgo on aswell
# Chris Fuller, March 2014


#import moduals
from atpy import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from os.path import join as pj

#switches
no_morph = 3 #number of morphology bins between 0 and 1
no_density = 5 #number of sigma and radial density bins between logmin and logmax
sfactor = 1.0
colours = ['r', 'g', 'b']
types = ['Early', 'Inter', 'Late']


morph= 'pS0' #columb morphology data is taken from
detect_col = 'F250_ngp_test' #'F250' 

#Inputs
#cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12.fits")
cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/test.fits")

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

################### end functions #####################

#create figure and subplots
fig, subs = plt.subplots(nrows=2, ncols=4, sharex=False, sharey=False, squeeze=False, figsize = (10.5,6.5),facecolor='w',edgecolor='w')

#density colums that we care about
cols = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']
xlabs = ['$Log_{10}(\Sigma 1)($Mpc$^{-2})$','$Log_{10}(\Sigma 5)($Mpc$^{-2})$','$Log_{10}(\Sigma 10) ($Mpc$^{-2})$', 'R / R$_{virial}$']

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
		morph_lower = bins_morph[j] 
		morph_upper = bins_morph[j+1]

		#print morph_lower, morph_upper
		#selected galaxies that are in morphological bin
		selection_detec = cat.where((cat[morph] > morph_lower) & (cat[morph] < morph_upper) & (np.nan_to_num(cat[detect_col]) != 0.0)) #total
		selection_total = cat.where((cat[morph] > morph_lower) & (cat[morph] < morph_upper)) #total

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
		subs[0,i].errorbar(bins_dens[:-1]+(bins_dens[1]-bins_dens[0])/2.0, detected / totals.astype(dtype=np.float), yerr=er(detected,totals.astype(dtype=np.float)), color = colours[j], label=types[j])

		#options 
		subs[0,i].yaxis.set_major_locator(MaxNLocator(4))
		subs[0,i].xaxis.set_major_locator(MaxNLocator(4))
		subs[0,i].set_ylim(0.0, 0.8)

		subs[0,i].tick_params(axis='x', labelbottom='off')

		if i != 0: 
			subs[0,i].tick_params(axis='y', labelleft='off')

		if i!=  len(cols)-1:
			subs[0,i].set_xlim(-1.4, 4.1)
			subs[0,i].invert_xaxis()

		#if i == len(cols)-1:
			#subs[1,i].set_xlim(0,20)
		#append temp_totals for next phase
		temp_totals.append(totals.astype(dtype=np.float))		

	#plot fraction of any given type dressler 1980 M-D relation
	#loop through each morphology and plot fraction of any morphological type on sub2
	for k in range(0,no_morph):
		focus = temp_totals[k]

		temp = np.zeros(no_density)
		#loop through temp_totals to create 
		for l in range(0,no_morph): temp += temp_totals[l]

		#plot dressler 1980 fraction of any given type 
		subs[1,i].errorbar(bins_dens[:-1]+(bins_dens[1]-bins_dens[0])/2.0, focus / temp, yerr=er(focus , temp), color = colours[k], label=types[k])
		
		subs[1,i].yaxis.set_major_locator(MaxNLocator(4))
		subs[1,i].xaxis.set_major_locator(MaxNLocator(4))
		subs[1,i].set_ylim(0.0, 0.95)

		if i ==0:
			subs[0,0].set_ylabel('Fraction Detected')
			subs[1,0].set_ylabel('Fraction of Population')

		subs[1,i].set_xlabel(xlabs[i])

		if i != 0: 
			subs[1,i].tick_params(axis='y', labelleft='off')
	
		if i!= len(cols)-1:
			subs[1,i].set_xlim(-1.4, 4.1)
			subs[1,i].invert_xaxis()
		#if i == len(cols)-1:
			#subs[1,i].set_xlim(0,20)

	#control what axis are hidden, and x and y lim

	#hide all x axies for sub[0,i]
	#subs.yaxis.set_major_locator(MaxNLocator(5))




"""
############################## add fornax and virgo #################################

#Inputs
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'

c2 = 'virgo_at_100mpc.fits'
c3 = 'fornax_at_100mpc.fits'

cats = [Table(pj(folder2,c2)), Table(pj(folder2,c3))] 

#sigma stats
colours2 = [ 'r', 'b' ]
namesy = [ 'early', 'late']
cluster = ['Virgo', 'Fornax']
linestypes = ['-.', '--']
for jj in range(len(cats)):
	#loop through density estimators
	cat = cats[jj]
	clust = cluster[jj]
	for i in range(0,len(cols)):
		#define cols 
		col = cols[i]

		#define density bins ****** these are created in log space ********
		bins_dens  = np.linspace(np.min(cats[0][col]),np.max(cats[0][col]),no_density+1)

		#define morphological bins  ****** these are created in log space ********
		bins_morph = [-3,2,20]


		#Add title to subplot
		#subs[0,i].text(0.05, 0.95, col, transform=subs[0,i].transAxes, fontsize=12, verticalalignment='top')

		#create empty array to hold fraction of morphological types
		temp_totals = []
		
		print col
		#loop through morphology
		for j in range(2):
			print j
			morph_lower = bins_morph[j] 
			morph_upper = bins_morph[j+1]

			#print morph_lower, morph_upper
			#selected galaxies that are in morphological bin
			selection_detec = cat.where((cat['goldmine'] > morph_lower) & (cat['goldmine'] < morph_upper) & (np.nan_to_num(cat[detect_col]) != 0.0)) #total
			selection_total = cat.where((cat['goldmine'] > morph_lower) & (cat['goldmine'] < morph_upper)) #total

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
			subs[0,i].errorbar(bins_dens[:-1]+(bins_dens[1]-bins_dens[0])/2.0, detected / totals.astype(dtype=np.float), yerr=er(detected,totals.astype(dtype=np.float)), ls= linestypes[jj], color = colours2[j], label= cluster[jj]+' '+namesy[j])
			temp_totals.append(totals.astype(dtype=np.float))	

			print  colours2[j], cluster[jj], namesy[j]
		#loop through each morphology and plot fraction of any morphological type on sub2
		for k in range(2):
			focus = temp_totals[k]

			temp = np.zeros(no_density)
			#loop through temp_totals to create 
			for l in range(2): temp += temp_totals[l]

			#plot dressler 1980 fraction of any given type 
			subs[1,i].errorbar(bins_dens[:-1]+(bins_dens[1]-bins_dens[0])/2.0, focus / temp, yerr=er(focus , temp), ls=linestypes[jj], color = colours2[k], label=cluster[jj]+' '+namesy[k])
			
			

	#hide all x axies for sub[0,i]
	#subs.yaxis.set_major_locator(MaxNLocator(5))

"""
subs[1,i].legend()
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Desktop/density_vs_radius+virgo+fornax.pdf')
plt.show()
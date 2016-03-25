#morph converter
#This program takes surface density and then looks at the fraction detected for a given morphological type
# Chris Fuller, Feb 2014


#import moduals
from atpy import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from lmfit import minimize, Parameters,report_fit

#switches
no_morph = 3 #number of morphology bins between 0 and 1
no_density = 4 #number of sigma and radial density bins between logmin and logmax
sfactor = 0.0
colours = ['b', 'g', 'r']
types = ['late','inter', 'early' ]


detect_col = 'F250' #'F250' 

#Inputs
cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12.fits")

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
#density colums that we care about

cols = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']

#switchs
N_start = 1
N_end = 10

Ns = np.arange(N_start,)
Ns= [1,5,10]

#quanties on x axcis
cols_pre  = 'SIGMA'


#create list to hold gradient and y-cepts
mm_list = []
dm_list = []

cc_list =[]
dc_list =[]


mm2_list = []
dm2_list = []


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

	mm_temp = []
	dm_temp = []

	mm2_temp = []
	dm2_temp = []

	#loop through morphology
	for j in range(3):
		
		#print morph_lower, morph_upper
		#selected galaxies that are in morphological bin
		selection_detec = cat.where((cat[types[j]] == 1) & (np.nan_to_num(cat[detect_col]) != 0.0)) #total
		selection_total = cat.where(cat[types[j]] == 1) #total

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
		# # # # #plot fraction detected # # # #
		x1 = bins_dens[:-1]+(bins_dens[1]-bins_dens[0])/2.0
		y1 = np.nan_to_num(detected / totals.astype(dtype=np.float))
		yerr1 = np.nan_to_num(er(detected,totals.astype(dtype=np.float)))

		yerr1[yerr1==0.0] = 0.5

		#fit straight line
		fit = lmfitter(x1 , y1, yerr1)
		m, dm = fit.params['m'].value, fit.params['m'].stderr
		c, dc = fit.params['c'].value, fit.params['c'].stderr
		x_line = x1
		y_line = model(fit.params, x_line)

		temp_totals.append(totals.astype(dtype=np.float))		

		
		mm_temp.append(m)
		dm_temp.append(dm)
	
	mm_list.append(mm_temp)
	dm_list.append(dm_temp)



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

		mm2_temp.append(m)
		dm2_temp.append(dm)
	
	mm2_list.append(mm2_temp)
	dm2_list.append(dm2_temp)
########################################################################################################################
########################################################################################################################


# Start Second Plot of grad and y-cept #


########################################################################################################################
########################################################################################################################

fig2, subs = plt.subplots(nrows=1, ncols=1, sharex=True, figsize = (4.,4.),facecolor='w',edgecolor='w')
early_mm = []
inter_mm = []
late_mm = []

early_dm = []
inter_dm = []
late_dm = []

early_mm2 = []
inter_mm2 = []
late_mm2 = []

early_dm2 = []
inter_dm2 = []
late_dm2 = [] 

#turn mm_col into mm_type
for ii in range(len(cols)):
	#create MM and DM for type
	early_mm.append(mm_list[ii][0])
	inter_mm.append(mm_list[ii][1])
	late_mm.append(mm_list[ii][2])

	mms = [early_mm, inter_mm, late_mm]


	early_dm.append(dm_list[ii][0])
	inter_dm.append(dm_list[ii][1])
	late_dm.append(dm_list[ii][2])

	dms = [early_dm, inter_dm, late_dm]

	#create MM and DM for type
	early_mm2.append(mm2_list[ii][0])
	inter_mm2.append(mm2_list[ii][1])
	late_mm2.append(mm2_list[ii][2])

	mms2 = [early_mm2, inter_mm2, late_mm2]


	early_dm2.append(dm2_list[ii][0])
	inter_dm2.append(dm2_list[ii][1])
	late_dm2.append(dm2_list[ii][2])

	dms2 = [early_dm2, inter_dm2, late_dm2]	
	
#loop through morphology
for j in range(3):
	xxx = [1,2,3,4]
	yyy =np.round( mms[j], decimals=3)
	ery =np.round( dms[j], decimals=3)

	yyy2 = np.round(mms2[j], decimals=3)
	ery2 = np.round(dms2[j], decimals=3)

	sub = subs

	print ''
	print types[j], '\\\\'
	print 'fraction detected','&', yyy[0],'$\pm$', ery[0], '&', yyy[1],'$\pm$', ery[1], '&', yyy[2],'$\pm$', ery[2], '&', yyy[3],'$\pm$', ery[3], '\\\\'
	print 'fraction population','&', yyy2[0],'$\pm$', ery2[0], '&', yyy2[1],'$\pm$', ery2[1], '&', yyy2[2],'$\pm$', ery2[2], '&', yyy2[3],'$\pm$', ery2[3], '\\\\'
	print ''
	#sub.errorbar(xxx, yyy, yerr=ery, color = colours[j], label=types[j], ls='-')
	#sub.scatter(xxx, yyy, color = colours[j], label=types[j])
	
	#sub.errorbar(xxx, yyy2, yerr=ery2, color = colours[j], label=types[j], ls='--')
	#sub.scatter(xxx, yyy2, color = colours[j], label=types[j])

sub.set_xticklabels(['', '1', '', '5', '','10','' ])
sub.set_xlim(0.5, 3.5)
#sub.set_ylim(-0.2, 0.1)


sub.set_ylabel('m')
sub.set_xlabel('$\Sigma_{N}$')
#plt.subplots_adjust(left=0.18, bottom=0.11, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
#plt.show()

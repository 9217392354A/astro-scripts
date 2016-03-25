

#program to plot FIR properties 
# Chris Fuller, April 2014

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
cat.add_column('DUST_STARS_BELL', cat.DMASS - cat.SMASS_BELL)
cat.add_column('all', cat.g)
cat.all = 1
#caculated extra columns
cat.add_column('D2G', - cat.HI_ALL2 +cat.DMASS)
cat.add_column('G2S', cat.HI_ALL2 - cat.SMASS)
cat.add_column('SFR2G', cat.SRF - cat.HI_ALL2)
cat.add_column('SFR2D', cat.SRF - cat.DMASS)
cat.add_column('colour', cat.g - cat.r)

#convert radius vir to log10
cat.RADIUS_VIR = np.log10(cat.RADIUS_VIR)

#select currentCat galaxies
firCat = cat.where(cat.DMASS_TYPE != 0)
gasCat = cat.where(np.nan_to_num(cat.HI_ALL2) > 0.0)
metalCat = cat.where(np.nan_to_num(cat.METAL) > 0.0)
gasfirCat = cat.where((cat.DMASS_TYPE != 0) & (np.nan_to_num(cat.HI_ALL2) > 0.0))

#switchs
N_den = 7 #number of density bins
N = 60
bin_type = 'fixed'
types = ['late', 'inter', 'early']#, 'all']
colours = ['b', 'g', 'r']#, 'k']


mms = {} #dic to hold grads
mes = {} # dic to hold grad errrors

#quanties on x axcis
xcols  = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']
#xcols =  ['SIGMA1_1000', 'SIGMA5_1000', 'SIGMA10_1000', 'RADIUS_VIR']
xlabs = ['$\log_{10}(\Sigma 1)($Mpc$^{-2})$','$\log_{10}(\Sigma 5)($Mpc$^{-2})$','$\log_{10}(\Sigma 10) ($Mpc$^{-2})$', '$\log_{10}$(R/R$_{Virial}$)']
xlims = [[-1.2, 4.1], [-1.2,3.3], [-1.2,3.3], [-1.7, 1.5]] 



if True:
	#quanties on the yaxis
	ycols = ['D2G', 'SFR2G', 'SFR2D']
	ylabs = ['$\log_{10} (M_{dust} / M_{star}$)', '$\log_{10} (M_{dust} / $M$_{\odot}$)']
	ylims = [[-4.8, -1.4 ]]*3



#create dict for lables
labels = {	'DMASS' 		: '$\log_{10} (M_{dust} / $M$_{\odot}$)', 
			'SMASS' 		: '$\log_{10} (M_{stars} / $M$_{\odot}$)',
			'DUST_STARS' 	: '$\log_{10} (M_{dust} / M_{star}$)',
			'SRF' 			: '$\log_{10}$(SFR) (M$_{\odot}$ yr$^{-1}$)',
			'sSFR'			: '$\log_{10}$(sSFR) (yr$^{-1}$)',
			'METAL'			: '12 + $\log_{10}$(O/H)',
			'D2G' 			: '$\log_{10} (M_{dust} / M_{gas}$)',
			'G2S' 			: '$\log_{10} (M_{gas} / M_{stars}$)',
			'colour'		: '$g-r$',
			'SFR2D'			: '$\log_{10} ( SFR/ M_{dust})$ (yr$^{-1}$)',
			'HI_ALL2'		: '$\log_{10} (M_{gas} / $M$_{\odot}$)',
			'SFR2G'			: '$\log_{10} ( SFR/ M_{gas})$ (yr$^{-1}$)'}

labels2 = {	'DMASS' 		: '$M_{dust}$', 
			'SMASS' 		: '$M_{stars}$',
			'DUST_STARS' 	: '$M_{dust} / M_{star}$',
			'SRF' 			: 'SFR',
			'sSFR'			: 'sSFR',
			'METAL'			: 'Z',
			'colour'		: 'g-r',
			'SFR2D'			:  'SFR / $M_{dust}$',
			'HI_ALL2' 		: '$M_{gas}$',
			'D2G'			: '$M_{dust} / M_{gas}$',
			'G2S'			: '$M_{gas} / M_{stars}$',
			'SFR2G'			: '$M_{SFR} / M_{gas}$'}

cats = {	'DMASS' 		: firCat, 
			'SMASS' 		: cat,
			'DUST_STARS' 	: firCat,
			'SRF' 			: cat,
			'sSFR'			: cat,
			'METAL'			: metalCat,
			'D2G' 			: gasfirCat,
			'G2S' 			: gasCat,
			'colour'		: cat,
			'SFR2D'			: firCat,
			'HI_ALL2'		: gasCat,
			'SFR2G'			: gasCat,
			'SFR2D'			: firCat }
	
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

 															# Functions #

# function that creates adaptive histogram 															
def adaptive_histogram(x, y, N):
	#intialise list to hold y values
	y_mean, x_mid, y_error, y_n, x_mean = [], [], [], [], []

	#sort arrays
	ii = np.argsort(x)

	#create new arrays sorted
	xx = x[ii]
	yy = y[ii]

	#empty list of bins
	bins = []

	#caculate how many bins needed
	number = np.int(np.ceil(len(xx)/ np.float(N)))

	print 'Number of bins: ', number
	#create intial lower and upper limits
	lower_i = 0
	upper_i = N

	#loop through and add N to each bin
	for i in range(number):
		#if not final bin
		if i != number -1:
			#create index upper and lower limits
			lower_i, upper_i = N*i, N*(i+1) - 1

		#if final bin
		else:
			lower_i, upper_i = N*i, len(ii) - 1 

		#physical bin widths
		lower, upper = xx[lower_i], xx[upper_i]

		#find width of bin
		width = upper-lower

		#select y values that fall in bins
		y_vals = yy[lower_i:upper_i+1]

		#mean, error in the mean, and midbin val
		y_n.append(np.float(len(y_vals)))
		y_mean.append(np.mean(y_vals))
		y_error.append(np.std(y_vals)/ np.sqrt(len(y_vals)))
		x_mid.append((lower+upper)/2.0)
		x_mean.append(np.mean(xx[lower_i:upper_i+1]))
	
	return y_n, y_mean, y_error, x_mid, x_mean


#function that bins up something and then finds errors that are error in mean
def mean_binup(x, y, bins):
	#intialise list to hold y values
	y_mean, x_mean, y_error, y_n = [], [], [], []

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

	if bin_type == 'adaptive':
		print 'adaptive method...'
		y_n, y_mean, y_error, x_mid, x_mean = adaptive_histogram(x, y, N)

	return np.array(x_mean), np.array(y_mean), np.array(y_error), np.array(y_n)

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

def fit_line1(x_data , y_data, y_error):

	#fit the data
	out = lmfitter(x_data , y_data, y_error)
	m, dm  = out.params['m'].value, out.params['m'].stderr

	#find redchi sqare
	return x_data , model(out.params, x_data), [m, dm], out.redchi

def fit_line2(x_data , y_data, y_error):

	#fit the data
	out = lmfitter2(x_data , y_data, y_error)
	#m, dm  = out.params['m'].value, out.params['m'].stderr
	return x_data , model2(out.params, x_data), [-100., -100.0], out.redchi

def residual2(params, x, y_data, y_error):
	y_model = model2(params, x)

	return (y_data-y_model)/y_error

def lmfitter2(x , y, y_error):
	params = Parameters()
	params.add('a', value=1., vary=True)
	params.add('b', value=1., vary=True)
	params.add('c', value=1., vary=True)


	# remove inf values in errors
	out = minimize(residual2, params, args=(x, y, y_error))
	#report_fit(params)
	return out

def model2(params, x):
	a = params['a'].value
	b = params['b'].value
	c = params['c'].value

	return a*x**2 + b*x + c

def print_grads(vals, err):
	# create final list that will hold vals
	#final = [[],[],[]]
	#loop through types
	labels = ['late&', 'intermediate&', 'early&']
	line =[]
	for mm in range(len(ycols)):
		print '\\\\'
		print '\\multicolumn{9}{l}{\\bf ' + labels2[ycols[mm]] + '} \\\\'
		#print '\n'*2
		final = [[],[],[]]
		for ii in range(len(types)):
			#loop through des
			line = ''
			line = [labels[ii]]
			value = ''
			count = 0
			for jj in range(len(xcols)):
				key = xcols[jj]+types[ii]+ycols[mm]

				#extract values from dics
				try:
					m, dm = vals[key], err[key]
					grad_diff = m/dm
					if abs(m/dm) > 3.0:
						quest = 'No'
						count +=1
					if abs(m/dm) < 3.0:
						quest = 'Yes'


					grad_diff = str(np.int(abs(np.round(m/dm, decimals=0))))
					value = str(np.round(m, decimals=2)) + '$\\pm$' + str(np.round(dm, decimals=2)) + '&' + quest + '(' +  grad_diff +'$\sigma$)' 



				except:
					value = '-' + '&' +'-'
	 			#final[ii].append(value)
	 			if jj != len(xcols)-1: line.append(value+'&')
	 			else: line.append(value + '\\\\')
	 		#print types[ii]
	 		print " ".join(line)
	 		line = ''

#lineFitter
# If given x,y and yerrors it will fit a straightline
# then if models has a poor fit ie p(chi2,dof) < tolerence
def lineFitter(x_data, y_data, y_error):
	tolerence = 0.1

	#fit linear
	out_lin = lmfitter(x_data , y_data, y_error)
	p_lin, chi_lin, nfree_lin = pStat(out_lin)
	info_1 = x_data , model(out_lin.params, x_data), [out_lin.params['m'].value, out_lin.params['m'].stderr], [p_lin, chi_lin, nfree_lin]

	#test if less than tolerence
	if p_lin < tolerence:
		fit = 'n=1'
		return info_1

	elif p_lin >= tolerence:
		fit = 'n=2'
		out_poly2 = lmfitter2(x_data , y_data, y_error)
		p_poly2, chi_poly2, nfree_poly2 = pStat(out_poly2)
		info_2 = x_data , model2(out_poly2.params, x_data), [-100., -100.], [p_poly2, chi_poly2, nfree_poly2]
		if p_poly2 < tolerence and nfree_poly2 > 2:
			return info_2

	return info_1

	#now test if p_poly2 is greater than p_lin

#for a given chisquare and nfree it cacualtes the pstat
def pStat(out):
	chi = out.redchi * out.nfree
	nfree = out.nfree
	pstat =  1 - scipy.stats.chi2.cdf(chi, nfree)
	return pstat, chi, nfree






##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 


	
#create figure and subplots
fig, subs = plt.subplots(nrows=len(ycols), ncols=len(xcols), sharex=False, sharey=False, squeeze=False, figsize = (8., 6.5), facecolor='w',edgecolor='w')

for j in range(len(ycols)):
	#select col
	col = ycols[j]

	currentCat = cats[col]
	

	#loop throught density tracers
	for i in range(len(xcols)):
		density = xcols[i]
		sub = subs[j,i]


		#create density bins
		#bins_density = np.linspace(np.floor(np.min(currentCat[density])), np.ceil(np.max(currentCat[density])), N_den+1) #automatic
		#bins_density = np.linspace(xlims[j][0], xlims[j][1], N_den+1)  #manual

		#loop through morphological types
		for k in range(len(types)):
			selection = currentCat.where(currentCat[types[k]] == 1) #select morphological type
			bins_density = np.linspace(np.floor(np.min(selection[density])), np.ceil(np.max(selection[density])), N_den+1) #automatic

			#extract raw data
			x_raw, y_raw = selection[density], selection[col]

			#find mean and errors for each point
			x, y, y_err, y_n = mean_binup(x_raw, y_raw, bins_density)

			#clean empty bins
			w_clean = y_n > 3.0

			x, y, y_err, y_n = x[w_clean], y[w_clean], y_err[w_clean], y_n[w_clean]


			sub.errorbar(x, y, yerr=y_err, c=colours[k], ls='')
			sub.scatter(x, y, marker='o', s=y_n*5, c=colours[k], edgecolor='None', alpha=0.4)

			#try fitting two models
			if len(x) > 3:

				#fitting function
				x_line, y_line, coeffs, chiStuff = lineFitter(x, y, y_err)

				#append values to dics
				mms[density+types[k]+col] = coeffs[0] #dic to hold grads
				mes[density+types[k]+col] = coeffs[1] #dic to hold errors

				sub.plot(x_line, y_line, c=colours[k], ls='--', lw=1)

			
		if i == 3: sub.axvline(x=np.log10(1.), c='k', ls='--')

		#plotting options
		sub.xaxis.set_major_locator(MaxNLocator(4))
		sub.yaxis.set_major_locator(MaxNLocator(4))

		#lable x on final y
		if j == len(ycols)-1:
			sub.set_xlabel(xlabs[i])

		#label y on final x
		if i == 0:
			sub.set_ylabel(labels[col])

		#if not on bottom switch off x axis
		if j != len(ycols)-1:
			sub.tick_params(axis = 'x', labelbottom = 'off')

		#if not on the left switch off y axis
		if i != 0:
			sub.tick_params(axis = 'y', labelleft = 'off')

		if i!=  len(xcols)-1:
			subs[j,i].set_xlim(-1., 3.3)
			#subs[0,i].invert_xaxis()

		if i == len(xcols)-1:
			subs[j,i].set_xlim(-1.5,1.2)

		sub.yaxis.set_label_coords(-0.25, 0.5)
		sub.xaxis.set_label_coords(0.5, -0.2)

print_grads(mms, mes)

##### plotting options
plt.subplots_adjust(left=0.11, bottom=0.1, right=0.97, top=0.98, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/ratios-density-v2.pdf')
plt.show()



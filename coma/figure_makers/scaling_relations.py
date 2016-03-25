#scaling relations plots
#Chris Fuller

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MaxNLocator
from scipy.stats import pearsonr 

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

#caculated extra columns
cat.add_column('G2S', cat.GMASS - cat.SMASS)
cat.add_column('G2D', cat.GMASS - cat.DMASS)
cat.add_column('SFR2D', cat.SRF - cat.DMASS)
cat.add_column('colour', cat.g - cat.r)

#select detected galaxies
detected = cat.where(cat.DMASS_TYPE != 0)




#find min dust
min_dust = min(detected.DMASS)

N = 8 #number of y bins

#columns to plot
plots = [	[['DMASS', 'SMASS'], ['SRF', 'DMASS']], 
			[['DUST_STARS', 'SMASS'], ['sSFR','SMASS']],
			[['DMASS','METAL'], ['DUST_STARS', 'colour']]]

plots = np.array(plots)

#create dict for lables
labels = {	'DMASS' 		: '$\log_{10} (M_{dust} / $M$_{\odot}$)', 
			'SMASS' 		: '$\log_{10} (M_{stars} / $M$_{\odot}$)',
			'DUST_STARS' 	: '$\log_{10} (M_{dust} / M_{star}$)',
			'SRF' 			: '$\log_{10}$(SFR) (M$_{\odot}$ yr$^{-1}$)',
			'sSFR'			: '$\log_{10}$(sSFR) (yr$^{-1}$)',
			'METAL'			: '$\log_{10}$(O/H) + 12',
			'colour'		: 'g-r'}

#morphological types
types = ['late', 'inter', 'early']
colours = ['b', 'g', 'r']

types_lab = ['late', 'intermediate', 'early']


coeffs_dict = {} #dic to hold coeffs

# extract i's and j's
j_range = plots.shape[0]
i_range = plots.shape[1]

#create a figure
#create figure and subplots
fig, subs = plt.subplots(nrows=j_range, ncols=i_range,
						figsize = (8., 9.5), facecolor='w',edgecolor='w')



mms = {} #dic to hold grads
mes = {} # dic to hold grad errrors
ppc_val = {} #dic to hold pcc values
ppc_pval = {} # dic to hold pcc p vales

#function that bins up something and then finds errors that are error in mean
def mean_binup(x, y, bins):
	#intialise list to hold y values
	y_mean, x_mean, y_error, y_n = [], [], [], []
	bin_type = 'fixed'
	if bin_type == 'fixed':
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

	#test correlations
	correlation = pearsonr(x_data, y_data)

	#find redchi sqare
	try: rchi = out.redchi	
	except: rchi = 0.0
	return x_data , model(out.params, x_data), [m, dm], rchi, correlation



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




#for a given chisquare and nfree it cacualtes the pstat
def pStat(out):
	chi = out.redchi * out.nfree
	nfree = out.nfree
	pstat =  1 - scipy.stats.chi2.cdf(chi, nfree)
	return pstat, chi, nfree

def printTable(dict):

	for j in range(len(ycols)):
		if len(ycols) != 1: print ycols[j]

		for ll in range(2):
			print '\\\\'
			print multiNotation(du[ll], 9) + '\\\\'
			printCoeff(dict, str(ll), ycols[j]) 

def printCoeff(dict, ll, ycol): #coeffs_dict[density+types[k]+col+ll] = coeffs #dic to hold grads
	for k in range(len(types)):
		line = types_lab[k] + '&'
		for i in range(len(xcols)):
			key = xcols[i] + types[k] + ycol + ll

			vals = dict[key]
			if i != len(xcols)-1: end ='&'
			else: end = '\\\\' 
			line += coeffOutput(vals) + end

		print line 
		#print types[k], types_lab[k]
		line = ''


def coeffOutput(vals):
	m = vals[0]
	dm = vals[1]

	sig = abs(m/dm)

	if sig > 3.:
		s = 'No'
	else:
		s = 'Yes'

	c=vals[2]
	dc=vals[3]

	#return ss(m) + '$\\pm$' + ss(dm) + '&' + ss(c) + '$\\pm$' + ss(dc) + '&' + s + '(' + ss(sig) + '$\\sigma$)' 
	return  s1(m) + '$\\pm$' + s1(dm) + '&' + s + '(' + sa(sig) + '$\\sigma$)' 
def s1(s):
	return str(np.round(s,decimals=2))

def sa(s):
	return str(np.int(np.round(s,decimals=0)))


def multiNotation(text,colN):
	return '\\multicolumn{'+ str(colN) +'}{l}{\\bf '+ text +'}'

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

#loop through i's and j's
for j in range(j_range):
	for i in range(i_range):
		#extract subplot
		sub = subs[j,i]

		#loop through each morphological type
		for k in range(len(types)):
			#select morphological type
			subCat = detected.where(detected[types[k]]==1)

			#extract x and y column headers
			x_col = plots[j,i,1]
			y_col = plots[j,i,0]

			#selected galaxies that are no null values
			subCat = subCat.where((subCat[x_col] > -100.0) & (subCat[y_col]> -100.0))

			#if single scaling relation desired selected true
			if False:
				#extract data
				x, y = subCat[x_col], subCat[y_col]
				#plot data
				sub.scatter(x, y, s=10, marker='o', c=colours[k], alpha=0.2, edgecolor='None')

				#mean bin
				bins= np.linspace(np.floor(np.min(subCat[x_col])), np.ceil(np.max(subCat[x_col])), N+1) #automatic
				x_mean, y_mean, y_error, y_n = mean_binup(x, y, bins)
				#clean empty bins
				w_clean = np.where((y_n > 3.0) & (y_error != 0.0))[0]

				x_mean, y_mean, y_error = x_mean[w_clean], y_mean[w_clean], y_error[w_clean]

				
				#plot line
				sub.errorbar(x_mean, y_mean, yerr=y_error, c=colours[k], ls='')
				sub.scatter(x_mean, y_mean, marker='o', s=20, c=colours[k])


				#fit data
				if len(y_mean) > 2:
					x_line , y_line, coeffs, rchi, corr = fit_line1(x_mean, y_mean, y_error)
					sub.plot(x_line, y_line, c=colours[k], ls='-')

				else:
					coeff = [-100, 100]
					corr = coeff				

				#append values to dics
				mms[x_col + y_col + types[k]] = coeffs[0] #dic to hold grads
				mes[x_col + y_col + types[k]] = coeffs[1] #dic to hold errors

				ppc_val[x_col + y_col + types[k]] = corr[0] #dic to hold pcc values
				ppc_pval[x_col + y_col + types[k]] = corr[1] # dic to hold pcc p vales

			#if cluster vs filament desired set true
			if True:
				cluster_subCat = subCat.where(subCat.RADIUS_VIR <= 1.0)
				filament_subCat = subCat.where(subCat.RADIUS_VIR > 1.0)

				sCats = [cluster_subCat, filament_subCat]

				#plot for first
				for l in range(2):
					subCat = sCats[l]

					#extract data
					x, y = subCat[x_col], subCat[y_col]

					bins= np.linspace(np.floor(np.min(subCat[x_col])), np.ceil(np.max(subCat[x_col])), N+1) #automatic
					x_mean, y_mean, y_error, y_n = mean_binup(x, y, bins)

					#clean empty bins
					w_clean = np.where((y_n > 3.0) & (y_error != 0.0))[0]
					x_mean, y_mean, y_error = x_mean[w_clean], y_mean[w_clean], y_error[w_clean]
									#fit data
					if len(y_mean) > 2:
						x_line , y_line, coeffs, rchi, corr = fit_line1(x_mean, y_mean, y_error)
						line = 1
					else:
						coeff = [-100, 100]
						corr = coeff
						line = 0	

					#append values to dics
					coeffs_dict[density+types[k]+col+str(ll)] = coeffs #dic to hold grads
					
					#plot line
					if l == 0: #cluster
						sub.errorbar(x_mean, y_mean, yerr=y_error, c=colours[k], ls='')
						sub.scatter(x_mean, y_mean, marker='o', s=60, facecolor= 'None', edgecolor = colours[k])
						if line == 1: sub.plot(x_line, y_line, c=colours[k], ls='--')

					if l == 1: #cluster
						sub.errorbar(x_mean, y_mean, yerr=y_error, c=colours[k], ls='')
						sub.scatter(x_mean, y_mean, marker='o', s=60, facecolor= colours[k], edgecolor = colours[k])
						if line == 1: sub.plot(x_line, y_line, c=colours[k], ls='-')

		#plot lims
		if x_col == 'DMASS':
			sub.axvline(x=min_dust, ls='--', c='k')

		if y_col == 'DMASS':
			sub.axhline(y=min_dust, ls='--', c='k')

		if y_col =='DUST_STARS' and x_col == 'SMASS':
			subsub = detected.where((subCat[x_col] > -100.0) & (subCat[y_col]> -100.0))
			#xxx = subsub[x_col]
			xxx = [8.0,11.0]
			yyy = min_dust - xxx
			sub.plot(xxx, yyy, c='k', ls='--')

		
		#label axis
		sub.set_xlabel(labels[x_col])
		sub.set_ylabel(labels[y_col])




printTable(coeffs_dict)

plt.subplots_adjust(left=0.11, bottom=0.08, right=0.97, top=0.97, wspace=0.35, hspace=0.35)
#fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/Scaling_relations.pdf')
plt.show()
#program to plot STELLAR MASS
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
cat.add_column('all', cat.g)
cat.all = 1

#convert radius vir to log10
cat.RADIUS_VIR = np.log10(cat.RADIUS_VIR)

#select detected galaxies
#cat = cat.where(np.nan_to_num(cat.K) != 0.0)
detected = cat.where(cat.DMASS_TYPE != 0)
undetected = cat.where(cat.DMASS_TYPE == 0)

cats = [detected, undetected]
du  = ['FIR-detected', 'FIR-undetected']
lisy = ['-', '--']
mark = ['*', '+']


#switchs
N_den = 6 #number of density bins
N = 60
bin_type = 'fixed'
types = ['late', 'inter', 'early', 'all']

types_lab = ['late', 'intermediate', 'early', 'all']
colours = ['b', 'g', 'r', 'k']


coeffs_dict = {} #dic to hold coeffs


#quanties on x axcis
xcols  = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']
xlabs = ['$\log_{10}(\Sigma 1)($Mpc$^{-2})$','$\log_{10}(\Sigma 5)($Mpc$^{-2})$','$\log_{10}(\Sigma 10) ($Mpc$^{-2})$', '$\log_{10}$(R/R$_{Virial}$)']
xlims = [[-1.2, 4.1], [-1.2,3.3], [-1.2,3.3], [-1.7, 1.5]] 



if True:
	#quanties on the yaxis
	ycols = ['SMASS']
	ylabs = ['$\log_{10} (M_{stars} / $M$_{\odot}$)']
	ylims = [[-4.8, -1.4 ]]
	texts = ['Brinchmann et al. (2004)', 'Bell et al. (2003)']
	

#create figure and subplots
fig, subs = plt.subplots(nrows=1, ncols=len(xcols), sharex=False, sharey=True, squeeze=False, figsize = (8., 3.), facecolor='w',edgecolor='w')

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

 															# Functions #

def print_test_K():
	print 'TEST K BAND'
	kcat = cat.where(np.nan_to_num(cat.K) != 0.0)
	for kk in range(len(types)):
		cluster = kcat.where((kcat[types[kk]]==1.0) & (kcat.RADIUS_VIR <= 0.0))
		filament = kcat.where((kcat[types[kk]]==1.0) & (kcat.RADIUS_VIR > 0.0))

		mean_cluster, err_cluster = np.round(np.mean(cluster.K), decimals=3), np.round(np.std(cluster.K)/np.sqrt(len(cluster)*1.0), decimals=3)
		mean_filament, err_filament = np.round(np.mean(filament.K), decimals=3), np.round(np.std(filament.K)/np.sqrt(len(filament)*1.0), decimals=3)


		print types[kk]
		print 'cluster: <K$_{cluster}$> = ', mean_cluster ,'$\pm$',err_cluster
		print 'filament: <K$_{filament}$> = ', mean_filament ,'$\pm$',err_filament
		print 'difference ($\sigma$): (<K$_{cluster}$> - <K$_{filament}$>) = ', (mean_cluster - mean_filament) / (err_cluster**2.0 + err_filament**2.0)**0.5

def print_test_sdss(band):
	print ''
	print 'TEST SDSS BAND', band
	kcat = cat
	for kk in range(len(types)):
		cluster = kcat.where((kcat[types[kk]]==1.0) & (kcat.RADIUS_VIR <= 0.0))
		filament = kcat.where((kcat[types[kk]]==1.0) & (kcat.RADIUS_VIR > 0.0))

		mean_cluster, err_cluster = np.round(np.mean(cluster[band]), decimals=3), np.round(np.std(cluster[band])/np.sqrt(len(cluster)*1.0), decimals=3)
		mean_filament, err_filament = np.round(np.mean(filament[band]), decimals=3), np.round(np.std(filament[band])/np.sqrt(len(filament)*1.0), decimals=3)


		print types[kk]
		print 'cluster: <K$_{cluster}$> = ', mean_cluster ,'$\pm$',err_cluster
		print 'filament: <K$_{filament}$> = ', mean_filament ,'$\pm$',err_filament
		print 'difference ($\sigma$): (<K$_{cluster}$> - <K$_{filament}$>) = ', (mean_cluster - mean_filament) / (err_cluster**2.0 + err_filament**2.0)**0.5
	print ''

def print_test_sdss2():

	print ' TEST SDSS 2 - MEAN ALL BANDS'
	mean_mag = cat.u_1
	sdss_bands = [ 'g', 'r', 'i', 'z']
	for band in sdss_bands:
		mean_mag += cat[band]

	mean_mag = mean_mag/5.0

	cat.add_column('MEAN_SDSS', mean_mag)

	kcat = cat
	for kk in range(len(types)):
		cluster = kcat.where((kcat[types[kk]]==1.0) & (kcat.RADIUS_VIR <= 0.0))
		filament = kcat.where((kcat[types[kk]]==1.0) & (kcat.RADIUS_VIR > 0.0))

		mean_cluster, err_cluster = np.round(np.mean(cluster['MEAN_SDSS']), decimals=3), np.round(np.std(cluster['MEAN_SDSS'])/np.sqrt(len(cluster)*1.0), decimals=3)
		mean_filament, err_filament = np.round(np.mean(filament['MEAN_SDSS']), decimals=3), np.round(np.std(filament['MEAN_SDSS'])/np.sqrt(len(filament)*1.0), decimals=3)


		print types[kk]
		print 'cluster: <K$_{cluster}$> = ', mean_cluster ,'$\pm$',err_cluster
		print 'filament: <K$_{filament}$> = ', mean_filament ,'$\pm$',err_filament
		print 'difference ($\sigma$): (<K$_{cluster}$> - <K$_{filament}$>) = ', (mean_cluster - mean_filament) / (err_cluster**2.0 + err_filament**2.0)**0.5
	print ''

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

	#find redchi sqare
	try: rchi = out.redchi	
	except: rchi = 0.0
	return x_data , model(out.params, x_data), [m, dm], rchi

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
	print '\\begin{table*}'
	print '\\begin{tabular}{ccccc}'
	print '\\hline'
	print '~& \\multicolumn{4}{c}{Gradient of straight line fit (m)} \\\\ '
	print 'Sample type & $\\Sigma_{1}$ & $\\Sigma_{5}$ & $\\Sigma_{10}$ & Radius \\\\'
	print '\\hline'
	print '\\\\'
	print '\\multicolumn{5}{l}{\\bf Detected}} \\\\'


	# create final list that will hold vals
	#loop through types
	for ii in range(len(types)):
		#loop through des
		for jj in range(len(xcols)):
			key = xcols[jj]+types[ii]+'0'

			#extract values from dics
			m, dm = vals[key], err[key]

			value = str(np.round(m, decimals=2)) + '$\\pm$' + str(np.round(dm, decimals=2))

 			final[ii].append(value)

 	for kk in range(3):
 		line = []
 		for pp in range(len(xcols)):
 			if pp != len(xcols) - 1: line.append(final[kk][pp]+'&')
 			else: line.append(final[kk][pp] + '\\\\')

 		#print types[kk]
 		print types_lab[kk] + " ".join(line)

	print '\\\\'
	print '\\multicolumn{5}{l}{\\bf Non-detected}} \\\\'
	# create final list that will hold vals
	final = [[],[],[]]
	#loop through types
	for ii in range(len(types)):
		#loop through des
		for jj in range(len(xcols)):
			key = xcols[jj]+types[ii]+'1'

			#extract values from dics
			m, dm, c, dc = vals[key][0], vals[key][1], vals[key][2], vals[key][3] 

			value = str(np.round(m, decimals=2)) + '$\\pm$' + str(np.round(dm, decimals=2))

 			final[ii].append(value)

 	for kk in range(3):
 		line = []
 		for pp in range(len(xcols)):
 			if pp != len(xcols) - 1: line.append(final[kk][pp]+'&')
 			else: line.append(final[kk][pp] + '\\\\')

 		#print types[kk]
 		print types_lab[kk] + " ".join(line)

 	print '\\hline'
 	print '\\end{tabular}'
 	print '\\caption{Above shows the gradients for the straight line fit parameters from Figure~\\ref{fig:stars-det_vs_non}.}'
 	print '\\label{tab:stars_det_non-grads}'
 	print '\\end{table*}'


#lineFitter
# If given x,y and yerrors it will fit a straightline
# then if models has a poor fit ie p(chi2,dof) < tolerence
def lineFitter(x_data, y_data, y_error):
	tolerence = 0.1

	#fit linear
	out_lin = lmfitter(x_data , y_data, y_error)
	p_lin, chi_lin, nfree_lin = pStat(out_lin)
	info_1 = x_data , model(out_lin.params, x_data), [out_lin.params['m'].value, out_lin.params['m'].stderr, out_lin.params['c'].value, out_lin.params['c'].stderr], [p_lin, chi_lin, nfree_lin]

	#test if less than tolerence
	if p_lin < tolerence:
		fit = 'n=1'
		return info_1

	elif p_lin >= tolerence:
		fit = 'n=2'
		out_poly2 = lmfitter2(x_data , y_data, y_error)
		p_poly2, chi_poly2, nfree_poly2 = pStat(out_poly2)
		info_2 = x_data , model2(out_poly2.params, x_data), [-100., -100.], [p_poly2, chi_poly2, nfree_poly2]
		if p_poly2 < tolerence and nfree_poly2 > 0:
			return info_2

	return info_1

	#now test if p_poly2 is greater than p_lin

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


	
#create figure and subplots
for j in range(len(ycols)):
	#select col
	col = ycols[j]

	
	print col

	#loop throught density tracers
	for i in range(len(xcols)):
		density = xcols[i]
		

		print density

		#create density bins
		#bins_density = np.linspace(np.floor(np.min(detected[density])), np.ceil(np.max(detected[density])), N_den+1) #automatic
		#bins_density = np.linspace(xlims[j][0], xlims[j][1], N_den+1)  #manual

		#loop through cats
		for ll in range(2):
			current_cat = cats[ll]
			#loop through morphological types
			for k in range(len(types)):
				sub = subs[0,i]
				selection = current_cat.where(current_cat[types[k]] == 1) #select morphological type
				bins_density = np.linspace(np.floor(np.min(selection[density])), np.ceil(np.max(selection[density])), N_den+1) #automatic

				#extract raw data
				x_raw, y_raw = selection[density], selection[col]

				#find mean and errors for each point
				x, y, y_err, y_n = mean_binup(x_raw, y_raw, bins_density)

				#clean empty bins
				w_clean = np.where((y_n > 3.0) & (y_err != 0.0))[0]

				x, y, y_err, y_n = x[w_clean], y[w_clean], y_err[w_clean], y_n[w_clean]

				#plot differntly depending on deteccted / undetected , ll = 0, 1 respectively
				if ll == 0:
					sub.errorbar(x, y, yerr=y_err, c=colours[k], ls='', alpha=0.35)
					sub.scatter(x, y, marker='*', s=80, c=colours[k], edgecolor=colours[k], alpha=0.35)
				else:
					sub.errorbar(x, y, yerr=y_err, c=colours[k], ls='', alpha=0.35)
					sub.scatter(x, y, marker='x', s=80, c=colours[k], alpha=0.35)

				#try fitting two models
				if len(x) > 3:

					#fitting function
					x_line, y_line, coeffs, chiStuff = lineFitter(x, y, y_err)

					#append values to dics
					sub.plot(x_line, y_line, c=colours[k], ls=lisy[ll], lw=1, alpha=1.)

				else:
					coeffs = [-100.0, -100.0]

				coeffs_dict[density+types[k]+col+str(ll)] = coeffs #dic to hold grads
				
				
				if i!= 3:
					sub.set_xlim(-1.4, 4.1)
					#subs[0,i].invert_xaxis()

				if i == 3:
					sub.set_xlim(-1.7,1.5)

				if i == 3: 
					sub.axvline(x=np.log10(1.), c='k', ls='--')



				#plotting options
				sub.xaxis.set_major_locator(MaxNLocator(4))
				sub.yaxis.set_major_locator(MaxNLocator(4))

				#lable x on final y
				if k == 2:
					sub.set_xlabel(xlabs[i])

				#label y on final x
				if i == 0:
					sub.set_ylabel(ylabs[0])

				#if not on bottom switch off x axis
				#if k != 2:
					#sub.tick_params(axis = 'x', labelbottom = 'off')

				#if not on the left switch off y axis
				if i != 0:
					sub.tick_params(axis = 'y', labelleft = 'off')

printTable(coeffs_dict)

sub.set_ylim(8.8, 11.2)
##### plotting options
plt.subplots_adjust(left=0.11, bottom=0.17, right=0.97, top=0.97, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/stars-density-v2.pdf')
plt.show()

printTable(coeffs_dict)

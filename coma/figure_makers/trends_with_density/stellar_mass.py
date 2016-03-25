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

#convert radius vir to log10
cat.RADIUS_VIR = np.log10(cat.RADIUS_VIR)

#select detected galaxies
#cat = cat.where(np.nan_to_num(cat.K) != 0.0)
detected = cat#.where(cat.DMASS_TYPE != 0)


#switchs
N_den = 7 #number of density bins
N = 60
bin_type = 'fixed'
types = ['early', 'inter', 'late']
colours = ['r', 'g', 'b']


mms = {} #dic to hold grads
mes = {} # dic to hold grad errrors

#quanties on x axcis
xcols  = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']
xlabs = ['$\log_{10}(\Sigma 1)($Mpc$^{-2})$','$\log_{10}(\Sigma 5)($Mpc$^{-2})$','$\log_{10}(\Sigma 10) ($Mpc$^{-2})$', '$\log_{10}$(R/R$_{Virial}$)']
xlims = [[-1.2, 4.1], [-1.2,3.3], [-1.2,3.3], [-1.7, 1.5]] 



if True:
	#quanties on the yaxis
	ycols = ['SMASS', 'SMASS_BELL']
	ylabs = ['$\log_{10} (M_{stars} / $M$_{\odot}$)', '$\log_{10} (M_{stars} / $M$_{\odot}$)'  ]
	ylims = [[-4.8, -1.4 ]]
	texts = ['Brinchmann et al. (2004)', 'Bell et al. (2003)']
	

#create figure and subplots
fig, subs = plt.subplots(nrows=len(ycols), ncols=len(xcols), sharex=False, sharey=False, squeeze=False, figsize = (8., 4.), facecolor='w',edgecolor='w')

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
	final = [[],[],[]]
	#loop through types
	for ii in range(len(types)):
		#loop through des
		for jj in range(len(xcols)):
			key = xcols[jj]+types[ii]

			#extract values from dics
			m, dm = vals[key], err[key]

			value = str(np.round(m, decimals=2)) + '$\\pm$' + str(np.round(dm, decimals=2))

 			final[ii].append(value)

 	for kk in range(3):
 		line = []
 		for pp in range(len(xcols)):
 			if pp != len(xcols) - 1: line.append(final[kk][pp]+'&')
 			else: line.append(final[kk][pp] + '\\\\')

 		print types[kk]
 		print " ".join(line)

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
		sub = subs[j,i]

		print density

		#create density bins
		#bins_density = np.linspace(np.floor(np.min(detected[density])), np.ceil(np.max(detected[density])), N_den+1) #automatic
		#bins_density = np.linspace(xlims[j][0], xlims[j][1], N_den+1)  #manual

		#loop through morphological types
		for k in range(len(types)):
			selection = detected.where(detected[types[k]] == 1) #select morphological type
			bins_density = np.linspace(np.floor(np.min(selection[density])), np.ceil(np.max(selection[density])), N_den+1) #automatic

			#extract raw data
			x_raw, y_raw = selection[density], selection[col]

			#find mean and errors for each point
			x, y, y_err, y_n = mean_binup(x_raw, y_raw, bins_density)

			#clean empty bins
			w_clean = y_n > 3.0

			x, y, y_err = x[w_clean], y[w_clean], y_err[w_clean] 

			sub.errorbar(x, y, yerr=y_err, c=colours[k], ls='')
			sub.scatter(x, y, marker='o', s=20, c=colours[k], edgecolor='None')

			#try fitting two models
			if len(y) > 3:
				x_line, y_line, coeffs, chiStuff = lineFitter(x, y, y_err)
				sub.plot(x_line, y_line, ls='--', color=colours[k])

			else:
				coeffs = [-999.0,-999.0]

			#append values to dics
			mms[density+types[k]] = coeffs[0] #dic to hold grads
			mes[density+types[k]] = coeffs[1] #dic to hold errors
			"""
			#linear
			x_line1, y_line1, coeffs1, redchi1 = fit_line1(x , y, y_err)

			#polynomial n=2
			try: x_line2, y_line2, coeffs2, redchi2 = fit_line2(x , y, y_err)

			except: redchi2 == 100000.0

			#test if redchi2 is less than redchi1
			if redchi2 < redchi1:
				x_line, y_line, coeffs, redchi = x_line2, y_line2, coeffs2, redchi2
				sub.plot(x_line1, y_line1, c=colours[k], ls='--', alpha=1.)

			if redchi2 > redchi1:
				x_line, y_line, coeffs, redchi = x_line1, y_line1, coeffs1, redchi1

			#if redchi2 > 2.0 or redchi2 < 0.7 :
				#x_line, y_line, coeffs, redchi = x_line1, y_line1, coeffs1, redchi1


			"""
			if i!= 3:
				subs[j,i].set_xlim(-1.4, 4.1)
				#subs[0,i].invert_xaxis()

			if i == 3:
				subs[j,i].set_xlim(-2.7,1.5)


			

			#sub.plot(x_line, y_line, c=colours[k], ls='-', lw=2)


			if True:
				print types[k]
				print x
				print y
				print y_err
				print y_n
				print w_clean
				
		if i == 3: sub.axvline(x=np.log10(1.), c='k', ls='--')



		#plotting options
		sub.xaxis.set_major_locator(MaxNLocator(4))
		sub.yaxis.set_major_locator(MaxNLocator(4))

		#lable x on final y
		if j == len(ycols)-1:
			sub.set_xlabel(xlabs[i])

		#label y on final x
		if i == 0:
			sub.set_ylabel(ylabs[j])

		#if not on bottom switch off x axis
		if j != len(ycols)-1:
			sub.tick_params(axis = 'x', labelbottom = 'off')

		#if not on the left switch off y axis
		if i != 0:
			sub.tick_params(axis = 'y', labelleft = 'off')


		font = {'family' : 'serif', 'weight' : 'normal', 'size'   : 8} 

		if i == 0: sub.text(0.05, 0.95, texts[j], transform=sub.transAxes, fontdict=font, verticalalignment='top')

	print texts[j]
	print_grads(mms, mes)

##### plotting options
plt.subplots_adjust(left=0.11, bottom=0.15, right=0.97, top=0.97, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/stars-density.pdf')
plt.show()

if False:
	print_grads(mms, mes)
	print_test_K()

	sdss_bands = ['u_1', 'g', 'r', 'i', 'z']

	for band in sdss_bands:
		print_test_sdss(band)

	print_test_sdss2()

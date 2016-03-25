#program to create mass and parameter fuctions
# Chris Fuller, June 2014

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
virgo = Table('/Users/chrisfuller/Dropbox/phd/herchel/virgo/virgo-all-data-v2.fits')
fornax = Table('/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/stellar-mass-fornax_final.fits')
#cat = cat.where(cat.F250 >= 50.0E-3)

D_coma = 100.0*3.0857E22



#volumes
volumes = [114.2, 2872.8] #volume Mpc using spere and cylender 

#cats  
cats = [cat.where(cat.RADIUS_VIR <= 1.0), cat.where(cat.RADIUS_VIR > 1.0)]
inputs = ['cat', 'cat', 'line_only', 'line_only']
catNames = ['Coma Cluster', 'Coma Filament', 'Field', 'Virgo Cluster (Fit Davies et al 2013)']

#plotting parameters
line_syles = ['-', '-', '-', '-', '-']
colours = ['r', 'b', 'Cyan', 'Orange']

N = 10 #number of y bins

#columns to plot
plots = [	['SRF'] ]

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
			'L250'			: '$\log_{10}$($L$) (W Hz$^{-1}$)'}

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
			'L250'			: 'FIR'}

select = {	'DMASS' 		:  0.0,
			'DTEMP'			:  0.0, 
			'SMASS' 		:  0.0,
			'DUST_STARS' 	: -7.0,
			'SRF' 			: -100.,
			'sSFR'			: -100.,
			'METAL'			: -100.,
			'colour'		: -100.,	
			'my_morph'		: -100.	}

#morphological types
types = ['late', 'inter', 'early']

types_lab = ['late', 'intermediate', 'early']


#set i and j range
i_range = plots.shape[1]
j_range = plots.shape[0]

#create figure and subplots
fig, subs = plt.subplots(nrows=j_range, ncols=i_range, 
						sharex=False, sharey=False, squeeze=False, 
						figsize = (8.0, 4.5), 
						facecolor='w',edgecolor='w')


# # # # Functions # # # #
#function to plot function
def functionPlotter(x, volume, bins, lineStyle, colour, subplot):
	#create histogram of data
	hist, bins = np.histogram(x, bins=bins)

	
	#cacuate function
	y  = hist/volume 

	#width of bin
	w = bins[1] - bins[0]

	#find x
	x = bins[:-1] + w*0.5

	#clean data
	clean = np.where(hist != 0)


	x = x[clean]
	y = y[clean]
	hist = hist[clean]

	y_err = np.sqrt(hist)/volume

	x = 10**x	

	#fit results
	x_line , y_line, coeffs, fit = fitSchechter(x , y, y_err)
	print x
	print y
	print y_err


	y_errorUp = abs(np.log10(y) - np.log10(y+y_err))
	y_errorDown = abs(np.log10(y) - np.log10(y-y_err))

	y_errors = np.append(y_errorDown.reshape(1,len(y_err)), y_errorUp.reshape(1,len(y_err)), axis=0)

	#plot schechter line
	subplot.plot(np.log10(x_line), np.log10(y_line), ls='-', c=colour)

	#plot results
	subplot.errorbar(np.log10(x), np.log10(y), yerr=y_errors, c=colour, ls='')
	#subplot.plot(np.log10(x), np.log10(y), ls=lineStyle, c=colour)
	subplot.scatter(np.log10(x), np.log10(y), marker='x', s=20, color=colour)

	printFit(fit.params)

#function to create bins
def binCreator(ts, para, N, meth):
	min_val = 0.0
	max_val = 0.0
	for i in range(len(ts)):
		t = ts[i]
		if len(t)<5.0: continue
		print meth
		if meth == 'FIR': t = t.where(t.D250 ==1)
		if i == 0: 
			min_val, max_val = np.min(t[para]) , np.max(t[para])
		else:
			min_temp, max_temp = np.min(t[para]) , np.max(t[para])
			if min_val > min_temp:
				min_val = min_temp
			if max_val < max_temp:
				max_val = max_temp

	return np.linspace(min_val,max_val, N+1)

def residual(params, x, y_data, y_error):
	return (y_data-model(params, x))/y_error

def lmfitter(x , y, y_error):
	params = Parameters()
	params.add('phiStar', value=0.3, vary=True)
	params.add('LStar', value=1., vary=True)
	params.add('alpha', value=-1.0, vary=True)

	# remove inf values in errors
	out = minimize(residual, params, args=(x, y, y_error))
	#report_fit(params)
	return out

def model(params, x):
	phiStar = params['phiStar'].value
	LStar = params['LStar'].value
	alpha = params['alpha'].value + 1

	L = (x/LStar)
	phi = phiStar*(L**alpha)*np.exp(-1.0*L)
	return phi

def myPrint(text,a,b):
	try:
		print text + ' ', np.round(a,decimals=3), '$\\pm$', np.round(b,decimals=3)

	except:
		print 'error ', text
def printFit(params):
	phiStar = params['phiStar'].value
	LStar = params['LStar'].value
	alpha = params['alpha'].value

	phiStarErr = params['phiStar'].stderr
	LStarErr = params['LStar'].stderr
	alphaErr = params['alpha'].stderr


	#L = (x/LStar)
	#phi = phiStar*(L**alpha)*np.exp(-1.0*L)

	try:
		myPrint( 'alpha:', alpha, alphaErr)
		myPrint( 'phiStar:', phiStar, phiStarErr)
		myPrint( 'LStar:', LStar, LStarErr)
	except:
		print 'error '
def fitSchechter(x_data , y_data, y_error):


	#fit the data
	out = lmfitter(x_data , y_data, y_error)
	params = out.params

	phiStar,  phiStarErr= params['phiStar'].value, params['phiStar'].stderr
	LStar, LStarErr = params['LStar'].value, params['LStar'].stderr
	alpha, alphaErr = params['alpha'].value, params['alpha'].stderr

	coeffs = [phiStar, phiStarErr, LStar, LStarErr, alpha, alphaErr]

	x_mod = np.linspace(x_data[0], x_data[-1], 200)
	
	return x_mod , model(out.params, x_mod), coeffs, out

def intTick(sub):
	import matplotlib.ticker as ticker

	start, end = sub.get_xlim()
	#sub.set_xlim(np.int(np.floor(start)), np.int(np.ceil(end)))
	#start, end = sub.get_xlim()
	sub.xaxis.set_ticks(np.arange(start, end, 1))
	#sub.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

	start, end = sub.get_ylim()
	sub.set_ylim(np.int(-4), np.int(np.ceil(1)))
	start, end = sub.get_ylim()
	sub.yaxis.set_ticks(np.arange(start, end, 1))

#loop throught plots
for j in range(j_range):
	for i in range(i_range):
		sub = subs[j,i]
		para = plots[j,i]



		if method[para] == 'FIR': bin_cat = cat.where(cat.D250 ==1)
		else: 
			bin_cat = cat 

		#create bins
		bins = binCreator(cats, para, N, method[para])
		#bins = np.linspace(min(bin_cat[para]),max(bin_cat[para]), N+1)


		#loop through cats
		for k in range(len(cats)):
			selection = cats[k]
			print ''
			print catNames[k]
			#check if line only or cat
			if inputs[k] == 'cat':
				if method[para] == 'FIR': selection = selection.where(selection.D250 ==1)
				#if method[para] == 'SF': selection = selection.where(selection.bptclass ==1)
				#for l in range(len(types)):
			#		subCat = selection.where(selection[types[l]]==1)
		#			functionPlotter(subCat[para], volumes[k], bins, line_syles[k], colours[l], sub)

				functionPlotter(selection[para], volumes[k], bins, line_syles[k], colours[k], sub)

			

		#ploting options
		sub.set_ylabel('$\log_{10}$(N)(Mpc$^{-3}$)')
		sub.set_xlabel(labels[para])
		#intTick(sub)


plt.subplots_adjust(left=0.11, bottom=0.13, right=0.97, top=0.97, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/sfr_Fucntion.pdf')
plt.show()



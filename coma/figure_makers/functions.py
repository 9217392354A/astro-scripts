#program to create functions


#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MaxNLocator
import scipy
#remove numpy runtime warings
np.seterr(invalid='ignore')


#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))
cat = cat.where(cat.SMASS > 6.8)


#volumes
volumes = [114.2, 2872.8] #volume Mpc using spere and cylender 
volumesGas = [114.2*0.5, 2872.8*0.3] #volume Mpc using spere and cylender 
#volumes  = [114.2, 38782.8]
#volumes= [1.0, 1.0]

#cats
cats = [cat.where(cat.RADIUS_VIR <= 1.0), cat.where(cat.RADIUS_VIR > 1.0)]
cats_names = ['Cluster&', 'Filament&']

#plotting parameters
line_syles = ['-', '--']
lines = [9.1, 9.3, 6.5]
limits = [[8,12], [8, 11], [6,9]]
Ns = [5,7,5]


#columns to plot
plots = [	['SMASS', 'HI_ALL2', 'DMASS']]#, 'SRF']] 

plots = np.array(plots)

guesses = [		[1.751, 4.0E8, -1.0],
				[0.6, 4.0E9,  -1.0],
				[0.6, 1.9E11,-1.],
				[0.6, 0.2,  -1.0]	]

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
def residual(params, x, y_data, y_error):
	return (y_data-model(params, x))/y_error

def lmfitter(x , y, y_error, guess):

	params = Parameters()
	params.add('phiStar', value=guess[0], vary=True)
	params.add('LStar', value=guess[1], vary=True)
	params.add('alpha', value=guess[2], vary=True)


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

def s(a,b):
	return s1(a) + '$\\pm$' + s1(b)

def s1(s):
	if s == 100.0 or s == np.float('nan'):
		return '-'
	else:
		if s > 0.001 and s < 0.001:
			string = ("{:.1e}".format(s)).split('e')
			return string[0] + '$10^{' + string[1]+ '}$' 
		else:
			return str(np.round(s,decimals=3))


def fitSchechter(x_data , y_data, y_error, guess,bins):
	#fit the data
	out = lmfitter(x_data , y_data, y_error, guess)
	params = out.params

	phiStar,  phiStarErr= params['phiStar'].value, params['phiStar'].stderr
	LStar, LStarErr = params['LStar'].value, params['LStar'].stderr
	alpha, alphaErr = params['alpha'].value, params['alpha'].stderr

	gamma2 = scipy.special.gamma(alpha+2)
	gamma1 = scipy.special.gamma(alpha+1)

	rho = phiStar*LStar*gamma2
	NN = phiStar*LStar*gamma1

	coeffs = [phiStar, LStar, alpha, rho, NN]
	sCoeffs = s(phiStar,phiStarErr)+ '&' + s(LStar/10.0**9,LStarErr/10.0**9) + '&' + s(alpha,alphaErr) + '&' + s1(rho/10.0**9) 

	x_mod = np.logspace(bins[0]-4, bins[-1]+4, 200)
	return x_mod , model(out.params, x_mod), sCoeffs, coeffs, out


#function to plot function
def functionPlotter(x, volume, bins, lineStyle, colour, subplot, guess, para):
	#create histogram of data
	hist, bins = np.histogram(x, bins=bins)

	#caculate errors
	hist_errors = np.sqrt(hist) 

	#cacuate function
	y  = hist /volume 
	y_err = hist_errors /volume

	#width of bin
	w = bins[1] - bins[0]

	#find x
	x = bins[:-1] + w*0.5

	#clean data
	clean = np.where(hist > 1)


	x = x[clean]
	y = y[clean]
	y_err = y_err[clean]



	if para == 'HI_ALL2' or para == 'SMASS':
		x_fit = x[1:]
		y_fit = y[1:]
		y_err_fit = y_err[1:]


		x_fitNOFIT = x[0]
		y_fitNOFIT = y[0]

		#print x_fitNOFIT
		#print x_fit

		subplot.scatter(x_fitNOFIT, np.log10(y_fitNOFIT), marker='s', s=40, facecolor='None', edgecolor = 'k')

	else:
		x_fit = x
		y_fit = y
		y_err_fit = y_err
	#fit results
	if len(y) > 3:
		x_line , y_line, scoeffs, coeffs, fit = fitSchechter(10**x_fit , y_fit, y_err_fit, guess,bins)

	y_log = np.log10(y)

	y_errorUp = abs(np.log10(y) - np.log10(y+y_err))
	y_errorDown = abs(np.log10(y) - np.log10(y-y_err))

	y_errors = np.append(y_errorDown.reshape(1,len(y_err)), y_errorUp.reshape(1,len(y_err)), axis=0)

	#plot results
	subplot.errorbar(x, y_log, yerr=y_errors,  c=colour, ls='')
	#subplot.plot(x, y_log, ls=lineStyle, c=colour)
	subplot.scatter(x, y_log, marker='o', s=10, facecolor=colour, edgecolor = colour)
	#plot schechter line
	if len(y) > 3:
		subplot.plot(np.log10(x_line), np.log10(y_line), ls='-', c=colour)

	#print'N: ', np.log10(coeffs[-2]*volume) 
	#print'density_N', coeffs[-2]
	return scoeffs, coeffs

#loop throught plots
for j in range(j_range):
	for i in range(i_range):
		sub = subs[j,i]
		para = plots[j,i]
		print para
		N = Ns[i]

		if method[para] == 'FIR': bin_cat = cat.where((cat.D250 ==1) & (cat.DMASS > 6.4))
		elif method[para] == 'gas': bin_cat = cat.where(np.nan_to_num(cat.HI_ALL2) > 8.0)
		elif para == 'METAL' : cat.where(np.nan_to_num(cat.METAL) > 0.0)
		#elif para == 'SRF': bin_cat = cat.where(cat.bptclass == 1)
		else: bin_cat = cat 

		#create bins
		bins = np.linspace(min(bin_cat[para]),max(bin_cat[para]), N+1)

		
		#loop through cats
		for k in range(len(cats)):
			selection = cats[k]
			if method[para] == 'FIR': selection = selection.where( (selection.D250 ==1) & (selection.DMASS > 6.4))
			elif para == 'HI_ALL2': selection = selection.where(np.nan_to_num(selection.HI_ALL2) > 5.0)
			elif para == 'SMASS': selection = selection.where(selection[para] > 5.0 )
			print len(selection)
			#elif para == 'SRF': selection = selection.where(selection.bptclass == 1)
			#for l in range(len(types)):
			#	subCat = selection.where(selection[types[l]]==1)
			#	functionPlotter(subCat[para], volumes[k], bins, line_syles[k], colours[l], sub, guesses[i])


			try:
				scoeffs, coeffs = functionPlotter(selection[para], volumes[k], bins, line_syles[k], colours_sam[k], sub, guesses[i], para)
				if k == 0: A = coeffs
				else: B = coeffs

				print cats_names[k] + scoeffs + '\\\\'

				C = np.array(A)/np.array(B)

			except:
				continue

		#print(C)
		print '\\\\'

		sub.yaxis.set_major_locator(MaxNLocator(5))
		sub.xaxis.set_major_locator(MaxNLocator(5))
		#ploting options
		sub.set_xlabel(labels[para])
		#sub.axvline(x=lines[i], ls='--', c='k')
		sub.set_xlim(limits[i][0], limits[i][1])


subs[0,0].set_ylabel('$\log_{10}$(N)(Mpc$^{-3}$)')
sub.set_ylim(-4, 1)

plt.subplots_adjust(left=0.1, bottom=0.15, right=0.97, top=0.97, wspace=0.15, hspace=0.3)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/mass_function.pdf')
plt.show()



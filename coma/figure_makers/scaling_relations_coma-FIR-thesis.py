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
np.seterr(all='ignore')
#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(pj(folder,fname))


#caculated extra columns
cat.add_column('D2G', - cat.HI_ALL2 +cat.DMASS)
cat.add_column('G2S', cat.HI_ALL2 - cat.SMASS)
cat.add_column('SFR2G', cat.SRF - cat.HI_ALL2)
cat.add_column('SFR2D', cat.SRF - cat.DMASS)
cat.add_column('colour', cat.g - cat.r)


#select currentCat galaxies
firCat = cat.where(cat.DMASS_TYPE != 0)
gasCat = cat.where(np.nan_to_num(cat.HI_ALL2) > 0.0)
metalCat = cat.where(np.nan_to_num(cat.METAL) > 0.0)
gasfirCat = cat.where((cat.DMASS_TYPE != 0) & (np.nan_to_num(cat.HI_ALL2) > 0.0))
sfCat = cat.where(cat.bptclass == 1.0)
sffirCat = cat.where((cat.DMASS_TYPE != 0) & (cat.bptclass == 1.0))

firsfmetalCat = cat.where((cat.DMASS_TYPE != 0) & (cat.bptclass == 1.0) & (np.nan_to_num(cat.METAL) > 0.0))

#find min dust
min_dust = min(firCat.DMASS)

N = 6 #number of y bins

#columns to plot
plots = [	[['DMASS','SMASS'],['DMASS','METAL'] ],  
			[['DUST_STARS','SMASS'],['', '']]]

plots = np.array(plots)



cats = [[firCat, firsfmetalCat], 
			[firCat, gasfirCat] ]


cats = np.array(cats)


#columns to plot
letters = [	['', ''], 
			['', ''],
			['e)', 'f)'],
			['g)', 'h)'],
			['i)', 'j)']]

letters = np.array(letters)

#create dict for lables
labels = {	'DMASS' 		: '$\log_{10} (M_{dust} / $M$_{\odot}$)', 
			'SMASS' 		: '$\log_{10} (M_{stars} / $M$_{\odot}$)',
			'DUST_STARS' 	: '$\log_{10} (M_{dust} / M_{star}$)',
			'SRF' 			: '$\log_{10}$(SFR) (M$_{\odot}$ yr$^{-1}$)',
			'sSFR'			: '$\log_{10}$(sSFR) (yr$^{-1}$)',
			'METAL'			: '12 + $\log_{10}$(O/H)',
			'colour'		: 'g-r',
			'SFR2D'			: '$\log_{10} ( SFR/ M_{dust})$ (yr$^{-1}$)',
			'HI_ALL2'		: '$\log_{10} (M_{gas} / $M$_{\odot}$)',
			'G2S'			: '$\log_{10} (M_{gas} / M_{star}$)'}

labels2 = {	'DMASS' 		: '$M_{dust}$', 
			'SMASS' 		: '$M_{stars}$',
			'DUST_STARS' 	: '$M_{dust} / M_{star}$',
			'SRF' 			: 'SFR',
			'sSFR'			: 'sSFR',
			'METAL'			: 'Z',
			'colour'		: 'g-r',
			'SFR2D'			:  'SFR / $M_{dust}$',
			'HI_ALL2' 		: '$M_{gas}$',
			'G2S'			: '$M_{gas} / M_{star}$'}



#morphological types
#types = ['late', 'inter', 'early']
#colours = ['b', 'g', 'r']

types = ['late','early']
colours = ['b', 'r']

types_lab = ['late', 'intermediate', 'early']


coeffs_dict = {} #dic to hold coeffs

# extract i's and j's
j_range = plots.shape[0]
i_range = plots.shape[1]

#create a figure
#create figure and subplots
fig, subs = plt.subplots(nrows=j_range, ncols=i_range,
						figsize = (8., 8.0), facecolor='w',edgecolor='w')



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

def residual(params, x, y_data):
	y_model = model(params, x)
	return y_data-y_model

def residual2(params, x, y_data, y_error):
	y_model = model(params, x)
	return (y_data-y_model) / y_error 

def lmfitter(x , y):
	intial_fit = np.polyfit(x, y, 1)
	params = Parameters()
	params.add('m', value=intial_fit[0], vary=True)
	params.add('c', value=intial_fit[1], vary=True)


	# remove inf values in errors
	out = minimize(residual, params, args=(x, y))
	#report_fit(params)
	return out

def lmfitter2(x , y, yerr):
	intial_fit = np.polyfit(x, y, 1)
	params = Parameters()
	params.add('m', value=intial_fit[0], vary=True)
	params.add('c', value=intial_fit[1], vary=True)


	# remove inf values in errors
	out = minimize(residual2, params, args=(x, y, yerr))
	#report_fit(params)
	return out

def model(params, x):
	m = params['m'].value
	c = params['c'].value

	return m*x + c

def fit_line1(x_data , y_data):

	#fit the data
	out = lmfitter(x_data , y_data)
	m, dm  = out.params['m'].value, out.params['m'].stderr
	c, dc  = out.params['c'].value, out.params['c'].stderr
	#test correlations
	correlation = pearsonr(x_data, y_data)

	return x_data , model(out.params, x_data), [m, dm, c, dc, correlation[0], correlation[1]]

def fit_line2(x_data , y_data, y_error, x, y):

	#fit the data
	out = lmfitter2(x_data , y_data, y_error)
	m, dm  = out.params['m'].value, out.params['m'].stderr
	c, dc  = out.params['c'].value, out.params['c'].stderr
	#test correlations
	correlation = pearsonr(x, y)
	spcc = sumPCC(correlation[0], correlation[1])

	return x_data , model(out.params, x_data), [m, dm, c, dc, correlation[0], correlation[1]], spcc

def sumPCC(PCC, pVal):
	#print '------------------ -------------------', PCC, pVal
	if pVal == np.float('nan') or pVal == 'nan':
		return 'None'

	if pVal > 0.05:
		return 'None'

	if PCC >= 0.75:
		return 'Strong'

	if PCC < 0.75 and PCC >= 0.5:
		return 'Moderate'

	if PCC < 0.5 and PCC >= 0.25:
		return 'Weak'

	if PCC < 0.25:
		return 'None'





def printShortTable(dict):

	#loop through i's and j's
	for j in range(j_range):
		for i in range(i_range):
			#loop through each morphological type
			x_col = plots[j,i,1]
			y_col = plots[j,i,0]
			if plots[j,i,0] == '': continue 
			print '\\\\'
			text = labels2[x_col] + ' vs. ' + labels2[y_col]
			multiNotation(text,10)
			
			for k in range(len(types)):
				#extract x and y column headers
				line = ''

				#extract coeeffs
				coeffs_0 = coeffs_dict[types[k]+x_col+y_col]
				m0, dm0, c0, dc0, pcc0, pstat0 = coeffs_0[0], coeffs_0[1], coeffs_0[2], coeffs_0[3], coeffs_0[4], coeffs_0[5]


				#test the pvalue
				sPCC0 = sumPCC(abs(pcc0), pstat0)


				#print diff_m, diff_c, diff


				#print sPCC0, sPCC1, diff

				try:
					line = types_lab[k] + '&' + sPCC0 +  '\\\\'
				except:
					line = types_lab[k] + '&' + '-' + '&' + '-'  + '&' +  '\\\\'
				print line


def printShortTable2(dict):

	#loop through i's and j's
	for j in range(j_range):
		for i in range(i_range):
			#loop through each morphological type
			x_col = plots[j,i,1]
			y_col = plots[j,i,0]
			if plots[j,i,0] == '': continue 
			text = labels2[x_col] + ' vs. ' + labels2[y_col]
			#multiNotation(text,10)

			lineFull = text + '&' 
			line = ''
			
			for k in range(len(types)):
				#extract x and y column headers
				line = ''

				#extract coeeffs
				coeffs_0 = coeffs_dict[types[k]+x_col+y_col]
				m0, dm0, c0, dc0, pcc0, pstat0 = coeffs_0[0], coeffs_0[1], coeffs_0[2], coeffs_0[3], coeffs_0[4], coeffs_0[5]


				#test the pvalue
				sPCC0 = sumPCC(abs(pcc0), pstat0)


				#print diff_m, diff_c, diff
				psformat(pcc0, pstat0)

				#print sPCC0, sPCC1, diff

				#try:
				line = psformat(pcc0, pstat0) #+  '\\\\'
				#except:
				#	line = types_lab[k] + '&' + '-' + '&' + '-'  + '&' +  #'\\\\'


				lineFull += line 
				if k!=2: lineFull += '&'
				else: lineFull += '\\\\'  
			print lineFull 

def psformat(a,b):
	A = str(np.round(a, decimals=2))

	if b < 0.01:
		string = ("{:.1e}".format(b)).split('e')
		B =  '$10^{' + string[1]+ '}$' 
	else:
		B = str(np.round(b,decimals=2))

	return A + '(' + B + ')'



def s1(s):
	if s == 100.0 or s == np.float('nan'):
		return '-'
	else:
		return str(np.round(s,decimals=2))

def s2(s):
	if s == 100.0 or s == np.float('nan'):
		return '-'
	else:
		if s < 0.001:
			string = ("{:.1e}".format(s)).split('e')
			return string[0] + '$10^{' + string[1]+ '}$' 
		else:
			return str(np.round(s,decimals=3))


def grad(m,dm):
	if m == 100.0: 
		return '-'
	else:
		return s1(m) + '$\pm$' + s1(dm)
def multiNotation(text,colN):
	print '\\multicolumn{'+ str(colN) +'}{l}{'+ text +'}' + '\\\\'

def clear_frame(ax=None):
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.itervalues():
        spine.set_visible(False)


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

#loop through i's and j's
for j in range(j_range):
	for i in range(i_range):
		#extract subplot
		try: 
			sub = subs[j,i]
		except:
			sub = subs[i]

		letter = letters[j,i]
		cat = cats[j,i]

		if plots[j,i,0] == '':
			clear_frame(sub)
			continue 

		sub.text(0.03, 0.97, letter, transform=sub.transAxes, fontsize=12, verticalalignment='top')
		#loop through each morphological type
		for k in range(len(types)):
			#select morphological type
			subCat = cat.where(cat[types[k]]==1)

			#extract x and y column headers
			x_col = plots[j,i,1]
			y_col = plots[j,i,0]


			#selected galaxies that are no null values
			#subCat = subCat.where((subCat[x_col] > -100.0) & (subCat[y_col]> -100.0))

			#if cluster vs filament desired set true
			if True:
				cluster_subCat = subCat.where(subCat.RADIUS_VIR <= 1.0)
				filament_subCat = subCat.where(subCat.RADIUS_VIR > 1.0)

				#sCats = [subCat]

				

				

				#extract data
				x, y = subCat[x_col], subCat[y_col]

				lines = []

				for N in range(1,8):

					bins= np.linspace(np.floor(np.min(subCat[x_col])), np.ceil(np.max(subCat[x_col])), N+1) #automatic
					x_mean, y_mean, y_error, y_n = mean_binup(x, y, bins)

					#clean empty bins
					w_clean = np.where((y_n > 3.0) & (y_error != 0.0))[0]
					x_mean, y_mean, y_error = x_mean[w_clean], y_mean[w_clean], y_error[w_clean]
					
					#fit data
					if len(y_mean) > 2:
						x_line , y_line, coeffs, spcc = fit_line2(x_mean, y_mean, y_error, x, y)
						line = 1
					else:
						correlation = pearsonr(x, y)
						coeffs = [100.0,100.0,100.0,100.0,correlation[0],correlation[1]]
						line = 0

					
					if False:
						try: 
							x_line , y_line, coeffs, spcc = fit_line1(x, y)
							line = 1
						except:
							line = 0
							coeffs = [100.0,100.0,100.0,100.0,100.0,100.0]

					lines.append(line)
				#select best N
				lines = np.array(lines)
				try:
					N = range(1,8)[np.where(lines==1)[0][-1]]
				except:
					line = 0
					N = 6

				#print 'N: ', N
				#print 'range(1,10) : ', range(1,10)
				#print 'Lines       : ', lines
				#print 'Where Result: ', range(1,10)[N[-1]]
				bins= np.linspace(np.floor(np.min(subCat[x_col])), np.ceil(np.max(subCat[x_col])), N+1) #automatic
				x_mean, y_mean, y_error, y_n = mean_binup(x, y, bins)

				#clean empty bins
				w_clean = np.where((y_n > 3.0) & (y_error != 0.0))[0]
				x_mean, y_mean, y_error = x_mean[w_clean], y_mean[w_clean], y_error[w_clean]
				
				#fit data
				if len(y_mean) > 2:
					x_line , y_line, coeffs, spcc = fit_line2(x_mean, y_mean, y_error, x, y)
					line = 1
				else:
					correlation = pearsonr(x, y)
					coeffs = [100.0,100.0,100.0,100.0,correlation[0],correlation[1]]
					line = 0

				if False:
					try: 
						x_line , y_line, coeffs, spcc = fit_line1(x, y)
						line = 1
					except:
						line = 0
						coeffs = [100.0,100.0,100.0,100.0,100.0,100.0]

				#append values to dics
				coeffs_dict[types[k]+x_col+y_col] = coeffs #dic to hold grads

				#if spcc == 'None': line = 0

				
				#sub.errorbar(x_mean, y_mean, yerr=y_error, c=colours[k], ls='')
				sub.scatter(x, y, marker='+', s=20, facecolor= 'None', edgecolor = colours[k], alpha=0.8)
				#sub.scatter(x_mean, y_mean, marker='*', s=y_n*5.0, facecolor='None' , edgecolor = colours[k])
				if line == 1 and k == 0: sub.plot(x_line, y_line, c=colours[k], ls='-', lw=2)


		#plot lims
		if x_col == 'DMASS':
			sub.axvline(x=min_dust, ls='--', c='k')

		if y_col == 'DMASS':
			sub.axhline(y=min_dust, ls='--', c='k')

		if y_col =='DUST_STARS' and x_col == 'SMASS':
			subsub = cat.where((subCat[x_col] > -100.0) & (subCat[y_col]> -100.0))
			#xxx = subsub[x_col]
			xxx = [8.0,11.0]
			yyy = min_dust - xxx
			sub.plot(xxx, yyy, c='k', ls='--')

		
		sub.yaxis.set_major_locator(MaxNLocator(4))
		sub.xaxis.set_major_locator(MaxNLocator(6))
		sub.set_xlabel(labels[x_col])
		sub.set_ylabel(labels[y_col])
		#sub.yaxis.set_label_coords(-0.14, 0.5)
		#sub.xaxis.set_label_coords(0.5, -0.2)



print ''*10
#printTable(coeffs_dict)
print ''*10

printShortTable(coeffs_dict)
print ''*10


printShortTable2(coeffs_dict)
plt.subplots_adjust(left=0.10, bottom=0.15, right=0.97, top=0.97, wspace=0.2, hspace=0.2)
fig.savefig('/Users/chrisfuller/Dropbox/Scaling_relations-COMA-FIR.pdf')
plt.show()
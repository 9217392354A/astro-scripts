#parameter average caculation for cluster and filament 
# Chris Fuller, 2014

import numpy as np
import atpy
from os.path import join as pj


#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" 
fname = 'coma_supercluster_cal12.fits' 

cat = atpy.Table(pj(folder,fname))
cat.add_column('colour', cat.g - cat.r)
cat.add_column('G2S', cat.HI_ALL2 - cat.SMASS)
cat.add_column('G2D', cat.HI_ALL2 - cat.DMASS)

#create dict for lables
labels = {	'DMASS' 		: '$\log_{10} (M_{dust} / $M$_{\odot}$)', 
			'SMASS' 		: '$\log_{10} (M_{stars} / $M$_{\odot}$)',
			'DUST_STARS' 	: '$\log_{10} (M_{dust} / M_{star}$)',
			'SRF' 			: '$\log_{10}$(SFR) (M$_{\odot}$ yr$^{-1}$)',
			'sSFR'			: '$\log_{10}$(sSFR) (yr$^{-1}$)',
			'METAL'			: '$\log_{10}$(O/H) + 12',
			'colour'		: 'g-r'}

types = ['early', 'inter', 'late']
type_names = ['early', 'intermediate', 'late']



params = [	'DMASS',
			'DTEMP', 
			'SMASS', 
			'SRF', 
			'DUST_STARS', 
			'sSFR',
			'METAL', 
			'colour',
			'my_morph',
			'G2S',
			'G2D',
			'HI_ALL2']




method = {	'DMASS' 		: 'FIR', 
			'SMASS' 		: 'non-FIR',
			'DTEMP'			: 'FIR',
			'DUST_STARS' 	: 'FIR',
			'SRF' 			: 'non-FIR',
			'sSFR'			: 'non-FIR',
			'METAL'			: 'non-FIR',
			'colour'		: 'non-FIR',
			'my_morph'		: 'non-FIR',
			'G2S'			: 'GAS-non-FIR',
			'G2D'			: 'GAS-FIR',
			'HI_ALL2'		: 'GAS-non-FIR'
					}

select = {	'DMASS' 		:  0.0,
			'DTEMP'			:  0.0, 
			'SMASS' 		:  0.0,
			'DUST_STARS' 	: -7.0,
			'SRF' 			: -100.,
			'sSFR'			: -100.,
			'METAL'			: -100.,
			'colour'		: -100.,	
			'my_morph'		: -100.	,
			'G2S'			: -100.,
			'G2D'			: -100.,	
			'HI_ALL2'		: 0.	}
# functions
def header(x):
	print '-'*50
	print '\n'
	print (parameter)

def footer():
	print '\n'
	print '-'*50

def meanError(x):
	return np.std(x)/np.sqrt(np.float(len(x)))

def mean2(a):
	m = np.mean(a)
	e = meanError(a)
	return m, e

def mean3(a,b):
	c = np.append(a, b)
	m = np.mean(c)
	e = meanError(c)
	return s([m,e])

def s(a):
	m = a[0]
	e = a[1]

	m1 = np.round(m, decimals=2)
	e1 = np.round(e, decimals=2)

	return str(m1) + '\\,$\\pm$\\,' + str(e1)

def diff(a,b):
	m1 = a[0]
	e1 = a[1]

	m2 = b[0]
	e2 = b[1]

	eTot = np.sqrt(e1**2 + e2**2)

	difference = np.round(abs(m2-m1)/eTot, decimals=2)
	return str(difference) + '$\\sigma$'

def st(val):
	return str(np.round(val,decimals=2)) 

def printReport(a,b):
	c = np.append(a, b)
	m = mean2(c)
	lower = min(c)
	upper = max(c)

	print 'range: ',st(lower), ' to ', st(upper)
	print 'mean: ', s(m)

def mean4type(c,para):
	a = c[para]
	res = mean2(a)
	return s(res)



for i in range(len(params)):
	parameter = params[i]
	header(parameter)

	#print 'WARNING '*10
	if method[parameter] == 'FIR':
		for k in range(len(types)):
			t = types[k]
			subCat = cat.where((cat[t]==1) & (cat.DMASS_TYPE != 0))
			print 'parameter mean:' + type_names[k] + ' ' + mean4type(subCat, parameter)
			cluster = subCat.where(subCat.RADIUS_VIR <= 1.0)
			filament = subCat.where(subCat.RADIUS_VIR > 1.0)

			
			clusterArray = cluster[parameter]
			filamentArray = filament[parameter]

			clusterMean = mean2(clusterArray)
			filamentMean = mean2(filamentArray)

			line = type_names[k] + ' Cluster: ' + s(clusterMean)
			line += ' filament: ' + s(filamentMean)
			line += '   difference: ' + diff(clusterMean, filamentMean)

			if False:
				print t
				printReport(clusterArray, filamentArray)
			else: 
				print line


	if method[parameter] == 'non-FIR':
		for k in range(len(types)):
			t = types[k]

			#main cat
			subCat = cat.where((cat[t]==1) & (cat[parameter] > select[parameter]))

			du = ['undetected: ','detected: ']

			
			print 'parameter mean:' + type_names[k] + ' ' + mean4type(subCat, parameter)
			for val in range(2):
			
				cluster = subCat.where((subCat.RADIUS_VIR <= 1.0) & (subCat.D250 == val))
				filament = subCat.where((subCat.RADIUS_VIR > 1.0) & (subCat.D250 == val))

				clusterArray = cluster[parameter]
				filamentArray = filament[parameter]

				clusterMean = mean2(clusterArray)
				filamentMean = mean2(filamentArray)

				line = du[val] + type_names[k] + ' Cluster: ' + s(clusterMean)
				line += ' filament: ' + s(filamentMean)
				line += '   difference: ' + diff(clusterMean, filamentMean)

				print line
				line = ''

				print 'overall mean:' + type_names[k] + ' ' +  mean3(filamentArray, clusterArray)
				print '\n'


	if method[parameter] == 'GAS':
		for k in range(len(types)):
			t = types[k]
			subCat = cat.where((cat[t]==1) & (cat.HI_ALL2 != 0.0))
			print 'parameter mean:' + type_names[k] + ' ' + mean4type(subCat, parameter)
			cluster = subCat.where(subCat.RADIUS_VIR <= 1.0)
			filament = subCat.where(subCat.RADIUS_VIR > 1.0)

			
			clusterArray = cluster[parameter]
			filamentArray = filament[parameter]

			clusterMean = mean2(clusterArray)
			filamentMean = mean2(filamentArray)

			line = type_names[k] + ' Cluster: ' + s(clusterMean)
			line += ' filament: ' + s(filamentMean)
			line += '   difference: ' + diff(clusterMean, filamentMean)

			if False:
				print t
				printReport(clusterArray, filamentArray)
			else: 
				print line

	if method[parameter] == 'GAS-FIR':
			for k in range(len(types)):
				t = types[k]
				subCat = cat.where((cat[t]==1) & (cat.HI_ALL2 != 0.0) & (cat.DMASS_TYPE != 0))
				print 'parameter mean:' + type_names[k] + ' ' + mean4type(subCat, parameter)
				cluster = subCat.where(subCat.RADIUS_VIR <= 1.0)
				filament = subCat.where(subCat.RADIUS_VIR > 1.0)

				
				clusterArray = cluster[parameter]
				filamentArray = filament[parameter]

				clusterMean = mean2(clusterArray)
				filamentMean = mean2(filamentArray)

				line = type_names[k] + ' Cluster: ' + s(clusterMean)
				line += ' filament: ' + s(filamentMean)
				line += '   difference: ' + diff(clusterMean, filamentMean)

				if False:
					print t
					printReport(clusterArray, filamentArray)
				else: 
					print line


	if method[parameter] == 'GAS-non-FIR':
			for k in range(len(types)):
				t = types[k]
				subCat = cat.where((cat[t]==1) & (cat.HI_ALL2 != 0.0) & (cat.DMASS_TYPE != 0))
				print 'parameter mean:' + type_names[k] + ' ' + mean4type(subCat, parameter)
				cluster = subCat.where(subCat.RADIUS_VIR <= 1.0)
				filament = subCat.where(subCat.RADIUS_VIR > 1.0)

				
				clusterArray = cluster[parameter]
				filamentArray = filament[parameter]

				clusterMean = mean2(clusterArray)
				filamentMean = mean2(filamentArray)

				line = type_names[k] + ' Cluster: ' + s(clusterMean)
				line += ' filament: ' + s(filamentMean)
				line += '   difference: ' + diff(clusterMean, filamentMean)

				if False:
					print t
					printReport(clusterArray, filamentArray)
				else: 
					print line
	
	if method[parameter] == 'GAS-non-FIR':
		for k in range(len(types)):
			t = types[k]

			#main cat
			subCat = cat.where((cat[t]==1) & (cat.HI_ALL2 != 0.0))
			print 'parameter mean:' + type_names[k] + ' ' + mean4type(subCat, parameter)
			du = ['undetected: ','detected: ']

			for val in range(2):
			
				cluster = subCat.where((subCat.RADIUS_VIR <= 1.0) & (subCat.D250 == val))
				filament = subCat.where((subCat.RADIUS_VIR > 1.0) & (subCat.D250 == val))

				clusterArray = cluster[parameter]
				filamentArray = filament[parameter]

				clusterMean = mean2(clusterArray)
				filamentMean = mean2(filamentArray)

				line = du[val] + type_names[k] + ' Cluster: ' + s(clusterMean)
				line += ' filament: ' + s(filamentMean)
				line += '   difference: ' + diff(clusterMean, filamentMean)

				print line
				line = ''

				print 'overall mean:' + type_names[k] + ' ' +  mean3(filamentArray, clusterArray)
				print '\n'


	footer()







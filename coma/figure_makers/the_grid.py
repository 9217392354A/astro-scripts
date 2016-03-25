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
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

conditions = 	{		'DMASS' 		: cat.DMASS_TYPE != 0, 
						'SMASS' 		: cat.DMASS_TYPE != 10,
						'SRF' 			: cat.bptclass == 1.0,
						'METAL'			: np.nan_to_num(cat.METAL) > 0.0,
						'colour'		: cat.DMASS_TYPE != 10,
						'HI_ALL2'		: np.nan_to_num(cat.HI_ALL2) > 0.0}


compList = ['SMASS', 'DMASS', 'SRF', 'HI_ALL2', 'METAL']


labels = {	'DMASS' 		: '$M_{dust}$', 
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
types = ['late', 'inter', 'early']

compLab = []
for x in compList: compLab.append(labels[x])

def pcctest(x_data , y_data):
	pcc,pval = pearsonr(x_data, y_data)
	return abs(pcc), pval, psformat(abs(pcc), pval)


def psformat(a,b):
	A = str(np.round(a, decimals=2))

	if b < 0.01:
		string = ("{:.1e}".format(b)).split('e')
		B =  '$10^{' + string[1]+ '}$' 
	else:
		B = str(np.round(b,decimals=2))

	return A + '(' + B + ')'

def printRow(line, list):
	for i in range(len(list)):
		if i != len(list) - 1: line += list[i] + '&'
		else: line +=  list[i] + '\\\\'

	print line




##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
pccs = np.zeros((len(compList), len(compList)))
pvalues = np.zeros((len(compList), len(compList)))
strings = np.chararray((len(compList), len(compList)), itemsize=40)

printRow('~&', compLab)

# Loop through comp list
for i in range(len(compList)):
	x_col = compList[i]
	line = compLab[i] + '&'
	for j in range(len(compList)):
		y_col = compList[j]
		subCat = cat.where((cat['early']==1) & (conditions[x_col]) & (conditions[y_col]))

		#extract data
		x_data = subCat[x_col]
		y_data = subCat[y_col]

		#test
		pcc, pval, string = pcctest(x_data , y_data)

		#put in arrays
		pccs[j,i] = pcc
		pvalues[j,i] = pval
		strings[j,i] = string
		if len(x_data) < 30: string = '-'

		if pcc == 1.0: string = '-'
		if string == 'nan(nan)': string = '-'

		line += string + '&'

	print line + '\\\\'











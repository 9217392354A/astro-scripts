# numbers of each sample
# Chris Fuller, July - 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj

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


#select currentCat galaxies
firCat = cat.where(cat.DMASS_TYPE != 0)
gasCat = cat.where(np.nan_to_num(cat.HI_ALL2) > 0.0)
metalCat = cat.where(np.nan_to_num(cat.METAL) > 0.0)
gasfirCat = cat.where((cat.DMASS_TYPE != 0) & (np.nan_to_num(cat.HI_ALL2) > 0.0))
sedCat = cat.where(cat.DMASS_TYPE == 2)

cats = [cat, metalCat, firCat, sedCat, gasCat, gasfirCat]
catNames = ['SDSS (full)', 'SDSS (emission line)', 'FIR (250\mic detected)', 'FIR (All bands detected - SED)', 'HI', 'HI (250\mic detected)']

types = ['late', 'inter', 'early']

def s(x):
	return str(np.int(x))

for i in range(len(cats)):
	c = cats[i]
	name = catNames[i]


	#find total
	total = len(c)
	nTypes = []



	for j in range(3):
		t = c.where(c[types[j]] == 1)
		nTypes.append(len(t))

	if i == 0:
		per = [100.0,100.0,100.0,100.0]
		nTypes.append(total)
		typeTotals = np.array(nTypes, dtype=np.float)
	else:
		nTypes.append(total)
		per = np.array(nTypes)*100.0/typeTotals

	line = name + ' & ' + str(nTypes[0]) + ' ('+ s(per[0]) + '\\,\\%) & ' + str(nTypes[1])+ ' ('+ s(per[1]) + '\\,\\%) & ' + str(nTypes[2]) + ' ('+ s(per[2]) + '\\,\\%) & ' + str(total) + ' ('+ s(per[3]) + '\\,\\%) \\\\'

	print line



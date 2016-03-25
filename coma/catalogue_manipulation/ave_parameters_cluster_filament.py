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
hCat = cat.where(np.nan_to_num(cat.H) != 0)
jCat = cat.where(np.nan_to_num(cat.J) != 0)
kCat = cat.where(np.nan_to_num(cat.K) != 0)
zCat = cat.where(np.nan_to_num(cat.z) != 0)

parameters = ['H','J','K', 'u_1', 'g', 'r','i', 'z','SMASS', 'SMASS']
names =      ['H','J','K', 'u'  , 'g', 'r','i', 'z','Stellar Mass(NIR) ', 'Stellar Mass (All)']
cats = [hCat, jCat, kCat, cat, cat, cat, cat, cat, kCat, cat]

types = ['late', 'inter', 'early']

def si(x):
	return str(np.int(x))

def sd(x):
	return str(np.round(x, decimals = 2))

for i in range(len(parameters)):
	name = names[i]
	line = name + ' & '

	for k in range(3):
		c = cats[i]

		c = c.where(c[types[k]]==1)

		#cluster
		cluster = c.where(c.RADIUS_VIR <= 1.0)
		filament = c.where(c.RADIUS_VIR > 1.0)


		#find total
		totalC = len(cluster)
		aveC = np.mean(cluster[parameters[i]]), np.std(cluster[parameters[i]])/np.sqrt(totalC*1.0) 

		totalF = len(filament)
		aveF = np.mean(filament[parameters[i]]), np.std(filament[parameters[i]])/np.sqrt(totalF*1.0) 

		diff = (aveC[0]-aveF[0])/np.sqrt(aveC[1]**2.0+aveF[1]**2.0) 


		if k != 2:
			line += sd(diff) + ' ('+ si(totalC+totalF) +') & ' 
		else:
			line += sd(diff) + ' ('+ si(totalC+totalF) +') \\\\ ' 

	print line









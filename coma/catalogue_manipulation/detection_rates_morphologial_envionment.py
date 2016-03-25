# Number of galaxies detected and undetected in the filament and cluster
# for each morphological type
# Chris Fuller, July - 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj

#Inputs
# input/output folder
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" 
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

#Types
types = ['late', 'inter', 'early']

#cats
cluster = cat.where(cat.RADIUS_VIR <= 1.0)
filament = cat.where(cat.RADIUS_VIR > 1.0)

cats = [cluster, filament]

def detRate(tab, m):
	t = tab.where(tab[morph] == 1)
	detected = len(t.where(t.D250 == 1))
	detError = np.sqrt(detected)

	total = len(t)*1.0
	per = detected*100.0 / total
	perError = detError*100.0 / total

	return detected, detError, per, perError

def diffCal(a, errA, b, errB):
	return abs(a-b)/np.sqrt(errA**2 + errB**2)

def subLine(stats):
	det, detError, per, perError = stats
	line = str(int(det)) + ' $\\pm$ ' +  str(int(detError)) 
	line += ' (' + str(int(per)) + ' $\\pm$ ' +  str(int(perError)) + ')'
	return line


#loop through types
for morph in types:

	#cluster and filament
	cStats = detRate(cluster, morph)
	fStats = detRate(filament, morph)

	#caculate difference
	diff = diffCal(cStats[2], cStats[3], fStats[2], fStats[3])

	#cluster stats
	cline = morph + '&' + subLine(cStats) + '&'
	fline = subLine(fStats) + '&' + str(np.round(diff, decimals=2)) + '\\\\'

	print cline + fline


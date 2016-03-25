# Program to make a figure of detection rate for cluster and filament
# parameters will be given as a list
# Chris Fuller, April 2014

#import
print 'importing modules...'
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
#remove numpy runtime warings
#np.seterr(invalid='ignore')

#Inputs
print 'reading in cats'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(pj(folder,cat_name))
cluster = cat.where(cat['RADIUS_VIR'] <=  1.0)
filament = cat.where(cat['RADIUS_VIR'] >  1.0)

bands = ["100","160","250","350","500"]


for band in bands:
	d_cluster = sum(cluster["D"+band])
	p_cluster = d_cluster * 100.0 / len(cluster)

	d_filament = sum(filament["D"+band])
	p_filament = d_filament * 100.0 / len(filament)

	print band, '&', d_cluster, '&', np.int(p_cluster), '&', d_filament, '&', np.int(p_filament), '\\\\'


print 'total cluster', len(cluster)
print 'total filament', len(filament)


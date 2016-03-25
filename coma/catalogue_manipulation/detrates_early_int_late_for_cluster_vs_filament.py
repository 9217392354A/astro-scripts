#program to caculate the number and persentage of galaxies detected in the filament for each morphologcail type
# Chris Fuller, April 2014


#import moduals
from atpy import Table
import numpy as np

types = ['late','inter', 'early' ]

#Inputs
cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12.fits")
cat = cat.where(cat.SMASS > 9.1)

for i in range(3):
	t = types[i]
	print t

	#sample
	sample_det, sample_tot = len(cat.where((cat[t]==1) & (cat.D250==1)))*1.0, len(cat.where(cat[t]==1))*1.0
	sample_per, sample_per_err = sample_det*100.0/sample_tot, np.sqrt(sample_det)*100.0/sample_tot

	sample =  str(int(sample_det)) + ' of '+ str(int(sample_tot)) +  ' (' + str(int(sample_per)) + '\\,$\\pm$\\,' + str(int(sample_per_err)) + '\\,\\%)'

	#cluster
	cluster_det, cluster_tot = len(cat.where((cat[t]==1) & (cat.D250==1) & (cat.RADIUS_VIR <= 1.0)))*1.0, len(cat.where((cat[t]==1) & (cat.RADIUS_VIR <= 1.0) ))*1.0
	cluster_per, cluster_per_err = cluster_det*100.0/cluster_tot, np.sqrt(cluster_det)*100.0/cluster_tot

	cluster =  str(int(cluster_det)) + ' of '+ str(int(cluster_tot)) +  ' (' + str(int(cluster_per)) + '\\,$\\pm$\\,' + str(int(cluster_per_err)) + '\\,\\%)'

	#filament
	filament_det, filament_tot = len(cat.where((cat[t]==1) & (cat.D250==1) & (cat.RADIUS_VIR > 1.0)))*1.0, len(cat.where((cat[t]==1) & (cat.RADIUS_VIR > 1.0) ))*1.0
	filament_per, filament_per_err = filament_det*100.0/filament_tot, np.sqrt(filament_det)*100.0/filament_tot

	filament =  str(int(filament_det)) + ' of '+ str(int(filament_tot)) +  ' (' + str(int(filament_per)) + '\\,$\\pm$\\,' + str(int(filament_per_err)) + '\\,\\%)'

	#difference
	sigma = np.sqrt(filament_per_err**2 + cluster_per_err**2)
	diff = abs(cluster_per - filament_per) / sigma


	print ''
	print 'SAMPLE:      ', sample
	print 'CLUSTER:     ',cluster
	print 'FILAMENT:    ',filament
	print 'DIFFERENCE:  ',diff
	print ''
	print ''
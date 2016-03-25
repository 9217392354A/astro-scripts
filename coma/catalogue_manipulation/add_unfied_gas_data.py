#program to plot STELLAR MASS
# Chris Fuller, April 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MaxNLocator
import scipy

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))



gas1 = np.nan_to_num(cat.GMASS2)
gas2 = np.nan_to_num(cat.logMsun)

gasTot = np.array([0.0]*len(cat))

#loop through col and select gas data
g1 = 0
g2 = 0
for i in range(len(cat)):
	if gas2[i] > 0.0:
		gasTot[i] = gas2[i]
		g2 += 1


	elif gas1[i] > 0.0:
		gasTot[i] = gas1[i]
		g1 += 1


print g1
print g2
cat.add_column('HI_ALL2', gasTot)
cat.write(pj(folder,'test1.fits'))



# Program to make a figure of detection rate vs any parameter of choice
# parameters will be given as a list
# Chris Fuller, March 2014

#import
print 'importing modules...'
from atpy import Table
import numpy as np
from os.path import join as pj

#Inputs
print 'reading in cats'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12.fits' #input name
t = Table(pj(folder,cat_name))


t.add_column('early', [0]*len(t), dtype=np.int16)
t.add_column('late', [0]*len(t), dtype=np.int16)
t.add_column('inter', [0]*len(t), dtype=np.int16)

t['early'][np.where(t.pE0 >= 0.8)[0]] = 1
t['late'][np.where(t.pS0 >= 0.8)[0]] = 1
t['inter'][np.where((t.pE0 < 0.8) & (t.pS0 < 0.8))[0]] = 1



for i in range(len(t)):
	total = t.early[i] + t.late[i] + t.inter[i]
	if total != 1: print 'error', i, t.early[i], t.late[i], t.inter[i], t.pE0[i], t.pS0[i], t.goldmine[i]


t.write(pj(folder, 'test.fits'), overwrite=True)


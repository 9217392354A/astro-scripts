#program to create dust-to-stars ratio
# Chris Fuller, April 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj



#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

w1 = np.where(cat.D250==1)[0]

new_col = np.array([0.0]*len(cat))

new_col[w1] = cat.DMASS[w1] - cat.SMASS[w1]

cat.add_column('DUST_STARS', new_col)
cat.write(pj(folder,'test-dust-v2.fits'))
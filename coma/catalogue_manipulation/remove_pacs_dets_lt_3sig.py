#program to remove pacs100 and pacs160 values that are above the 3sigma noise limit
# Chris Fuller, Oct 2014

import numpy as np
from atpy import Table
from os.path import join as pj

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

#noise limits
pacs100lim = (24.0*10**-3)*3.0
pacs160lim = (43.8*10**-3)*3.0


#run through each band and select all galaxies needing removal
rem_100 = np.where((cat.F100 < pacs100lim) & (cat.F100 > 0.0))[0]
rem_160 = np.where((cat.F160 < pacs160lim) & (cat.F160 > 0.0))[0]
sed_change_type = np.where((cat.F160 < pacs160lim) & (cat.F100 < pacs100lim) & ((cat.F100 > 0.0) | (cat.F160 > 0.0)))[0]


#edit it table


#100
cat.F100[rem_100] = 0.0
cat.SN100[rem_100] = 0.0
cat.R100[rem_100] = 0.0
cat.D100[rem_100] = 0
cat.DMASS_TYPE[rem_100] = 1
cat.DMASS[rem_100]  = cat.DMASS_250[rem_100]
cat.DUST_STARS[rem_100] = cat.DMASS[rem_100] - cat.SMASS[rem_100]


#100
cat.F160[rem_160] = 0.0
cat.SN160[rem_160] = 0.0
cat.R160[rem_160] = 0.0
cat.D160[rem_160] = 0
cat.DMASS_TYPE[rem_160] = 1
cat.DMASS[rem_160]  = cat.DMASS_250[rem_160]
cat.DUST_STARS[rem_160] = cat.DMASS[rem_160] - cat.SMASS[rem_160]

#change dmass_type
#cat.DMASS_TYPE[sed_change_type] = 1
#cat.DMASS[sed_change_type]  = cat.DMASS_250[sed_change_type]
#cat.DUST_STARS[sed_change_type] = cat.DMASS[sed_change_type] - cat.SMASS[sed_change_type]


cat.write(pj(folder,"coma_supercluster_cal12_pacscorrected.fits"), overwrite=True)


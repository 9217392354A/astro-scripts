#scaling relations plots
#Chris Fuller

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(pj(folder,fname))

w_metal = np.where(np.nan_to_num(cat.METAL) > 0.0) 
w_gas_metal = np.where((np.nan_to_num(cat.METAL) > 0.0) & (np.nan_to_num(cat.HI_ALL2) > 0.0))

tempGas = np.nan_to_num(np.nan_to_num(cat.HI_ALL2) / np.nan_to_num(cat.HI_ALL2))
tempMetal = np.nan_to_num(np.nan_to_num(cat.METAL) / np.nan_to_num(cat.METAL))


zGas = 29.2 * 10.0**(cat.METAL- 12.0)*tempMetal

MzGas =  (zGas * 10**np.nan_to_num(cat.HI_ALL2)) * tempGas * tempMetal

Mmetals = 1.2 * MzGas + 10**cat.DMASS  * tempGas * tempMetal

Mgas = 2.5*10**np.nan_to_num(cat.HI_ALL2) * tempGas * tempMetal


Ztot = Mmetals / Mgas


"""
MzGas[w_gas_metal] = np.log10(MzGas[w_gas_metal])
Mmetals[w_gas_metal] = np.log10(Mmetals[w_gas_metal])
Mgas[w_gas_metal] = np.log10(Mgas[w_gas_metal])


#caculated extra columns
cat.add_column('zGas', zGas)
cat.add_column('MZGAS', MzGas)
cat.add_column('Mmetals', Mmetals)
cat.add_column('Ztot', Ztot)
cat.add_column('Mgastot', Mgas)

#cat.write(pj(folder,'test-dust-v2.fits'), overwrite=True)
"""

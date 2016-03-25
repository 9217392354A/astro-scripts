# Program to Name galaxies

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
cat = Table(pj(folder,fname)


cat.add_column('OBJECT_NAME_PAPER', Mgas)

#cat.write(pj(folder,'test-dust-v2.fits'), overwrite=True)


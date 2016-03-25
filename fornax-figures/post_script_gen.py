#program to create ps of all galaxies, split the figures up by morphological type
# Chris Fuller, September 2013

#import modules
from os.path import join as pj
import numpy as np
import aplpy as ap
from atpy import Table
from os import mkdir, listdir
import matplotlib.pyplot as plt

###################### inputs ######################
cat = Table(pj("/Users/chrisfuller/Dropbox/phd/plots/upper_limits/","fornax.fits"))
folder_fir = "/Users/chrisfuller/dropbox/phd/herchel/fornax/galaxies/"

#################### functions ######################

#################### main ###########################

morph = [0,1] # this controls the selection of what morphological class we are using

#make selection of galaxies based on morph
w = np.where((cat.goldmine >= morph[0]) & (cat.goldmine <= morph[1]) & (np.nan_to_num(cat.F250) != 0.0))[0]
print 'number of galaxies ', len(w)

dx = 1.0 / np.ceil(np.sqrt(len(w)))
dy = 1.0 / np.ceil(np.sqrt(len(w)))


XX = np.arange(0.,1.,dx)
YY = np.arange(1.,0.,-1*dy)

#create fig
fig = plt.figure(figsize = (8.5,8.5),facecolor='w',edgecolor='w')

#set counter to 0
i = 0

#first loop through XX
for Y in YY:
	for X in XX:
		if i >= len(w): continue
		obj =  cat['OBJECT_1'][w[i]]

		fits = ap.FITSFigure(folder_fir +  str(obj) +'/' + str(obj) + "-PSWmap-rawmap.fits", figure=fig,subplot=[X,Y,dx,dy] )
		fits.show_grayscale()
		#increse counter by 1
		i += 1


fig.savefig("/Users/chrisfuller/Desktop/early.eps")
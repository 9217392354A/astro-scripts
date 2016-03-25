# program to plot a single optical image and then overplot fir contors

import aplpy as ap
from os.path import join as pj
from matplotlib.pyplot import figure,show
#from astropy.io import 

folder = '/Users/chrisfuller/Dropbox/phd/plots/FCC215/'

help('modules')

# galaxy info
ra = 54.657083333333325
dec = -35.75749999999999
a = 0.19009469085522004/60.0
b = 0.19009469085522004/60.0
pa = 0.

fig = figure(figsize = (4.5,4.5),facecolor='w',edgecolor='w')
f = ap.FITSFigure(pj(folder,"215-opt-smit.fits"),figure=fig)

#f.recenter(ra,dec, a*3.0)
f.show_grayscale()

fig.savefig('/Users/chrisfuller/Desktop/fcc215.pdf')


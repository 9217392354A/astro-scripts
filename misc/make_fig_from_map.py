#program to make figure for paper using ngp maps, by Chris Fuller, Jan 2013

import aplpy
import numpy as np
import atpy as at
from os.path import join as pj
import matplotlib.pyplot as plt

#folder = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/"
#cat = at.Table(pj(folder,"HeFoCS-fluxes-260313-mybgsub-v2-handmeasure-bad-detections-removed-morphology.fits"),type='fits')

folder = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/Catalogues/"
cat = at.Table(pj(folder,"fcc.fits"),type='fits')

print 'reading in fits file'
#fig = aplpy.FITSFigure('/home/herdata/spx6cff/maps/NGP-entire-PLWmap-mosaic-20120512.fits')
fig = aplpy.FITSFigure('/Users/chrisfuller/Documents/maps/fornax/HeFoCS-All-PSWmap-mosaic_MS-DR35.fits', figsize=(9.5,9.5))
#fig.show_grayscale()
fig.show_colorscale(cmap='copper')
#fig.recenter(54.6289,-35.4545,4.0)


x = cat.FCC_RA2000
y = cat.FCC_DEC2000

#comax = np.array([194.903])
#comay = np.array([27.903])
#radius = np.array([1.67])
 

fig.show_markers(x, y,c='#90EE90',marker='+', s=80)
#fig.show_circles(comax,comay,radius)
fig.add_scalebar(3.0*1.0)
fig.scalebar.set_label('1 Mpc')
#fig.show_grid()
fig.set_tick_labels_format(xformat='hh:mm',yformat='dd:mm')
#fig.set_tick_labels_style('latex')
#fig.set_labels_latex(True)


#fig.save('/Users/chrisfuller/Dropbox/phd/papers/fornax/fornax.pdf')#, dpi=600) 
fig.save('/Users/chrisfuller/Desktop/pdfs/fornax.pdf')#, dpi=600) 









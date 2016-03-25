#program to make figure out of ngp
# Chris Fuller, April 2014

import atpy as at
import aplpy as ap
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder

c1 = 'coma_supercluster_cal12.fits' #input name
cat = at.Table(pj(folder,c1))
alfaCat = cat.where(np.nan_to_num(cat.GMASS) > 0.0)
litgCat = cat.where((np.nan_to_num(cat.HI_ALL2) > 0.0) & (np.nan_to_num(cat.GMASS) == 0.0))

ra_all = cat.GRA2000
dec_all = cat.GDEC2000

ra_litg = litgCat.GRA2000
dec_litg = litgCat.GDEC2000

ra_alfa = alfaCat.GRA2000
dec_alfa = alfaCat.GDEC2000



imroot = '/Users/chrisfuller/Dropbox/NGP-entire-PLWmap-mosaic-20131121.fits'

#create figure
fig = plt.figure(figsize = (8.0,8.0),facecolor='w',edgecolor='w') 

f = ap.FITSFigure(imroot, hdu=1, figure=fig, subplot=(1, 1, 1), downsample=1)
f.show_colorscale(cmap='gist_heat')
#f.show_(ra, dec, edgecolor='white')
f.add_scalebar(1.67*2)
f.scalebar.set_label('2 Mpc')
f.scalebar.set_color('black')
f.scalebar.set_font(size='medium', weight='heavy', stretch='normal', family='sans-serif', style='normal', variant='normal')

#f.show_grid()
f.set_tick_labels_format(xformat='hh:mm',yformat='dd:mm')
f.show_markers(ra_all, dec_all, c='red',marker='x', edgecolor = 'r')
#f.show_markers(ra_litg, dec_litg,c='green',marker='s', edgecolor = 'green', facecolor = 'None')
#f.show_markers(ra_alfa, dec_alfa,c='blue',marker='s', edgecolor = 'blue', facecolor = 'None')




f.show_circles(194.9531, 27.9807, 1.67, edgecolor='black', lw=1)
#f.show_rectangles(199.5, 28.0, 15.0, 0.02, edgecolor='blue', lw=1)

#f.add_beam(major=beam, minor=beam,angle=0.0)
#f.show_beam(major=beam, minor=beam,angle=0.0,fc='white')
#f.axis_labels.set_font(size=8)
#f.tick_labels.set_font(size=6)
f.set_theme('publication')

plt.subplots_adjust(left=0.12, bottom=0.01, right=0.96, top=0.99, hspace=0.0, wspace=0.0)



#save and show
fig.savefig(pj('/Users/chrisfuller/Dropbox/phd/thesis/main/chapter-coma/','ngp-big-ccconly.pdf'))
plt.show()

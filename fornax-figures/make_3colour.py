# program to make nice optical figure
# Chris Fuller, August 2013

import aplpy
from os.path import join as pj
from atpy import Table
from os import chdir

#paths
optical_path = '/Users/chrisfuller/Documents/maps/fornax/sss/'
fir_path     = '/Users/chrisfuller/Documents/maps/fornax/'
folder       = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs"

#change to sss dir
chdir(optical_path)

#files
optical_files = ['fornax_R.fits', 'fornax_J.fits', 'fornax_I.fits']
spire_file	 = 'HeFoCS-All-PLWmap-mosaic_MS-DR35.fits'
pacs_file	 = 'F-R1R2R3R4-pacs160-scanamorphos.fits'

#import optical catalouge
cat = Table(pj(folder,"HeFoCS-final-morph.fits"),type='fits')

print 'reading in fits file'
#fig = aplpy.FITSFigure(pj(optical_path, optical_files[0]), downsample=300)
#fig.show_grayscale()
#fig.recenter(54.6289,-35.4545,4.0)
fig = aplpy.FITSFigure(pj(fir_path, spire_file), hdu=1)
fig.show_contour(pj(fir_path, spire_file), hdu=1, levels=[0.001])


x = cat.GRA2000
y = cat.GDEC2000
a = cat.FULLMAJAX
b = cat.FULLMINAX
pa = cat.PA

fig.show_markers(x, y,c='r',marker='o')
#fig.show_ellipses(x, y, a/60., b/60., angle=pa)
fig.add_scalebar(1.0/0.3)
fig.scalebar.set_label('1 Mpc')
fig.show_grid()
fig.set_tick_labels_format(xformat='hh:mm',yformat='dd:mm')


fig.save('/Users/chrisfuller/Desktop/fornax.png')#, dpi=600) 

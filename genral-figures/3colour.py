#Chris Fuller, May 21st 2013
#Figure maker

import aplpy as ap
from os.path import join as pj
import atpy as at
import os

folder = "/Users/chrisfuller/Documents/maps/fornax/sss"

cat = at.Table(pj(folder,"cat.fits"),type='fits')


os.chdir(folder)


#ap.make_rgb_cube(['fornax_R.fits', 'fornax_J.fits','fornax_I.fits'], 'fornax_cube.fits')

# Make an RGB image
#ap.make_rgb_image('fornax_cube.fits', 'fornax_rgb.png')

# Plot the RGB image using the 2d image to indicate the projection
#f = ap.FITSFigure('fornax_cube_2d.fits')
#f.show_rgb('fornax_rgb.png')

# save
#f.save(pj(folder,"output.png"))

ra = cat.GRA2000
dec = cat.GDEC2000


 
f = ap.FITSFigure('fornax_J.fits',downsample=10)
f.show_grayscale()
#f.show_contour('PLW-edgemask-mask.fits', levels=[1.0],color='k')
f.show_markers(ra,dec)
f.add_scalebar(0.3331)
f.scalebar.set_label('100 kpc')
f.scalebar.set_color('white')
f.show_grid()
f.axis_labels.set_font(size=14)
f.tick_labels.set_font(size=12)
f.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')

f.save('fornax.pdf')
f.save('fornax.png')
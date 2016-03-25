# fig gen for paper, create figure with 3 galaxies 1 good 1 blended 1 background 
# Chris Fuller May 2013
import os
import numpy as np  
from os.path import join as pj
import aplpy
from atpy import Table
import matplotlib.pyplot as plt
###################### inputs ######################
folder = "/Users/chrisfuller/Dropbox/phd/plots/fig_gen/"

print "reading in cat"
cat = Table(pj(folder,"cat.dat"),type='ascii')
os.chdir(folder)

#################### functions ###################

#function to plot FIR contours
def plotta(name,texts,y):
    #find optical parameters
    i = np.where(cat.OBJECT == name)[0][0]
    a,b,pa = cat.FULLMAJAX[i],cat.FULLMINAX[i], cat.PA[i]
    ra,dec = cat.GRA2000[i], cat.GDEC2000[i]
    beam = 18.2/3600.0
    name =int(name)
    #make plot
    f = aplpy.FITSFigure(folder+str(name)+".fits",figure=fig,subplot=[xL,y,dx,dy])
    print folder + str(name) + ".fits"
    f.recenter(ra,dec,beam*3.0)
    f.show_grayscale()
    f.show_ellipses(ra, dec, a/60.0, b/60.0, angle=pa-90, edgecolor='white', lw=1, alpha=0.3)
    f.show_contour(str(name)+"-PSWmap-PSF-convolvedmap.fits", levels=50, filled=True, cmap='cool') 
    f.add_scalebar(1.0/60.0)
    f.scalebar.set_label('1 arcmin')
    f.scalebar.set_color('white')
    #f.show_grid()
    f.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    #f.show_markers(ra, dec,c='r',marker='x')
    #f.add_beam(major=beam, minor=beam,angle=0.0)
    f.show_beam(major=beam, minor=beam,angle=0.0,fc='white')
    f.axis_labels.set_font(size=8)
    f.tick_labels.set_font(size=6)
    
    # add text 
    plt.figtext(xR, y+0.23,texts, size=12, weight="semibold", color='white')
    
 
#################### control ######################
#create figure
fig = plt.figure(figsize = (3.6,9.5),facecolor='w',edgecolor='w')

#width and height of subplots
dx = 0.70
dy = 0.25            
#cord for subplots
xL = 0.22
xR = 0.2215

y1 = 0.70
y2 = 0.4
y3 = 0.1


plotta(117.,"(a) Good",y1)
plotta(135.,"(b) Background",y2)
plotta(136.,"(c) Blended",y3)



#plt.savefig("/Users/chrisfuller/Dropbox/phd/papers/fornax/good_back_blend.pdf")
plt.savefig('/Users/chrisfuller/Desktop/pdfs/good_back_blend.pdf')
plt.show()

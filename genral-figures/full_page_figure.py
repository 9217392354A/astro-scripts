# program to make a full page figure of all detections
# Chris Fuller, 9th December 2013

import numpy as np
import aplpy as ap
import atpy as at
from os.path import join as pj
import matplotlib.pyplot as plt

# inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(folder+fname,type='fits')

cat = cat.where(cat.DMASS_TYPE != 0)


######################## Functions ###########################
# plot image
def plotta(im,xx,yy,j):
    #make plot
    f = ap.FITSFigure(im,figure=fig,subplot=[xx,yy-dy,dx,dy])
    

    f.show_colorscale(cmap = 'hot', vmin = -1.245E-02, vmax  = 2.405E-02 )
    f.recenter(ra, dec, re_cen)

    f.show_ellipses(ra, dec, a/60.0, b/60.0, angle=pa-90, edgecolor='grey', lw=2)
    #f.add_scalebar(1.0/60.0)
    #f.scalebar.set_label('1 arcmin')
    #f.scalebar.set_color('white')

    #f.show_grid()
    #f.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    #f.show_markers(ra, dec, c='red',marker='x')
    #f.add_beam(major=beam, minor=beam,angle=0.0)
    f.show_beam(major=beam, minor=beam,angle=0.0,fc='white')
    f.axis_labels.hide()
    f.tick_labels.hide()
    
    # add text 
    plt.figtext(xx+0.012 , yy-dy+0.14, name, size=20, weight="bold", color='white')
    
######################## End Functions #######################

#sort by object
#cat.sort('OBJECT')

#find all detected galaxies
flux = cat['F250']
OBJECT = cat['PAPER_NAME']
detected = np.where(np.nan_to_num(flux)>0.0)[0]
ext = cat.EXTENDEDNESS

print 'number of detected: ', len(detected)

# x and y
XX = 5
YY = 6

#size of image
dx = 1.0 / XX
dy = 1.0 / YY


#create figure
fig = plt.figure(figsize = (8.0 , 11.5), facecolor='w', edgecolor='w')

# start i
i = 0

#set cdelt
cdelt = 6.0/3600.

#loop through all the detected galaxies
for y in np.arange(1.0, 0.0, -dy):
	for x in np.arange(0.0, 1.0, dx):
		#try:
		if True:
			if i >= len(detected): continue
			#folder where the FIR images are stored
			root = pj(folder,str(OBJECT[i]))
			

			#find optical parameters
			a,b,pa = cat.FULLMAJAX[i],cat.FULLMINAX[i], cat.PA[i]
			ra,dec = cat.GRA2000[i], cat.GDEC2000[i]
		 	beam = 20./3600.0
		 	name = cat.OBJECT[i]
		 	D_FIR  = cat.R250[i]

		 	e = b/a

		 	if name == 121: D_FIR = 0.7

		 	print OBJECT[i],flux[i], ext[i][2], x, y
			#if image is extended
			if ext[i][2] == "E":
				re_cen = np.sqrt(((D_FIR * 2.0) /60.0)**2 + ((0.15*5.0)/60.0)**2)
				#find FIR image
				image = pj(root, str(OBJECT[i])+ "-PSWmap-rawmap.fits")


			#if image is a poinst source
			elif ext[i][2] == "P":
				re_cen = (0.15*9.0) /60.0
				image = pj(root, str(OBJECT[i])+ "-PSWmap-PSF-convolvedmap.fits")

			else: 
				print 'error image not found'

			#plot images
			plotta(image,x,y,i)


			#update counter
			i += 1

		#except: print i
#save figure
#plt.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/fornax','detections.pdf'),dpi=600)
plt.savefig('/Users/chrisfuller/Desktop/detections.pdf')
plt.show()



#Program created to overlay dust contors onto SDSS images.
#Chris Fuller, Nov 2012

import numpy as np
import aplpy as ap
from os.path import join as pj
import os
import matplotlib.pyplot as plt

b = 18.2/3600.0 



jpegfolder = "/home/fossil/spx6cff/coma/dust_contours/jpeg-files/"
fitsfolder = "/home/fossil/spx6cff/coma/dust_contours//fits-files-convolved/"
outfolder = "/home/fossil/spx6cff/coma/dust_contours/outfiles/"

fitsFiles = os.listdir(fitsfolder)
jpegFiles = os.listdir(jpegfolder)

fitsFiles = [x for x in fitsFiles if x[-4:] == "fits"]
jpegFiles = [x for x in jpegFiles if x[-3:] == "jpg"]

fits = np.sort(fitsFiles)
jpeg = np.sort(jpegFiles)
#fig = plt.figure(figsize=(10,9))
#f1 = ap.FITSFigure(pj('CCC2040-PSWmap-rawmap.fits'), hdu=0, figure=fig)
#f1.show_rgb("CCC2040-SDSS.jpg")
#f1.show_contour('CCC2040-PSWmap-rawmap.fits')
#outName = "plot.png"
#f1.save(outName)

if len(fitsFiles) != len(jpegFiles):
    raise

os.chdir(outfolder)
for i in range(0,len(fitsFiles)):
    outName = fits[i].split('-')[0]+"-dust_contour_plot.png"
    fig = plt.figure(figsize=(10,9))
    f1 = ap.FITSFigure(pj(fitsfolder,fits[i]), hdu=0, figure=fig)
    f1.show_rgb(pj(jpegfolder,jpeg[i]))
    f1.show_contour(pj(fitsfolder,fits[i]), levels=15)
    f1.save(outName)
    print 'saving...',outName, ' :', np.around(float(i)*100/float(len(fitsFiles))), '%'
    
print 'finished'

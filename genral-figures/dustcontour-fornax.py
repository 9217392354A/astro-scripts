#program to make dust contours for fornax
#designed to be used with superCOSMOS-part2 products
#Chris Fuller, Feb 2013

import os
import numpy as np
from os.path import join as pj
import aplpy as ap
import atpy as at

###################### inputs ######################
rootSuper = "/Users/chrisfuller/Desktop/current/superCOSMOS/fits/"
rootFIR = "/Users/chrisfuller/Desktop/current/flux/galaxies/"
outfolder = "/Users/chrisfuller/Desktop/current/dust_contors/"
catfolder = "/Users/chrisfuller/Desktop/current/flux/"
PGC =  "/Users/chrisfuller/Desktop/current/"
try: os.mkdir(outfolder)
except: print "outfolder aready exists..."

#bands = ["PLW","PMW","PSW","pacs160","pacs100"] #insert bands as needed
bands = ["PSW"]

print "reading in cat"
cat = at.Table(pj(catfolder,"fcc_input_ready.fits"),type='fits')
pgc = at.Table(pj(PGC,"FCC-PGC-matched-fulloptical.fits"))
##################### functions #####################
#find if there is a PSF file or just and raw map and return the src fro it
def source(folder,band):
    print folder
    PSF = [s for s in os.listdir(folder) if band+'map-PSF' in s]
    EXT = [s for s in os.listdir(folder) if band+'map-rawmap' in s]
    if len(PSF) == 1.0:
        x = pj(folder,PSF[0])
    elif len(EXT) == 1.0: 
        #print EXT[0]
        x = pj(folder,EXT[0])
    else: 
        print 'no file found in ',folder," .....",band
    return x
    
def dustcontor(FIR, SUPER, outfolder, obj,band):
    os.chdir(outfolder)
    outName = obj+"-"+band+"-dust_contour.eps"
    #fig = plt.figure(figsize=(10,9))
    f1 = ap.FITSFigure(SUPER)
    f1.show_grayscale()
    f1.show_contour(FIR, levels=75)
    #x,y =findcen(cat,fits[i].split('-')[0])
    x,y,a,ba,pa = findop(cat,i)
    #find pgc stuff
    xPGC,yPGC,t,sep = findopPGC(pgc,i)
    f1.add_scalebar(1.0/60.0)
    f1.scalebar.set_label('1 arcmin')
    f1.show_grid()
    f1.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    f1.set_tick_labels_style('latex')
    #f1.set_labels_latex(True)
    f1.show_markers(x, y,c='r',marker='x')
    f1.show_markers(xPGC, yPGC, c='g',marker='o')
    f1.show_ellipses(x, y, a/60.0, ba/60.0, angle=pa, edgecolor='red')
    #f1.add_beam()
    #f1.beam.set_major(0.03) # degrees
    #f1.beam.set_minor(0.02) # degrees
    #f1.beam.set_angle(45.) # degrees
    #f1.beam.show()
    f1.save(outName)
    
def findop(cat,i):
    x = float(cat.GRA2000[i])
    y = float(cat.GDEC2000[i])
    a = float(cat.FULLMAJAX[i])
    b = float(cat.FULLMINAX[i])
    pa = float(cat.PA[i])
    return x,y,a,b,pa
    
def findopPGC(cat,i):
    x = float(cat.PGC_RAJ2000[i])
    y = float(cat.PGC_DECJ2000[i])
    t = str(cat.MType_FCC[i])
    sep = float(cat.Separation[i])
    return x,y,t,sep
    
    
##################### main program #################
#make a list of root folders in FIR
#foldersFIR = os.listdir(rootFIR)
#superFiles = os.listdir(rootSuper)

#loop through each band and make dust contour maps
for band in bands:
    print "starting band ",band
    #loop through each galaxy
    for i in range(0,len(cat)):
        #sourceFIR
        gal = str(cat.OBJECT[i])
        galFolder = pj(rootFIR,gal)
        #try:
        if gal =="121":
            srcFIR = source(galFolder,band)
            srcSUPER= pj(rootSuper,gal+'.fits')
            dustcontor(srcFIR,srcSUPER,outfolder,gal,band)
        #except:
         #   print ".............galaxy not measured................"
          #  print ">>>>>>>>>>>>>>>>>", gal,"<<<<<<<<<<<<<<<<<<<<<<<"
           # print ">>>>>>>>>>>>>>>>>",band,"<<<<<<<<<<<<<<<<<<<<<<<"
            
            
    




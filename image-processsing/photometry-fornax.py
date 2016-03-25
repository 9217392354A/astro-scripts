#program to make post script for each dectected galaxy
#Chris Fuller, Feb 2013

import os
import numpy as np  
from os.path import join as pj
import aplpy as ap
import atpy as at
import matplotlib.pyplot as plt
from scipy import pi,sin,cos
import math as m
import pyfits as pf
from pygame import mixer
mixer.init() #you must initialize the mixer
alert=mixer.Sound('/Users/chrisfuller/Dropbox/personal/sounds/R2D2e.wav')


#switches
galaxy = '69'
stype  = 'EXTENDED'
bands = ['PLWmap','PMWmap','PSWmap','pacs160','pacs100']
mapfile=['HeFoCS-All-PLWmap-mosaic_MS-DR35.fits','HeFoCS-All-PMWmap-mosaic_MS-DR35.fits','HeFoCS-All-PSWmap-mosaic_MS-DR35.fits','F-R1R2R3R4-pacs160-scanamorphos.fits','F-R1R2R3R4-pacs100-scanamorphos.fits']

#FWHM of different beams
fwhms = [36.0/60.,24.5/60.,18.2/60.,13.4/60.,9.4/60.]
#Calibration uncertainties
calnoise = [0.07,0.07,0.07,0.12,0.12]
#set the beam area
beams=[1587./3600.,751./3600.,423./3600,206./3600.,100./3600.]# beam area in sq arcmin PACS maps are already in Jy/px

###################### inputs ######################
rootFIR = "/Users/chrisfuller/Desktop/fornax/galaxies/"
catfolder = "/Users/chrisfuller/Desktop/fornax/"
mapdir = '/Users/chrisfuller/Desktop/fornax/'

print "reading in cat"
cat = at.Table(pj(catfolder,"HeFoCS-fluxes-120313_re-run.fits"),type='fits')
##################### functions #####################
#find if there is a PSF file or just and raw map and return the src fro it
def source(folder,band,t):
    PSF = [s for s in os.listdir(folder) if band+'-PSF' in s]
    EXT = [s for s in os.listdir(folder) if band+'-rawmap' in s]
    
    #if len(PSF) == 1.0:
    if t =="EXTENDED":
        x = pj(folder,EXT[0])
    #elif len(EXT) == 1.0: 
        #print EXT[0]
    elif t == "POINT":
        x = pj(folder,PSF[0])
    else: 
        print 'no file found in ',folder," .....",band
    return x

#convole beam with optical radius to set the size of the diagnostic frame
def con(opa,band):
    #find band
    if "PLW" in band: beam = 35.0/60.0
    if "PMW" in band: beam = 25.0/60.0
    if "PSW" in band: beam = 18.0/60.0
    if "160" in band: beam = 12.0/60.0
    if "100" in band: beam = 8.0/60.0
    
    r = np.sqrt(opa**2 + beam**2)
    return r
    
def findop(f,i):
    x = float(f.GRA2000[i])
    y = float(f.GDEC2000[i])
    a = float(f.FULLMAJAX[i])
    b = float(f.FULLMINAX[i])
    pa = float(f.PA[i])
    morph = f.MTYPE_PGC[i]
    altname = f.NAMES[i]
    return x,y,a,b,pa
    
#creates an ellipse frame from which to selected pixels
def dE(x,y,a,b,PA):
    the = -PA+90.0
    cos_a,sin_a=cos(the*pi/180.0),sin(the*pi/180.0)
    Xs = (x*cos_a - y*sin_a)**2.0 
    Ys = (x*sin_a + y*cos_a)**2.0 
    r = np.sqrt((Xs/a**2) + (Ys/b**2)) 
    return r
    
#creates an ellipse frame from which to selected pixels    
def ellipse(image,cdelt,a,b,PA):
    c = cdelt*60.0
    #find central pixel and shape
    x0 = np.array(image.shape[0])/2.0
    y0 = np.array(image.shape[1])/2.0
    xMax = x0*2
    yMax = y0*2
    x = np.arange(0-x0, xMax-x0)
    y = np.arange(0-y0, yMax-x0)
    #pixel cordinates 
    xx, yy = np.meshgrid(x, y)
    #creat image of ellipses for appature photometry
    new_image= dE(xx,yy,a,b,PA)*c
    return new_image
 


 

##################### main program #################
#loop through each band and make dust contour maps
#for j in range(0,5):
for j in range(0,1):
    band = bands[j]
    fBigmap = mapfile[j]
    #import big map
    if j<3: #SPIRE
        print "SPIRE MAP: ", band, " FILENAME: ", fBigmap
        hdulist = pf.open(pj(mapdir,fBigmap))
        bigmap = hdulist[0].data
        bighdr = hdulist[0].header
        cdelt = float(bighdr['CDELT2'])
    if j>=3: #PACS
        print "PACS MAP: ", band, " FILENAME: ", fBigmap
        bigmap = hdulist[0].data
        bighdr = hdulist[0].header
        cdelt = float(bighdr['CDELT2'])
    
    print "starting band ",band
    #loop through each galaxy
    for i in range(0,len(cat)):
        gal = str(cat.OBJECT[i])
        galFolder = pj(rootFIR,gal)
        
        #test if detected        
        if gal == galaxy:
            print "Starting FCC", gal
            os.chdir(galFolder)
            #optical dust gal parameters            
            x,y,a,b,pa = findop(cat,i)
            
            print x,y
            print a,b,pa

            FIRraw = source(galFolder,band,stype)
            #hdulist = pf.open(FIRraw)
            #RAWimage = hdulist.data
            #RAWhdr = hdulist.he
            #flux = photometry()
            #FIRellipse = pj(galFolder,gal+"-"+band+"-ellipses.fits")




print "finished!"
#alert.play()

            
    




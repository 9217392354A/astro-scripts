#program to make post script for each dectected galaxy
#Chris Fuller, Feb 2013
# v2

import os
import numpy as np  
from os.path import join as pj
import aplpy as ap
#import atpy as at
import matplotlib.pyplot as plt
import pyfits as pf
#from pygame import mixer#
#mixer.init() #you must initialize the mixer
#alert=mixer.Sound('/Users/chrisfuller/Dropbox/personal/sounds/R2D2e.wav')

###################### inputs ######################
rooty = "/Users/chrisfuller/Dropbox/phd/herchel/coma/"
rootFIR = "/Users/chrisfuller/Documents/maps/coma/galaxies-mybgsub"
outfolder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/all-detections/"
catfolder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/'

try: os.mkdir(outfolder)
except: print "outfolder aready exists..."


bands = ["PSWmap"]

print "reading in cat"

#cat = at.Table(pj(catfolder,"HeFoCS-fluxes-120313_re-run.fits"),type='fits')
hdulist = pf.open(pj(catfolder,"cluster+filament-mybg-130114-3.1arcsec.fits"))
cat = hdulist[1].data
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
    
def background(folder,band):
    print folder
    PSF = [s for s in os.listdir(folder) if band+'-PSF' in s]
    EXT = [s for s in os.listdir(folder) if band+'-rawmap' in s]
    if len(PSF) == 1.0:
        x = pj(folder,PSF[0])
    elif len(EXT) == 1.0: 
        #print EXT[0]
        x = pj(folder,EXT[0])
    else: 
        print 'no file found in ',folder," .....",band
    return x
    
def findop(f,i):
    x = float(f.field('GRA2000')[i])
    y = float(f.field('GDEC2000')[i])
    a = float(f.field('FULLMAJAX')[i])
    b = float(f.field('FULLMINAX')[i])
    pa = float(f.field('PA')[i])

    return x,y,a,b,pa

#finds FIR flux and other suff
def galParamFIR(cat,i,band):  
    #find band
    #print cat.EXTENDEDNESS[i]
    if "PLW" in band: 
        flux = cat.field('F500')[i]
        error = cat.field('E500')[i]
        SN = cat.field('SN500')[i]
        t = cat.field('EXTENDEDNESS')[i][0]
        r = cat.field('R500')[i]
    if "PMW" in band: 
        flux = cat.field('F350')[i]
        error = cat.field('E350')[i]
        SN = cat.field('SN350')[i]
        t = cat.field('EXTENDEDNESS')[i][1]
        r = cat.field('R350')[i]
    if "PSW" in band: 
        flux = cat.field('F250')[i]
        error = cat.field('E250')[i]
        SN = cat.field('SN250')[i]
        t = cat.field('EXTENDEDNESS')[i][2]
        r = cat.field('R250')[i]
    if "160" in band: 
        flux = cat.field('F160')[i]
        error = cat.field('E160')[i]
        SN = cat.field('SN160')[i]
        t = cat.field('EXTENDEDNESS')[i][3]
        r = cat.field('R160')[i]
    if "100" in band: 
        flux = cat.field('F100')[i]
        error = cat.field('E100')[i]
        SN = cat.field('SN100')[i]
        t = cat.field('EXTENDEDNESS')[i][4]
        r = cat.field('R100')[i]
        
    if t== "E": t="EXTENDED"
    if t== "P": t="POINT"
    
    #if flux not present then return a string of undetected
    if str(flux) == "": return "not included at this band", "NO FLUX", ""
    
    #test if detected
    if float(flux) == 0.0 or float(SN) == 0.0:
        return "Total flux = "+str(round(float(flux),3))+"+/-"+str(round(float(error),3))+" Jy,  Total SN = " +str(round(float(SN),1)), "UPPERLIMIT",t,r
        
    else: return "Total flux = "+str(round(float(flux),3))+"+/-"+str(round(float(error),3))+" Jy,  Total SN = " +str(round(float(SN),1)), "DETECTED",t,r
    
#convole beam with optical radius to set the size of the diagnostic frame

def rad(opa,band):
    #find band
    if "PLW" in band: beam = 35.0/60.0
    if "PMW" in band: beam = 25.0/60.0
    if "PSW" in band: beam = 18.0/60.0
    if "160" in band: beam = 12.0/60.0
    if "100" in band: beam = 8.0/60.0
    
    r = np.sqrt(opa**2 + beam**2)
    if r < 0.5: return 1.0/60.0
    if r > 2.5: return 2.5/60.0
    return r/60.0
    
def beamy(band):
    #find band
    if "PLW" in band: beam = 36.0/3600.0
    if "PMW" in band: beam = 24.5/3600.0
    if "PSW" in band: beam = 18.2/3600.0
    if "160" in band: beam = 13.4/3600.0
    if "100" in band: beam = 9.4/3600.0
    
    return beam
    
def makeCumlative(x):
	temp = 0
	result = []
	for val in x:
		temp = temp + val
		result.append(temp)
	return result
 
 #find if detected
def dtest(x):
    try: 
        x =float(x)
        if x == 0.0: return 0
        else: return 1
    except: return 0
##################### main program #################
#make a list of root folders in FIR
#foldersFIR = os.listdir(rootFIR)
#superFiles = os.listdir(rootSuper)

#loop through each band and make dust contour maps
for band in bands:
    print "starting band ",band
    #loop through each galaxy
    for i in range(0,len(cat)):
        gal = str(cat.field('OBJECT')[i])
        galFolder = pj(rootFIR,gal)
        
        #test if detected
        #
        
        #if gal == "67": #continue
        if dtest(cat.field('F250')[i]) == 0: continue
        #try:
        if True:
            print "Starting ", band," CCC", gal
            os.chdir(galFolder)
            sFIR,detype,stype,Rdust = galParamFIR(cat,i,band)
            #optical dust gal parameters            
            x,y,a,ba,pa = findop(cat,i)
            
            #comment out below only use for hand fitting ellipse
            
            #a = a + 0.25
            #ba = a + 0.25
            #x = 55.03212
            #y=-35.62766
            
            pa = pa+90.0
            print a,ba,pa
            print x,y
            #sourceFIR


            #create figure
            fig = plt.figure(figsize = (8,9.5),facecolor='w',edgecolor='w')
            
            FIRraw = source(galFolder,band,stype)
            FIRback = pj(galFolder,gal+"-"+band+"-bgmap.fits")
            FIRbacksub=pj(galFolder,gal+"-"+band+"-bgsubmap.fits")
            FIRellipse = pj(galFolder,gal+"-"+band+"-ellipses.fits")

            #graphs
            #aper noise profile         
            apNoise = np.loadtxt(pj(galFolder,gal+"-"+band+"-aperfit.csv"))
            fluxes  = np.loadtxt(pj(galFolder,gal+"-"+band+"-profiles.csv"),skiprows=1,delimiter=",")
            
            #split fluxes into various intresting componants to plot
            opRadius = fluxes[:,6]
            surBright = fluxes[:,4]
            radSN =np.abs(fluxes[:,0]/fluxes[:,1])
            fluxCum = fluxes[:,2]
            
            X = np.arange(0.0,10000.0)
            Y = np.arange(0.0,10000.0)
            
            dx = 0.40
            dy = 0.24

            #cord for subplots
            xL = 0.07
            xR = 0.57
            
            y1 = 0.72
            y2 = 0.385
            y3 = 0.05
            

            y5 = 0.02
            
            upFac=1.05
            yt1=y1+dy*upFac
            yt2=y2+dy*upFac
            
            
            
            #text size
            fsize = 8
            space = 3.0/60.0
            


            
            Rdust = float(Rdust)
            
            
            #make title        
            plt.figtext(xL,yt1-0.08,gal, size=14,weight="bold")
            plt.figtext(xL,yt1-0.10,sFIR,size=8,weight="semibold")
            plt.figtext(xL,yt1-0.12,stype,size=8,weight="semibold")

            #make titles for aplpy images
            plt.figtext(xR*1.14,yt1,"FIR sub-image",size=10,weight="roman")
            
            
            
            
            ######## make dust/optical plot #########
            ######## top left ################
            
            

            ######## make raw map plot #########
            ######## top right ################
            
            f2 = ap.FITSFigure(FIRraw,figure=fig,subplot=[xR,y1,dx,dy] )
            f2.show_contour(FIRellipse,levels=[Rdust],colour='green', smooth=1)
            f2.recenter(x,y, rad(a,band))
            f2.show_colorscale(exponent="log")
            f2.axis_labels.set_font(size=8)
            f2.tick_labels.set_font(size=6)
            #f2.set_tick_xspacing(space*0.2)
            #f2.set_tick_yspacing(space*0.2)
            f2.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
            f2.show_markers(x, y,c='r',marker='x')
            #f2.show_markers(xPGC, yPGC, c='g',marker='o')
            f2.add_scalebar(1.0/60.0)
            f2.scalebar.set_label('1 arcmin')
            f2.show_ellipses(x, y, a/60.0, ba/60.0, angle=pa, edgecolor='red', lw=1)
            f2.add_beam(major=beamy(band), minor=beamy(band),angle=0.0)
            #f2.show_beam(major=beamy(band), minor=beamy(band),angle=0.0,fc='red',alpha=0.25)
            
            
            #f3.show_beam(major=beamy(band), minor=beamy(band),angle=0.0,fc='red',alpha=0.25)
            ######## make background subtracted map plot #########
            ######## middle right ################            

            
            ######## surface brightness profile #########
            ######## middle lower left ################  
            f5 = plt.axes([xL,y2,dx,dy])
            f5.set_title("Surface brightness profile",size=10)
            f5.set_xlabel("Radius (arcmin)",size=8)
            f5.set_ylabel("Intensity (Jy/beam)",size=8)
            for tick in f5.xaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            for tick in f5.yaxis.get_major_ticks():
                tick.label.set_fontsize(6) 
            f5.semilogy()
            f5.plot(opRadius,surBright,'k-')
            f5.axvline(x=Rdust, color='k', ls='--')
            f5.axvline(x=a, color='r', ls='--')
            
            ######## noise app #########
            ######## middle lower right ################  
            f6 = plt.axes([xR,y2,dx,dy])
            f6.set_title("Noise-aperture dependency",size=10)
            f6.set_xlabel("Radius (arcmin)",size=8)
            f6.set_ylabel("Noise (Jy)",size=8)
            for tick in f6.xaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            for tick in f6.yaxis.get_major_ticks():
                tick.label.set_fontsize(6) 
            f6.plot(opRadius,apNoise[:len(opRadius)],'k-')
            f6.axvline(x=Rdust, color='k', ls='--')
            f6.axvline(x=a, color='r', ls='--')
                        
            
            ######## SN profile #########
            ######## bottom left ################  
            f7 = plt.axes([xL,y3,dx,dy])
            f7.set_title("Intensity S/N profile",size=10)
            f7.set_xlabel("Radius (arcmin)",size=8)
            f7.set_ylabel("S/N",size=8)
            for tick in f7.xaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            for tick in f7.yaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            f7.plot(opRadius,radSN,'k-')
            f7.axvline(x=Rdust, color='k', ls='--')
            f7.axvline(x=a, color='r', ls='--')
            
            
            ######## intensity cumlative #########
            ######## bottom right ################  
            f8 = plt.axes([xR,y3,dx,dy])
            f8.set_title("Cumulative intensity profile",size=10)
            f8.set_xlabel("Radius (arcmin)",size=8)
            f8.set_ylabel("Intensity (Jy)",size=8)
            for tick in f8.xaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            for tick in f8.yaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            f8.plot(opRadius,makeCumlative(surBright),'k-')
            f8.axvline(x=Rdust, color='k', ls='--')
            f8.axvline(x=a, color='r', ls='--')
            
            #plt.show()
            plt.savefig(pj(outfolder,gal+"-"+band+"-postscript.eps"))
            plt.savefig(gal+"-"+band+"-postscript.eps")
            #plt.show()
            plt.clf()

   			         
        #else: continue
        #except:
        #    print ".............galaxy not measured................"
        #    print ">>>>>>>>>>>>>>>>>", gal,"<<<<<<<<<<<<<<<<<<<<<<<"
        #    print ">>>>>>>>>>>>>>>>>",band,"<<<<<<<<<<<<<<<<<<<<<<<"


print "finished!"
#alert.play()

            
    




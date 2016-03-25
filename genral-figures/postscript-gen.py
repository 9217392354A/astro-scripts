#program to make post script for each dectected galaxy
#Chris Fuller, Feb 2013

import os
import numpy as np  
from os.path import join as pj
import aplpy as ap
import atpy as at
import matplotlib.pyplot as plt
from pygame import mixer
mixer.init() #you must initialize the mixer
alert=mixer.Sound('/Users/chrisfuller/Dropbox/personal/sounds/R2D2e.wav')


#switches
donly  = 1 # set to 1 if only wanting detections made
dustOn = 1 # set to 1 if dust contors to be created also
###################### inputs ######################
rootSuper = "/Users/chrisfuller/Desktop/fornax/superCOSMOS/superCOSMOS-pgc/fits/"
rootFIR = "/Users/chrisfuller/Desktop/fornax/correct_pa/galaxies/"
#rootFIR = "/Volumes/niven.astro.cf.ac.uk/galaxies-initialmask"
outfolder = "/Users/chrisfuller/Desktop/all-detections/"
catfolder = "/Users/chrisfuller/Desktop/fornax/"
#catfolder = "/Volumes/niven.astro.cf.ac.uk/"
try: os.mkdir(outfolder)
except: print "outfolder aready exists..."

#bands = ["PLWmap","PMWmap","PSWmap","pacs160","pacs100"] #insert bands as needed
bands = ["PSWmap"]
#bands = ["pacs100"]
#bands = ["PLWmap","PMWmap","pacs160","pacs100"]
print "reading in cat"

cat = at.Table(pj(catfolder,"HeFoCS-fluxes-120313_re-run.fits"),type='fits')

cat
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
    x = float(f.GRA2000[i])
    y = float(f.GDEC2000[i])
    a = float(f.FULLMAJAX[i])
    b = float(f.FULLMINAX[i])
    pa = float(f.PA[i])
    morph = f.MTYPE_PGC[i]
    altname = f.NAMES[i]
    return x,y,a,b,pa,morph,altname

#finds FIR flux and other suff
def galParamFIR(cat,i,band):  
    #find band
    #print cat.EXTENDEDNESS[i]
    if "PLW" in band: 
        flux = cat.F500[i]
        error = cat.E500[i]
        SN = cat.SN500[i]
        t = cat.EXTENDEDNESS[i][0]
        r = cat.R500[i]
    if "PMW" in band: 
        flux = cat.F350[i]
        error = cat.E350[i]
        SN = cat.SN350[i]
        t = cat.EXTENDEDNESS[i][1]
        r = cat.R350[i]
    if "PSW" in band: 
        flux = cat.F250[i]
        error = cat.E250[i]
        SN = cat.SN250[i]
        t = cat.EXTENDEDNESS[i][2]
        r = cat.R250[i]
    if "160" in band: 
        flux = cat.F160[i]
        error = cat.E160[i]
        SN = cat.SN160[i]
        t = cat.EXTENDEDNESS[i][3]
        r = cat.R160[i]
    if "100" in band: 
        flux = cat.F100[i]
        error = cat.E100[i]
        SN = cat.SN100[i]
        t = cat.EXTENDEDNESS[i][4]
        r = cat.R100[i]
        
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
        gal = str(cat.OBJECT[i])
        galFolder = pj(rootFIR,gal)
        
        #test if detected
        #
        
        if gal != "67": continue
        #if dtest(cat.F250[i]) == 0 and donly==1: continue
        #try:
        if False:
            print "Starting ", band," FCC", gal
            os.chdir(galFolder)
            sFIR,detype,stype,Rdust = galParamFIR(cat,i,band)
            #optical dust gal parameters            
            x,y,a,ba,pa,morph,altname = findop(cat,i)
            
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
            fig = plt.figure(figsize = (8,11.5),facecolor='w',edgecolor='w')
            
            FIRraw = source(galFolder,band,stype)
            FIRback = pj(galFolder,gal+"-"+band+"-bgmap.fits")
            FIRbacksub=pj(galFolder,gal+"-"+band+"-bgsubmap.fits")
            FIRellipse = pj(galFolder,gal+"-"+band+"-ellipses.fits")
            SUPER= pj(rootSuper,gal+'.fits')

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
            
            dx = 0.22
            dy = 0.14            
            #cord for subplots
            xL = 0.20
            xR = 0.55
            
            y1 = 0.7
            y2 = 0.5
            y3 = 0.3
            y4 = 0.1
            y5 = 0.02
            
            upFac=1.05
            yt1=y1+dy*upFac
            yt2=y2+dy*upFac
            
            
            
            #text size
            fsize = 8
            space = 3.0/60.0
            


            
            Rdust = float(Rdust)
            
            
            #make title        
            plt.figtext(0.15,0.95,"FCC"+gal+" "+detype+" - "+stype, size=18,weight="heavy")
            plt.figtext(0.15,0.925,sFIR,size=14,weight="bold")
            plt.figtext(0.15,0.91,"Morphological classification: "+morph, size=8,weight="semibold")
            plt.figtext(0.15,0.898,altname, size=8,weight = 'semibold')
            
            #make titles for aplpy images
            plt.figtext(xL*1.2,yt1,"SUPERcosmos r-band",size=10,weight="roman")
            plt.figtext(xR*1.14,yt1,"Raw image",size=10,weight="roman")
            plt.figtext(xL*1.29,yt2,"Background fit",size=10,weight="roman")
            plt.figtext(xR*0.97,yt2,"Background-subtracted image",size=10,weight="roman")
            
            
            
            
            
            ######## make dust/optical plot #########
            ######## top left ################
            #f1 = plt.subplot(6,4,1)
            f1 = ap.FITSFigure(SUPER,figure=fig,subplot=[xL,y1,dx,dy] )
            f1.show_contour(FIRellipse,levels=[Rdust],colour='green', smooth=1)
            f1.show_contour(FIRellipse,levels=[a],colour='red',smooth=1)
            f1.recenter(x,y, rad(a,band))
            f1.show_grayscale(vmin=4000.0,vmax =10000.0)
            f1.add_scalebar(1.0/60.0)
            f1.scalebar.set_label('1 arcmin')
            f1.show_grid()
            f1.axis_labels.set_font(size=8)
            f1.tick_labels.set_font(size=6)
            f1.set_tick_xspacing(space)
            f1.set_tick_yspacing(space)
            f1.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
            f1.show_markers(x, y,c='r',marker='x')
            #f1.show_markers(xPGC, yPGC, c='g',marker='o')
            f1.show_ellipses(x, y, a/60.0, ba/60.0, angle=pa, edgecolor='red', lw=1)
            

            ######## make raw map plot #########
            ######## top right ################
            
            f2 = ap.FITSFigure(FIRraw,figure=fig,subplot=[xR,y1,dx,dy] )
            f2.show_contour(FIRellipse,levels=[Rdust],colour='green', smooth=1)
            f2.recenter(x,y, rad(a,band))
            f2.show_colorscale(exponent="log")
            f2.axis_labels.set_font(size=8)
            f2.tick_labels.set_font(size=6)
            f2.set_tick_xspacing(space*2.0)
            f2.set_tick_yspacing(space*2.0)
            f2.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
            f2.show_markers(x, y,c='r',marker='x')
            #f2.show_markers(xPGC, yPGC, c='g',marker='o')
            f2.add_scalebar(1.0/60.0)
            f2.scalebar.set_label('1 arcmin')
            f2.show_ellipses(x, y, a/60.0, ba/60.0, angle=pa, edgecolor='red', lw=1)
            f2.add_beam(major=beamy(band), minor=beamy(band),angle=0.0)
            f2.show_beam(major=beamy(band), minor=beamy(band),angle=0.0,fc='red',alpha=0.25)
            
            ######## make background map plot #########
            ######## middle left ################            
            f3 = ap.FITSFigure(FIRback,figure=fig,subplot=[xL,y2,dx,dy] )
            f3.show_contour(FIRellipse,levels=[Rdust],colour='green', smooth=1)
            f3.show_colorscale()
            f3.axis_labels.set_font(size=8)
            f3.tick_labels.set_font(size=6)
            f3.set_tick_xspacing(space*2.0)
            f3.set_tick_yspacing(space*2.0)
            f3.add_scalebar(1.0/60.0)
            f3.scalebar.set_label('1 arcmin')
            f3.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
            f3.add_beam(major=beamy(band), minor=beamy(band),angle=0.0)
            f3.show_beam(major=beamy(band), minor=beamy(band),angle=0.0,fc='red',alpha=0.25)
            ######## make background subtracted map plot #########
            ######## middle right ################            
            f4 = ap.FITSFigure(FIRbacksub,figure=fig,subplot=[xR,y2,dx,dy] )
            f4.show_contour(FIRellipse,levels=[Rdust],colour='green', smooth=1)
            f4.show_colorscale()
            f4.axis_labels.set_font(size=8)
            f4.tick_labels.set_font(size=6)
            f4.set_tick_xspacing(space*2.0)
            f4.set_tick_yspacing(space*2.0)
            f4.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
            f4.add_scalebar(1.0/60.0)
            f4.scalebar.set_label('1 arcmin')
            f4.show_ellipses(x, y, a/60.0, ba/60.0, angle=pa, edgecolor='red', lw=1)
            f4.add_beam(major=beamy(band), minor=beamy(band),angle=0.0)
            f4.show_beam(major=beamy(band), minor=beamy(band),angle=0.0,fc='red',alpha=0.25)
            
            ######## surface brightness profile #########
            ######## middle lower left ################  
            f5 = plt.axes([xL,y3,dx,dy])
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
            f6 = plt.axes([xR,y3,dx,dy])
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
            f7 = plt.axes([xL,y4,dx,dy])
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
            f8 = plt.axes([xR,y4,dx,dy])
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
            plt.clf()

            
            #################### Diagnostics image - Dust contors ect
            
            if dustOn ==1:
                #create figure
                fig2 = plt.figure(figsize = (8,9),facecolor='w',edgecolor='w')
                
                #make title        
                plt.figtext(0.05,0.95,"FCC"+gal+" "+detype+" - "+stype, size=18, weight="heavy")
                plt.figtext(0.05,0.91,sFIR,size=14,weight="bold")
                plt.figtext(0.05,0.88,"Morphology = "+morph, size=14,weight = "bold")
                #plt.figtext(0.05,0.85, location_s(cat,pgc,i), size=8,weight = "semibold")
                plt.figtext(0.05,0.85,altname, size=8,weight = 'semibold')
            
                        ######## make dust/optical plot #########
                        ######## top left ################
                #f1 = plt.subplot(6,4,1)
                fa = ap.FITSFigure(SUPER,figure=fig2,subplot=[0.1,0.20,0.8,0.6])
                fa.show_contour(FIRellipse,levels=[Rdust],colour='green', smooth=1)
                fa.recenter(x,y, rad(a,band))
                fa.show_grayscale(vmin=4000.0,vmax =10000.0)
                fa.show_contour(FIRraw,levels=25) # for optical
                fa.add_scalebar(1.0/60.0)
                fa.scalebar.set_label('1 arcmin')
                fa.show_grid()

                #f1.set_tick_xspacing(space)
                #f1.set_tick_yspacing(space)
                fa.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
                fa.show_markers(x, y,c='r',marker='o')
                #fa.show_markers(xPGC, yPGC, c='g',marker='o')
                fa.show_ellipses(x, y, a/60.0, ba/60.0, angle=pa, edgecolor='red', lw=1)
                fa.add_beam(major=beamy(band), minor=beamy(band),angle=0.0)
                fa.show_beam(major=beamy(band), minor=beamy(band),angle=0.0,fc='red',alpha=0.25)
                fa.axis_labels.set_font(size=8)
                fa.tick_labels.set_font(size=6)
                plt.savefig(pj(outfolder,gal+"-"+band+"-dustcontour.eps"))
                plt.clf() 
                #plt.show()            
   			         
        #else: continue
        #except:
        #    print ".............galaxy not measured................"
        #    print ">>>>>>>>>>>>>>>>>", gal,"<<<<<<<<<<<<<<<<<<<<<<<"
        #    print ">>>>>>>>>>>>>>>>>",band,"<<<<<<<<<<<<<<<<<<<<<<<"


print "finished!"
alert.play()

            
    




#detection rate caculator, Chris Fuller March 2013

import os
import numpy as np  
from os.path import join as pj
import atpy as at
import pylab as pl

###################### inputs ######################
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/'
#figFolder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/figures/"

#folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/catalogue_creation_phase/'
figFolder = '/Users/chrisfuller/Desktop/'

#figFolder = "/Users/chrisfuller/Desktop/final_cat/"
saveFig = False #set to true if fig to be outwritten
print "reading in cats"
cats = [s for s in os.listdir(folder) if ".fits" in s]
bands = ["F500","F350","F250","F160","F100"]
#bands = ["F350"]

for catName in cats:
    #create empty lists to hold detection numbers
    t= []
    d = []
    p = []
    master = []
    cat = at.Table(pj(folder,catName),type='fits')
    print " "*70
    print "-"*70
    print "Detection rates: " + catName
    print " band  | Detected | Not Detected | Not Measured | detection rate "
    print "-"*70
    for band in bands:
        notDet = 0
        Det = 0
        notMeasured = 0

        for i in range(0,len(cat)):
            #if cat.HEVICS_PLW[i] == 1 and band == "F100": cat[band][i]= -1.0
            #if cat.HEVICS_PLW[i] == 1 and band == "F160": cat[band][i]= -1.0            
            #get flux
            try: 
                flux = float(cat[band][i])
                if flux == 0.0: notDet += 1# not detected
                elif flux > 0.0: Det +=1
                else: notDet +=1                    # not det detected
            except:
                notDet += 1                 # not measured
        total = Det + notDet
        try:   percentD = np.round((float(Det)/float(total))*100.0, decimals=2)
        except:percentD = 0.0
        
        #append lists
        t.append(total)
        d.append(Det)
        p.append(percentD)
        print str(band)+"um |    "+str(Det) + "    |      " + str(notDet)+"     |      " + str(notMeasured)+"       |        " + str(percentD)+"%  "

    tot = pl.bar(np.arange(1,len(bands)+1),t,align = 'center',color ="grey",width=0.6,log=True, alpha = 1.0) #totals
    det = pl.bar(np.arange(1,len(bands)+1),d,align = 'center',color ="red",width=0.6,log=True, alpha = 1.0) #detected
    pl.xticks(np.arange(1,len(bands)+1),bands)
    
    for j in range(0,len(bands)): 
        pl.text(j+1,100,str(np.round(p[j],decimals=1))+"%",horizontalalignment='center')
        pl.text(j+1,d[j]+2,str(d[j]),horizontalalignment='center')
        pl.text(j+1,t[j]+10,str(t[j]),horizontalalignment='center')
        
    pl.xlabel("Herchel Band")
    pl.ylabel("Galaxies")
    pl.legend([tot,det],["Total","Detected"])
    if saveFig: 
        pl.savefig(pj(figFolder,"detection-rates-"+catName.split(".")[0]+".eps"))
        pl.show()    
        
    
    
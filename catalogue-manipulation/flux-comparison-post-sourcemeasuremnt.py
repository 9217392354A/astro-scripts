#flux comparsion, Chris Fuller March 2013

import os
import numpy as np  
from os.path import join as pj
import atpy as at
import matplotlib.pyplot as plt
#import scipy as sp

###################### inputs ######################
folder =  "/Users/chrisfuller/Desktop/fornax/"

print "reading in cats"
cats = [s for s in os.listdir(folder) if "all_bands-a0p5" in s]
bands = ["F500","F350","F250","F160","F100"]
#bands = ["F350"]

catA = at.Table(pj(folder, "HeFoCS-fluxes-220313-mybgsub.fits"),type='fits')
catB = at.Table(pj(folder, "HeFoCS-fluxes-220313.fits"),type='fits')

for band in bands:
    fluxA = catA[band]
    fluxB = catB[band]
    newFA = []
    newFB = []
    for i in range(0,len(catA)):
        try: 
            fA = float(fluxA[i])
            newFA.append(float(fluxA[i]))
        except: continue
        try: 
            fB = float(fluxB[i])
            newFB.append(float(fluxB[i]))
        except: continue
    # Find and plot 1st order line of best fit 
    A = np.array(newFA,dtype=float)
    B = np.array(newFB,dtype=float)
    coeff = np.polyfit( A, B, 1 ) 
    p = np.poly1d( coeff ) 
    x = np.arange( 0, np.around(max(A)), 0.01 ) 
    plt.plot( x, p(x), label=band + " m = "+ str(np.round(coeff[0],decimals=3)) + " c = "+ str(np.round(coeff[1],decimals=3)) ) 
    plt.plot(A,B,'x', label =band)
#plt.plot(0.03,0.03,"o", label = "30mJy ~ 3 sigma noise level")    
# Add titles/labels and a legend to the graph 
plt.xlabel( "Flux using initial masks (Jy)") 
plt.ylabel( 'Flux using newmasks (Jy)' ) 
plt.legend( loc='best' ) 
plt.grid()
plt.loglog()
plt.show()

#plt.hist(np.log10(A),bins=10, log=True)
#plt.hist(np.log10(B),bins=10, log=True)

plt.show()
    

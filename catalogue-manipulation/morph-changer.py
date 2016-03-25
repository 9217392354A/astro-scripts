#program to change morpholigies into 4 bin method
# Chris Fuller,Mar 2013 

import numpy as np  
from os.path import join as pj
import atpy as at
import pylab as pl
################ inputs ################
figFolder = "/Users/chrisfuller/Dropbox/phd/herchel/fonax/figures/"
catFolder = "/Users/chrisfuller/Dropbox/phd/herchel/fonax/catalogues/"
cat = at.Table(pj(catFolder,"HeFoCS-FLUXES-190313-fcc-matched.fits"),type='fits')
outname = "HeFOCS-FLUXES-190313-morphology.fits"

m4 = ["dE","dS","Im","d:","SBmIII","SmIV"]
m3 = ["Sc","Sd","SBbc","SBbd","SBcdIII"]
m2 = ["Sa","Sb","SBa","SBbc"]
m1 = ["E","S0","SB0"]

############# Functions ###############
#check if any of the sub - strings are in the string
def sis(Slist,s):
    for substring in Slist:
        if substring in s: return True
def dtest(flux):
    try: 
        flux = float(flux)
        if flux == 0.0: return 0
        else: return 1 
    except:
        print "error with flux: ", flux
        return 0
############# Control ################## 
#initalise counters
t1,t2,t3,t4 =0,0,0,0
d1,d2,d3,d4 =0,0,0,0
#add empty col to hold new morphology
cat.add_empty_column('MORPH', np.int16)   
#loop through cat and add new morphology to new col
for i in range(0,len(cat)):
    morph = cat.MType[i]
    if   sis(m4,morph): 
        cat.MORPH[i] = 4
        t4 += 1
        d4 += dtest(cat.F250[i])
    elif sis(m3,morph): 
        cat.MORPH[i] = 3
        t3 += 1
        d3 += dtest(cat.F250[i])
    elif sis(m2,morph): 
        cat.MORPH[i] = 2
        t2 += 1
        d2 += dtest(cat.F250[i])
    elif sis(m1,morph): 
        cat.MORPH[i] = 1
        t1 += 1
        d1 += dtest(cat.F250[i])
    else: print "error.......................", morph 

#write out cat
try:cat.write(pj(catFolder,outname))
except: print "cat already outwritten force overwrite?"
######## plot ############


tot=pl.bar([1,2,3,4],[t1,t2,t3,t4],align = 'center',color ="grey",width=0.6,log=True, alpha = 1.0) #totals
det=pl.bar([1,2,3,4],[d1,d2,d3,d4],align = 'center',color ="red",width=0.6,log=True, alpha = 0.8) #detections

pl.xticks([1,2,3,4],('E/S0','Sa/Sb','Sc/Sd','dwarf/irr/late'))
pl.text(1,0.5,str(np.round(d1*100.0/t1,decimals=1))+"%",horizontalalignment='center')
pl.text(2,0.5,str(np.round(d2*100.0/t2,decimals=1))+"%",horizontalalignment='center')
pl.text(3,0.5,str(np.round(d3*100.0/t3,decimals=1))+"%",horizontalalignment='center')
pl.text(4,0.5,str(np.round(d4*100.0/t4,decimals=1))+"%",horizontalalignment='center')

#number detected
pl.text(1,d1,str(d1),horizontalalignment='center')
pl.text(2,d2,str(d2),horizontalalignment='center')
pl.text(3,d3,str(d3),horizontalalignment='center')
pl.text(4,d4,str(d4),horizontalalignment='center')

#number in sample
pl.text(1,t1,str(t1),horizontalalignment='center')
pl.text(2,t2,str(t2),horizontalalignment='center')
pl.text(3,t3,str(t3),horizontalalignment='center')
pl.text(4,t4,str(t4),horizontalalignment='center')

pl.xlabel("Morphological Type")
pl.ylabel("Galaxies")
pl.legend([tot,det],["Total","Detected"])
pl.savefig(pj(figFolder,"morphology-detections.eps"))
pl.show()


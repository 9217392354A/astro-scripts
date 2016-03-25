# morphology plotter

#program to change morpholigies into 4 bin method
# Chris Fuller,Mar 2013 
#version 2

import numpy as np  
from os.path import join as pj
import atpy as at
import pylab as pl

################ inputs ################
figFolder = "/Users/chrisfuller/Desktop/working/"
catFolder = "/Users/chrisfuller/Desktop/working/"
catA = at.Table(pj(catFolder,"morphology-added_virgo.fits"),type='fits')
catB = at.Table(pj(catFolder,"morphology-added_fornax.fits"),type='fits')
cats = [catA,catB]
names= ['Virgo','Fornax']

mList = [1,2,3,4]
mLables=['E/S0','Sa/Sb','Sc','/irr/late/dwarf']


############# Functions ###############
#check if any of the sub - strings are in the string
def myMatch(substring, string):
    if substring == string: return True
    else: return False     
    
def stringy(x,h):
    return str(np.int(np.round(x,decimals=h)))

def whichBin(binlist, Mtype):
    #which bin?
    #mtype == morph
    #binlist == mList
    Bnum = -1
    for BIN in binlist:
        Bnum +=1
        for inBin in BIN:
            #does our galaxy match this bin?
            if myMatch(inBin,Mtype):
                return Bnum
    print "No match found....",Mtype
    return 0
            
            
def dtest(flux):
    try: 
        flux = float(flux)
        if flux == 0.0: return 0
        elif flux > 0.0: return 1
        else: return 0
    except:
        print "error with flux: ", flux
        return 0
        
def pCal(x,y):
   try:   percentD = np.round((float(x)/float(y))*100.0, decimals=1)
   except:percentD = 0.0
   return percentD
   
def e2c():
    try:
        input("press enter to continue")
    except: pass
def zoomorph(c,n):
    if c.SPIRAL[n] == 1:      return 1
    elif c.ELLIPTICAL[n] ==1: return 2
    else: return 3

############# Control ################## 


al = [np.arange(1,len(mList)+1)*2-1, np.arange(1,len(mList)+1)*2]
colours = [["0.8",'r'],['0.8','b']]
for j in range(0,len(cats)):
    #if j == 1: continue
    cat = cats[j]
    p = np.zeros(shape=(len(mList)))    
    t = np.zeros(shape=(len(mList)),dtype=int)
    d = np.zeros(shape=(len(mList)),dtype=int)
    for i in range(0,len(cat)):
        morph = cat.MORPH[i]
        Gflux = cat.F250[i]
        binNo = int(morph)-1
        t[binNo] += 1
        d[binNo] += dtest(Gflux)
        cat.MORPH[i] = binNo + 1
        #print mList[binNo], morph
        #e2c() #uncomment to check each as the role through


    # if you want it to be normed by total detections hit true
    daaa =d
    taaa = t
    error = daaa/np.sqrt(daaa*1.0)
    if True: 
        norm = np.float(np.sum(t))
        t = t / norm
        d = d / norm
    percent = d*100.0/t
    
    pl.bar(al[j],t,align = 'center',width=0.8,log=False, alpha = 1.0, color =colours[j][0]) #totals
    pl.bar(al[j],d,align = 'center',width=0.8,log=False, alpha = 1.0, label=names[j], color=colours[j][1]) #detections
    for k in range(0,len(mList)): 
        pl.text(al[j][k],d[k] + 0.005,stringy(percent[k],0)+ "%",horizontalalignment='center')
        pl.text(al[j][k], 0.002,stringy(taaa[k],0),horizontalalignment='center',color ='white')

pl.xticks(al[0]+0.5,mLables)

pl.ylim(ymax=0.2)
pl.xlabel("Morphological Type")
pl.ylabel("Normalised Number of Galaxies")
pl.legend(loc=2)
pl.savefig(pj(figFolder,"morphology-detections-virgo-fornax-norm-0.2.eps"))
pl.show()



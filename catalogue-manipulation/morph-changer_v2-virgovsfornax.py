#program to change morpholigies into 4 bin method
# Chris Fuller,Mar 2013 
#version 2

import numpy as np  
from os.path import join as pj
import atpy as at
import pylab as pl
import re

################ inputs ################
figFolder = "/Users/chrisfuller/Desktop/det_fig"
catFolder = "/Users/chrisfuller/Desktop/det_fig"
catName= "cat.fits"
cat = at.Table(pj(catFolder,catName),type='fits')
#cat = at.Table(pj(catFolder,"virgo.fits"),type='fits')

outname = "morphology-added_"+catName


m1 = ["E","S0", "SB0"]
m2 = ["Sa","Sb","SBa"]
m3 = ["Sc","SBbc"]
m4 = ["Im","Sd","SBmIII","SmIV","dE","dS","d:","SBbd","SBcdIII"]

mList = [m1,m2,m3,m4]
mLables=['E/S0','Sa/Sb','Sc','/irr/late/dwarf']


############# Functions ###############
#check if any of the sub - strings are in the string
def myMatch(substring, string):
    result = re.match(substring,string)
    if result == None: return False #no match found
    else: return True      
    
    
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
        else: return 1 
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

############# Control ################## 
p = np.zeros(shape=(len(mList)))    
t = np.zeros(shape=(len(mList)),dtype=int)
d = np.zeros(shape=(len(mList)),dtype=int)
#add empty col to hold new morphology
cat.add_empty_column('MORPH', str) 
#loop through whole cat:

for i in range(0,len(cat)):
    morph = str(cat.MType[i])
    Gflux = cat.F250[i]
    binNo = whichBin(mList,morph)
    t[binNo] += 1
    d[binNo] += dtest(Gflux)
    cat.MORPH[i] = binNo + 1
    #print mList[binNo], morph
    #e2c() #uncomment to check each as the role through



#caculate %
for k in range(0,len(mList)): p[k] = pCal(d[k],t[k])

#write out cat
cat.write(pj(catFolder,outname),overwrite=True)
#except: print "cat already outwritten force overwrite?"
######## plot ############


tot=pl.bar(np.arange(1,len(mList)+1),t,align = 'center',color ="grey",width=0.6,log=True, alpha = 1.0) #totals
det=pl.bar(np.arange(1,len(mList)+1),d,align = 'center',color ="red",width=0.6,log=True, alpha = 0.8) #detections

pl.xticks(np.arange(1,len(mList)+1),mLables)

for j in range(0,len(mList)): 
        pl.text(j+1,1.20,str(np.round(p[j],decimals=1))+"%",horizontalalignment='center')
        pl.text(j+1,d[j],str(d[j]),horizontalalignment='center')
        pl.text(j+1,t[j],str(t[j]),horizontalalignment='center')

pl.xlabel("Morphological Type")
pl.ylabel("Log Number of Galaxies")
pl.legend([tot,det],["Total","Detected"])
pl.savefig(pj(figFolder,"morphology-detections-"+catName.split(".")[0]+".eps"))
pl.show()


# program to make detection rates of an input catalogue and put them into a figure for a 
# paper
# Chris Fuller, May 2013

import numpy as np  
from os.path import join as pj
import atpy as at
import matplotlib.pyplot as plt
import pylab as pl



############### inputs #######################
folder = "/Users/chrisfuller/Dropbox/phd/plots/det_fig"
fornax = at.Table(pj(folder,"morphology-added_cat.fits"))
virgo = at.Table(pj(folder,"virgo_fluxes.fits"))
#bands = ["F500","F350","F250","F160","F100"]
bands =["F100","F160","F250","F350","F500"]
#mLables=['','1','','2','','3','','4']
mList=['E/S0','Sa/Sb/Sc','/irr/late/dwarf']
############# functions ######################
#create a col of optical mags for the detected galaxies
def rtn_det(f,op):
    return op[np.where(f>0.0)]
    
#find the morphologies and how many are in each bin
def w(cat,x):
    return len(np.where(cat[mor]==x)[0])
    
def dw(cat,m,f):
    return len(np.where((cat[mor]==m) & (f>0.0))[0])
    
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
    
############ control ########################

#create morphology plot

#work out each bin for virgo and fornax

flux_fornax , flux_virgo = np.nan_to_num(fornax["F250"]) , np.nan_to_num(fornax["F250"])









totals = [w("1"),w("2")+w("3"),w("4")]
det = [dw("1",flux),dw("2",flux)+dw("3",flux),dw("4",flux)]



"""
pl.bar(np.arange(1,4),np.log10(totals),align = 'center',color ="black",width=1.0, alpha = 1.0) #totals
pl.bar(np.arange(1,4),np.log10(det),   align = 'center',color ="purple",width=1.0,alpha = 1.0) #det
#pl.semilogy()
#f2.tick_params(axis='x', labelbottom='off')
#f2.set_ylim(1,380)   
pl.ylabel("N")  
pl.xticks(np.arange(1,len(mList)+1),mList)
#pl.yticks(np.arange(1,4),['0','1','10','100'])

d = np.array(det)
t = np.array(totals)
p = d*100.0/t
for j in range(0,len(mList)): 
        pl.text(j+1,0.20,str(np.round(p[j],decimals=1))+"%",horizontalalignment='center',color = 'white')
        pl.text(j+1,d[j],str(d[j]),horizontalalignment='center',color = 'white')
        pl.text(j+1,t[j],str(t[j]),horizontalalignment='center')

#pl.ylim(ymin=1.0)
pl.xlabel("Morphological Type")
#f2.tick_params(axis='both',labelleft='on', labelbottom='on')
#.savefig(pj(folder,"detection_rates_morph.eps"))
pl.show()
"""
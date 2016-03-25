# flux histogram creator, opens cataloges all droped in to folder and plots them as a hisogram
# but changing them all to the same distance
# Chris Fuller, April 13


import numpy as np  
from os.path import join as pj
import atpy as at
import pylab as pl

folder = "/Users/chrisfuller/Desktop/working/"
#import cataloges
coma =   at.Table(pj(folder,"coma.fits"),type='fits')
fornax = at.Table(pj(folder,"fornax.fits"),type='fits')
virgo =  at.Table(pj(folder,"virgo.fits"),type='fits')
filament =  at.Table(pj(folder,"filament.fits"),type='fits')

#names = ['coma','fornax','virgo']
#cats = [coma,fornax,virgo]

names = ['virgo','coma','fornax', 'filament']
cats = [virgo,coma,fornax,filament]

distance = [16.0,100.0,19.0,100.0]
fluxD = 16.0
####### functions ########
def Dtest(x):
    if x == "" or x == " ": return False
    elif float(x) == 0.0: return False
    elif float(x) > 0.0: return True
    else: return False


def emt2zero(f):
    for i in range(0,len(f)):
        if Dtest(f[i]): continue
        else: f[i] = 0.0
    return f
    
def rmempt(f):
    newF = []
    for i in range(0,len(f)):
        if Dtest(f[i]):
            newF.append(f[i])
    return np.array(newF)
            
######## control ##########

fluxes = []

for j in range(0,len(cats)):
    #if names[j]!= 'fornax': continue
    cat = cats[j]
    flux =[]
    #flux = emt2zero(cat.F250)
    flux =rmempt(cat.F250)
    
    #scale to distance of fluxD
    sfactor = (distance[j]/fluxD)**2
    
    newflux = flux*sfactor#*1000.0
    fluxes.append(newflux)

    pl.hist(newflux, bins=np.logspace(-2.0, 2.2, 15),normed=True,label=names[j],alpha=1.0,log=True,histtype='step',align='mid',lw=4.0)
    print names[j],' mean: ',np.mean(newflux),' min: ',np.min(newflux),' max: ',np.max(newflux)
pl.gca().set_xscale("log")
pl.legend()
pl.grid()
pl.xlabel("250um flux at the distance of Virgo, (Jy)")
pl.ylabel("Normalised count")

pl.savefig(pj(folder,"fluxes_at_virgo_all.eps"),dpi=600)
pl.show()
    
    
    
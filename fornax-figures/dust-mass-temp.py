#make chi squared plots
# Chris Fuller, May 2013

import pylab as pl 
from atpy import Table
#from os.path import join as pj
from os import chdir
import numpy as np
import math as m
import matplotlib.pyplot as plt

#inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/SED-fits/"
chdir(folder)
cat = Table('cat.fits',type='fits')



######
chi = cat.CHISQ
mass = cat.MASS
temp = cat.T
mtype = cat.Morphology
######

### clean ###
clean = np.where(chi<15.0)

chi = chi[clean]
mass = mass[clean]
temp = temp[clean]
mtype = mtype[clean]

############# functions #############
def chi2pdf(xdata,k,norm):
    y=[]
    for x in xdata:
        y.append(norm*(1.0 / (2.0*m.gamma(k/2)) * (x/2.0)**(k/2.0-1.0) * m.exp(-x/2.0)))
    return y

########### control ##################


#sep into early and late
early = np.where(mtype=="E")
late  = np.where(mtype=="L")

if True and False:
    #chi pdf
    xchi=np.arange(0,np.ceil(np.max(chi)),0.01)
    ychi = chi2pdf(xchi,3.0,16.0)
    
    b=np.arange(0,np.ceil(np.max(chi)),0.57)
    #chisq plot
    pl.hist(chi[early],bins=b,histtype='step',color='red',lw=4.0)
    pl.hist(chi[late],bins=b, histtype='step',color='blue',lw=4.0)
    pl.plot(xchi,ychi,'k')
    pl.xlabel('$\chi^{2}$')
    pl.ylabel('N')
    pl.savefig("chisq_dis.pdf")
    
    pl.show()

#mass and temp
bm=np.arange(np.floor(np.min(mass)),np.ceil(np.max(mass)),0.5)

fig = plt.figure(figsize = (4.0,7.5),facecolor='w',edgecolor='w')

f1 = plt.subplot(2,1,1)
f2 = plt.subplot(2,1,2)

#mass 
f1.hist(mass[early],bins=bm,histtype='step',color='red',lw=2.0)
f1.hist(mass[late],bins=bm,histtype='step',color='blue',lw=2.0)
f1.axvline(np.mean(mass[early]),color='red',ls='--')
f1.axvline(np.mean(mass[late]),color='blue',ls='--')
f1.set_xlim(4,10)
f1.set_ylim(0,5)
f1.set_xlabel('log($M_{Dust}$/$M_{\odot}$)' )
f1.set_ylabel('N')

bt = np.arange(np.floor(np.min(temp)),np.ceil(np.max(temp)),3.0)
#mass 
f2.hist(temp[early],bins=bt,histtype='step',color='red',lw=2.0)
f2.hist(temp[late],bins=bt,histtype='step',color='blue',lw=2.0)
f2.axvline(np.mean(temp[early]),color='red',ls='--')
f2.axvline(np.mean(temp[late]),color='blue',ls='--')
f2.set_xlim(6,26)
f2.set_ylim(0,4)
f2.set_xlabel('Dust Temp. (K))' )
f2.set_ylabel('N')
ya = f2.get_yaxis()
ya.set_major_locator(plt.MaxNLocator(integer=True))

plt.subplots_adjust(hspace=0.3,wspace=0.0)
fig.savefig('dust_mass_temp.eps')

print 'temp range early: ',np.min(temp[early]),np.max(temp[early]),np.mean(temp[early])
print 'temp range late: ',np.min(temp[late]),np.max(temp[late]),np.mean(temp[late])

print 'mass range early: ',np.min(mass[early]),np.max(mass[early]),np.mean(mass[early])
print 'mass range late: ',np.min(mass[late]),np.max(mass[late]),np.mean(mass[late])
fig.show()


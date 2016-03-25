# making number counts function with the aim of showing no excess galaxies not found with the optical
# catalogue
# Chris Fuller, May 2013

import numpy as np
from atpy import Table
import matplotlib.pyplot as plt
from os.path import join as pj
import pylab as pl


folder = '/Users/chrisfuller/Dropbox/phd/plots/SoureCounts/tables'

object_size = '5' #pixels
threshold = '1p5' #sextractor detection threshold

#import cats
#hevics = Table(pj(folder,'Hevics_os' + object_size + '_' + threshold + '.FIT'), type='fits')
#hefocs = Table(pj(folder,'Hefocs_os' + object_size + '_' + threshold + '.FIT'), type='fits')
#ngp    = Table(pj(folder,'ngp_os'    + object_size + '_' + threshold + '.FIT'), type='fits')

blank  = Table(pj(folder, 'ngp_nw_1p24v2.FIT'), type='fits')
hefocs = Table(pj(folder, '1.24_fornax.fits'))


#import cluster cats
#virgo = Table(pj(folder,'virgo.fits'))
fornax = Table(pj(folder, 'fornax.fits'))

#cats holds all the ...cats?
cats = [blank, hefocs, fornax]
name = ['blank', 'hefocs map', 'fornax cluster']

############################# functions ###################################################
#removes nan
def remove_nan(x):
    #w1 = np.where(np.nan_to_num(x) == 0)[0]
    #return x[w1]
    return np.nan_to_num(x)
############################# control  ####################################################
#turn cats into fluxes (Jy)
X_jy = [ 25.0/423.0, 36.0/423.0]

#X_jy = 25.0 /  (423.0)

X_sr = 1.0 / 3282.80635 #converts sq degree's to str^-1

bins_a = np.arange(np.log10(15*10**-3),2,0.25)

lins=['-','-','--'] 
cols = ['k','r','aqua']


fig = plt.figure(figsize = (8.,8.),facecolor='w',edgecolor='w')
f1 = plt.subplot(1,1,1)
#f2 = plt.subplot(2,1,2)
# loop through cats and create new list of fluxes
for i in range(0, len(cats)):
    try:
        flux = cats[i]['FLUX_AUTO'] * X_jy[i]
        #clean all that are less than 12
        Num = 12.0 
        w_size = np.where(cats[i]['ISOAREA_IMAGE'] > Num)
        flux = flux[w_size]
    except:
        flux = cats[i]['F250']
    flux = remove_nan(flux)
    flux = np.log10(flux)
    
    #create histogram
    n, bins = np.histogram(flux, bins=bins_a)
    #n, bins, patches = pl.hist(flux, bins=bins_a, histtype = 'step', alpha=0.0)

    if i == 0: 
        n_sr = n / (61.06*0.00030461741978670857)
        print name[i]
    else: 
        n_sr = n  / (20.0*0.00030461741978670857)

    #plt.plot(bins[:len(bins)-1], n_sr, label=name[i])

    #er = 1.0/np.log10(np.sqrt(n_sr)) 
    
    er = np.sqrt( n ) / n    
    er = np.nan_to_num(er)
    print name[i]
    print 'n'
    print n
    print 'error'
    print er
    print 'n_sr'
    print n_sr
    print 'err n_sr'
    print n_sr*er    
    
    
    
    
    
    x=bins[:len(bins)-1] #+ (bins[1]-bins[0])*0.5
    f1.errorbar(x, n_sr, yerr=(n_sr*er), label=name[i],ls=lins[i], color=cols[i])
    f1.axvline(x=np.log10(15.*10**-3), ls ='--', color='k')
    f1.set_ylabel("$Log_{10}(N)$ ($Sr^{-1}$)")
    f1.semilogy()
    f1.set_yticklabels([1,2,3,4,5,6,7])
    f1.set_xlim(-2.,1)


#plt.legend()

#create risidual plot
#error = np.sqrt(er_blank**2 + er_cluster**2 )#+ er_for**2)
#y = ( n_cluster- n_blank)# - n_for)# / n_cluster

#f2.errorbar(x,np.log10(y), yerr=np.log10(error), ls='-',color='black')
#f2.axhline(0.0,ls='--',color='black')

#f2.set_ylim(2.5,-2.5)
#f1.tick_params(axis='x', labelbottom='off')
#f2.set_xlim(-2,1)
#f2.set_ylabel('Residual %')
f1.set_xlabel("$Log_{10}(Flux)$ (Jy)")
plt.subplots_adjust(bottom=0.07, left=0.08, right=0.98, top=0.98)

fig.savefig(pj('/Users/chrisfuller/Desktop/','numbercounts.pdf'))
plt.show()




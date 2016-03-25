#program to genorate stellar masses from optical colours and 2mass magnitudes. For
# galaxies that only have a bt mag then we will fit a function to the 
# other galaxies giving a fuction for btmag to stellar mass.

# Chris Fuller August  2013

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/plots/stellar_mass/"
cat = Table(pj(folder,"stellar_mass.fits"))

log_mass = cat['SMASS'] 
morphology = cat['MTYPE_PGC_2'] 

dwarf = np.where((morphology >= -3) & (morphology <= -1))[0]
early = np.where((morphology >= 0) & (morphology <= 1))[0]
late  = np.where((morphology >= 2) & (morphology <= 6))[0]
later = np.where((morphology >= 7) & (morphology <= 20))[0]

dwarf_detected = np.where((morphology >= -3) & (morphology <= -1) & (np.nan_to_num(cat['F250'])>0))[0]
early_detected = np.where((morphology >= 0) & (morphology <= 1) & (np.nan_to_num(cat['F250'])>0))[0]
late_detected = np.where((morphology >= 2) & (morphology <= 6) & (np.nan_to_num(cat['F250'])>0))[0]
later_detected = np.where((morphology >= 7) & (morphology <= 20) & (np.nan_to_num(cat['F250'])>0))[0]

types = [dwarf, early, late, later]
detected = [dwarf_detected, early_detected, late_detected, later_detected]
types_names = ['Dwarf', 'Early', 'Late', 'Later']

figM = plt.figure(figsize=(4.5, 8.5), facecolor='w', edgecolor='k')

width = 0.5

binsa = np.arange(6.,12.0,width)
for i in range(0,len(types)):
    subfi = plt.subplot(4,1,i+1)
    hist, bin_edges = np.histogram(log_mass[types[i]], bins=binsa) 
    hist_d, bin_edges = np.histogram(log_mass[detected[i]], bins=binsa) 
    bin_centers = bin_edges[:-1]+ ((bin_edges[1] - bin_edges[0]) / 2.0)
    #plot undetected
    subfi.bar(bin_centers,hist, width=width,align ='center', color='black')
    subfi.bar(bin_centers,hist_d, width=width,align ='center', color='aqua')
    subfi.set_xlim(5.,12.)
    #subfi.set_ylim(0,3)
        
    #plot lines of limit

    subfi.axvline(x=5.6, color='k', ls='--')
    subfi.axvline(x=6.6, color='k', ls='--')
    subfi.axvline(x=7.6, color='k', ls='--')
    subfi.axvline(x=8.6, color='k', ls='--')
    subfi.axvline(x=9.6, color='k', ls='--')
    subfi.axvline(x=10.6, color='k', ls='--')

    if i == 2:
        subfi.text(5.6,8,'10$^1$')
        subfi.text(6.6,8,'10$^2$')
        subfi.text(7.6,8,'10$^3$')
        subfi.text(8.6,8,'10$^4$')
        subfi.text(9.6,8,'10$^5$')
        subfi.text(10.6,8,'10$^6$')

    if i != 0: 
        subfi.set_ylim(0.,12.9)
        hight = 12*0.95
    if i == 0: hight = 60*0.90
    if i != 3: 
        subfi.set_xticklabels([])
        
        
    if i == 3: subfi.set_xlabel('$\log_{10}$$(M_{star} / M_{\odot})$')
    subfi.set_ylabel('N')

    subfi.text(5.3,hight,types_names[i],fontsize=16, weight = 'semibold' )


plt.subplots_adjust(left=0.125,bottom=0.075, top = 0.985, hspace=0.0,wspace=0.0)
#figM.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/fornax','smass-detected.pdf'))
plt.show()


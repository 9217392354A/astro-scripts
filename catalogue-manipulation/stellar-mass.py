#program to genorate stellar masses from optical colours and 2mass magnitudes. For
# galaxies that only have a bt mag then we will fit a function to the 
# other galaxies giving a fuction for btmag to stellar mass.

# Chris Fuller 15 July 2013

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/"
cat = Table(pj(folder,"stellar-mass-fornax_v2.fits"))

###################### Functions ############################################################

# bell et al 2003
# log10 (M/L) - a + (b x color)
def stellar_mass(band, col, sun_band):
    #constants used below    
    a = -0.206
    b =  0.135
    
    #distance modules of fornax
    dist_mod = 31.24
    
    #convernt to abs mag
    Band = band - dist_mod
    
    #convert into solar luminositys
    L = -0.4 * (Band + sun_band)
    
    #find conversion factor for each galaxies based on a,b, and B_V
    ML_con = a + (b * col)
    
    
    MASS = ML_con + L     
    #for i in range(0,len(L)): print MASS[i], ML_con[i] , L[i], Band[i]
    return MASS


######################################################################################





#extract only galaxies with colours and K band magnitude
colors =  cat["B-V"]
K = cat["K_2MASS"]
BTmag = cat["BTmag"]

#find where K and colour are true
w = np.where((np.nan_to_num(K) != 0.0 ) & (np.nan_to_num(colors) != 0.0 ))[0]

#Convert to stellar mass
log_mass = stellar_mass(K[w], colors[w], -3.28)


# fit BTmag to stellar_mass and then find all galaxies that don't have a K band flux and convert
fit = np.polyfit(BTmag[w], log_mass, 1)
BT2Stellar = np.poly1d(fit)

print "m" + str(fit)


#set stellar mass coloumn with all galaxies that have mass's
cat['SMASS'] = 1.0
cat['SMASS'][w] = log_mass

#extract all galaxies without colours and K band magnites
if True:
    for i in range(0, len(cat)):
        #test to see if there is already a stellar mass there
        if cat.SMASS[i] == 1.0: cat.SMASS[i] = BT2Stellar(cat.BTmag[i])
        else:
            print 'mass already present from K band and colour'
        
#load complete cols into mass and BT
mass = cat.SMASS
BT = cat.BTmag

cat.write(pj(folder,'stellar-mass-fornax_final.fits'),overwrite=True)

figM = plt.figure(figsize=(4.5, 4.5), facecolor='w', edgecolor='k')
subfi = plt.subplot(1,1,1)
#plot
subfi.plot(BT, BT2Stellar(BT), 'k')
#subfi.scatter(BT, mass, s=30, c='r', marker='+')
subfi.scatter(BT[w], mass[w], s=30, c='b', marker='+')
#subfi.xaxis_inverted()
subfi.set_ylabel('$\log_{10}$$(M_{star} / M_{\odot})$')
subfi.set_xlabel('$m_{BT}$')
subfi.set_xlim(9.8 ,15.8)
subfi.set_ylim(8.4 ,11.5)
plt.subplots_adjust(left =0.175, bottom= 0.125, hspace=0.0,wspace=0.0)
figM.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/fornax','BT2Stellarmass.pdf'))
plt.show()


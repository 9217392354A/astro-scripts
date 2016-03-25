# making sed fit output
# Chris Fuller, May 2013

import numpy as np
from atpy import Table
import matplotlib.pyplot as plt
from astropy import units as u
from scipy.optimize import fmin



bands = ['100', '160', '250', '350', '500']

#inverse_bands =  np.array([500.,350.,250.,160.,100.], dtype=np.float)

#physical constants
h = 6.626e-34
c = 3.0e+8
k = 1.38e-23

########## functions ###############
#this gives a flux in jy for a wavelength, this uses a fixed beta
#wavelength = microns, distance = Mpc, mass = log10 solar mass
def modb(mass,temp,distance,wavelength):
    m = (10**mass)*2E30 #converts mass into kg
    d = distance*3.08567758E22 #converts Mpc to m
    wav = wavelength*1E-6 #convert micron to m
    
    k_l = kappa(wav)
    b_l = planck(wav,temp) 
    
    return (k_l*m*b_l / d**2)
    
    
    
#finds kappa
def kappa(ll):
    return 0.192*((350.0E-6/ll)**2.0)

def planck(wav, T):
   aa = 2.0*h*c**2
   bb = h*c/(wav*k*T)
   return aa / ( np.expm1(bb) * (wav**5) )


def bbmin(M,*args):
        #ajust for different distance temp ect
        T = 20.
        D = 17.5
        F = args[0]
        new = cjy*modb(M,T,D,250.0) 
        return (new - F)**2.

def flux2mass(flux): 
    # first make a guess out the mass
    return fmin(bbmin, 8., args=(flux,10.), disp=0)
    
    
   

   
########## control #################

#Convert into jys
#convert s into jy
lm3 = u.W*u.m**-2*u.m**-1
cjy = lm3.to(u.Jy, equivalencies=u.spectral_density(u.um,250.0))


flux_list = np.arange(0.010, 3.0, 0.010)
mass_list = []

for s in flux_list:
    mass_list.append(flux2mass(s))

#convert into numpy array
mass_list =np.array(mass_list, dtype=np.float)
flux_list = np.log10(flux_list)


#now fit a straight line to the fit
fit = np.polyfit(flux_list, mass_list,1)
print fit
plt.plot(flux_list, mass_list)
plt.show()











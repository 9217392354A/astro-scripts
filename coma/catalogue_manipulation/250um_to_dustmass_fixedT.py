#program that takes 250um flux assumes a dust temp of 20K and then asigns dust mass based on this
#Chris Fuller Feb 2014

import numpy as np
from atpy import Table
from astropy import units as u
from os.path import join as pj

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
input_cat = 'coma_supercluster_cal12.fits' #input name
output_name= 'coma_supercluster_cal12-test.fits'

#read in catalogue
cat = Table(pj(folder,input_cat))

#column headers
flux_col = 'F250'

redshift = 0.023100
wav = 250.0 #um
T_dust = 20 #K
distance = 100.0 #Mpc

#all my galaxies are at a fixed distance and thus readshift so the program will
#create arrays of these, athought they could be fed into the fuction for each 
# individual galaxy

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#function that takes a flux at a single wavelenght, and returns a dust mass   #
# it requires a redshift, distance, flux, and wavelength_observed			  #
# Genrally T = 20K assumed for normal galaxy								  #
#input wavelength = microns, distance = Mpc, temp = K 						  #
#returns mass = log10 solar MASS                    						  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def flux_to_dust(flux, z, temp, distance, wavelength_observed):
	#convert lamda obs into lamda emmit as well as convert micron to m
	wavelength = (wavelength_observed * 1E-6)   / (z + 1.0)
	d = distance*3.08567758E22 #converts Mpc to m

	#find conversion factor for jy to si
	lm3 = u.W*u.m**-2*u.m**-1
	jy_to_si = u.Jy.to(lm3, equivalencies=u.spectral_density(u.m,wavelength))

	#convert flux in Jy to Si Units
	s = flux*jy_to_si

	#find kappa note beta is fixed at 2
	k_l = kappa(wavelength)

	#find value of plack equation
	b_l = planck(wavelength,temp) 

	#caculate dust mass in kg
	mass = (s * d**2) / (k_l * b_l)

	#convert to log10(Dust/Msol)
	m = np.log10(mass/2E30)

	return m
    
#finds kappa
def kappa(ll):
    return 0.192*((350.0E-6/ll)**2.0)

def planck(wav, T):
	#physical constants
	h = 6.626e-34
	c = 3.0e+8
	k = 1.38e-23

	aa = 2.0*h*c**2
	bb = h*c/(wav*k*T)
	return aa / ( np.expm1(bb) * (wav**5) )
   
########## control #################
"""
#find all detected galaxies
w_detected = np.where(np.nan_to_num(cat[flux_col]) != 0.0)[0]

#get cols ready
flux = cat[flux_col][w_detected]
z = np.array([redshift]*len(flux), dtype=np.float)
wav = np.array([wav]*len(flux), dtype=np.float)
dist = np.array([distance]*len(flux), dtype=np.float)
temp = np.array([T_dust]*len(flux), dtype=np.float)

#caculated dust mass
mass = flux_to_dust(flux, z, temp, dist, wav)

#create new empty column in table
cat.add_empty_column('DMASS_250', np.float, unit='Log10(Msun)', null='', description='Dust Mass caculated from 250um flux of detected galaxies')

#add data to table 
cat['DMASS_250'][w_detected] = mass

#save table
cat.write('/Users/chrisfuller/Desktop/'+ output_name)
"""
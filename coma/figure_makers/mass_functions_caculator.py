# This program is designed to caculate the error in a mass functions cacaluted mass density

from scipy.special import gamma as Gamma
from numpy import round, log10, sqrt

if False:
	###### Coma ######
	#coma stellar mass
	Stellar_phi, Stellar_dphi 			= 2.5, 		0.08
	Stellar_Mstar, Stellar_dMstar		= 89.0E9, 	3.0E9
	Stellar_alpha, Stellar_dalpha		= -1.0,  	0.01

	#coma gas mass
	Gas_phi, Gas_dphi 					= 0.3, 		0.2
	Gas_Mstar, Gas_dMstar				= 4.8E9, 	3.2E9
	Gas_alpha, Gas_dalpha				= -0.7,  	0.4

	#coma stellar mass
	Dust_phi, Dust_dphi 				= 0.4, 		0.2
	Dust_Mstar, Dust_dMstar				= 0.06E9, 	0.03E9
	Dust_alpha, Dust_dalpha				= -0.8,  	0.2


if False:
	###### Filament #########
	#Filament stellar mass
	Stellar_phi, Stellar_dphi 			= 0.05, 	0.01
	Stellar_Mstar, Stellar_dMstar		= 175.0E9, 	48.0E9
	Stellar_alpha, Stellar_dalpha		= -1.2,  	0.1

	#Filament gas mass
	Gas_phi, Gas_dphi 					= 0.06, 	0.01
	Gas_Mstar, Gas_dMstar				= 3.9E9, 	1.5E9
	Gas_alpha, Gas_dalpha				= -0.08,  	0.5

	#Filament stellar mass
	Dust_phi, Dust_dphi 				= 0.04,		0.02
	Dust_Mstar, Dust_dMstar				= 0.1E9, 	0.04E9
	Dust_alpha, Dust_dalpha				= -1.0,  	0.2

if False:
	###### Virgo #####
	#Virgo stellar mass
	Stellar_phi, Stellar_dphi 			= 0.3, 		0.1
	Stellar_Mstar, Stellar_dMstar		= 192.0E9, 	117.0E9
	Stellar_alpha, Stellar_dalpha		= -1.2,  	0.1

	#Virgo gas mass
	Gas_phi, Gas_dphi 					= 0.6, 		0.3
	Gas_Mstar, Gas_dMstar				= 4.5E9, 	1.6E9
	Gas_alpha, Gas_dalpha				= -1.0,  	0.2

	#Virgo stellar mass
	Dust_phi, Dust_dphi 				= 0.7,		0.1
	Dust_Mstar, Dust_dMstar				= 0.06E9, 	0.01E9
	Dust_alpha, Dust_dalpha				= -0.9,  	0.1

if True:
	###### Field #####
	#Field stellar mass
	Stellar_phi, Stellar_dphi 			= 0.002, 		0.0001
	Stellar_Mstar, Stellar_dMstar		= 100.0E9, 	20.0E9
	Stellar_alpha, Stellar_dalpha		= -1.2,  	0.1

	#Field gas mass
	Gas_phi, Gas_dphi 					= 0.009, 	0.001
	Gas_Mstar, Gas_dMstar				= 5.0E9, 	0.1E9
	Gas_alpha, Gas_dalpha				= -1.50,  	0.05

	#Field stellar mass
	Dust_phi, Dust_dphi 				= 0.006,		0.001
	Dust_Mstar, Dust_dMstar				= 4.0E7, 	0.4E7
	Dust_alpha, Dust_dalpha				= -1.0,  	0.2


rho_g_s_virgo = [0.040263, 0.037849]
rho_d_s_virgo = [0.000596, 0.000439]

rho_g_s_coma = [0.0058083497, 0.0055145699]
rho_d_s_coma = [9.90384e-05, 7.04471e-05]

rho_g_s_filament = [0.040263, 0.037849]
rho_d_s_filament = [0.000596, 0.000439]

rho_g_field = [0.251188643150958, 0.0]
rho_d_s_field = [0.000630957344480193, 0.0000]

#function to caculate mass density
# equation rho = phi*Mstar*Gamma(alpha+2)
def massDensity(phi, Mstar, alpha):
	return phi*Mstar*Gamma(alpha+2)

def r3(x):
	return str(x) 

#function to cacualate error in massDensity
def errorDensity(typeof, phi, dphi, Mstar, dMstar, alpha, dalpha):
	val = massDensity(phi, Mstar, alpha)


	#find error in gamma
	gamma = Gamma(2.0 + alpha)
	lowerGamma = abs(gamma - Gamma(2.0 + alpha - dalpha))
	upperGamma = abs(gamma - Gamma(2.0 + alpha + dalpha))
	dGamma = 0.5*(lowerGamma + upperGamma)

	#print gamma, lowerGamma, upperGamma, dGamma
	#print gamma/dGamma

	dVal = val*sqrt((dphi/phi)**2 + (dMstar/Mstar)**2 + (dGamma/gamma)**2)

	#print and return result
	print typeof + ':  rho =   ' + r3(val/1.0E9 ) + '$\\pm$' + r3(dVal/1.0E9)

	return val, dVal

def densityRatio(a, b):
	return log10(a/b)
	
#function to caculate error in log10(rho_1 /rho_steller)
def errorDensityRatio(rhoX, rhoStars):
	
	Z = rhoX[0] / rhoStars[0]
	dZ = Z*sqrt((rhoX[1] /rhoX[0] )**2 + (rhoStars[1]/rhoStars[0])**2)



	print round(Z,decimals=10), round(dZ,decimals=10)
	
	Zmax = log10(Z + dZ)
	Zmin = log10(Z - dZ)
	
	Z = log10(Z)

	upper = abs(Z - Zmax)
	mina = abs(Z - Zmin)


	print ''
	#print Z, upper, mina

	print str(round(Z, decimals=1)) + '$^{' + str(round(upper, decimals=1)) + '}_{' + str(round(mina, decimals=1)) + '}$'
	




rho_Stellar = errorDensity('Stellar', Stellar_phi, Stellar_dphi, Stellar_Mstar, Stellar_dMstar, Stellar_alpha, Stellar_dalpha)
rho_Gas = errorDensity('Gas', Gas_phi, Gas_dphi, Gas_Mstar, Gas_dMstar, Gas_alpha, Gas_dalpha)
rho_Dust = errorDensity('Dust', Dust_phi, Dust_dphi, Dust_Mstar, Dust_dMstar, Dust_alpha, Dust_dalpha)

errorDensityRatio(rho_Gas,rho_Stellar)
errorDensityRatio(rho_Dust,rho_Stellar)






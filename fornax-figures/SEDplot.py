# making sed fit output
# Chris Fuller, May 2013

import numpy as np
from atpy import Table
import matplotlib.pyplot as plt

folder = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/'
cat = Table(folder+"sed.fits",type='fits')
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
   

   
########## control #################
#create fig to hold plots
fig = plt.figure(figsize = (8,10.5),facecolor='w',edgecolor='w')

#number of subplots in the x and y
xx = 5
yy = 5

#plot number start and end
start = 0
end = len(cat)

#create x values
x = np.arange(90.0,800.0,0.5)    
#convert s into jy
from astropy import units as u
lm3 = u.W*u.m**-2*u.m**-1
cjy = lm3.to(u.Jy, equivalencies=u.spectral_density(u.um,x))



#loop through predefined galaxies
for i in range(start,end):
    lx = np.array(bands, dtype=np.float)
    #create subplot
    f = plt.subplot(xx,yy,i+1)

    #get name, chisq, mass and temp
    chisq = cat.CHISQ[i]
    name = cat.OBJECT[i]
    mass = cat.MASS[i]
    temp = cat.T[i]
    
    s = modb(mass,temp, 17.5,x)
    
    
    s_jy = s*cjy
    
    #get a list of fluxes and errors
    flux =[]
    error = []
    for band in bands: 
        flux.append(cat['F'+band][i])
        error.append(cat['E'+band][i])
    
    #turn list into floats
    flux =  np.array(flux , dtype=np.float)
    error = np.array(error, dtype=np.float)
    #f.scatter(lx,flux,s=10, c='k', marker='x')
    #where flux not equal to 0
    
    w1 = np.where(flux != 0.0)[0]
    error =  error[w1]
    flux = flux[w1]
    lx = lx[w1]
    
    if name == 97: name = "97*" 


    f.plot(x,s_jy)
    f.errorbar(lx, flux, yerr=error, xerr=None, c='k', marker='o',ms=2,ls='none')
    
    #plot settings
    f.set_xlim(50.0,1000.0)
    f.set_ylim(0.0011,8000.0)
    f.loglog()
    f.text(70.0,1000.0,"FCC"+str(name))
    f.text(70.0,100.0, "$\chi^{2}$ = "+str(np.round(chisq,decimals=2)), fontsize=10)
    if i != 20: 
        f.tick_params(axis='both',labelleft='off', labelbottom='off')

    if i == 20: 
        f.set_xlabel('Wavelength ($\mu m$)')
        f.set_ylabel('Flux (Jy)')
        f.set_xticks([100., 160., 250., 350., 500.])
        f.set_xticklabels([100,160,250,350,500])
        for tick in f.yaxis.get_major_ticks():
                tick.label.set_fontsize(8) 
                # specify integer or one of preset strings, e.g.
                #tick.label.set_fontsize('x-small') 
                #tick.label.set_rotation(-45)
        for tick in f.xaxis.get_major_ticks():
                tick.label.set_fontsize(8) 
                # specify integer or one of preset strings, e.g.
                #tick.label.set_fontsize('x-small') 
                tick.label.set_rotation(-90)
    #f.text(100,100,i)
plt.subplots_adjust(left=0.07, bottom=0.05, right=0.98, top=0.992, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/fornax/SED_fits.pdf')
#fig.savefig('/Users/chrisfuller/Desktop/pdfs/SED_fits.pdf')
plt.show()



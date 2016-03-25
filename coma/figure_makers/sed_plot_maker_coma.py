# making sed fit output for Coma Galaxies for thesis
# Chris Fuller, May 2013

import numpy as np
from atpy import Table
import matplotlib.pyplot as plt

folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(folder+fname,type='fits')

cat = cat.where(cat.DMASS_TYPE == 2)

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


def printLongFigure(j,fname):

    #if first print caption
    if j == 0:
        print '\\begin{figure}'
        print '\\centering'
        print '\\includegraphics[width=\\linewidth]{'+ fname + '}'
        print '\\caption{First caption}'
        print '\\label{fig:cont}'
        print '\\end{figure}'
    else:
        print '\\begin{figure}'
        print '\\ContinuedFloat '
        print '\\centering'
        print '\\includegraphics[width=\\linewidth]{'+ fname + '}'
        print '\\label{fig:cont}'
        print '\\end{figure}'



   
########## control #################

#create x values
x = np.arange(90.0,800.0,0.5)

#convert s into jy
from astropy import units as u
lm3 = u.W*u.m**-2*u.m**-1
cjy = lm3.to(u.Jy, equivalencies=u.spectral_density(u.um,x))

#number of subplots in the x and y
xx = 5
yy = 5

#plot number start and end
start = 0
N = len(cat)


#caculate number of plots per figure
N_per_fig = xx * yy

#number of frames
N_frames = np.int(np.ceil(np.float(N)/N_per_fig*1.0) )

#loop through each frame
for j in range(N_frames):



    #create fig to hold plots
    fig = plt.figure(figsize = (8,10.5),facecolor='w',edgecolor='w')

    #define end
    end = start + N_per_fig

    if  N - start < N_per_fig: end = N

    #print 'starting frame...', j+1, ' of ', N_frames

    count = 0

    #loop through predefined galaxies
    for i in range(start,end):
        lx = np.array(bands, dtype=np.float)
        #create subplot
        f = plt.subplot(xx,yy,count+1)
        count += 1

        #get name, chisq, mass and temp
        chisq = cat.CHISQ[i]
        name = cat.PAPER_NAME[i]
        mass = cat.DMASS[i]
        temp = cat.DTEMP[i]
        
        s = modb(mass,temp, 100.,x)
        
        
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
        

        f.plot(x,s_jy)
        f.errorbar(lx, flux, yerr=error, xerr=None, c='k', marker='o',ms=2,ls='none')
        
        #plot settings
        f.set_xlim(50.0,1000.0)
        f.set_ylim(0.0011,8000.0)
        f.loglog()
        f.text(70.0,1000.0,str(name))
        f.text(70.0,100.0, "$\chi^{2}$ = "+str(np.round(chisq,decimals=2)), fontsize=10)
        if count != 21: 
            f.tick_params(axis='both',labelleft='off', labelbottom='off')

        if count == 21: 
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

    # ajust and save
    plt.subplots_adjust(left=0.07, bottom=0.05, right=0.98, top=0.992, wspace=0.0, hspace=0.0)
    #fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/fornax/SED_fits.pdf')
    fig.savefig('/Users/chrisfuller/Desktop/SED_fits_' + str(j) + '.pdf')
    #plt.show()
    
    printLongFigure(j,'chapter-coma/SED_fits_' + str(j) + '.pdf')
    start = end

#CFC772

#progam to model chemical evolution
#scaling relations plots
#Chris Fuller

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit

#Inputs
folder = "/Users/chrisfuller/Documents/phd/herchel/coma/final_outputs/"  # input/output folder
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
fname2= 'test-dust-v2.fits'
cat = Table(pj(folder,fname2))

################ Functions ######################

#function for finding errors in a straight line of two arras x and y
# assuming least squares and equal error on each point
def fitLine(x_data,y_data):
    x_data = np.array(x_data, dtype=np.float)
    y_data = np.array(y_data, dtype=np.float)
    
    n = len(x_data)
    D = np.sum(x_data**2) - 1./n * np.sum(x_data)**2
    x_bar = np.mean(x_data)
    p_coeff, residuals, _, _, _ = np.polyfit(x_data, y_data, 1, full=True)
    dm_squared = 1./(n-2)*residuals/D
    dc_squared = 1./(n-2)*(D/n + x_bar**2)*residuals/D
    
    return p_coeff[0], p_coeff[1], np.sqrt(dm_squared[0]),np.sqrt(dc_squared[0])

def residual(params, x, y_data):
    a = params['m'].value
    y_model = a*x
    return np.sqrt(abs(y_data-y_model))

def lmfitter(x, y):
    params = Parameters()
    params.add('m', value=0.01, vary=True)

    out = minimize(residual, params, args=(x, y))
    report_fit(params)
    return out.params['m'].value, 0.0, out.params['m'].stderr, 0.0

def model(params, x):
    a = params['m'].value
    return a*x

##################################################

#select galaxies that have all the metallicities we need
zcat = cat.where((np.nan_to_num(cat.METAL) > 0.0) & (np.nan_to_num(cat.HI_ALL2) > 0.0)  & (cat.DMASS_TYPE != 0))

#caculate the fraction of metals in the dust and in the gas
print np.mean(10**zcat.DMASS/ 10**zcat.MZGAS ) 
print np.std(10**zcat.DMASS/ 10**zcat.MZGAS ) / np.sqrt(len(zcat)*1.0)

#caculate total mass
total = np.log10(10**zcat.Mgastot + 10**zcat.SMASS + 10**zcat.DMASS)

fig, subs = plt.subplots(nrows=1, ncols=1, figsize = (4., 4.), facecolor='w',edgecolor='w')
cluster = np.where(zcat.RADIUS_VIR <= 1.0)
filament = np.where(zcat.RADIUS_VIR > 1.0)
################################################################################

# xxx and yyy
xx = 10**zcat.SMASS / 10**total
yyy = zcat.Ztot 



ff = 1 - xx
xxx = np.log(1/ff)

#fit straight line using np
ecoeffs = fitLine(xxx,yyy)
ecoeffsCluster  = fitLine(xxx[cluster],yyy[cluster])
ecoeffsFilament = fitLine(xxx[filament],yyy[filament])

#fit straight line using np
coeffs = lmfitter(xxx,yyy)
coeffsCluster  = lmfitter(xxx[cluster],yyy[cluster])
coeffsFilament = lmfitter(xxx[filament],yyy[filament])

#draw fit
xxx_line = np.linspace(0,4,100)

sub1 = subs#[0]
sub1.plot(xxx[cluster],yyy[cluster],'ro')
sub1.plot(xxx[filament],yyy[filament],'bo')
sub1.plot(xxx_line, ecoeffs[0]*xxx_line+ ecoeffs[1], '-k')
sub1.plot(xxx_line, ecoeffsCluster[0]*xxx_line + ecoeffs[1], '--r')
sub1.plot(xxx_line, ecoeffsFilament[0]*xxx_line + ecoeffs[1], '--b')
#sub1.plot(xxx_line, 0.01*xxx_line, ls='--', color='green')


sub1.set_ylim(0.0, 0.048)
sub1.set_xlim(0.0, 2.75)
sub1.set_xlabel('$ln(1/f)$')
sub1.set_ylabel('Metallicity $(Z_{Total})$')



################################################################################


xx = 10**zcat.SMASS / 10**total
yy = zcat.DMASS - total

x_line = np.linspace(0.0,1.0, 100)

f = (1 - x_line) 
y_line = np.log10(0.37*ecoeffs[0]*(f)*np.log(1/f)) 
y_lineC = np.log10(0.37*ecoeffsCluster[0]*(f)*np.log(1/f)) 
y_lineF = np.log10(0.37*ecoeffsFilament[0]*(f)*np.log(1/f)) 

y_lineMod = np.log10(0.37*0.01*(f)*np.log(1/f)) 

print 'Peff', ecoeffs[0], '\\pm', ecoeffs[2]
print 'Peff(Filament)', ecoeffsFilament[0] , '\\pm', ecoeffsFilament[2]
print 'Peff(Cluster)', ecoeffsCluster[0], '\\pm', ecoeffsCluster[2]



#sub2 = subs[1]

#sub2.plot(xx[cluster],yy[cluster],'ro')
#sub2.plot(xx[filament],yy[filament],'bo')
#sub2.plot(x_line, y_line, '-k')
#sub2.plot(x_line, y_lineC, '--r')
#sub2.plot(x_line, y_lineF, '--b') 
#sub2.plot(x_line, y_lineMod, ls= '--', color='green') 


#sub2.set_xlabel('$1-f$')
#sub2.set_ylabel('$\log_{10}(M_{Dust}/M_{Total})$')

plt.subplots_adjust(left=0.2, bottom=0.14, right=0.97, top=0.97, wspace=0.2, hspace=0.2)
fig.savefig('/Users/chrisfuller/Desktop/chemical_evolution.pdf')
plt.show()



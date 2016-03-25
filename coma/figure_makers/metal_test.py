#progam to plot metallicity
#Chris Fuller

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
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


fig, subs = plt.subplots(nrows=1, ncols=1, figsize = (8., 8.), facecolor='w',edgecolor='w')

zcat = cat.where(np.nan_to_num(cat.METAL) > 0.0)

cluster = np.where(zcat.RADIUS_VIR <= 1.0)
filament = np.where(zcat.RADIUS_VIR > 1.0)
################################################################################
M, Z = zcat.SMASS, zcat.METAL

coeffs = fitLine(M, Z)

M_line = np.linspace(8., 11.5, 1000)

subs.plot(M_line, M_line*coeffs[0] + coeffs[1], '-k')
subs.plot(M_line, M_line*0.2 + 6.66, '-g')
subs.plot(M_line, - 0.0803*M_line**2 + 1.847*M_line - 1.492, '-r')
subs.plot( M, Z, 'k+')


subs.set_xlabel('$\log_{10}(M_{Star}/M_{\odot})$')
subs.set_ylabel('$12 + \log_{10}(O/H)$')

plt.subplots_adjust(left=0.11, bottom=0.07, right=0.97, top=0.97, wspace=0.2, hspace=0.2)
#fig.savefig('/Users/chrisfuller/Dropbox/chemical_evolution.pdf')
plt.show()



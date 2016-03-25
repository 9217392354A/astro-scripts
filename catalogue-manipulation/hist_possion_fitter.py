#program to fit a possion function to sdss_x_ngp and mysdss_x_npg
#need to check if is nearest ney
import math
import numpy as np
from scipy.misc import factorial
from os.path import join as pj
from atpy import Table
import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import fmin



#i/o
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/positions' 
mycat = Table(pj(folder,'ngp_x_coma+fila-sdss.fits'),type='fits')
fullcat = Table(pj(folder,'ngp_x_full-sdss.fits'),type='fits')

#######################################Functions###################################
#Program for fitting a gaussian distribution 
def possion(p,*args):    
    xDat = args[0]
    yDat = args[1]
    yError=args[2]
     
    lamda = p[0]

    ci = 0.0
    for i in range(0,len(xDat)):
        fx = ((lamda**xDat[i])*np.exp(-lamda)) / factorial(xDat[i])
        
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci

def possion_ci(lamda, xDat, yDat, yError):
    ci = 0.0
    for i in range(0,len(xDat)):
        fx = ((lamda**xDat[i])*np.exp(-lamda)) / factorial(xDat[i])
        
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci

# Program for genorating a gaussian distribution         
def genPossion(fit, x):
	lamda = fit[0]
	return ((lamda**x)*np.exp(-lamda)) / factorial(x)

def mod(x):
	return np.sqrt(x**2.0)

def	itFit(p, *args):

	sep_my = args[0]

	sigma = mod(p[0])
	b = mod(p[1])
	b = 35

	w = np.where(sep_my<sigma*2.0)
	x = sep_my[w]
	bins_mean, n = np.histogram(x, bins=mod(b),density=True)
	fit = fmin(possion,[1.7],args=(n[:-1],bins_mean,np.sqrt(bins_mean)))
	print fit, sigma, b
	ci = possion_ci(fit,n[:-1],bins_mean,np.sqrt(bins_mean))
	return ci

############ main ###########

sep_my = mycat['Separation']
sep_full =  fullcat['Separation']

sep_my = sep_my[np.where(sep_my<10.0)]

bn= 35.0
w = np.where(sep_my<100.0)
x = sep_my[w]
bins_mean, n = np.histogram(x, bins=bn,density=True)
fit = fmin(possion,[1.7],args=(n[:-1],bins_mean,np.sqrt(bins_mean)))

print fit
#fit = fmin(itFit,[10,20],args=(sep_my,sep_my))
fit = fmin(possion,[1.7],args=(n[:-1],bins_mean,np.sqrt(bins_mean)))

xx= np.arange(0,10,0.00001)
y = genPossion(fit,xx)

print 'P(0.05) = ', xx[np.where(y <=0.05)[0][0]], 'P(0.10) = ', xx[np.where(y <=0.10)[0][0]], 'P(0.20) = ', xx[np.where(y <=0.20)[0][0]]  


bins_mean, n = np.histogram(sep_my, bins=bn,density=True)

print 'lamda = ', fit[0], ' 3sigma = ', fit[0]*3.0, '5sigma = ', fit[0]*5.0

#y = fitfunc(p1, x)

fig = plt.figure(figsize = (4.6,4.6),facecolor='w',edgecolor='w')
sub = fig.add_subplot(1,1,1)

sub.bar(n[:-1],bins_mean,width=max(sep_my)/len(n),log=False, alpha = 1.0, fc ='None',hatch='/') 
sub.plot(xx,y, 'r-')
#sub.axvline(x=fit[0]*3.0, ls = '-', c='k')
sub.set_ylabel('Probability Density')
sub.set_xlabel('Angular Separation, (arcsec)')

plt.subplots_adjust(left=0.15, bottom=0.12, right=0.97, top=0.98)
fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/angular_separation.pdf', dpi=600.0)
#plt.show()

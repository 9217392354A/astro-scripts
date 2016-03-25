#Measures veloctiy dispersion
#Chris Fuller, 6 Dec 2012

import numpy as np
import math
import os
from os.path import join as pj
from scipy.optimize import fmin
import matplotlib.pyplot as plt 
from matplotlib import rc

#key varaiables
Rv = 1.68*60.0
vBin = 100.0
Vmax = 970000.22211218
Vmin = -4268.75271268
#files/folders
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/old_stuff/Catalogues/"
outfolder="/Users/chrisfuller/Desktop/"
catfile = 'vel-rad.csv'

#read in catalogues
print 'reading in catalogue'
cat = np.loadtxt(pj(folder,catfile), dtype=float,  skiprows=1, delimiter=",", unpack=False)
catHead = np.loadtxt(pj(folder,catfile), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(cat[0,:]))

#######################################Functions###################################
#Program for fitting a gaussian distribution 
def gauss(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError=args[2]
       
    hmax = x[0]
    hmin = x[1]
    mean =x[2]
    sigma = x[3]

    ci = 0.0
    for i in range(0,len(xDat)):
        A = (1.0/np.sqrt(2*math.pi*sigma**2))
        B = -((xDat[i]-mean)**2)/(2.0*sigma**2)
        fx = (hmax*(A*math.exp(B)))+hmin
        
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci           

# Program for genorating a gaussian distribution         
def genGauss(hmax, hmin, mean, sigma, x):
    y = [] 
    for i in range(0,len(x)):
        A = (1.0/np.sqrt(2*math.pi*sigma**2.0))
        B = -((x[i]-mean)**2.0)/(2.0*sigma**2.0)
        temp = (hmax*A*math.exp(B))
        temp += hmin
        y.append(temp)
    return y

##########################Main program###################

#select galaxies within the virial radius
Rcol = np.where(catHead == "RADIUS")[1][0]
Vcol = np.where(catHead == "VELOCITY")[1][0]
R1 = np.where(cat[:,Rcol] < Rv)[0]
vel = cat[R1,Vcol]
V1 = np.where((vel < Vmax) & (vel > Vmin))[0]
vel = vel[V1]

radius = cat[:,Rcol]
velocity = cat[:,Vcol]


################Peform vel cut and fitting################
v1 = min(vel)
v2 = 0.0
vMid = []
noGal=[]
while v2 < max(vel)  : 
    v2 = v1 + vBin
    vMid = np.append(vMid,(v2 - vBin/2.0))
    temp = []
    for gal in vel:
        if v1 <= gal <= v2:
            temp.append(gal)
    noGal = np.append(noGal,(len(temp) + 1))    
    v1 = v2
#Guess 
x = vMid
sigma = 905.24489991
mean = 6984.48741243
hmin = 0.81587744303
hmax = 42446.2407087
y=[]
y = genGauss(hmax, hmin, mean, sigma, x)
guess = [hmax,hmin,mean,sigma]
###############Fitting data with Gaussian######################
#fit = fmin(gauss,guess,args=(vMid,noGal,np.sqrt(noGal)))


fit = [hmax,hmin,mean,sigma]

fity = genGauss(fit[0],fit[1],fit[2],fit[3],x)

print 'hmax = ', fit[0]
print 'hmin = ', fit[1]
print 'mean = ', fit[2]
print 'sigma = ', fit[3]

#plt.figure(1)
sigmatext = ('$\sigma$ = ' + str(int(fit[3])) + ' km s$^{-1}$')
meantext = ('$\mu$ = ' + str(int(fit[2])) + ' km s$^{-1}$')
backgroundtext = ('background = ' + str(int(fit[1])) + 'galaxy / 200km s$^{-1}$')
bar=plt.bar(vMid, noGal, vBin, color='white',align = 'center', linewidth = 1.0, alpha = 1.0)
#plt.plot(x,y, '--')
ga=plt.plot(x,fity, 'k', linewidth = 1)
#plt.text(10000, 50, sigmatext , fontdict=None, fontsize=18)
#plt.text(10000, 60, meantext , fontdict=None, fontsize=18)
#plt.text(8100, 64, backgroundtext , fontdict=None)
plt.xlabel('Velocity, km s$^{-1}$ ', fontsize=10)
plt.ylabel('Number of Galaxies', fontsize=10)
plt.xlim(xmin=4000, xmax=10000)
plt.ylim(ymax=(max(noGal)+max(noGal)/10.0))
plt.tick_params(axis='both', labelsize = 'small')
#plt.title('Coma Cluster Velocity Distribution', fontsize=20)


#mean1 = fit[2]
#sigma1 = fit[3]

sigma1 = 905.244899
mean1 = 6984.4874

#a = plt.plot([(sigma1*2 + mean1),(sigma1*2 + mean1)], [0,max(noGal)], '--b')
#plt.plot([(mean1 - sigma1*2),(mean1 - sigma1*2)], [0,max(noGal)], '--b')
b = plt.plot([(sigma*3 + mean),(sigma*3 + mean)], [0,(max(noGal)+max(noGal)/10.0)], '--k')
plt.plot([(mean - sigma*3),(mean - sigma*3)], [0,(max(noGal)+max(noGal)/10.0)], '--k')
#plt.legend([a,b],['2$\sigma$','3$\sigma$'])

print '3sigma min max', (sigma*3 + mean) , (mean - sigma*3)
#plt.show()
plt.savefig(pj(outfolder,"vel.pdf"))

#plt.figure(2)

rect = plt.Rectangle((0,mean - 3*sigma),Rv, sigma*6, color="purple", alpha=0.4)
plt.gca().add_patch(rect)

rect = plt.Rectangle((Rv,mean - 3*sigma),600-Rv, sigma*6, color="grey", alpha=0.4)
plt.gca().add_patch(rect)



plt.scatter(radius,velocity, s=5, c='black' )
plt.xlim(xmin=0, xmax=600)
plt.ylim(ymin=3000, ymax=11000)
plt.ylabel('Velocity, km s$^{-1}$ ', fontsize=10)
plt.xlabel('Distance to Cluster Centre, arcsec', fontsize=10)
plt.tick_params(axis='both', labelsize = 'small')
plt.plot([0,10000.0],[(sigma*3 + mean),(sigma*3 + mean)], '--k')
plt.plot([0,10000.0],[(mean - sigma*3),(mean - sigma*3)], '--k')
plt.plot([Rv,Rv],[0,100000000],'--k')

plt.savefig(pj(outfolder,"rad-vel.pdf"))
#plt.show()


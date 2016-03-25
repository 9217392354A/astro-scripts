#Program to fit a gausian to a velocity profile, 
# Chris Fuller March 2012
from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt 
################user varaibles##################################
catalogefile = "coma_ngp_tight.csv"
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/old_stuff/Catalogues/"
vBin = 200.0

###################Initilising Variables########################
noGal =[]
vMid = [] 
####################read in catalogue############################

cat = loadtxt(pj(folder,catalogefile), dtype=float,  skiprows=1, delimiter=",", unpack=False)
vel = cat[:,8]

########################Functions##################################
def gauss(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError=args[2]
       
    hmax = x[0]
    hmin = x[1]
    mean =x[2]
    sigma = x[3]
    A = (1.0/sqrt(2*math.pi*sigma**2))
    B = -((xDat-mean)**2)/(2.0*sigma**2)
    fx = (hmax*(A*math.exp(B)))+hmin
    
    ci = ((yDat - fx)**2)/ (yError**2)
    return ci           
        
def genGauss(hmax, hmin, mean, sigma, x):
    y = [] 
    for i in range(0,len(x)):
        A = (1.0/sqrt(2*math.pi*sigma**2.0))
        B = -((x[i]-mean)**2.0)/(2.0*sigma**2.0)
        temp = (hmax*A*math.exp(B))
        temp += hmin
        y.append(temp)
    return y

#bins up all the galaxies in the veloicty col    
def bineybinbin(vel,vBin):
    v1 = min(vel)
    v2 = 0.0    
    while v2 < max(vel):
        v2 = v1 + vBin
        vMid.append(v2 - vBin/2.0)
        temp = []
        for gal in vel:
            if v1 <= gal <= v2:
                temp.append(gal)
        if len(temp) == 0:
            temp.append(1.0)
        noGal.append(len(temp))    
        v1 = v2
    return noGal,vMid
################Peform vel cut and fitting################
binOutput = bineybinbin(vel,vBin)
noGal = binOutput[0]
vMid = binOutput[1]
#Guess 
x = vMid
sigma = 950.0
mean = 7000.0
hmin = 3.0
hmax = 160000.0
y=[]
y = genGauss(hmax, hmin, mean, sigma, x)
guess = [hmax,hmin,mean,sigma]
###############Fitting data with Gaussian######################
fit = fmin(gauss,guess,args=(vMid,noGal,sqrt(noGal)))
print fit 
fity = genGauss(fit[0],fit[1],fit[2],fit[3],x)

oldSigma = sigma
newSigma = 0.0
i=0
while abs(newSigma-oldSigma <= 1.0):
    i +=1
    #bin everything
    binOutput = bineybinbin(vel,vBin)
    noGal = binOutput[0]
    vMid = binOutput[1]
    fit = fmin(gauss,guess,args=(vMid,noGal,sqrt(noGal)))
    print "iter = ", i, "  $\sigma$ = ", fit[3], "  mean = ", fit[2]
    print "diff= " , abs(newSigma-oldSigma), "len new vel" 
    newSigma = fit[3]
    
    

###################Plotting Final Product################################
sigmatext = ('$\sigma$ = ' + str(int(fit[3])) + ' km s$^{-1}$')
meantext = ('$\mu$ = ' + str(int(fit[2])) + ' km s$^{-1}$')
backgroundtext = ('background = ' + str(int(fit[1])) + 'galaxy / 200km s$^{-1}$')
bar=plt.bar(vMid, noGal, vBin, color='red',align = 'center', alpha = 0.3)
#plt.plot(x,y, '--')
ga=plt.plot(x,fity, '--' )
plt.text(8100, 80, sigmatext , fontdict=None)
plt.text(8100, 72, meantext , fontdict=None)
#plt.text(8100, 64, backgroundtext , fontdict=None)
plt.xlabel('Velocity, km s$^{-1}$ ', fontsize=12)
plt.ylabel('NN', fontsize=12)
plt.title('Coma Cluster Velocity Distribution', fontsize=20)
plt.show()

###create bins and work out how many galaxies fall into that bin###

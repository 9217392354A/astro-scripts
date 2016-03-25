#Program to fit ak gausian to a velocity profile, and plot other stuff! 
# Chris Fuller March 2012
#v1.0

from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt 
################user varaibles##################################

catalogefile = "coma_ngp_tight.csv"
folder = "/Users/chrisfuller/Dropbox/coma/Catalogues"
vBin = 200.0

velLower = 0000.0
velUpper = 12000.0
radLower = 0.0
radUpper = 1.8

r = 0.2

##################Initilising Variables########################
noGal =[]
vMid = [] 

####################read in catalogue############################

cat = loadtxt(pj(folder,catalogefile), dtype=float,  skiprows=1, delimiter=",", unpack=False)
velInitial = cat[:,8]
radInitial = cat[:,7]


########################Functions##################################
def select(cat,col,lower,upper):
    data = cat[:,col]
    newcat = []
    for i in range(0,len(cat)):
        if lower <= data[i] <= upper:
            newcat.append([cat[i,0],cat[i,1],cat[i,2],cat[i,3],cat[i,4],cat[i,5],cat[i,6],cat[i,7],cat[i,8],cat[i,9]])
        
    newcat = array(newcat)            
    return newcat         

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
        A = (1.0/sqrt(2*math.pi*sigma**2))
        B = -((xDat[i]-mean)**2)/(2.0*sigma**2)
        fx = (hmax*(A*math.exp(B)))+hmin
        
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci           

# Program for genorating a gaussian distribution         
def genGauss(hmax, hmin, mean, sigma, x):
    y = [] 
    for i in range(0,len(x)):
        A = (1.0/sqrt(2*math.pi*sigma**2.0))
        B = -((x[i]-mean)**2.0)/(2.0*sigma**2.0)
        temp = (hmax*A*math.exp(B))
        temp += hmin
        y.append(temp)
    return y
    
#????
def radBin(r,rad):
    counts = []
    rMid = []
    r2 = 0.0
    r1 = 0.0
    while r2 < max(rad):
        r2 = r1 + r
        temp =[]
        for j in rad:
            if r1 <= j <= r2:
                temp.append(j)        
        area = math.pi*((r2**2)-(r1**2))
        counts.append((len(temp)/area))
        rMid.append((r1+r2)/2.0) 
        r1 = r2
    return rMid,counts
  
#function for fitting an exponential curve  
def curve(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError=args[2]
    A = x[0]
    B = x[1]
    C = x[2]
    ci = 0.0
    for i in range(0,len(xDat)):
        fx = A*math.exp(-xDat[i]/B)+ C
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci
 
# Function for genorating a exponential curve    
def genCurve(x,xDat):
    y =[]
    A = x[0]
    B = x[1]
    C = x[2]
    
    for i in range(0,len(xDat)):
        y.append(A*math.exp(-xDat[i]/B)+ C)
    return y
           
    
############Creating initial data sets##########

ben = select(cat,9,velLower,velUpper)
alfie = select(ben,8,radLower,radUpper)

vel = alfie[:,9]
rad = alfie[:,8]


guessrad = [1200.0,4.0,10.0]
radcounts = radBin(r,rad)
radfit = fmin(curve,guessrad, args =[radcounts[0],radcounts[1],sqrt(radcounts[1])])


yC = genCurve(radfit, radcounts[0])

plt.plot(radcounts[0],yC)
plt.bar(radcounts[0],radcounts[1],width = r,color='red',align = 'center', alpha = 0.5, log='true')
plt.xlabel('Radius from cluster center, r /deg', fontsize=12)
plt.ylabel('Galaxies/degree$^2$', fontsize=12)
plt.title('Coma Cluster angular distribution', fontsize=20)
plt.grid()
plt.show()    
 
################Peform vel cut and fitting################
v1 = min(vel)
v2 = 0.0
while v2 < max(vel)  : 
    v2 = v1 + vBin
    vMid.append(v2 - vBin/2.0)
    temp = []
    for gal in alfie:
        if v1 <= gal[9]<= v2:
            temp.append(gal[9])
    noGal.append(len(temp) + 1)    
    v1 = v2
#Guess 
x = vMid
sigma = 750.0
mean = 7100.0
hmin = 3.0
hmax = 160000.0
y=[]
y = genGauss(hmax, hmin, mean, sigma, x)
guess = [hmax,hmin,mean,sigma]
###############Fitting data with Gaussian######################
fit = fmin(gauss,guess,args=(vMid,noGal,sqrt(noGal)))



fity = genGauss(fit[0],fit[1],fit[2],fit[3],x)
print 'hmax = ', fit[0]
print 'hmin = ', fit[1]
print 'mean = ', fit[2]
print 'sigma = ', fit[3]
print 'r = ', radUpper

sigmatext = ('$\sigma$ = ' + str(int(fit[3])) + ' km s$^{-1}$')
meantext = ('$\mu$ = ' + str(int(fit[2])) + ' km s$^{-1}$')
backgroundtext = ('background = ' + str(int(fit[1])) + 'galaxy / 200km s$^{-1}$')
bar=plt.bar(vMid, noGal, vBin, color='white',align = 'center', linewidth = 2.0, alpha = 1.0)
#plt.plot(x,y, '--')
ga=plt.plot(x,fity, 'k', linewidth = 2)
#plt.text(10000, 50, sigmatext , fontdict=None, fontsize=18)
#plt.text(10000, 60, meantext , fontdict=None, fontsize=18)
#plt.text(8100, 64, backgroundtext , fontdict=None)
plt.xlabel('Velocity, km s$^{-1}$ ', fontsize=18)
plt.ylabel('Number of Galaxies', fontsize=18)
plt.xlim(xmin=min(vel))
plt.ylim(ymax=(max(noGal)+max(noGal)/10.0))
plt.tick_params(axis='both', labelsize = 'large')
plt.title('Coma Cluster Velocity Distribution', fontsize=20)

#mean1 = fit[2]
#sigma1 = fit[3]

sigma1 = 922.03
mean1 = 6983.02

#a = plt.plot([(sigma1*2 + mean1),(sigma1*2 + mean1)], [0,max(noGal)], '--b')
#plt.plot([(mean1 - sigma1*2),(mean1 - sigma1*2)], [0,max(noGal)], '--b')
b = plt.plot([(sigma1*3 + mean1),(sigma1*3 + mean1)], [0,(max(noGal)+max(noGal)/10.0)], '--k')
plt.plot([(mean1 - sigma1*3),(mean1 - sigma1*3)], [0,(max(noGal)+max(noGal)/10.0)], '--k')
#plt.legend([a,b],['2$\sigma$','3$\sigma$'])

print '2sigma min max', (sigma1*2 + mean1) , (mean1 - sigma1*2)

plt.show()

###plot


plt.scatter(vel,rad,s=5, marker ='o', c='k' )
b = plt.plot([(sigma1*3 + mean1),(sigma1*3 + mean1)], [0,max(rad)], '--k', linewidth =2)
plt.plot([(mean1 - sigma1*3),(mean1 - sigma1*3)], [0,max(rad)], '--k', linewidth =2)
plt.plot([0,max(vel)],[1.68,1.68], '--k', linewidth =2)
#plt.legend([a,b],['2$\sigma$','3$\sigma$'])
plt.xlabel('Velocity, km s$^{-1}$ ', fontsize=18)
plt.ylabel('Radius, degree', fontsize=18)
plt.title('Coma Cluster Velocity/Radial Distribution', fontsize=20)
plt.ylim(0,(max(rad)))
plt.xlim(min(vel), max(vel))
plt.show()
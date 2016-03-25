#contor mapping
# Chris Fuller July 2012
from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import time
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import * 
################user varaibles##################################
catFile = "coma_lines_test.csv"
folder = "/Users/chrisfuller/Dropbox/coma/Catalogues/current/"
#switch: mean = 1, median = 2, number = 3 sum = 4
switch = 2
delta = 3.0
#plotting options
conPlot = 1
imPlot = 0
####################read in catalogue############################
cat = loadtxt(pj(folder,catFile), dtype=float,  skiprows=1, delimiter=",", unpack=False)
catHead = loadtxt(pj(folder,catFile), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(cat[0,:]))
########################Functions##################################
#Function to create bins and bin up, for whatever columb that is desired with control over the bin size and 
#switchs for mean and median
def contorBin(cat,catHead,col, switch, delta):
    RAcol = where(catHead == "ra")[1][0]
    DECcol = where(catHead == "dec")[1][0]
    #create the bins in x and y / ra/dec
    x = np.arange(floor(min(cat[:,RAcol])), ceil(max(cat[:,RAcol])), delta)
    y = np.arange(floor(min(cat[:,DECcol])), ceil(max(cat[:,DECcol])), delta)
    X, Y = np.meshgrid((x + (delta/2.0)), (y + (delta/2.0)))
    Z = np.zeros((len(y),len(x)))
    #cycle through x values first
    starttime = time.time()
    for i in range(0,len(x)-1):
        for j in range(0, len(y)-1):
            Z[j,i] = contorBina(cat,catHead,col, switch,i,j,x,y)
        if i == 0:
            print "estimated time for", switch, " = ", ((time.time()-starttime)*len(x))/60.0
    return X,Y,Z

def contorBina(cat,catHead,col, switch,i,j,x,y):
    RAcol = where(catHead == "ra")[1][0]
    DECcol = where(catHead == "dec")[1][0] 
    temp = []                  
    for gal in cat:
        if x[i] <= gal[RAcol] <= x[i+1] and y[j] <= gal[DECcol] <= y[j +1]:
            temp = np.append(temp,float(gal[col]))
    if switch == 1 and len(temp) > 0:
        temp = np.mean(temp)
    elif switch == 2 and len(temp) > 0:
        temp = np.median(temp)
    elif switch == 3:
        temp = len(temp)
    elif switch ==4 and len(temp) > 0:
        temp = sum(temp)
    else:
        temp = 0.0
    return temp

#Function that selected in velocity and lines of constant ra and dec
def selecta(cat,col,lower,upper):
    newCat = []
    count = 0
    for row in cat:
        if lower <= row[col] <= upper and count == 0:
            newCat = row
            newCat = newCat.reshape(1,len(row)) 
            count += 1
        elif lower <=row[col] < upper:
            newCat = append(newCat,row.reshape(1,len(row)), axis=0)
    return newCat

############################## Main Program ###################

X,Y,Z = contorBin(cat,catHead,where(catHead == "SFRtot")[1][0], switch, delta)

Xden,Yden,Zden = contorBin(cat,catHead,where(catHead == "SFRtot")[1][0], 3, delta)

Zden = Zden / (delta**2)
RAcol = where(catHead == "ra")[1][0]
DECcol = where(catHead == "dec")[1][0]

RA = cat[:,RAcol]
DEC = cat[:,DECcol]
#pcolor(X, Y, Z)
#colorbar()

#imshow(Z, interpolation='bicubic')
#figure(1)
#grid(True)
#colorbar()

if conPlot == 1:
    #V = log10(arange(10,ceil(max(Z.flat)),10))
    matplotlib.rcParams['xtick.direction'] = 'out'  
    matplotlib.rcParams['ytick.direction'] = 'out'
    plt.figure()
    CS = plt.contourf(X,Y,Z,linestyles ="dotted", alpha = 0.3,nchunk=1)
    colorbar()
    Den = plt.contour(Xden,Yden,Zden )
    Pos= plt.plot(RA,DEC, 'ko', markersize=1, alpha=0.2)
    colorbar()
    plt.title("Contour Map SFR vs Local Dencity")
    plt.xlabel("RA", fontsize=18)
    plt.ylabel("DEC", fontsize=18)
    plt.xlim(185,215)
    plt.ylim(15, 40)
    plt.gca().invert_xaxis()
    plt.show()

if imPlot == 1:
    plt.imshow(Z, interpolation='nearest')
    plt.show()

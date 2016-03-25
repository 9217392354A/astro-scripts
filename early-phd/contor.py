# Chris Fuller May 2012
from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt 

################user varaibles##################################
folder = "/Users/chrisfuller/Dropbox/coma/Catalogues"
catFile = "coma_ngp_tight.csv"
###################Initilising Variables########################
R=[]
####################read in catalogue###########################
cat = loadtxt(pj(folder,catFile), dtype=float,  skiprows=1, delimiter=",", unpack=False)
catHead = loadtxt(pj(folder,catFile), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(cat[0,:]))
#read in colour cut line in the form [m,c]
cut = array(open(pj(folder,"redblue-divider.csv")).read().split(","),dtype=float)

########################Functions###############################
def genStraight(x,xDat):
    y =[]
    M = x[0]
    C = x[1]    
    for i in range(0,len(xDat)):
        y.append(M*xDat[i] + C)
    return y 

def Zfind(x,y):
    xBins=[]
    yBins=[]
    xBins = range(min(x)*1000.0, max(x)*1000.0, (-min(x)+max(x))/50.0)/1000.0
    yBins = range(min(y)*1000.0, max(y)*1000.0, (-min(y)+max(y))/50.0)/1000.0
    Z
    for i in range(0,max(xBins) - 1):
        tempX =[]
        tempX = selecta(xbins[i],xbin[i+1],x)
        colY = []
        for j in range(0,yBins):
            tempY = []
            tempY = selecta(ybins[j],ybin[j+1],tempX)
            colY = append(colY, len(tempY))
    return Z
    
def selecta(lower,upper,array):
    temp = []
    for element in array:
        if lower <= element <= upper:
            temp = append(temp,element)
        else:
            print 'warining selecta'
    return temp

        
    
###################Main Program##################################
redMr = []
redC = []
blueMr = []
blueC = []
col = []
raRed= []
decRed = []
raBlue = []
decBlue = []
redVel =[]
blueVel =[]
redRad = []
blueRad = []
for row in cat:
    if 1==1:
        r = row[5]
        g = row[4]
        ra = row[1]
        dec = row[2]
        col = r-g
        if col  < r*cut[0] + cut[1]:
            redC = append(redC, col)
            redMr = append(redMr, (r))
            raRed = append(raRed, ra)
            decRed = append(decRed, dec)
            redVel = append(redVel,row[8])
            redRad = append(redRad,row[12])
        if col >  r*cut[0] + cut[1]:
            blueC = append(blueC, col)
            blueMr = append(blueMr, (r))
            raBlue = append(raBlue, ra)
            decBlue = append(decBlue, dec)
            blueVel = append(blueVel,row[8])
            blueRad = append(blueRad,row[12])

R = float(len(blueMr))/float(len(redMr))
print "blue / red = ", R, " Filament"

plt.contour(raRed,decRed)
plt.show()
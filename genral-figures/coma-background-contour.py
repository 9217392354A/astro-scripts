#program to plot contours of surface density of the background to the coma cluster
#Chris Fuller feb 2013

import numpy as np 
import atpy as at
from os.path import join as pj
import time
import datetime
import matplotlib.pyplot as plt


#user variables
delta = 0.5 #pixel size degree's


#folders and files
fileName = "virgo_spx6cff.fit"
folder = "/Users/chrisfuller/Desktop/cluster_dust/"



print "reading in cat"
cat = at.Table(pj(folder,fileName),type='fits')

folder = "/Users/chrisfuller/Desktop/"

#import data
ra = cat.ra
dec = cat.dec

g=cat.g
r=cat.r

#clean data
w1 = np.where((g> -100.0) & (r > -100.0) & (r < 40.00))
g = g[w1]
r = r[w1]
c= g-r

ra=ra[w1]
dec=dec[w1]


##########################################################   Functions   ##########################################################
#work out each pixel colour and galaxy density
##### Functions ######
def contorBinDencity(ra, dec, delta, X,Y,x,y,fa):
    #create the bins in x and y / ra/dec
    #x = np.arange(np.floor(np.min(ra)), np.ceil(np.max(ra)), delta)
    #y = np.arange(np.floor(np.min(dec)), np.ceil(np.max(dec)), delta)
    #X, Y = np.meshgrid((x + (delta/2.0)), (y + (delta/2.0)))
    Z = np.zeros((len(y),len(x)))
    den = np.array([],dtype=float)
    #cycle through x values first
    starttime = time.time()
    for i in range(0,len(x)):
        for j in range(0, len(y)):
            Z[j,i] = len(whatsInTheBin(ra, dec, x[i], y[j],delta/2.0,fa))
            den = np.append(den,Z[j,i])
        if i == 0:
            print "estimated time for surface density bin task"," = ", str(datetime.timedelta(seconds=((time.time()-starttime)*len(x))))
    return Z,den
    
def contorBinMean(ra, dec, delta, c,fa):
    #create the bins in x and y / ra/dec
    x = np.arange(np.floor(np.min(ra)), np.ceil(np.max(ra)), delta)
    y = np.arange(np.floor(np.min(dec)), np.ceil(np.max(dec)), delta)
    X, Y = np.meshgrid((x + (delta/2.0)), (y + (delta/2.0)))
    col = np.array([],dtype=float)
    Z = np.zeros((len(y),len(x)))
    #cycle through x values first
    starttime = time.time()
    for i in range(0,len(x)):
        for j in range(0, len(y)):
            Z[j,i] = np.mean(c[whatsInTheBin(ra, dec, x[i], y[j],delta/2.0,fa)])
            col = np.append(col,Z[j,i])
        if i == 0:
            print "estimated time for mean bin task"," = ", str(datetime.timedelta(seconds=((time.time()-starttime)*len(x))))
    return X,Y,Z,x,y,col

def whatsInTheBin(ra, dec, x, y, delta,fa):
    delta = delta*fa
    return np.where((ra> (x-delta)) & (ra < (x +delta)) & (dec > (y -delta)) & (dec < (y +delta)))[0]



########################################## main program ############################################################

delt = [1.0]#pixel size degree's
#factor = np.arange(1.0,10.0,0.25)
factor = [10.0]
for delta in delt:
    for fa in factor:
        print delta
        X,Y,Zbackground,x,y, colour = contorBinMean(ra, dec, delta, c,fa)
        # plotplot
        plt.figure(1)
        csb =  plt.contour(X,Y,Zbackground)
        plt.xlabel("ra")
        plt.ylabel("dec")
        plt.title("delta ="+str(delta)+" faintest magnitude = " + str(max(r)))
        #plt.savefig(pj(folder,"contour-both-delta"+str(delta)+ "pixel_factor-"+str(fa)+".eps"))        
        #plt.clf()        
        
        plt.figure(2)
        plt.imshow(Zbackground, interpolation= "bicubic")
        plt.xlabel("ra")
        plt.ylabel("dec")
        plt.title("delta ="+str(delta)+" faintest magnitude = " + str(round(max(c))))
        #plt.savefig(pj(folder,"background-delta"+str(delta)+"pixel_factor-"+str(fa)+".eps"))
        #plt.clf()
              
plt.show()


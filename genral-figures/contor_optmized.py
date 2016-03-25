#program to plot contours of surface density of the background to the coma cluster
#Chris Fuller feb 2013

import numpy as np 
import atpy as at
from os.path import join as pj
import matplotlib.pyplot as plt
import multiprocessing as m

#user variables
delta = 5.0 #pixel spacing degree's
fa = delta #pixel size degrees

#folders and files
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/reddening-background-coma/"

print "reading in cat"
cat = at.Table(pj(folder,"npg_big_spx6cff-comaRemoved.fit"),type='fits')

#import data
ra,dec = cat.ra,cat.dec
g, r= cat.g, cat.r

#clean data
w1 = np.where((g> -100.0) & (r > -100.0) & (r < 40.00))
c = g[w1] - r[w1]
ra,dec=ra[w1], dec[w1]

##### Functions ######

def testy(x):
	y = x*100.0
	return y
    
def contorBinMean(ra, dec, delta, c,fa):
    #create the bins in x and y / ra/dec to make meshgrid
    x = np.arange(np.floor(np.min(ra)), np.ceil(np.max(ra)), delta)
    y = np.arange(np.floor(np.min(dec)), np.ceil(np.max(dec)), delta)
    X, Y = np.meshgrid((x + (delta/2.0)), (y + (delta/2.0)))
    
    #create array to hold final values
    Z = np.zeros((len(y),len(x)))
    
    #cycle through x values first
    
    p = m.Pool(processes=2)
    r_x = range(0,len(x))
    r_y = range(0,len(y))
    
    print p.map(testy,r_x)
    
    for i in range(0,len(x)):
        for j in range(0, len(y)):
            Z[j,i] = whatsInTheBin(ra, dec, x[i], y[j],delta/2.0,fa,c)
    return X,Y,Z,x,y

#histogram2d 

def whatsInTheBin(ra, dec, x, y, delta,fa,c):
    delta = delta*fa
    w = np.where((ra> (x-delta)) & (ra < (x +delta)) & (dec > (y -delta)) & (dec < (y +delta)))[0]
    return np.mean(c[w])



##################### main program ##################

X,Y,Zbackground,x,y = contorBinMean(ra, dec, delta, c,fa)

plt.contour(X,Y,Zbackground)
plt.show()      
            
      


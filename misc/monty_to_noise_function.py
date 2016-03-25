#program written by Chris Fuller May 2012
#program that works out th standerd deviation of an image, using a clipping method
#files must be saved in the same dir with name mc-10-250.txt 


#import stuff
from numpy import *
import numpy
import scipy
import math
import sys
from scipy.optimize import fmin
import matplotlib.pyplot as plt 
import os
from os.path import join as pj

#File stuff
cat = "bigcoma.csv"
catfolder = "/Users/chrisfuller/Dropbox/coma/Catalogues"
catout ="comaTEST.csv" 
valuesFile = "valueserror.txt"
folder = "/Users/chrisfuller/Dropbox/coma/flux/montycarlo2/"
files = os.listdir(folder)
m160 = [x for x in files if x[-7:] == "160.txt"]
m250 = [x for x in files if x[-7:] == "250.txt"]
m350 = [x for x in files if x[-7:] == "350.txt"]
m500 = [x for x in files if x[-7:] == "500.txt"]

######################Functions####################################### 
def funtoflux(filename, band):
    if int(band[:3]) == 250:
        x = 59.1
    elif int(band[:3]) == 350:
        x = 133.0
    elif int(band[:3]) == 500:
        x = 63.0
    elif int(band[:3]) == 160 or int(band[:3]) == 70 or int(band[:3]) == 100:
        x = 1000.0
    else:
        raise 'no band set'
    output = ["flux " + str(band) +"um" , "pixels"]
    output = array(output,ndmin=2)  
    output = output.reshape(1,2)
    funtoolsfile = open(pj(folder, filename), 'r')
    count = 0
    for line in funtoolsfile.readlines():
        count += 1
        try: 
            if line[3] == "1" and count > 16 and line[2] == " ":
                firstline = count -1
            else:
                continue
        except:
            continue
    count = 0    
    funtoolsfile.close()
    funtoolsfile = open(pj(folder, filename), 'r')    
    for row in funtoolsfile.readlines():
        count += 1
        info = row.split()
        try:        
            if count > firstline:
                selection = array((str(x*float(info[1])), str(info[2])))
                selection = selection.reshape(1,len(selection)) 
                output = append(output,selection, axis=0)
        except:
            continue 
    return output
    
def noisefinder(flux):
    a = flux[1:,0]
    b = array(a,dtype=float)
    sigma = std(b)
    newsigma = sigma*2.0
    count = 0
    while newsigma/sigma > 0.999:    
        sigma = std(b)
        meanF =[0.0]
        meanF = mean(b)
        selection = where((b < (meanF + 3*sigma))&(b > (meanF - 3*sigma)))
        b = b[selection]
        count += 1
        newsigma = std(b)
    return newsigma

# for fitting a second order polynomial
def curvefitter(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError = 1.0
    A = x[0]
    B = x[1]
    C = x[2]
    ci = 0.0
    for i in range(0,len(xDat)):
        fx = A*xDat[i]**2 + B*xDat[i] + C
        ci += ((yDat[i] - fx)**2)/ (yError**2)
    return ci
 # for genorating a second order polynomial     
def curveGen(x,xval):
     xGen = range(0,int(max(xval)))
     y = []
     for i in range(0,len(xGen)):
         fx = []
         fx = fx = x[0]*xGen[i]**2 + x[1]*xGen[i] + x[2]
         y = append(y,fx)
     return xGen,y


#######################################################################################

noise250 = []
numPix250 =[]
noise350 = []
numPix350 = []
noise500 = []
numPix500 = []
noise160 = []
numPix160 = []
#find the noise


for y in m250:
    flux = []
    noise =[]
    band = y.split("-")
    flux = funtoflux(y,band[2])
    noise = noisefinder(flux)
    numPix = median(double(flux[1:,1]))
    noise250 = append(noise250, noise)
    numPix250 = append(numPix250, numPix)    

for y in m350:
    flux = []
    noise =[]
    band = y.split("-")
    flux = funtoflux(y,band[2])
    noise = noisefinder(flux)
    numPix = median(double(flux[1:,1]))
    noise350 = append(noise350, noise)
    numPix350 = append(numPix350, numPix)
    
for y in m500:
    flux = []
    noise =[]
    band = y.split("-")
    flux = funtoflux(y,band[2])
    noise = noisefinder(flux)
    numPix = median(double(flux[1:,1]))
    noise500 = append(noise500, noise)
    numPix500 = append(numPix500, numPix)
    
for y in m160:
    flux = []
    noise =[]
    band = y.split("-")
    flux = funtoflux(y,band[2])
    noise = noisefinder(flux)
    numPix = median(double(flux[1:,1]))
    noise160 = append(noise160, noise)
    numPix160 = append(numPix160, numPix)

noise250.sort()
numPix250.sort()
noise350.sort()
numPix350.sort()
noise500.sort()
numPix500.sort()
noise160.sort()
numPix160.sort()

x = [1,1,1]


fit160 = fmin(curvefitter,x,args=(numPix160,noise160))
curve160 = curveGen(fit160,numPix160)
fit250 = fmin(curvefitter,x,args=(numPix250,noise250))
curve250 = curveGen(fit250,numPix250)
fit350 = fmin(curvefitter,x,args=(numPix350,noise350))
curve350 = curveGen(fit350,numPix350)
fit500 = fmin(curvefitter,x,args=(numPix500,noise500))
curve500 = curveGen(fit500,numPix500)

table =[]
table = append(fit160.reshape(3,1),fit250.reshape(3,1), axis=1)
table = append(table, fit350.reshape(3,1), axis=1)
table = append(table,fit500.reshape(3,1),axis=1)



savetxt(pj(folder, valuesFile), table, delimiter=",")
plt.plot(numPix160, noise160, 'bx', label='PACS 160')
plt.plot(numPix250, noise250, 'rx', label='SPIRE250')
plt.plot(numPix350, noise350, 'gx', label='SPIRE350')
plt.plot(numPix500, noise500, 'cx', label='SPIRE500')
plt.plot(curve160[0],curve160[1],'b', label='PACS 160')
plt.plot(curve250[0],curve250[1],'r', label='SPIRE 250')
plt.plot(curve350[0],curve350[1],'g', label='SPIRE 350')
plt.plot(curve500[0],curve500[1],'c', label='SPIRE 500')
plt.xlabel('Pixels per aperture', fontsize=18)
plt.ylabel('Standard deviation', fontsize=18)
#plt.title('Sigma vs Pixels for ATLAS NGP', fontsize=20)
#plt.ylim(ymax=1000)
#plt.xlim(xmax=1000)
plt.loglog()
plt.legend(loc=4)
plt.show()





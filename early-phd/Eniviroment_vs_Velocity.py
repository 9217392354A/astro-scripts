#program created to plot veldispersion of the coma filiament vs starformation, or other traces of 
#enviromental process's
#Chris Fuller September 2012

from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import time
import matplotlib
import matplotlib.pyplot as plt
from scipy import fftpack
import pylab as py
import radialprofile

catFile = "tight_GalParam_030812_1129_spx6cff_0.csv"
folder = "/Users/chrisfuller/Dropbox/coma/AuxData/SFR/"
catFile = "Wall_galparam_dr9_spx6cff.csv"

#tasks localDen means surface plotting everything first if 0 then it will be drawn 
#from the saved cataloge
surfaceDen = 1
#switch: mean = 1, median = 2, sum = 3
switch = 2
#Bin size degree
delta = 3.0

#Make initial selection from catalogue **Selection not working yet requries code from catChanger**
#clean just selected from log sfr total -10<sfrtotal<10 remove if flag's to be included
#Selection if  selection = 0 then selection bypassed
selection = 1
clean = 1
vel_lower = 4000.0
vel_upper = 10000.0
ra_lower = 160.0
ra_upper = 209.0
dec_lower = 10.0
dec_upper = 36.0

#plotting options
conPlot = 0
imPlot = 0
faPlot =0
simPlot = 0
#if fourrier analisis hit 1
fa = 0

####################read in catalogue############################
print 'Reading in Cataloge'
cat = loadtxt(pj(folder,catFile), dtype=float,  skiprows=1, delimiter=",", unpack=False)
catHead = loadtxt(pj(folder,catFile), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(cat[0,:]))
starttime = time.time()
########################Functions##################################
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

#finding surface dencity of various parameters
def surfaceDencity(cat,catHead,col, switch,i,j,x,y):
    RAcol = where(catHead == "ra")[1][0]
    DECcol = where(catHead == "dec")[1][0] 
    temp = []   
    n =[]               
    for gal in cat:
        if x[i] <= gal[RAcol] <= x[i+1] and y[j] <= gal[DECcol] <= y[j +1]:
            temp = np.append(temp,float(gal[col]))
    n = len(temp)
    if switch == 1 and len(temp) > 0:
        temp = np.mean(temp)
    elif switch == 2 and len(temp) > 0:
        temp = np.median(temp)
    elif switch ==3 and len(temp) > 0:
        temp = sum(temp)
    elif switch == 10 and len(temp) > 3:
        temp = std(temp)
    else:
        temp = 0.0
    return temp,n

#Function to create bins and bin up, for whatever columb that is desired with control over the bin size and 
#switchs for mean and median
def pixelBin(cat,catHead,col, switch, delta):
    RAcol = where(catHead == "ra")[1][0]
    DECcol = where(catHead == "dec")[1][0]
    #create the bins in x and y / ra/dec
    x = np.arange(floor(min(cat[:,RAcol])), ceil(max(cat[:,RAcol])), delta)
    y = np.arange(floor(min(cat[:,DECcol])), ceil(max(cat[:,DECcol])), delta)
    X, Y = np.meshgrid((x + (delta/2.0)), (y + (delta/2.0)))
    Z = np.zeros((len(y),len(x)))
    Zden = np.zeros((len(y),len(x)))
    #cycle through x values first
    starttime = time.time()
    for i in range(0,len(x)-1):
        for j in range(0, len(y)-1):
            Z[j,i], Zden[j,i] = surfaceDencity(cat,catHead,col, switch,i,j,x,y)
        if i == 0:
            print "estimated time for", switch, " = ", ((time.time()-starttime)*len(x))/60.0
    return X,Y,Z,Zden
#save's any array to text comma delimited file
def savetext(folder,filename,array):
    outfile = open(pj(folder, filename), 'w')
    for row in array:
        count = 0
        for element in row:
            count += 1
            if count < len(row):
                x = str(element) + ","
            else:
                x = str(element)
            outfile.write(x)
        outfile.write("\n")
    outfile.close()
#function to print out a report on key values of selected columbs    
def reportA(A,s):
    x=A
    print "--------------------------------------------------------------------------------"
    print "Report ", s
    print "mean = ", mean(x), " median = ", median(x), " min = ", amin(x) , " max = ", amax(x)
    
def lineFit(x,*args):    
    xDat = args[0]
    yDat = args[1]
    yError=args[2]
       
    m = x[0]
    c = x[1]

    ci = 0.0
    for i in range(0,len(xDat)):
        fx = m*xDat[i] + c
        ci += ((yDat[i] - fx)**2)/ (yError[i]**2)
    return ci  

# Program for genorating a linar distribution         
def genLine(m,c, x):
    y = [] 
    for i in range(0,len(x)):
        y.append(x[i]*m + c)
    return y

##################################### Main Program #################################
#remove flagged values
if clean == 1:
    print 'cleaning data'
    cat = selecta(cat,where(catHead == "log_SFR_total")[1][0], -10.0,10.0)
#Makes an initial selection insuring that no edge effects scew local dencity calculation
if selection == 1:
    print 'making initial selections'
    cat = selecta(cat,where(catHead == "velocity")[1][0],vel_lower,vel_upper)
#find the local dencity using various parameters
if surfaceDen == 1:
    print 'making maps'
    col = where(catHead == "log_SFR_total")[1][0]
    X,Y,Z,Zden = pixelBin(cat,catHead,col, switch, delta)
    Zden = Zden / (delta**2)
    col = where(catHead == "velocity")[1][0]
    X,Y,Zvel, n = pixelBin(cat,catHead,col, 10, delta)
#report on Z and Zden
reportA(Z,"SFR")
reportA(Zden, "Dencity")
#turn from log sfr to sfr copying itself into a save arraray
Z_log =  Z
Z  = np.power(10,Z)
Z = clip(Z,0,2)
rows,cols = Z.shape
Z[:rows-1,:cols-1]


if fa == 1:
    # First for SFR
    # Take the fourier transform of the image.
    F1 = fftpack.fft2(Z)
    # Now shift the quadrants around so that low spatial frequencies are in
    # the center of the 2D fourier transformed image.
    F2 = fftpack.fftshift(F1)
    # Calculate a 2D power spectrum
    psd2D = np.abs( F2 )**2
    # Calculate the azimuthally averaged 1D power spectrum
    psd1D = radialprofile.azimuthalAverage(psd2D)
    
    #Secondly for Dencity
    # Take the fourier transform of the image.
    F1den = fftpack.fft2(Zvel)
    # Now shift the quadrants around so that low spatial frequencies are in
    # the center of the 2D fourier transformed image.
    F2den = fftpack.fftshift(F1vel)
    # Calculate a 2D power spectrum
    psd2Dden = np.abs( F2vel )**2
    # Calculate the azimuthally averaged 1D power spectrum
    psd1Dden = radialprofile.azimuthalAverage(psd2Dden)
    
    print 'creating cross power spectrum'    
    G1 = F2
    G2con = np.conjugate(F2vel)
    r = (G1*G2con) 
    R = r/abs(r)
    xps2D = np.real(R)**2
    xps = radialprofile.azimuthalAverage(xps2D)
    
guess = [10.0,10.0]
m,c = fmin(lineFit,guess, args =[Z.flatten(),Zvel.flatten(),sqrt(n.flatten())])
x = arange(min(Z.flatten()),max(Z.flatten()), 0.001)
y = genLine(m,c, x)

    
################################Plotting#############################################
RAcol = where(catHead == "ra")[1][0]
DECcol = where(catHead == "dec")[1][0]
RA = cat[:,RAcol]
DEC = cat[:,DECcol]

if simPlot == 1:
    plt.figure(1)
    plt.subplot(221)
    plt.imshow(Z, interpolation = 'nearest')
    plt.gca().invert_xaxis()
    plt.title('SFR Surface Dencity Image')
    
    plt.subplot(222)
    plt.imshow(Zvel, interpolation = 'nearest')
    plt.gca().invert_xaxis()
    plt.title('Galaxy Surface Dencity Image')

    plt.subplot(223)
    plt.semilogy( psd1D )
    plt.title('SFR Surface Dencity Image')
    
    plt.subplot(224)
    plt.semilogy( psd1Dden )
    plt.title('Galaxy Surface Dencity Image')
    #plt.tight_layout()
    
    plt.figure(2)
    plt.clf()
    plt.semilogy(xps)
    plt.title('cross power specturm')
    plt.show()

delta = str(delta)
savetext(folder,("SFR_Surface_Dencity_"+delta+"_binsize.csv"),Z)
savetext(folder,"Galaxy_Surface_Dencity"+delta+"_binsize.csv",Zden)    

plt.plot(Z.flatten(), Zvel.flatten(), 'kx')
plt.plt(x,y)
plt.show()

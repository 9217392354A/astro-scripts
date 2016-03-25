#Image/map maker
#Program created by Chris Fuller Aug2012 to create maps of starformation and local dencity

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
################################################################################
catFile = "tight_GalParam_030812_1129_spx6cff_0.csv"
folder = "/Users/chrisfuller/Dropbox/coma/AuxData/SFR/"
catFile = "Wall_galparam_dr9_spx6cff.csv"
#folder = "/export/home/spx6cff/coma/AuxData/SFR/"
#catFile= "tight_GalParam_030812_1129_spx6cff_0.csv"
#tasks localDen means surface plotting everything first if 0 then it will be drawn 
#from the saved cataloge
surfaceDen = 1
#switch: mean = 1, median = 2, sum = 3
switch = 1
#Bin size degree
delta = 2.0
#Make initial selection from catalogue **Selection not working yet requries code from catChanger**
#clean just selected from log sfr total -10<sfrtotal<10 remove if flag's to be included
#Selection if  selection = 0 then selection bypassed
selection = 1
clean = 1
vel_lower = 4000.0
vel_upper = 9000.0
ra_lower = 160.0
ra_upper = 209.0
dec_lower = 10.0
dec_upper = 36.0

#plotting options
PresPlot = 1
conPlot = 0
imPlot = 0
faPlot =0
simPlot = 1
#if fourrier analisis hit 1
fa = 1
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
    elif switch == 3 and len(temp) > 0:
        temp = sum(temp)
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
    
##################################### Main Program #################################
#remove flagged values
if clean == 1:
    print 'cleaning data'
    cat = selecta(cat,where(catHead == "log_SFR_total")[1][0], -10.0,10.0)
#Makes an initial selection insuring that no edge effects scew local dencity calculation
if selection == 1:
    print 'making initial selections'
    cat = selecta(cat,where(catHead == "velocity")[1][0],vel_lower-1000.0,vel_upper +1000.0)
    cat = selecta(cat,where(catHead == "ra")[1][0],ra_lower - 5.0,ra_upper + 5.0)
    cat = selecta(cat,where(catHead == "dec")[1][0],dec_lower - 5.0 ,dec_upper + 5.0)
#find the local dencity using various parameters
if surfaceDen == 1:
    print 'making maps'
    col = where(catHead == "log_SFR_total")[1][0]
    X,Y,Z,Zden = pixelBin(cat,catHead,col, switch, delta)
    Zden = Zden / (delta**2)
#report on Z and Zden
reportA(Z,"SFR")
reportA(Zden, "Dencity")
if selection == 1:
    cat = selecta(cat,where(catHead == "velocity")[1][0],vel_lower,vel_upper)
    cat = selecta(cat,where(catHead == "ra")[1][0],ra_lower,ra_upper)
    cat = selecta(cat,where(catHead == "dec")[1][0],dec_lower,dec_upper)
    
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
    F1den = fftpack.fft2(Zden)
    # Now shift the quadrants around so that low spatial frequencies are in
    # the center of the 2D fourier transformed image.
    F2den = fftpack.fftshift(F1den)
    # Calculate a 2D power spectrum
    psd2Dden = np.abs( F2den )**2
    # Calculate the azimuthally averaged 1D power spectrum
    psd1Dden = radialprofile.azimuthalAverage(psd2Dden)
    
    print 'creating cross power spectrum'    
    G1 = F2
    G2con = np.conjugate(F2den)
    r = (G1*G2con) 
    R = r/abs(r)
    xps2D = np.abs(R)**2
    xps = radialprofile.azimuthalAverage(xps2D)
    
################################Plotting#############################################
RAcol = where(catHead == "ra")[1][0]
DECcol = where(catHead == "dec")[1][0]
RA = cat[:,RAcol]
DEC = cat[:,DECcol]

if PresPlot == 1:
    V = [0.0,4.0,7.5,10.0,15.0,20.0,30.0,40.0,50.0,60.0]
    matplotlib.rcParams['xtick.direction'] = 'out'  
    matplotlib.rcParams['ytick.direction'] = 'out'
    plt.figure()
    #CS = plt.contourf(X,Y,Z,3,linestyles ="dotted", alpha = 0.3,nchunk=1)
    #plt.colorbar()
    Den = plt.contourf(X,Y,Zden,V,nchunk=0)
    Pos= plt.plot(RA,DEC, 'ko', markersize=1, alpha=0.2)
    #plt.colorbar()
    #plt.title("Local Dencity")
    plt.xlabel("RA", fontsize=18)
    plt.ylabel("DEC", fontsize=18)
    plt.xlim(190,210)
    plt.ylim(21, 37)
    plt.gca().invert_xaxis()
    plt.show()
    

if conPlot == 1:
    #V = log10(arange(10,ceil(max(Z.flat)),10))
    matplotlib.rcParams['xtick.direction'] = 'out'  
    matplotlib.rcParams['ytick.direction'] = 'out'
    plt.figure()
    CS = plt.contourf(X,Y,Z,3,linestyles ="dotted", alpha = 0.3,nchunk=1)
    plt.colorbar()
    Den = plt.contour(X,Y,Zden,5)
    Pos= plt.plot(RA,DEC, 'ko', markersize=1, alpha=0.2)
    plt.colorbar()
    plt.title("Contour Map SFR vs Local Dencity")
    plt.xlabel("RA", fontsize=18)
    plt.ylabel("DEC", fontsize=18)
    plt.xlim(160,210)
    plt.ylim(15, 37)
    plt.gca().invert_xaxis()
    plt.show()

if imPlot == 1:
    plt.imshow(Z, interpolation = 'bicubic')
    plt.title("Galaxy Surface Dencity")
    plt.gca().invert_xaxis()
    plt.show()
    plt.imshow(Zden, interpolation = 'bicubic')
    plt.title("SFR Surface Dencity")
    plt.gca().invert_xaxis()
    plt.show() 

    # Now plot up both
if faPlot == 1:
    # First for SFR
    py.figure(1)
    py.clf()
    py.imshow( Z, cmap=py.cm.Greys)
    plt.gca().invert_xaxis()
    py.title('SFR Surface Dencity Image')
    
    py.figure(2)
    py.clf()
    py.imshow( np.log10(psd2D))
    py.title('2D FFT of SFR Surface Dencity')
    
    py.figure(3)
    py.clf()
    py.semilogy( psd1D )
    py.xlabel('Spatial Frequency')
    py.ylabel('Power Spectrum')
    py.title('Azimultly averaged FFT of SFR Surface Dencity')
    
    #Secondly for Dencity
    py.figure(4)
    py.clf()
    py.imshow( Zden, cmap=py.cm.Greys)
    plt.gca().invert_xaxis()
    py.title('Galaxy Surface Dencity Image')
    
    py.figure(5)
    py.clf()
    py.imshow( np.log10(psd2Dden))
    py.title('2D FFT of Galaxy Surface Dencity')
    
    py.figure(6)
    py.clf()
    py.semilogy( psd1Dden )
    py.title('Azimultly averaged FFT of Galaxy Surface Dencity')
    py.xlabel('Spatial Frequency')
    py.ylabel('Power Spectrum')
    
    #show all fa plots
    py.show()

if simPlot == 1:
    plt.figure(1)
    plt.subplot(221)
    plt.imshow(Z, interpolation = 'bicubic')
    plt.gca().invert_xaxis()
    plt.title('SFR Surface Dencity Image')
    
    plt.subplot(222)
    plt.imshow(Zden, interpolation = 'bicubic')
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


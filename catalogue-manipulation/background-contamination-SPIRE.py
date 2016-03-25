#Program written by Chris Fuller January 2013
#program to find out contimanation in SPIRE images from rejected and accepted sources

import numpy as np
from os.path import join as pj
import atpy as at

#file stuff
folder = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/"

badgal = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/SPIRE-contamination/bad-gals.txt'
outfile= 'contamination-SPIRE-fornax.txt'

#read in catalogues
#print 'reading in catalogue'

badgals = np.loadtxt(badgal, dtype=str)
cat = at.Table(pj(folder,"final_outputs/HeFoCS-new_spire.fits"),type='fits')

#######################################Functions###################################
def binup(n,bins):
    n = np.array(n)
    N = []
    for b in bins:
        temp = n[np.where((n >= b[0]) & (n <= b[1]))]        
        N = np.append(N,len(temp))
    return N
    
def printResults(band,bins,Nr,Na,Nu, cont, Nc):
    print "-"*60
    print band, " um "
    print "-"*60
    print "Flux bin (mJy), Na  ,  Nr ,   Nu  , N ,  cont %, Nc  "
    for i in range(0,3):
        print bins[i], ", ", Na[i], ", ",Nr[i], ", ",Nu[i], ", ",Na[i]+Nr[i]+Nu[i], cont[i],", " ,Nc[i]
        
def contamination(beam,Nr,Na,Nu,Np):
    Nc = []
    cont = []
    for i in range(0,len(Nr)):
        #total number of sources in bin
        total1 = 242.0 #Nr[i]+Na[i] #+ total
        #total area of all appertures
        totalArea = total1 * beam
        #append number of cont sources by multiping by Np the numberof sources per sq degree
        Nc.append(int(np.round(Np[i]*totalArea, decimals=0)))
        cont.append((np.round((Np[i]*totalArea)*100.0/total, decimals=2)))
    return np.array(Nc), np.array(cont)

def writeheader(folder,filename):
    outfile = open(pj(folder, filename), 'w')
    #print header to table 
    print '\\begin{table}'
    print '\\begin{tabular}{cccccc} '
    print '\\hline'
    print 'Flux bin & Na & Nr & Np & contamination & Nc \\\\ '
    print '(mJy) & ~ & ~ & (deg$^{-2}$) & \\% & ~ \\\\ \\hline '
    outfile.close()

def s(x):
    return str(x)

def interger(x):
    try:
        x = int(x)
    except:
        x = "+++"
        print 'errrrrrror'
    return x
    
    
def saveResult(band,bins,Nr,Na,Np, Nu, cont, Nc, folder,filename):
    outfile = open(pj(folder, filename), 'w')
    print '~ & ~ & ~ & ~ & ~ & ~ \\\\'
    print '\\textbf{' + s(band) +'$\\mu$m} & ~ & ~ & ~ & ~ & ~ \\\\ '
    for i in range(0,len(Nr)):
        if i == len(Nr) -1:
            x = s(int(bins[i][0])) + '+ &' + s(int(Na[i])) + '&' + s(int(Nr[i])) + '&' + s(int(Np[i])) + '&' + s(cont[i]) + '&' + s(int(Nc[i])) + '\\\\ '
        else:
            x = s(int(bins[i][0])) + '-' + s(interger(bins[i][1])) + '&' + s(int(Na[i])) + '&' + s(int(Nr[i])) + '&' + s(int(Np[i])) + '&' + s(cont[i]) + '&' + s(int(Nc[i])) + '\\\\'
        print x
    outfile.close()


def writefooter(folder,filename):
    outfile = open(pj(folder, filename), 'w')
    print '\\hline '
    print '\\end{tabular} '
    print '\\end{table} '
    outfile.close()
    
    
###################################Main Program####################################

#bands = [500,350,250]
bands = [250,350,500]

i=0
writeheader(folder,outfile)
for band in bands:  
    total = 0
    na = []
    nr = [] 
    nu = []
    #bins are different for 250
    if band == 250:
        #bin widths
        bins = [[15.0,20.0],[20.0,45.0],[45.0,100.0],[100.0,1000000000.0]]
        #bins = [[15.0,60.0],[30.0,100.0],[100.0,1000.0]]
        #sources per square degree        
        #Np = np.array([1694.0, 1824.0, 313.0, 25.0]) #glen et al
        Np = np.array([260.8,578.75,59.5, ])
        #beam area        
        beam = np.pi*(17.6/3600.0)**2
    elif band == 350:
        bins = [[15.0,20.0],[20.0,45.0],[45.0,100.0],[100.0,1000000000.0]]
        #bins = [[15.0,60.0],[30.0,100.0],[100.0,1000.0]]
        #beam = np.pi*(24.0/3600.0)**2
        #Np = np.array([566.0,1209.0,154.0,73.0])#glen et al
        Np = np.array([79.65,340.6,37.9])
    elif band == 500:
        bins = [[15.0,20.0],[20.0,45.0],[45.0,100.0],[100.0,1000000000.0]]
        #bins = [[15.0,60.0],[30.0,100.0],[100.0,1000.0]]
        beam = np.pi*(36.0/3600.0)**2
        #Np = np.array([185.0,614.0,19.0,1.0])#glen et al
        Np = np.array([197.2,143.95,10.0])
    #no cycle through each galaxies only intrested in the point population
    for j in range(0,len(cat)):
        obj = str(cat.OBJECT[j])
        stype = cat.EXTENDEDNESS[j]
        flux = np.nan_to_num(cat['F'+str(band)][j])
        sn = np.nan_to_num(cat['SN'+str(band)][j])
        if stype[i] == "P":
            total += 1 
            #for galaxies that are rejected
            if (len(np.where(obj == badgals)[0]) == 1) and (float(sn) > 3.0):
                nr.append(float(flux)*1000.0)
            #for galaxies that were accepted
            elif (len(np.where(obj == badgals)[0]) == 0) and (float(sn) > 3.0):
                na.append(float(flux)*1000.0)
            else:
                nu.append(float(np.sqrt(flux*flux))*1000.0)
    #add 1 to the count        
    i+=1
    
    #now bin up and print/saveout out results
    Nr = binup(nr,bins)
    Na = binup(na,bins)
    Nu = binup(nu,bins)
    Nc, cont = contamination(beam,Nr,Na,Nu,Np)
    #print nu
    #print and save results to text file and to terminal
    #printResults(band,bins,Nr,Na,Nu, cont, Nc)
    saveResult(band,bins,Nr,Na,Np,Nu, cont, Nc, folder,outfile)
        
writefooter(folder,outfile)  
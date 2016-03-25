#Program created by Chris Fuller May 2012
#Program to take fun tools ds9 output and convert into flux created for ATLAS NGP
#Will need to be addapted for other maps
#v1.0 (9/5/12)

#import stuff
from numpy import *
import numpy
import scipy
import math
import sys
import os
from os.path import join as pj

#File stuff
cat = "coma.csv"
catout ="coma-Fixed-3sigma.csv" 
outFolder = "/Users/chrisfuller/Dropbox/coma/flux/finaloutputs/"
folder = "/Users/chrisfuller/Dropbox/coma/flux/"
files = os.listdir(folder)

cat = loadtxt(pj(folder,cat), dtype=str, delimiter=",")
# find files for respective jobbies
montycarloFiles = [x for x in files if x[-6:] == "mc.txt"]
apertures = [x for x in files if x[-5:] == "s.txt"]
backgrounds = [x for x in files if x[-5:] == "d.txt"]

########################################################################################################
#Functions

#This function takes the raw funtools output and then turn it into a flux
def funtoflux(filename, band):
    if int(band) == 250:
        x = 59.1
    elif int(band) == 350:
        x = 133.0
    elif int(band) == 500:
        x = 63.0
    elif int(band) == 160 or int(band) == 70 or int(band) == 100:
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
    fivesigma = 5*newsigma
    return fivesigma
        
def backgroundsubtraction(flux,backgrounds,n,band):
    table = array([(flux[0,0] + " backgroundsubtracted"),flux[0,0],"background", "threshold", "signel/noise",("detected"+band)],ndmin=2)
    if band == "160":
        error = n[0]
    elif band == "250":
        error = n[1]
    elif band == "350":
        error = n[2]
    elif band == "500":
        error = n[3]
    else: 
        raise 'no noise information present'
    count = 0  
    detected = 0
    for i in range(1,len(flux)):
        d =[]
        count += 1
        rawflux = float(flux[i,0])
        back = (float(backgrounds[i,0])/float(backgrounds[i,1]))*float(flux[i,1])
        fluxSub = rawflux - back
        fluxnoise = sqrt(abs(rawflux))*error        
        backnoise = sqrt(abs(back))*error
        totalnoise = sqrt(fluxnoise**2 + backnoise**2)
        SN = fluxSub / error
        if fluxSub >= 5.0*error:
            detected += 1
            d = 1
        else:
            d = 0
        selection = array([fluxSub, rawflux,back, error*3.0, SN, d],ndmin=2)
        table = append(table,selection, axis=0)
    print band," detections = ",detected, " above ", error*3.0
    return table 

def savetext(folder,filename,array):
    outfile = open(pj(folder, filename), 'w')
    for row in array:
        count = 0
        for element in row:
            count += 1
            if count < len(row):
                x = element + ","
            else:
                x = element
            outfile.write(x)
        outfile.write("\n")
    outfile.close()
            
    
#############################################################################       
 
noise160 = []           
noise250 = []
noise350 = []
noise500 = []
#find the noise
for y in montycarloFiles:
    flux = []
    band = y.split("-")
    flux = funtoflux(y,band[1])
    noise = noisefinder(flux)
    print "band ", band[1], "detection threshold 5sigma ", noise
    if band[1] == "160":
        noise160 = append(noise160,noise)  
    elif band[1] == "250":
        noise250 = append(noise250,noise)
    elif band[1] == "350":
        noise350 = append(noise350,noise)
    elif band[1] == "500":
        noise500 = append(noise500,noise)
    else:
        raise 'no noise information present'

n160 = mean(noise160)    
n250 = mean(noise250)
n350 = mean(noise350)
n500 = mean(noise500)
no = array([n160/5.0, n250/5.0, n350/5.0, n500/5.0],dtype=float)       

#background subtraction
f160 = [x for x in files if x[5:11] == "160-ap"]
b160 = [x for x in files if x[5:11] == "160-ba"]
f250 = [x for x in files if x[5:11] == "250-ap"]
b250 = [x for x in files if x[5:11] == "250-ba"]
f350 = [x for x in files if x[5:11] == "350-ap"]
b350 = [x for x in files if x[5:11] == "350-ba"]
f500 = [x for x in files if x[5:11] == "500-ap"]
b500 = [x for x in files if x[5:11] == "500-ba"]

fluxes = [f160,f250, f350, f500]
backs = [b160, b250, b350, b500]

for i in range(0,len(fluxes)):
    selection = []
    backa =[]
    back = []
    info = str(fluxes[i])
    b = info.split("-")
    band = b[1]
    f = str(fluxes[i])
    flux = funtoflux(f.split("\'")[1], band)
    backa = str(backs[i])
    back = funtoflux(backa.split("\'")[1], band)
    selection = backgroundsubtraction(flux,back,no,band)
    cat = append(cat,selection, axis=1) 

savetext(outFolder,catout,cat)  
print 'finished'







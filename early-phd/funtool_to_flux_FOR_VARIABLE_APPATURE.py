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
cat = "bigcoma.csv"
catfolder = "/Users/chrisfuller/Dropbox/coma/Catalogues"
catout ="coma-VARIABLE-3sigma-TESTONLY.csv" 
outFolder = "/Users/chrisfuller/Dropbox/coma/flux/finaloutputs/"
folder = "/Users/chrisfuller/Dropbox/coma/variable_app_circle/"
files = os.listdir(folder)
errorValues = loadtxt("/Users/chrisfuller/Dropbox/coma/flux/montycarlo2/valueserror.txt",delimiter=",")
cat = loadtxt(pj(catfolder,cat), dtype=str, delimiter=",")
# find files for respective jobbies
apertures = [x for x in files if x[3:] == "app"]
backgrounds = [x for x in files if x[3:] == "bac"]

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


#taking values for flux in jy or mjy and then taking away th ebackground and working out s/n
# this program requires an array to be defined globally call errorVal that has all the values 
# of sigma for combinded signel to noise.
def backgroundsubtraction(flux,backgrounds,band):
    table = array([(flux[0,0] + " backgroundsubtracted"),flux[0,0],"background", "detection threshold",("detected"+band)],ndmin=2)
    if band == "160":
        errorVal = errorValues[:,0]     
    elif band == "250":
        errorVal = errorValues[:,1]  
    elif band == "350":
        errorVal = errorValues[:,2]  
    elif band == "500":
        errorVal = errorValues[:,3]  
    else:
        raise 'no noise information present'
    count = 0  
    detected = 0
    for i in range(1,len(flux)):
        d =[]
        count += 1
        rawflux = float(flux[i,0])
        try:
            back = (float(backgrounds[i,0])/float(backgrounds[i,1]))*float(flux[i,1])
        except:
            back = 0.0
        fluxSub = rawflux - back
        sigmaFlux = errorVal[0]*float(flux[i,1])**2 + errorVal[1]*float(flux[i,1]) + errorVal[2]
        if fluxSub >= 5.0*sigmaFlux:
            detected += 1
            d = 1
        else:
            d = 0
        selection = array([fluxSub, rawflux, back, (5.0*sigmaFlux), d],ndmin=2)
        table = append(table,selection, axis=0)
    print band," detections = ",detected
    return table 

# for saving an array to a text file .csv
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

#background subtraction
f160 = [x for x in files if x[:7] == "app-160"]
b160 = [x for x in files if x[:8] == "back-160"]
f250 = [x for x in files if x[:7] == "app-250"]
b250 = [x for x in files if x[:8] == "back-250"]
f350 = [x for x in files if x[:7] == "app-350"]
b350 = [x for x in files if x[:8] == "back-350"]
f500 = [x for x in files if x[:7] == "app-500"]
b500 = [x for x in files if x[:8] == "back-500"]

fluxes = [f160, f250, f350, f500]
backs = [b160, b250, b350, b500]

for i in range(0,len(fluxes)):
    selection = []
    backa =[]
    back = []
    info = str(fluxes[i])
    b = info.split("-")
    band = b[1]
    band = band.split(".")[0]
    f = str(fluxes[i])
    flux = funtoflux(f.split("\'")[1], band)
    backa = str(backs[i])
    back = funtoflux(backa.split("\'")[1], band)
    selection = backgroundsubtraction(flux,back,band)
    cat = append(cat,selection, axis=1)

savetext(outFolder,catout,cat)  
print 'finished'







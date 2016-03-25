#Program created by Chris Fuller to take fluxes from robbies source extraction software compare between point and extended 
#Oct 2012

import numpy as np
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt

catFile1 = "ngp-flux-point.csv"
catFile2 = "ngp-flux-extended.csv"
folder = "/Users/chrisfuller/Dropbox/coma/Souce_Measurement/"
outName = "ngp_flux_exteneded_and_point.csv"
#######Functions###########
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
####################read in catalogue############################
#Extended
catE = np.loadtxt(pj(folder,catFile2), dtype=str,  skiprows=1, delimiter=",", unpack=False)
catHead = np.loadtxt(pj(folder,catFile1), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(catE[0,:]))
#Point
catP = np.loadtxt(pj(folder,catFile1), dtype=str,  skiprows=1, delimiter=",", unpack=False)
newcatHead = catHead

x=1.0/60.0
#wave bands to complete
bands = ["500","350","250"]
beam = [36.3*x,24.9*x,18.2*x]

j = -1
#loop over all bands
for band in bands:
    print "starting band ", band,"um"
    j += 1
    best_flux_col = []
    ptoe = 0
    etop = 0
    #loop through all values in table    
    for i in range(0,len(catP)):
        best_flux = 0.0
        #find the columns for Flux, Error, and Extendedness
        Fcol = np.where(catHead == "F"+band)[1][0]
        Ecol = np.where(catHead == "E"+band)[1][0]
        Tcol = np.where(catHead == "EXTENDEDNESS")[1][0]
                
        # First Check if one is a point and one is exteneded, if both are the same check the values are the same then insert as best_flux
        if catP[i,Tcol][j] == catE[i,Tcol][j]:
            best_flux = catP[i,Fcol]
            pORe = catP[i,Tcol][j]
            if (float(catE[i,Fcol]) - float(catE[i,Fcol])) > 0.001: 
                print "error: different flux with same method, CCC",i,"  in band  ", band, "um"
                raise 
        # Now test if they are different if they are work out which is larger and change best flux according as well as saying what type of source it is
        # and if they are different by greater than the error        
        elif catP[i,Tcol][j] != catE[i,Tcol][j] and (float(catP[i,Fcol])-float(catE[i,Fcol]))**2 >= (float(catP[i,Ecol]) + float(catE[i,Ecol]))**2:
            if float(catP[i,Fcol]) > float(catE[i,Fcol]):
                best_flux = catP[i,Fcol]
                pORe = "P"
                #print "Galaxy ", catP[i,np.where(catHead == "OBJECT")[1][0]], "extended = ", catE[i,Fcol], "to point = ", catP[i,Fcol], "radius = ", float(catE[i,np.where(catHead == ("R"+band))[1][0]])/beam[j]
                etop += 1 
                
            elif float(catE[i,Fcol]) > float(catP[i,Fcol]):
                best_flux = catE[i,Fcol]
                pORe = "E"
                #print "Galaxy ", catP[i,np.where(catHead == "OBJECT")[1][0]], "point = ", catP[i,Fcol]," to extended = ", catE[i,Fcol], "radius = ", float(catE[i,np.where(catHead == ("R"+band))[1][0]])/beam[j]
                ptoe += 1
        else:
            #print "error, no value?", "Galaxy ", catP[i,np.where(catHead == "OBJECT")[1][0]], "point = ", catP[i,Fcol]," extended = ", catE[i,Fcol], "radius = ", float(catE[i,np.where(catHead == ("R"+band))[1][0]])/beam[j]
            #print "Point, ", catP[i,Tcol][j], "Extended ", catE[i,Tcol][j]
            best_flux =  catE[i,Fcol]
        best_flux_col = np.append(best_flux_col, best_flux)
    
    newcatHead = np.append(newcatHead, ("best_F"+str(band)))
    catP = np.append(catP, best_flux_col.reshape(len(best_flux_col),1), axis=1)
    print "number of points to extended = ", ptoe, " number of extended to point = ", etop
    
newcatHead = newcatHead.reshape(1,len(newcatHead))
newcat= np.append(newcatHead, catP, axis=0)
best_fluxs =  newcat[:,-3:]
objects = newcat[:,:2]
new = np.append(objects, best_fluxs,axis=1)
best_fluxs = new[:,1:]
savetext(folder,outName, newcat)
savetext(folder,"best_fluxs.csv",best_fluxs)
print "FINISHED!"


#Chris Fuller, 10 Dec 2012
# Program for making coma tables into the right format


import numpy as np
import math
import os
from os.path import join as pj

#key variables
bands = ["500_BEST","350_BEST","250_BEST"]

#files/folders
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs"
clusterfile = 'coma_cluster-full-101212.csv'
filamentfile = 'filament-full-101212.csv'
outfileC = 'paper-only-table-CCC.txt'
root = "/Users/chrisfuller/Dropbox/phd/herchel/coma/source-measurement-work/"


#read in catalogues
print 'reading in catalogue'
clusterCat = np.loadtxt(pj(folder,clusterfile), dtype=str,  skiprows=1, delimiter=",", unpack=False)
filamentCat = np.loadtxt(pj(folder,filamentfile), dtype=str,  skiprows=1, delimiter=",", unpack=False)
catHead = np.loadtxt(pj(folder,clusterfile), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(clusterCat[0,:]))

#######################################Functions###################################
def createObject(cat, catHead, prefix):
    col = np.where(catHead == "OBJECT")[1][0]
    oldcol = cat[:,col]
    newcol = []
    for i in range(0,len(oldcol)):
        newObj =[]
        newObj = prefix + str(i+1)
        newcol = np.append(newcol, newObj)
    newcol = np.append("Object",newcol)
    return newcol.reshape(len(newcol), 1)

def fluxCols(cat, catHead, bands, badfiles):
    #for each band loop through and return flux(error) or upper limit
    obcol = np.where(catHead == "OBJECT")[1][0]

    for band in bands:
        Fcol = np.where(catHead == 'F'+band)[1][0]
        Ecol = np.where(catHead == 'E'+band)[1][0]
        
        #now loop through each line of the cat working out the string value 
        newcol = []
        for i in range(0,len(cat)):
            flux = cat[i,Fcol]
            if flux == '0.0'or flux =='':
                newflux = "$<$("+cat[i,Ecol][:5]+")"
            else:
                newflux = flux[:5] +"("+cat[i,Ecol][:5]+")"
                
            #check if it was one of th fluxes excluded from cat by-eye if so *
            if len(np.where(badfiles == cat[i,obcol])[0]) == 1:
                newflux = newflux + '*'
            newcol= np.append(newcol, newflux)
        colname = "S_{" + band[:3]+"}($\sigma$S_{" + band[:3]+"})"   
        newcol = np.append(colname,newcol)
        if band == bands[0]:
            fluxes = newcol.reshape(len(newcol), 1)
        else:
            fluxes = np.append(fluxes,newcol.reshape(len(newcol), 1),axis=1)
    return fluxes

def positionCols(cat,catHead):
    raCol = np.where(catHead == 'GRA2000')[1][0]
    decCol = raCol + 1
    RA = np.append("RA",cat[:,raCol])
    DEC = np.append("Dec",cat[:,decCol])
    return np.append(RA.reshape(len(RA),1),DEC.reshape(len(RA),1),axis=1)

def typeCol(cat,catHead):
    ecol = np.where(catHead == 'ELLIPTICAL')[1][0]
    scol = np.where(catHead == 'SPIRAL')[1][0]
    typecol = []
    for i in range(0,len(cat)):
        if cat[i,scol] == "1":
            t = 'S'
        elif cat[i,ecol] == "1":
            t = 'E'
        else:
            t = "-"
        typecol = np.append(typecol, t)
    typecol = np.append("Type", typecol)
    return typecol.reshape(len(typecol),1)

def colMerge(cols):
    a = []
    for col in cols:
        if col.flat[0] == cols[0][0]:
            a = col
        else:
            a = np.append(a,col,axis=1)
    return a
    
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
    
def saveLaytex(folder,filename,array):
    outfile = open(pj(folder, filename), 'w')
    for row in array:
        outfile.write("        ")
        for element in row:
            if element != row[-1]:
                x = element + "&"
            else:
                if row[0] == array[0,0]:
                    x = element + "\\\\ \hline"
                else:
                    x = element +"\\\\"
            outfile.write(x)
        outfile.write('\n')
    outfile.close()
    
def file2name(filename):
    return filename.split("-")[0]

def filelist2list(filelist):
    new = []
    for filename in filelist:
        new = np.append(new,file2name(filename))
    return new    
        
###################################Main Program####################################

#create list of good and bad detections 
checkedfiles = os.listdir(pj(root,"checked"))
badfiles = os.listdir(pj(root,"faluse-detections"))

#remove any other file types from list
checkedfiles = [x for x in checkedfiles if x[-2:] == "ps"]
badfiles     = [x for x in badfiles if x[-2:] == "ps"] 

#change filename to gal name: 
checkedfiles = filelist2list(checkedfiles)
badfiles = filelist2list(badfiles)

#create new catlogue to write out!
cols =  createObject(clusterCat, catHead, "CCC"),positionCols(clusterCat,catHead),typeCol(clusterCat, catHead), fluxCols(clusterCat, catHead, bands, badfiles)
clusterCat = colMerge(cols)

#create new catlogue to write out!
fcols =  createObject(filamentCat, catHead, "FGC"),positionCols(filamentCat,catHead),typeCol(filamentCat, catHead), fluxCols(filamentCat, catHead, bands, badfiles)
filamentCat = colMerge(fcols)



saveLaytex(folder,outfileC,clusterCat)
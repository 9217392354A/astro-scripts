#final catalogue maker for the source extraction programs. 
#Chris Fuller, 4 Dec 2012

import numpy as np
import os
from os.path import join as pj

#files/folders
folder = "/Users/chrisfuller/Desktop/checked/"
pointfolder = '/Users/chrisfuller/Desktop/force_as_point/'
catfile = 'NGP-flux-best-afterpointforce.csv'
pointcat = '/Users/chrisfuller/Dropbox/coma/Souce_Measurement/ngp-flux-point.csv'
outname = 'NGP-fluxes-final-041212.csv'
####################################Functions#############################################
def data2zero(row,catHead):
    bands = [500,350,250]
    ext = np.where(catHead == "EXTENDEDNESS")[1][0]
    for band in bands:
        col = np.where(catHead == str("F"+str(band)+"_BEST"))[1][0]
        row[col] = ""
        col = np.where(catHead == str("R"+str(band)+"_BEST"))[1][0]
        row[col] = ''
    row[ext] = "PPP"
    return row
    
def data2point(row,i,cat,catHead):
    bands = [500,350,250]
    ext = np.where(catHead == "EXTENDEDNESS")[1][0]
    row = cat[i,:]
    for band in bands:
        col = np.where(catHead == str("R"+str(band)+"_BEST"))[1][0]
        row[col] = ''
    return row
    
    #row[ext] = "PPP"
    
    
 # for saving an array to a text file .csv
def savetext(folder,header,filename,array):
    outfile = open(pj(folder, filename), 'w')
    array = np.append(header, array,axis=0)
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
    
        
###########################Main Program######################################################
    
#read in catalogues
print 'reading in catalogue'
cat = np.loadtxt(pj(folder,catfile), dtype=str,  skiprows=1, delimiter=",", unpack=False)
catHead = np.loadtxt(pj(folder,catfile), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(cat[0,:]))
newcat = cat
pointcat = np.loadtxt(pointcat, dtype=str,  skiprows=1, delimiter=",", unpack=False)

#create list of PSW checked galaxies. 
d250=[]
listd250 = np.sort([x for x in os.listdir(folder) if x[-2:] == "ps"])
for obj in listd250: 
    d250 = np.append(d250, obj.split('-')[0])
    
#create list of point sources. 
points=[]
listpoints = np.sort([x for x in os.listdir(pointfolder) if x[-2:] == "ps"])
for obj in listpoints: 
    points = np.append(points, obj.split('-')[0])

name = np.where(catHead == "OBJECT")[1][0]
flux = np.where(catHead == "OBJECT")[1][0]


i = -1
for row in cat:
    i += 1
    if len(np.where(d250 == row[name])[0]) > 0:
        continue
    elif len(np.where(points == row[name])[0]) > 0:
        newcat[i,:] = data2point(row,i,pointcat, catHead)
        
    else:
        newcat[i,:] = data2zero(row,catHead)


savetext(folder,catHead,outname, newcat)
        
        
        #len(np.where(points == row[name])[0]) > 0:
        
        
        
    






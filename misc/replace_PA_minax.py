#code for replacing PA and minaxis
#Chris Fuller, Jan 2013

import numpy as np
from os.path import join as pj

#files and folders
#files/folders
folder = "/Users/chrisfuller/Desktop/"
clusterfile = 'fcc_insurvey_input_ready.csv'
filename = "fcc_insurvey_opticaladded.csv"
#import cat
print 'reading in catalogue'
cat = np.loadtxt(pj(folder,clusterfile), dtype=str,  skiprows=1, delimiter=",", unpack=False)
catHead = np.loadtxt(pj(folder,clusterfile), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(cat[0,:]))


########### Functions ###############
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
#######################################
#cols
pa = np.where(catHead == "PA")[1][0]
a = np.where(catHead == "FULLMAJAX")[1][0]
b = np.where(catHead == "FULLMINAX")[1][0]
col = np.where(catHead == "OBJECT")[1][0]
for i in range(0,len(cat)):
    maja = cat[i,a]
    mina = cat[i,b]
    pan = cat[i,pa]
    gal = cat[i,col]
    
    if pan == "":
        cat[i,pa] = "0.0"
    if mina == "":
        cat[i,b] = cat[i,a]
    
array = np.append(catHead,cat,axis=0)
savetext(folder,filename,array)

        
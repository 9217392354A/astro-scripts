#program to gen a list for supercomos extraction
#Chris Fuller Feb 2012
import os
import numpy as np
from os.path import join as pj
#import atpy as at #cat = at.Table(pj(folder,fileName),type='fits')

##################################### varaibles ##############
#len of square to withdraw from supercosmos
X = 5.0 #arcmin

####################################### inputs ##############
folder = "/Users/chrisfuller/Desktop/superCOSMOS"
fileName = "fcc_all.csv"

#gen cat from text file with a header
cat = np.loadtxt(pj(folder,fileName), dtype=str,  skiprows=1, delimiter=",", unpack=False)
catHead = np.loadtxt(pj(folder,fileName), dtype=str,  skiprows=0, delimiter=",", unpack=False)[0,:].reshape(1,len(cat[0,:]))
###################################### Functions #############
 # for saving an array to a text file 
def savetext(folder,filename,array,c):
    outfile = open(pj(folder, filename), 'w')
    for row in array:
        count = 0
        for element in row:
            count += 1
            if count < len(row):
                x = element + c
            else:
                x = element
            outfile.write(x)
        outfile.write("\n")
    outfile.close()
    
#finds the number possible for a given x size
def possNum(x):
    return np.floor(1000.0/(x**2))
    
#finds the number folders needed
def dSec(x):
    return np.ceil(len(cat)/possNum(x))
    
############################## Main Program ################
# work out how many file groups needed
N = dSec(X)
perFol = possNum(X) 
print 'number of groups needed to download supercosmos data...', N

#asign each galaxy a number from 1 to N refering to its designation
sec = 1
fol = np.array([], dtype=int)
for i in range(0,len(cat)):
    #only increass the dif
    if i % perFol == 0.0 and i != 0:
        sec += 1 
    fol = np.append(fol,sec)
fol = fol.reshape(len(fol),1)
cat = np.append(cat,fol, axis=1)

#now run through each galaxy and turn ra and dec into supercosmos friendly formats
for i in range(0,len(cat)):
    cat[i,1] = cat[i,1].replace(":"," ")
    cat[i,2] = cat[i,2].replace(":"," ")
        
#loop through each section writing a list to disk in each folder to upload
for i in range(1,int(N)+1):
    f = pj(folder, "super-" + str(i))
    try:
        print 'making dir...', f
        os.mkdir(f)
    except:
        print 'dir already exits'
    
    #select new array
    loc = np.where(cat[:,3] == str(i))[0]
    #array for upload
    upload = cat[loc,1:3]
    #array for next program recording all needed info
    record = cat[loc,:]
    
    #save each array to the respective folder
    savetext(f,"upload-" + str(i) + ".txt", upload, " ")
    savetext(f,"record-" + str(i) + ".txt", record, " ")
    
#save entire cat to desktop
savetext(folder, "superCOSMOS.csv", cat, ",")

print "program finished successfully..."
    
    
        



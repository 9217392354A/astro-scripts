#creates extra coulmbs for catalogues
# Chris Fuller July 2012
from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt 
import time
################user varaibles##################################
catFile = "Huge_GalParam_030812_1114_spx6cff.csv"
folder = "/Users/chrisfuller/Dropbox/coma/AuxData/SFR/"

#Tasks - set to 1 to includ task
loctask = 1
SFRtask = 0

#Selection if  selection = 0 then selection bypassed
selection = 0
vel_lower = 4000.0
vel_upper = 9000.0
ra_lower = 165.0
ra_upper = 220.0
dec_lower = 10.0
dec_upper = 46.0
##################Initilising Variables########################

####################read in catalogue############################
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

#Function that takes a standard recipie for Halpa and converts it into SFR, z=redshift, d=distance modulus
def SFRcalc(Ha_ew, Hb_ew, rPetro, rFiber, Ha_cont, Hb_cont, z):
    # find each line flux with steller abobsobtion flux's
    # NOTE: we can inprove this by making a correction that vaires with type
    d = 102.0 * 3.08568025E22
    Sa = ((Ha_ew + 1.3)/Ha_ew)*Ha_cont
    Sb = ((Hb_ew + 1.3)/Hb_ew)*Ha_cont
    #Aperture correction 
    A = 10**(-0.4*(rPetro-rFiber))
    #Dust correction
    D = ((Sa/Sb)/2.86)**2.114
    #Flux with Corrections
    F= Sa*A*D
    #Lumiocity of line
    L= F*4*pi*d**2
    #Find SFR using standard recipie
    SFR = L / 1.27E34*1e18
    return SFR

#Function taking SFRcalc and looping over it for all rows of catalogue
def SFR(cat,catHead):
    #find columbs with respective headings
    Ha = where(catHead == "Ha_ew")[1][0]  
    Hb = where(catHead == "Hb_ew")[1][0]
    Ha_con = where(catHead == "Ha_continuum")[1][0]
    Hb_con = where(catHead == "Hb_continuum")[1][0]
    rpet = where(catHead == "petroMag_r")[1][0]
    rfib= where(catHead == "fiberMag_r")[1][0]
    z = where(catHead == "velocity")[1][0]
    Ha_err = where(catHead == "Ha_ewErr")[1][0]
    Hb_err = where(catHead == "Hb_ewErr")[1][0]
    SFRall = []
    # Loop over either cat calling SFRcalc to make caculations replacing with 0 if less than 5 sigma
    #detection and replaceing with -9999 if poor spectrum
    i=0
    for row in cat:

        if row[Ha] == -9999.0 or row[Hb] == -9999.0:
            print "Warning Flag: poor spectrum"
            SFRgal = 0.0
        elif (row[Ha]/row[Ha_err]) < 5.0 or (row[Hb]/row[Hb_err]) < 5.0:
            SFRgal = 0.0
            i += 1
            
        else:
            SFRgal = SFRcalc(row[Ha],row[Hb],row[rpet],row[rfib],row[Ha_con],row[Hb_con],row[z])
        SFRall = append(SFRall, SFRgal)
        
    print i
    return SFRall
    
    
#Function to work out local dencity from the distance to the Nth nearest
#meber above Mr >= -20 or at the distance of coma mr > 15
#Obj is a row in cat
def localDencity(cat,catHead):
    #Create a list of objects brighter than mr>= 15 within the vel range
    # of obj and within 0.5 degree    
    # decreass threshold to increass speed.         
    threshold = 5.0
    # conversion between deg and Mpc this is fixed at the distance of the x-ray center but a more complicated
    # system could be implimented baced on redshift at some future point.
    conversion = 1.798
    Local = array([0,0,0]).reshape(1,3)
    RAcol = where(catHead == "ra")[1][0]
    Deccol = where(catHead == "dec")[1][0]
    velCol = where(catHead == "velocity")[1][0]
    mrCol = where(catHead == "r")[1][0]
    count = 0
    starttime= time.time()
    counter = 0 
    for Obj in cat:
        RaObj = Obj[RAcol]
        DecObj = Obj[Deccol]
        VelObj = Obj[velCol]
        n = []
        counter += 1
        if counter == 100:
            deltaT = time.time() - starttime
            estimate = (len(cat)/100.0)*deltaT
            print "estimated time for localDencity calculation = ", round(estimate/60.0), "minutes"
        if counter%100 ==0:
            print round((float(counter)/float(len(cat)))*100.0, 2),"% complete"       
        for gal in cat:
            RaGal = gal[RAcol]
            DecGal = gal[Deccol]
            VelGal = gal[velCol]
            mr = gal[mrCol]
            radius = sqrt(((RaObj - RaGal)*math.cos(DecObj))**2 + (DecObj - DecGal)**2)
            if radius <= threshold and mr >= 15.0 and (VelObj - 500.0) <= VelGal <= (VelObj + 500.0) and radius != 0.0:
                n = append(n, radius)
        if len(n) < 10:
            count += 1
            n = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        n = sort(n)   
        D3 = n[2]*conversion
        D5 = n[4]*conversion
        D10 = n[5]*conversion
        SIGMA3 = 3 / math.pi*D3**2
        SIGMA5 = 5 / math.pi*D5**2
        SIGMA10= 10 / math.pi*D10**2
        temp = []
        temp = [SIGMA3, SIGMA5, SIGMA10]
        temp = array(temp).reshape(1,3)
        Local = append(Local,temp, axis=0)
    Local = delete(Local,0,0)
    print "number of galaxies with less than 10 galaxies within ", threshold, "degs = ", count
    return Local
    
#save's any array to text comma delimited file
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
#function to print out a report on key values of selected columbs    
def report(cat, col, catHead):
    x=cat[:,col]
    print "--------------------------------------------------------------------------------"
    print "Report ", catHead[0,col]
    print "mean = ", mean(x), " median = ", median(x), " min = ", min(x) , " max = ", max(x)
##################################### Main Program #################################
#Makes an initial selection insuring that no edge effects scew local dencity calculation
if selection == 1
    cat = selecta(cat,where(catHead == "velocity")[1][0],vel_lower-1000.0,vel_upper +1000.0)
    cat = selecta(cat,where(catHead == "ra")[1][0],ra_lower - 5.0,ra_upper + 5.0)
    cat = selecta(cat,where(catHead == "dec")[1][0],dec_lower - 5.0 ,dec_upper + 5.0)
#if set to true works out local dencity using closest member calculations and SFR from Ha line
if loctask == 1:
    print "Commencing Local Dencity Task"
    local = localDencity(cat,catHead)
    newCols = ["SIGMA3","SIGMA5","SIGMA10"]
    catHead = append(catHead,newCols)
    catHead = catHead.reshape(1,len(catHead))
    cat = append(cat,local, axis=1)
if SFRtask == 1:
    print "Commencing SFR Task"
    SFRGalaxies = SFR(cat,catHead).reshape(len(cat),1)
    newCol = ["SFR_Ha"]
    catHead = append(catHead,newCol)
    catHead = catHead.reshape(1,len(catHead))
    cat = append(cat,SFRGalaxies, axis=1)
outCatFull = append(catHead, cat, axis=0)
#final selection 
cat = selecta(cat,where(catHead == "velocity")[1][0],vel_lower,vel_upper)
cat = selecta(cat,where(catHead == "ra")[1][0],ra_lower,ra_upper)
cat = selecta(cat,where(catHead == "dec")[1][0],dec_lower,dec_upper)
outCat = append(catHead, cat, axis=0)

#write both full and selected catalogues to disk
print "writing file to disk"
savetext(folder,"SFR_Local_Dencity_selected.csv", outCat)
savetext(folder,"SFR_Local_Dencity_full.csv", outCatFull)

#report to uses vital satistics for columbs in local dencity
report(cat,where(catHead == "SIGMA3")[1][0],catHead)
report(cat,where(catHead == "SIGMA5")[1][0],catHead)
report(cat,where(catHead == "SIGMA10")[1][0],catHead)    
report(cat,where(catHead == "log_SFR_total")[1][0],catHead) 
report(cat,where(catHead == "reliable")[1][0],catHead) 
#report time talken
timeTaken= time.time() - starttime
print "finished; task took ", round((time.time() - starttime)/60.0), "minutes"    
        
        
        
        
    
    
    
    
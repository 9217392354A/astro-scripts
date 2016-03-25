# program reads in a fits table and either a list or a folder and then sets all galxies
# not detected at 250 or declared as faluse detections and removes them and sets the col
# remove_250 ect to 1 so that later a latex table can be created easily

# Chris Fuller, Mar 2013

import os
import numpy as np  
from os.path import join as pj
import atpy as at

################ inputs ################
badFolder = "/Users/chrisfuller/Desktop/output/bad/"
catFolder = "/Users/chrisfuller/Desktop/"

cat = at.Table(pj(catFolder,"HeFoCS-fluxes-260313-mybgsub-v2-handmeasure.fits"),type='fits')

bands = ['500','350','250','160','100']
################ functions #############

# takes a list of galaxies and turns it into a list of names removing duplicates removing duplicates
def folder2list(folder):
    galaxies = [s for s in os.listdir(folder) if "postscript" in s]
    l =[]
    for gal in galaxies:
        l.append(str(gal.split("-")[0]))
    l = np.array(l)
    l.sort()
    return l

#find if detected at 250
def detectionTest(x):
    try: y = float(x)
    except: return 0
            
    if y  == 0.0: return 0
    elif y > 0.0: return 1
    else: raise 'problem with flux = ',y 
        
# find if on the bad list
def badTest(x,blist):
    
    test = len(np.where(blist == str(x))[0])
    if test==0: return 0
    elif test ==1: return 1
    else: raise 'problem with gal name in bad list'


# removes bad galaxies from catalogue and then returns cat with new columbs remove_250    
def removeBadgals(t, badG):
    # test if table already has out column col  
    for band in bands:
        test = [s for s in t.columns if "remove-"+band in s]
        if len(test) == 0:
            t.add_empty_column('remove-'+band, np.int16)
    #loop through entire table
    for i in range(0,len(t)):
        dtest = detectionTest(t.F250[i]) #test if detected at 250
        btest = badTest(t.OBJECT[i],badG)         #test if in bad gal list
        
        if t.HEVICS_PLW[i] == 0: continue #not in either pacs or spire
        
    
        if   dtest == 1 and btest == 0: continue # not in bad list and detected at 250 no action required
        elif dtest == 0 or btest == 1: # not detected at 250 or on bad list so set all to zero
            for band in bands:
                #test if it was detected anyway
                band_dtest = detectionTest(t['F'+band][i])
                
                if band_dtest == 1: 
                    t['remove-'+band][i] = 1 #was detected in band and had been removed
                    t['E'+band][i] = t['E'+band][i]*3.0 
                    print t['remove-'+band][i], " flux: ", t['F'+band][i], " band: ",band
                
                if t.HEVICS_PLW[i] == 3 and (band=="160" or band =="100"): val = 0.0
                else: val = 0.0
                #set all bands to zero unless pacs
                t['F'+band][i] = val                
                t['R'+band][i] = val
                t['SN'+band][i] = val 
    return t
    
################# main ################
bad = folder2list(badFolder)
new_cat = removeBadgals(cat, bad)
################# outputs #############

new_cat.write(pj(catFolder,"HeFoCS-fluxes-260313-mybgsub-v2-handmeasure-bad-detections-removed.fits"), overwrite=True)
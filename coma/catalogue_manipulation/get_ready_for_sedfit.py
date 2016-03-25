#get table ready for matts prog
#chris fuller, april 2013

import os
import numpy as np  
from os.path import join as pj
from atpy import Table

#get table

###################### inputs ######################
#cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12.fits")
cat = Table("/Users/chrisfuller/Dropbox/IC4051.fits")

#variables
bands = ['PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']
instrument = ['PACS', 'PACS', 'SPIRE', 'SPIRE', 'SPIRE']
distance = 100.0
pacs_cal = 0.0
spire_cal= 0.095


# write out 
outFile = open('/Users/chrisfuller/Desktop/sed_fit.csv', 'w')

#create header
hdr = "OBJECT, z, z_error, distance, distance_error, "

#add flux headers to hdr master
for band in bands:
    if band != bands[-1]: hdr += band + '_w, ' + band + ', '  + band +'_e, ' #with final comma 
    else: hdr += band + '_w, ' + band + ', ' + band +'_e ' #without comma

outFile.write( hdr + " \n")


#now add data to line and write out row wise
for i in range(0, len(cat)):
    #first add obj to line
    line = str(cat.OBJECT[i]) + ', '

    #now add fixed data
    line += '0.0, 0.0, ' + str(distance) + ', 0.0, '

    #now cycle through fluxes and enter values into file
    for band in bands:

        #figure out final charater bassed on newline and comma
        if band != bands[-1]: final = ', '
        else: final =' \n'

        #set sigma cal
        if band[1] == 'S': sigma_cal = spire_cal
        else: sigma_cal = pacs_cal

        my_band = band[-3:]

        #add wavelenght col
        line += my_band + '.0, ' 

        #undetected 
        if np.nan_to_num(cat['F' + my_band][i]) == 0.0: 
            line += '-, -' + final
        
        #detected
        elif np.nan_to_num(cat['F' + my_band][i]) > 0.0:

            #remove cal error qudratically
            error = np.sqrt( (cat['E' + my_band][i]**2)  - ((cat['F' + my_band][i]*sigma_cal)**2))
            if str(error) == 'nan': error = cat['E' + my_band][i]

            #add flux and error to line
            line += str(cat['F' + my_band][i]) + ', ' + str(error) + final

        else:
            raise 'error with flux reconition'


    #write line to file
    outFile.write(line) 

#close outfile 
outFile.close()

print 'program complete'



# Program to grab SDSS images for the HRS matched to size of SPIRE PLW
# written 29th May 2012

# import modules
import numpy as np
import os
import sys
from os.path import join as pj
from atpy import Table

print 'reading in cats'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,cat_name))

cat = cat.where((cat.goldmine < 2) & (cat.goldmine > -10) & (cat.late == 1))

npix=400
pixScale = 0.4

# outfolder
outFolder = '/Users/chrisfuller/Desktop/'

os.chdir(outFolder)

for i in range(len(cat)):
    name = cat.OBJECT[i]
    coord = cat.GRA2000[i], cat.GDEC2000[i]

 
    outFile = name  + "-SDSS.jpg"


    
    command = 'wget -A.jpg "http://skyservice.pha.jhu.edu/DR8/ImgCutout/getjpeg.aspx?ra={0:5f}%20%20&dec={1:5f}%20%20%20&scale={2:5f}&width={3:}&height={4:}&opt=&query=" -O '\
              .format(coord[0], coord[1], pixScale, npix, npix) + outFile
    os.system(command)
    print command
print "Program Finished Successfully"
# Program to grab SDSS images for the HRS matched to size of SPIRE PLW
# written 29th May 2012

# import modules
import pyfits
import numpy
import os
import sys
from os.path import join as pj

### select inputs

# fits image folder
fitsFolder = '/Users/chrisfuller/Desktop/fits-files/'

# outfolder
outFolder = '/Users/chrisfuller/Desktop/jpeg-files/'

try:
    os.mkdir(outFolder)
except:
    print "dir already created"

# select quality factor (i.e., pix size SDSS compared to PSW)
# e.g. a factor of 6 gives 1'' pixels
factor = 12

# Website sets a maximum (enter here)
maximum = 2048

#######################################################################################################

# search FITS folder
fitsFiles = os.listdir(fitsFolder)

# change directory to output folder
os.chdir(outFolder)

# loop over each file -  only download if PSW map exists
for file in fitsFiles:
    # check if fits file
    if file[-5:] != ".fits":
        continue
    
    # check if PSW image
    if file.count("PSW") == 0:
        continue
    elif file.count("PSW") > 1:
        raise "Problem with file name"
    
    #  load fits file to get needed info
    fits = pyfits.open(pj(fitsFolder, file))
    
    # get header for image
    header = fits[0].header
    
    # calculate dimensions and position of requested SDSS image
    # correct crval to middle position
    midPix = [(header["naxis1"] + 1) / 2.0, (header["naxis2"] + 1) / 2.0]
    try:
        cdelt1 = header["cdelt1"]
        cdelt2 = header["cdelt2"]
    except:
        raise "You need to Program harder!!"
    deltaPix = [midPix[0] - header["crpix1"], midPix[1] - header["crpix2"]]
    coord = [header["crval1"] + deltaPix[0] * cdelt1, header["crval2"] + deltaPix[1] * cdelt2]
    
    # calculate dimensions of image and pixel scale
    npix = [header["naxis1"]*factor, header["naxis1"]*factor]
    pixScale = cdelt2*3600.0/factor
    if npix[0] > maximum:
        pixScale = pixScale * npix[0] / maximum
        npix = [maximum,maximum]

    # write command
    splitName = file.split("-")
    outFile = splitName[0]    
    #outFile = splitName[0] + "-" + splitName[1]
    #if splitName[2][0] != "P":
    #    outFile = outFile + "-" + splitName[2]
    outFile = outFile  + "-SDSS.jpg"
    
    command = 'wget -A.jpg "http://skyservice.pha.jhu.edu/DR8/ImgCutout/getjpeg.aspx?ra={0:5f}%20%20&dec={1:5f}%20%20%20&scale={2:5f}&width={3:}&height={4:}&opt=&query=" -O '\
              .format(coord[0], coord[1], pixScale, npix[0], npix[1]) + outFile
    os.system(command)
    print command
print "Program Finished Successfully"
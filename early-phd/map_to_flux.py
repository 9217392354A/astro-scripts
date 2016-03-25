#Program created by Chris Fuller to test a function for extracting flux's from a fits file using appature photomotry 

#import stuff
from numpy import *
import numpy
import scipy
import math
import sys
import os
from os.path import join as pj

#File stuff
cat = "bigcoma.csv"
catfolder = "/Users/chrisfuller/Dropbox/coma/Catalogues"
catout ="comaTEST.csv" 
folder = "/Users/chrisfuller/Dropbox/coma/flux2/"
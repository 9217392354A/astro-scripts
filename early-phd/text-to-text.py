import sys
import os
from os.path import join as pj
from numpy import *
import numpy as np

folder = "/Users/chrisfuller/Dropbox/coma/Catalogues/Radio/"
catFile = "alfalfa-pre.csv"
cat = loadtxt(pj(folder,catFile), dtype=str,  skiprows=0, delimiter=",", unpack=False)

for row in cat:
    info = row.split(',')
    
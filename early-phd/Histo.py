# -*- coding: utf-8 -*-
"""
Created on Thu May 31 18:20:02 2012

@author: chrisfuller
"""

# Chris Fuller May 2012
from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt 

################user varaibles##################################
folder = "/Users/chrisfuller/Dropbox/HeViCS/"
V250 = loadtxt(pj(folder,"V250.csv"), dtype=float, delimiter=",", unpack=False)
C250 = loadtxt(pj(folder,"C250.csv"), dtype=float, delimiter=",", unpack=False)

plt.hist(250),bins=10,align='mid')
plt.hist(V250),bins=10,align='mid')
plt.show()


    

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 08:36:44 2013

Xsquare = ((ra[i] - centre[0])*math.cos(dec[i] / 180.0 * math.pi)*math.cos(PA) + (dec[i] - centre[1])*math.sin(PA))**2
Ysquare = (-(ra[i] - centre[0])*math.cos(dec[i] / 180.0 * math.pi)*math.sin(PA) + (dec[i] - centre[1])*math.cos(PA))**2
radius[i] = math.sqrt(Xsquare + Ysquare / math.cos(inclin)**2.0)
radius[i] = 2.0 * dist * 1000.0 * math.tan(radius[i]*math.pi/(180.0*2.0))

@author: chrisfuller
"""

import numpy as np
from scipy import pi,sin,cos
import matplotlib.pyplot as plt
import math as m



image = np.zeros(shape=(4147,4117))

a=4
b=3
pa = 240
cdelt = 5.0/3600.0
    
def dEllipse(x,y,a,b,PA):
    the = -PA+90.0
    cos_a,sin_a=cos(the*pi/180.0),sin(the*pi/180.0)
    Xs = (x*cos_a - y*sin_a)**2.0 
    Ys = (x*sin_a + y*cos_a)**2.0 
    r = np.sqrt((Xs/a**2) + (Ys/b**2)) 
    return r
def photo(image,cdelt,a,b,PA):
    c = cdelt*60.0
    #find central pixel and shape
    x0 = np.array(image.shape[0])/2.0
    y0 = np.array(image.shape[1])/2.0
    xMax = x0*2
    yMax = y0*2
    x = np.arange(0-x0, xMax-x0)
    y = np.arange(0-y0, yMax-x0)
    #pixel cordinates 
    xx, yy = np.meshgrid(x, y)
    #creat image of ellipses for appature photometry
    new_image= dEllipse(xx,yy,a,b,PA)*c
    return new_image
 



i = photo(image,cdelt,a,b,pa)
plt.imshow(i)
plt.show()
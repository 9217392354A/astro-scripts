# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 18:02:21 2012

@author: chrisfuller
"""

v1 = min(vel)
v2 = 0.0

while v2 < max(vel):
    v2 = v1 + vBin
    vMid.append(v2 - vBin/2.0)
    temp = []
    for gal in cat:
        if v1 <= gal[8]<= v2:
            temp.append(gal[8])
    if len(temp) == 0:
        temp.append(1.0)
    noGal.append(len(temp))    
    v1 = v2
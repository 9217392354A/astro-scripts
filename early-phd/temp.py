# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 12:09:30 2012

@author: spx6cff
"""

###create bins 
#caculate outer radius r2

r1 = 0.0
r2 = r
i=0
binNum = 0
while binNum < bins:
    if binNum == 0:
        r1 = 0.0   
    else:
        r2 = r1 + r
        
    
    
    temp = []
    for i in rad:
        if r1 <= i <= r2:
            temp.append(i)    
    countsArea.append(len(temp))    
    
    rMidArea = (r1+r2)/2.0 
    rListArea.append(rMid)
    
    print binNum,r1,r2, (r2-r1),len(temp)
    
    
    if binNum ==0:
        r1 = r
    else:
        r1 = r2  
    binNum += 1
    
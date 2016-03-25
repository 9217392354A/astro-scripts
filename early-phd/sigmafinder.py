#Program created by Chris Fuller
#program takes a table and works out its velocity and angular dispersion


# import stuff
from numpy import *
import numpy as np
from scipy.optimize import fmin
import sys
import os
from os.path import join as pj
import matplotlib.pyplot as plt 


#user varaibles
catalogefile = "roughcomasdsscat.csv"
folder = "/home/fossil/spx6cff/coma/tables/calatogues"
r = 0.8
rArea = 0.1
bins = 20
width = 5.0/60.0
distance = 8.0
vBin = 100.0


#Center of X-ray emssion ie center of mass
centerRA = 194.833
centerDEC = 27.903

rList = []
counts =[]
rListArea = []
countsArea =[]
widthR = []
widthA =[]
noGal =[]
vMid = []
#############################################################
#Functions
#results = fmin(graph,p0,args=(x,y,yerror)))
def graph(x,*args):
    xDat = args[0]
    yDat = args[1]
    yError=args[2]

    ci = 0.0
    for i in range(0,len(xDat)):

        step = ((yDat[i]-x[0]**xDat[i]*math.exp(-x[0])/math.factorial(int(xDat[i]))**2)/yError[i]**2)
        ci += step
    return ci

    

#Caculated each value sepratly
#read in catalogue 

cat = loadtxt(pj(folder,catalogefile), dtype=float,  skiprows=1, delimiter=",", unpack=False)

vel = cat[:,8]
rad = cat[:,7]


###create bins and work out how many galaxies fall into that bin###

##########################Fixed Area!###############from numpy import *


r1 = 0.0
r2 = r

binNum = 0
while r2 < distance:
    if binNum == 0:
        r1 = 0.0   
    else:
        r2 = ((r1**2) + (r**2))**0.5
    area = math.pi*((r2**2)-(r1**2))
    temp = []
    for i in rad:
        if r1 <= i <= r2:
            temp.append(i)    
    counts.append((len(temp)/area))    
    rMid = (r1+r2)/2.0 
    rList.append(rMid)
    widthA.append(r2-r1)
    
    print binNum,r1,r2, (r2-r1),len(temp),'area = ', area
    if binNum ==0:
        r1 = r
    else:
        r1 = r2  
    binNum += 1
    

##same but;
###########################Fixed Radius then counts/area######################

r1 = 0.0
r2 = rArea

binNum = 0
while r2 < distance:
    if binNum == 0:
        r1 = 0.0   
    else:
        r2 = r1 + rArea
        
    area = math.pi*((r2**2)-(r1**2))
    
    
    temp = []
    for j in rad:
        if r1 <= j <= r2:
            temp.append(j)
    if len(temp) == 0:
        continue
            
    countsArea.append((len(temp)/area))    
    
    
    
    rMidArea = (r1+r2)/2.0 
    rListArea.append(rMidArea)
    
    print binNum,r1,r2, (r2-r1),len(temp),'area = ', area
    
    
    if binNum ==0:
        r1 = rArea
    else:
        r1 = r2  
    binNum += 1
    
##############Fitting ###############
p0 = [1.0]
#results = fmin(graph,p0,args=(rList,counts,sqrt(counts)))
#results2 = fmin(graph,p0,args=(rListArea,countsArea,sqrt(countsArea)))

#print 'fixed area', results
#print 'fixed radius', results2
#plt.subplot(111)


#rects2 = plt.bar(rListArea, countsArea, widthA, color='blue',align = 'center',log='true')
#plt.subplot(111)
#rects1 = plt.bar(rList, counts, r, color='red',align = 'center',log='true')
#plt.title('Coma cluster angluar dispersion, fixed area vs fixed radius/area ')
#plt.legend((rects1[0], rects2[0]), ('Fixed area    lamda = ',round(results[0],3)), ('Fixed radius / area      lamda = ',round(results2[0],3)))
#plt.show()

scat1 = plt.scatter(rList, (counts), s=20, c='b', marker='o')
plt.subplot(111)
scat2 = plt.scatter(rListArea, (countsArea), s=20, c='r', marker='o')

plt.xlabel('Radius from cluster center, r /deg', fontsize=12)
plt.ylabel('Galaxies/degree^2', fontsize=12)
plt.title('Coma Cluster angular distribution', fontsize=20)
plt.xlim(0, distance) 
plt.ylim( ymin =0 )
plt.legend()
plt.grid(True)
plt.semilogy()
#plot(xRad, yRad), 'go-', label='line 1', linewidth=2)



plt.show()

   
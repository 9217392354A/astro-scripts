#galaxy evolution simulation, Chris Fuller April 2013

import os
import numpy as np  
from os.path import join as pj
import atpy as at
import pylab as pl
import random as ra

####### initial conditions
Ngal = 10 # Number of galaxies
box= 10 # Size of box Mpc
steps = 10 #numer of steps
years = 1E6

# convert box into m
box = box*1E6*1.6E15
####### functions ###########
#function to assign galaxyies initial locations in xyz
def intialCon(N,length):
    x,y,z=[],[],[]
    for i in range(0,N):
        x.append((ra.random()-0.5)*2*length)
        y.append((ra.random()-0.5)*2*length)        
        z.append((ra.random()-0.5)*2*length)
    x = np.array(x,dtype=float)
    y = np.array(y,dtype=float)
    z = np.array(z,dtype=float)
    emptThree = np.zeros((len(x),3), dtype=float)
    locs = np.column_stack((x,y,z))
    return locs ,emptThree  
     
def d(x1,x2):
    return (x2 - x1)
# distance finder
def Cal(locs,u,t):
    G = 6.7E-7 #grav constant
    newgalLoc = locs
    v = u 
    for i in range(0,len(locs)):
        galLoc = locs[i]
        
        mgal1 = 1E10
        fx = []
        fy = []
        fz = []

        for j in range(0,len(locs)):

            if j == i: continue #if same as gal continue
            mgal2 = 1E10
            
            #find distance to galaxy
            dx = d(galLoc[0],locs[j][0])
            dy = d(galLoc[1],locs[j][1])
            dz = d(galLoc[2],locs[j][2])
            
            #find force      
            fx.append(G*mgal1*mgal2/dx**2)
            fy.append(G*mgal1*mgal2/dy**2)
            fz.append(G*mgal1*mgal2/dz**2)
            
        # net force
        fx = np.sum(fx)
        fy = np.sum(fy)
        fz = np.sum(fz)
        
        # accel
        ax = fx/mgal1
        ay = fy/mgal1
        az = fz/mgal1
        
        # vel
        vx = u[0] + ax*t
        vy = u[1] + ay*t
        vz = u[2] + az*t
        
        v[j,0] = vx
        v[j,1] = vy
        v[j,2] = vz
        
        # new coords
        newgalLoc[j,0] = galLoc[0] + vx*t
        newgalLoc[j,1] = galLoc[1] + vy*t
        newgalLoc[j,2] = galLoc[2] + vz*t
    
    return newgalLoc,v
            
            
            
    
            

####### control #############

#first create initial list of positions
locations, velocity  = intialCon(Ngal,box)
 
# find time step

tstep = years/steps
#now run through time steps
for i in range(0,steps):
    #first caculate the velocity of each galaxy and its new position
    
    new_locations,new_velocity = Cal(locations,velocity,tstep)
    print "step: ",i , " veldis:", np.std(new_velocity[2])


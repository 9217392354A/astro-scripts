#flux_verification, Chris Fuller May 2013

import numpy as np  
from os.path import join as pj
import atpy as at
#from astropy.table import Table
import matplotlib.pyplot as plt


###################### inputs ######################
folder = "/Users/chrisfuller/Dropbox/phd/plots/flux_verifacation/"
#folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/flux_comparsion_deep_pacs/"
#fname = "deep_pacs_fluxes_130513.fits"
fname = 'flux_verification.fits'

print "reading in cat"
cat = at.Table(pj(folder,fname),type="fits")

#################### functions ###################
#plotta takes each band arguments and then what you want to plot against it
#then finds out how may values were detected in each band and plots it up
def plotta(band,litfluxs,dX,dY,x,y,tune):
    print 'starting band: ', band
    #extract hefocs flux and create list of litfluxs
    hefocs = np.nan_to_num(cat['F'+str(band)])
    listfluxs = []
    for name in litfluxs: listfluxs.append(np.nan_to_num(cat[str(name)]))
    
    # create plot axes
    f = plt.axes([x,y,dX,dY])
    f.set_title(band+"$\mu m$",size=10)
    f.xaxis.set_visible(False)
    f.set_ylabel("HefoCS Flux (Jy)",size=8)
    for tick in f.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)
    for tick in f.yaxis.get_major_ticks():
        tick.label.set_fontsize(6) 
    f.loglog()
    
    # create risidual axes
    r = plt.axes([x,y-(dY/y)+tune,dX,dY])
    r.set_ylabel("Difference (%)",size=8)
    r.set_xlabel("Literature Flux (Jy)",size=8)
    for tick in r.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)
    for tick in r.yaxis.get_major_ticks():
        tick.label.set_fontsize(6) 
    r.semilogx()
    r.grid(True,axis='y')
    
    #now plot hefocs against all the lit fluxes
    for litflux in listfluxs:
        litflux_c,hefocs_c = cleaner(litflux,hefocs)
        fit_x,fit_y, ris_y, ris_x = fitta(litflux_c,hefocs_c)
        f.plot(fit_x,fit_y,'k')
        f.scatter(litflux_c,hefocs_c,s=6, marker='o', color='k')
        r.scatter(ris_x,ris_y,s=6,marker='o', color='k')
        r.axhline(0.0, color='k')
        #r.autoscale(enable=True, axis='both', tight=True)
        #f.autoscale(enable=True, axis='both', tight=True)
        r.set_ylim( -5.0, 5.0 )
        #r.set_xlim( np.min(litflux_c)*0.8, np.max(litflux_c)*1.2)
        #f.set_xlim( np.min(litflux_c)*0.8, np.max(litflux_c)*1.2)
        #f.set_ylim( np.min(hefocs_c)*0.8, np.max(hefocs_c)*1.2)
                
        
        #r.set_yticks([-0.2,-0.1,0.0,0.1,0.2])
        
#fitting rotuine and caculates risiduals        
def fitta(x_array,y_array):
    
    # Find and plot 1st order line of best fit 
    coeff = np.polyfit( x_array, y_array, 1 ) 
    p = np.poly1d( coeff ) 
    x = np.logspace(np.log10(np.min(x_array)), np.log10(np.max(x_array)))
    y_mod = p(x)
    
    #find risidual
    r = (y_array - p(x_array))/y_array
    return x,y_mod,r,x_array
        
def cleaner(x_u,y_u):
    x_c, y_c= [],[]
    for i in range(0,len(x_u)):
        #check if both are real
        if x_u[i] == 0.0 or y_u[i] == 0.0: continue
        else:
            x_c.append(x_u[i])
            y_c.append(y_u[i])

    return x_c,y_c

    
def plotta_sub(band,litfluxs,i):
    #print 'starting band: ', band
    #extract hefocs flux and create list of litfluxs
    hefocs = np.nan_to_num(cat['F'+str(band)])
    listfluxs = []
    for name in litfluxs: listfluxs.append(np.nan_to_num(cat[str(name)]))
 
    
    #create subplots
    f = plt.subplot(yy,xx,ff[i]) #this is the subplot that holds both figures    
    r = plt.subplot(yy,xx,rr[i])
    
    # create plot axes
    #f.set_title(band+"$\mu m$",size=10)
    f.text(0.21,8.5 ,band+"$\mu m$",size=10 )
    for tick in f.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in f.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)     
    f.set_xlim(0.19,15)
    f.set_ylim(0.19,15)
    r.set_xlim(0.19,15)
    f.loglog()
    # create risidual axes
    for tick in r.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in r.yaxis.get_major_ticks():
        tick.label.set_fontsize(10) 
    r.semilogx()
    r.grid(True,axis='y')
    #f.yaxis.set_visible(False)
    
    #now plot hefocs against all the lit fluxes
    for litflux in listfluxs:
        litflux_c,hefocs_c = cleaner(litflux,hefocs)
        fit_x,fit_y, ris_y, ris_x = fitta(litflux_c,hefocs_c)
        f.plot(fit_x,fit_y,'k')
        f.scatter(litflux_c,hefocs_c,s=6, marker='o', color='k')
        r.scatter(ris_x,ris_y,s=6,marker='o', color='k')
        r.axhline(0.0, color='k')        
        r.set_ylim( -0.5, 0.5 )
        
        #errors 
        dm,dc = errorInLine(litflux_c,hefocs_c)
        
        print 'band: ', band, 'banderrors: dm=', np.round(dm,decimals=3), ', dc=', np.round(dc, decimals=3)         
        

        #r.autoscale(enable=True, axis='both', tight=True)
        #f.autoscale(enable=True, axis='both', tight=True)
        #r.set_xlim( np.min(litflux_c)*0.8, np.max(litflux_c)*1.2)
        #f.set_xlim( np.min(litflux_c)*0.8, np.max(litflux_c)*1.2)
        #f.set_ylim( np.min(hefocs_c)*0.8, np.max(hefocs_c)*1.2)
                
        
        #r.set_yticks([-0.2,-0.1,0.0,0.1,0.2])
        #bands = ["100","160","250","350","500"]

    if i == 2: #label 
        r.set_ylabel("Difference (%)",size=10)
        r.set_xlabel("Literature Flux (Jy)",size=10)
        f.set_ylabel("HefoCS Flux (Jy)",size=10)
        
    #if i == 4:
    #    r.set_xlabel("Literature Flux (Jy)",size=10)
    
    if i == 0 or i == 1 or i == 3:
        for tick in r.xaxis.get_major_ticks(): tick.label.set_visible(False)
        for tick in f.xaxis.get_major_ticks(): tick.label.set_visible(False)
        
    if i == 3 or i == 4:
        for tick in r.yaxis.get_major_ticks(): tick.label.set_visible(False)         
        for tick in f.yaxis.get_major_ticks(): tick.label.set_visible(False)   
    #if i == 3 or i == 4: #don't label
    #    f.yaxis.set_visible(False)
    #    r.yaxis.set_visible(False)
    #print '--max-- lit:', np.max(litflux_c)," hefocs: ",np.max(hefocs_c)
    #print '--min-- lit:', np.min(litflux_c)," hefocs: ",np.min(hefocs_c)



#function for finding errors in a straight line of two arras x and y
# assuming least squares and equal error on each point
def errorInLine(x_data,y_data):
    x_data = np.array(x_data, dtype=np.float)
    y_data = np.array(y_data, dtype=np.float)
    
    n = len(x_data)
    D = np.sum(x_data**2) - 1./n * np.sum(x_data)**2
    x_bar = np.mean(x_data)
    p_coeff, residuals, _, _, _ = np.polyfit(x_data, y_data, 1, full=True)
    dm_squared = 1./(n-2)*residuals/D
    dc_squared = 1./(n-2)*(D/n + x_bar**2)*residuals/D
    
    return np.sqrt(dm_squared[0]),np.sqrt(dc_squared[0])
    
    
#################### control #####################
#create figure
fig = plt.figure(figsize = (4.,9.6),facecolor='w',edgecolor='w')

#width and height of subplots
dx = 0.35
dy = 0.10            
#cord for subplots
xL = 0.1
xR = 0.6

tuneA = 0.0230
tuneB = 0.097
tuneC = 0.398

y1 = 0.8
y2 = 0.5
y3 = 0.2

#BGS
#plotta("100",["BGS_F100"],dx,dy,xL,y1,tuneA)
#plotta("160",["BGS_F160"],dx,dy,xR,y1,tuneA)
#plotta("250",["BGS_F250"],dx,dy,xL,y2,tuneB)
#plotta("350",["BGS_F350"],dx,dy,xR,y2,tuneB)
#plotta("500",["BGS_F500"],dx,dy,xL,y3,tuneC)

xx = 2
yy = 6 
ff = [1,5,9,2,6]
rr = [3,7,11,4,8]

bands = ["100","160","250","350","500"]
xxvals = [1,3,5,7,9]
for i in range(0,len(bands)): 
    plotta_sub(bands[i], ["BGS_F"+bands[i]], i) 

#plotta_sub("100",["BGS_F100"],1)
#plotta_sub("160",["BGS_F160"],2)
#plotta_sub("250",["BGS_F250"],3)
#plotta_sub("350",["BGS_F350"],4)
#plotta_sub("500",["BGS_F500"],5)

#plotta("100",["DEEP100"],0.8,0.2,0.1,0.75,0.067)
#plotta("160",["DEEP160"],0.8,0.2,0.1,0.25,0.6)


plt.subplots_adjust(left=0.14, bottom=0.06, right=0.98, top=0.992, wspace=0.0, hspace=0.0)
#plt.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/fornax',"flux_verfication.pdf"))
plt.savefig('/Users/chrisfuller/Desktop/flux_verfication.pdf')
plt.show()
# fig gen for paper, create figure that is for a sigle ps output drop all of that fits files and the cat
#into the folder and let it work its magic. This dose'nt inlcude that bg fit 
# Chris Fuller May 2013
import os
import numpy as np  
from os.path import join as pj
import aplpy as ap
from atpy import Table
import matplotlib.pyplot as plt
###################### inputs ######################
folder = "/Users/chrisfuller/Dropbox/phd/plots/fig_gen/"

print "reading in cat"
cat = Table(pj(folder,"cat.text"),type="ascii")
os.chdir(folder)

#################### functions ###################

#function to plot FIR contours
def plottabw(name,fname,texts,x,y):
    print name
    #find optical parameters
    i = np.where(cat.OBJECT == name)[0][0]
    a,b,pa = cat.FULLMAJAX[i]* 1.1,cat.FULLMINAX[i], cat.PA[i]
    rDust = cat.R250[i] 
    ra,dec = cat.GRA2000[i], cat.GDEC2000[i]
    beam = 18.2/3600.0
    #make plot
    f = ap.FITSFigure(fname,figure=fig,subplot=[x,y,dx,dy])
    f.show_contour(str(name)+"-PSWmap-ellipses.fits", levels=[rDust], smooth=1,colors='black')
    f.recenter(ra,dec,beam*7.0)
    f.show_colorscale(cmap='gray_r')
    f.add_scalebar(1.0/60.0)
    f.scalebar.set_label('1 arcmin')
    f.scalebar.set_color('black')
    #f.show_grid()
    f.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    #f.show_markers(ra, dec,c='r',marker='x')
    f.show_ellipses(ra, dec, a/60.0, b/60.0, angle=pa-90, edgecolor='red', lw=1)
    f.axis_labels.set_font(size=8)
    f.tick_labels.set_font(size=6)
    
    # add text 
    
    plt.figtext(x+tune, y+0.23,texts, size=12, weight="semibold", color='black')
    
    #function to plot FIR contours
def plotta(name,fname,texts,x,y):
    print name
    #find optical parameters
    i = np.where(cat.OBJECT == name)[0][0]
    a,b,pa = cat.FULLMAJAX[i]* 1.1,cat.FULLMINAX[i], cat.PA[i]
    ra,dec = cat.GRA2000[i], cat.GDEC2000[i]
    beam = 18.2/3600.0
    rDust = cat.R250[i] 
    #make plot
    f = ap.FITSFigure(fname,figure=fig,subplot=[x,y,dx,dy])
    f.show_contour(str(name)+"-PSWmap-ellipses.fits",levels=[rDust], smooth=1, colors='black')
    f.recenter(ra,dec,beam*7.0)
    f.show_colorscale()
    #f.show_contour(str(name)+"-PSWmap-PSF-convolvedmap.fits",levels=15,colour='green', alpha=0.5) 
    f.add_scalebar(1.0/60.0)
    f.scalebar.set_label('1 arcmin')
    f.scalebar.set_color('black')
    #f.show_grid()
    f.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    #f.show_markers(ra, dec,c='r',marker='x')
    f.show_ellipses(ra, dec, a/60.0, b/60.0, angle=pa-90, edgecolor='red', lw=1)
    f.add_beam(major=beam, minor=beam,angle=0.0)
    f.show_beam(major=beam, minor=beam,angle=0.0,fc='white')
    f.axis_labels.set_font(size=8)
    f.tick_labels.set_font(size=6)
    
    # add text 
    plt.figtext(x+tune, y+0.23,texts, size=12, weight="semibold", color='black')
    
def makeCumlative(x):
	temp = 0
	result = []
	for val in x:
		temp = temp + val
		result.append(temp)
	return result
    
#################### control ######################
#create figure
fig = plt.figure(figsize = (8.0,8.5),facecolor='w',edgecolor='w')

#width and height of subplots
dx = 0.30
dy = 0.25            
#cord for subplots
xL = 0.15
tune= 0.0025
xR = 0.60

y1 = 0.7
y2 = 0.38
y3 = 0.05


name = 312
i = np.where(cat.OBJECT == name)[0][0]
Rdust = cat.R250[i]
a,b,pa = cat.FULLMAJAX[i],cat.FULLMINAX[i], cat.PA[i]


#read in graphs
apNoise = np.loadtxt("312-PSWmap-aperfit.csv")
fluxes  = np.loadtxt("312-PSWmap-profiles.csv",skiprows=1,delimiter=",")

#split fluxes into various intresting componants to plot
opRadius = fluxes[:,6]
surBright = fluxes[:,4]
radSN =np.abs(fluxes[:,0]/fluxes[:,1])
fluxCum = fluxes[:,2]

#fits images
plottabw(312,"312.fits","(a)",xL,y1)
plotta(312,"312-PSWmap-rawmap.fits","(b)",xR,y1)


#ajust a
a = a * 0.8


# surface brightness
f = plt.axes([xL,y2,dx,dy])
f.set_title("Surface brightness profile",size=10)
f.set_xlabel("Radius (arcmin)",size=8)
f.set_ylabel("Intensity (Jy/beam)",size=8)
for tick in f.xaxis.get_major_ticks():
    tick.label.set_fontsize(6)
for tick in f.yaxis.get_major_ticks():
    tick.label.set_fontsize(6) 
f.semilogy()
f.plot(opRadius,surBright,'k-')
f.axvline(x=Rdust, color='k', ls='--')
f.axvline(x=a, color='r', ls='--')
plt.figtext(xL+tune, y2+0.23,"(c)", size=12, weight="semibold")




######## noise app #########
######## middle lower right ################  
f = plt.axes([xR,y2,dx,dy])
f.set_title("Noise-aperture dependency",size=10)
f.set_xlabel("Radius (arcmin)",size=8)
f.set_ylabel("Noise (Jy)",size=8)
for tick in f.xaxis.get_major_ticks():
    tick.label.set_fontsize(6)
for tick in f.yaxis.get_major_ticks():
    tick.label.set_fontsize(6) 
f.plot(opRadius,apNoise[:len(opRadius)],'k-')
f.axvline(x=Rdust, color='k', ls='--')
f.axvline(x=a, color='r', ls='--')
plt.figtext(xR+tune, y2+0.23,"(d)", size=12, weight="semibold")

######## SN profile #########
######## bottom left ################  
f = plt.axes([xL,y3,dx,dy])
f.set_title("Intensity S/N profile",size=10)
f.set_xlabel("Radius (arcmin)",size=8)
f.set_ylabel("S/N",size=8)
for tick in f.xaxis.get_major_ticks():
    tick.label.set_fontsize(6)
for tick in f.yaxis.get_major_ticks():
    tick.label.set_fontsize(6)
f.plot(opRadius,radSN,'k-')
f.axvline(x=Rdust, color='k', ls='--')
f.axvline(x=a, color='r', ls='--')
plt.figtext(xL+tune, y3+0.23,"(e)", size=12, weight="semibold")

######## intensity cumlative #########
######## bottom right ################  
f = plt.axes([xR,y3,dx,dy])
f.set_title("Cumulative intensity profile",size=10)
f.set_xlabel("Radius (arcmin)",size=8)
f.set_ylabel("Intensity (Jy)",size=8)
for tick in f.xaxis.get_major_ticks():
    tick.label.set_fontsize(6)
for tick in f.yaxis.get_major_ticks():
    tick.label.set_fontsize(6)
f.plot(opRadius,makeCumlative(surBright),'k-')
f.axvline(x=Rdust, color='k', ls='--')
f.axvline(x=a, color='r', ls='--')
plt.figtext(xR+tune, y3+0.23,"(f)", size=12, weight="semibold")


plt.savefig("312.eps")
plt.show()

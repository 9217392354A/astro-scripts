# fig gen for paper, create figure that is for a sigle ps output drop all of that fits files and the cat
#into the folder and let it work its magic. This dose'nt inlcude that bg fit 
# Chris Fuller March 2014
import os
import numpy as np  
from os.path import join as pj
import aplpy as ap
from atpy import Table
import matplotlib.pyplot as plt
###################### inputs ######################
gal_name = 'CCC2099'

folder_fig = '/Users/chrisfuller/Dropbox/phd/herchel/coma/fig_gen/' + gal_name + '/'
folder_cat = '/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/'

print "reading in cat"
cat = Table(pj(folder_cat,"coma_supercluster_cal12.fits"))

#select galaxy
galaxy = cat.where(cat.OBJECT == gal_name)

#change dir to folder directory
os.chdir(folder_fig)

#################### functions ###################

#function to plot FIR contours
def plottabw(name,fname,texts,x,y):
    print name
    #find optical parameters
    a,b,pa = galaxy.FULLMAJAX[0],galaxy.FULLMINAX[0], galaxy.PA[0]
    #pa = pa + 90.
    rDust = galaxy.R250[0]*1.3 
    ra,dec = galaxy.GRA2000[0], galaxy.GDEC2000[0]
    beam = 18.2/3600.0
    #make plot
    f = ap.FITSFigure(fname,figure=fig,subplot=[x,y,dx,dy])
    #f.show_contour(str(name)+"-PSWmap-ellipses.fits", levels=[rDust], smooth=1,colors='black')
    #f.show_grid()
    f.recenter(ra,dec,beam*4.)
    f.show_grayscale(vmin=-9.709e-02, vmax=2.549e-01)
    f.add_scalebar(0.5/60.0)
    f.scalebar.set_label('30 arcsec')
    f.scalebar.set_color('black')
    #f.show_grid()
    f.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    #f.show_markers(ra, dec,c='r',marker='x')
    f.show_ellipses(ra, dec, a*1.15/60.0, a*1.15/60.0, angle=pa, edgecolor='red', lw=1)
    f.show_ellipses(ra, dec, rDust/60.0, rDust/60.0, angle=pa, edgecolor='black', lw=1)
    f.axis_labels.set_font(size=8)
    f.tick_labels.set_font(size=6)
    #f.set_theme('publication')
    # add text 
    
    plt.figtext(x+tune, y+0.23,texts, size=12, weight="semibold", color='black')
    
    #function to plot FIR contours
def plotta(name,fname,texts,x,y):
    print name
    #find optical parameters
    a,b,pa = galaxy.FULLMAJAX[0],galaxy.FULLMINAX[0], galaxy.PA[0]
    ra,dec = galaxy.GRA2000[0], galaxy.GDEC2000[0]
    #pa = pa + 90.
    beam = 18.2/3600.0
    rDust = galaxy.R250[0]*1.4
    #make plot
    f = ap.FITSFigure(fname,figure=fig,subplot=[x,y,dx,dy])
    #f.show_contour(str(name)+"-PSWmap-ellipses.fits",levels=[rDust], smooth=1, colors='black')
    f.recenter(ra,dec,beam*4.)
    f.show_colorscale(vmin=-6.02e-02, vmax=20.664e-02, cmap='gist_earth')
    #f.show_contour(fname,levels=150, filled=True, cmap='cool') 
    f.add_scalebar(0.5/60.0)
    f.scalebar.set_label('30 arcsec')
    f.scalebar.set_color('black')
    #f.show_grid()
    f.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
    #f.show_markers(ra, dec,c='r',marker='x')
    f.show_ellipses(ra, dec, a*1.15/60.0, a*1.15/60.0, angle=pa, edgecolor='red', lw=1)
    f.show_ellipses(ra, dec, rDust/60.0, rDust/60.0, angle=pa, edgecolor='black', lw=1)
    f.add_beam(major=beam, minor=beam,angle=0.0)
    f.show_beam(major=beam, minor=beam,angle=0.0,fc='white')
    f.axis_labels.set_font(size=8)
    f.tick_labels.set_font(size=6)
    #f.set_theme('publication')
    
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


Rdust = galaxy.R250[0]
a,b,pa = galaxy.FULLMAJAX[0],galaxy.FULLMINAX[0], galaxy.PA[0]


#read in graphs
apNoise = np.loadtxt(gal_name + "-PSWmap-aperfit.csv")
fluxes  = np.loadtxt(gal_name + "-PSWmap-profiles.csv",skiprows=1,delimiter=",")

#split fluxes into various intresting componants to plot
opRadius = fluxes[:,6]
surBright = fluxes[:,4]
radSN =np.abs(fluxes[:,0]/fluxes[:,1])
fluxCum = fluxes[:,2]

#fits images
plottabw(gal_name,gal_name+".fits","(a)",xL,y1)
plotta(gal_name, gal_name + "-PSWmap-rawmap.fits","(b)",xR,y1)


#ajust a
#a = a * 0.8


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


plt.savefig(gal_name + ".pdf")
plt.show()

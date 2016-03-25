#program to make figure for paper using ngp maps, by Chris Fuller, Jan 2013

import aplpy
import numpy as np
import aplpy as ap
import atpy as at
from os.path import join as pj
import matplotlib.pyplot as plt

#folder = "/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/"
#cat = at.Table(pj(folder,"HeFoCS-fluxes-260313-mybgsub-v2-handmeasure-bad-detections-removed-morphology.fits"),type='fits')

folder = "/Users/chrisfuller/Desktop/working"
coma = at.Table(pj(folder,"ccc.fits"),type='fits')
fornax = at.Table(pj(folder,"fcc.fits"),type='fits')
virgo = at.Table(pj(folder,"vcc.fits"),type='fits')
coma_cluster = at.Table(pj(folder,"coma-cluster.fits"),type='fits')

fname_coma = pj(folder,"ngp-plw.fits")
fname_fornax = pj(folder,"fornax-plw.fits")
fname_virgo = pj(folder,"virgo-plw.fits")
dx = 0.35
dy = 0.35            
            #cord for subplots
xL = 0.10
xR = 0.60
            
y1 = 0.5
y2 = 0.1

#create figure
fig = plt.figure(figsize = (10,10),facecolor='w',edgecolor='w')

#coma
c = ap.FITSFigure(fname_coma,figure=fig,subplot=[xL,y1,dx,dy],downsample=5)
c.show_grayscale()
c.show_circles(194.87844322,27.88420422,1.67)
c.show_markers(coma.RA, coma.DEC,c='g',marker='x')
c.show_markers(coma_cluster.RA, coma_cluster.DEC,c='r',marker='x')
c.add_scalebar(0.56*1.0)
c.scalebar.set_label('1 Mpc')
#c.show_circles(comax,comay,1.67)
c.axis_labels.set_font(size=8)
c.tick_labels.set_font(size=6)

#Virgo
v = ap.FITSFigure(fname_virgo,figure=fig,subplot=[xR,y1,dx,dy],downsample=5)
v.show_grayscale()
v.recenter(185.7,9.8,3.5*3.5)
v.show_markers(virgo.RA, virgo.DEC,c='r',marker='x')
v.add_scalebar(3.5*1.0)
v.scalebar.set_label('1 Mpc') #1.7
#v.show_circles(194.87844322,27.88420422,1.7*3.5)
v.axis_labels.set_font(size=8)
v.tick_labels.set_font(size=6)

#Fornax
f = ap.FITSFigure(fname_fornax,figure=fig,subplot=[xL,y2,dx,dy],downsample=5)
f.show_grayscale()
f.recenter(54.6289,-35.4545,3.5*3.0)
f.show_markers(fornax.RA, fornax.DEC,c='r',marker='x')
f.add_scalebar(3.0*1.0)
f.scalebar.set_label('1 Mpc')
f.axis_labels.set_font(size=8)
f.tick_labels.set_font(size=6)

#coma cluster
co = ap.FITSFigure(fname_coma,figure=fig,subplot=[xR,y2,dx,dy],downsample=5)
co.show_grayscale()
co.recenter(194.87844322,27.88420422,0.56*3.5)
co.show_markers(coma_cluster.RA, coma_cluster.DEC,c='r',marker='x')
co.add_scalebar(0.56*1.0)
co.show_circles(194.87844322,27.88420422,1.67)
co.scalebar.set_label('1 Mpc')
co.axis_labels.set_font(size=8)
co.tick_labels.set_font(size=6)

fig.savefig(pj(folder,"figure-all-maps.eps"))










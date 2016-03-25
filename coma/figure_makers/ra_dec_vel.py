#make ra dec plots
# Chris Fuller, March 2014



#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D

#Inputs
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'
c1 = 'big_ngp.fits' #input name

#read in cat
cat = Table(pj(folder2,c1))

cat = cat.where((cat.velocity < 11000.0) & (cat.velocity > 3000.0))
cat = cat.where((cat.dec > 15.) & (cat.dec < 35.))
cat = cat.where((cat.ra < 212.) & (cat.ra > 170.))

early = cat.where(cat.elliptical == 1)
late =  cat.where(cat.spiral == 1)
inter=  cat.where(cat.uncertain == 1)

tps = [inter, late, early]
colours = ['g', 'b', 'r']



fig1 = plt.figure(figsize = (10.5,6.5), facecolor='w',edgecolor='w')
for i in range(len(tps)):
	t = tps[i]
	sub = fig1.add_subplot(111)
	sub.scatter(t.ra, t.dec, marker='o', s=25, color = colours[i], alpha=0.6)
	sub.autoscale_view(tight=True, scalex=True, scaley=True)

sub.invert_xaxis()

sub.set_xlabel('RA (J2000)')
sub.set_ylabel('Dec (J2000)')

plt.subplots_adjust(left=0.1, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
fig1.savefig('/Users/chrisfuller/Desktop/ra_vs_dec.pdf')



fig2 = plt.figure(figsize = (10.5,6.5), facecolor='w',edgecolor='w')
for i in range(len(tps)):
	t = tps[i]
	sub = fig2.add_subplot(111)
	sub.scatter(t.ra, t.velocity, marker='o', s=25, color = colours[i], alpha=0.6)
	sub.autoscale_view(tight=True, scalex=True, scaley=True)

sub.invert_xaxis()

sub.set_xlabel('RA (J2000)')
sub.set_ylabel('Velocity ($km s^{-1}$)')

plt.subplots_adjust(left=0.1, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
fig2.savefig('/Users/chrisfuller/Desktop/ra_vs_vel.pdf')






fig3 = plt.figure(figsize = (10.5,6.5), facecolor='w',edgecolor='w')
sub = fig3.add_subplot(111, projection='3d')
for i in range(len(tps)):
	t = tps[i]
	sub.scatter(t.ra, t.dec, t.velocity, marker='o', s=5, color = colours[i], alpha=0.99)
	sub.autoscale_view(tight=True, scalex=True, scaley=True)

sub.invert_xaxis()

sub.set_xlabel('RA (J2000)')
sub.set_ylabel('Dec (J2000)')
sub.set_zlabel('Velocity ($km s^{-1}$)')

plt.subplots_adjust(left=0.1, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
fig3.savefig('/Users/chrisfuller/Desktop/the_great_wall_3d.pdf')


plt.show()
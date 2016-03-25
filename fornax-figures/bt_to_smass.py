#program to genorate stellar masses from optical colours and 2mass magnitudes. For
# galaxies that only have a bt mag then we will fit a function to the 
# other galaxies giving a fuction for btmag to stellar mass.

# Chris Fuller August  2013

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/plots/upper_limits/"
cat = Table(pj(folder,"fornax.fits"))

log_mass = cat['MSTAR'] 
bt = cat['BTmag'] 

# # # # # # # # # # # # # # # # # # # # Main Program # # # # # # # # # # # # # # # # # # # # # # # 
#create fig
fig= plt.figure(figsize = (8.,8.),facecolor='w',edgecolor='w')

#sub fig
p1 = plt.subplot(111)
p1.plot(bt, log_mass,'+b')


p1.set_xlabel('$Log_{10}(m_{bt})$')
p1.set_ylabel('$Log_{10}(M_{star}/M_{\odot})$')
#p1.set_ylim(4,9)
#p1.set_xlim(-2.5,2.5)   
# old fit  [ 0.75986014  6.60244077]
# new fit  [ 0.78038038  6.47259655]
# new fit  [ 0.78937274  6.48619262]
#produce_fit(p2,virgo.F250,virgo.TEMP,'blue')
#produce_fit(p2,fornax.F250,fornax.T,'red')
#line_only(p2, virgo.F250,virgo.TEMP,fornax.F250,fornax.T)







#plt.subplots_adjust(left=0.125,bottom=0.075, top = 0.985, hspace=0.0,wspace=0.0)
#figM.savefig(pj('/Users/chrisfuller/Dropbox/phd/papers/fornax','smass-detected.pdf'))
#plt.show()


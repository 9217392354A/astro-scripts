# program to make detection rates of an input catalogue and put them into a figure for a 
# paper
# Chris Fuller, May 2013

import numpy as np  
from os.path import join as pj
import atpy as at
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib.ticker import MaxNLocator


############### inputs #######################
folder = "/Users/chrisfuller/Dropbox/phd/plots/det_fig"
cat = at.Table(pj(folder,"morphology-added_cat.fits"))

opt ="BTmag" #name of optical magnitude column
mor ="MORPH"
#bands = ["F500","F350","F250","F160","F100"]
bands =["F100","F160","F250","F350","F500"]
#mLables=['','1','','2','','3','','4']
mList=['E/S0','Sa/Sb/Sc','/irr/late/dwarf']
############# functions ######################
#create a col of optical mags for the detected galaxies
def rtn_det(f,op):
    f = np.nan_to_num(f)
    return op[np.where(f!=0.0)]
    
#function creats the bins in the optimal spacing
def binman(x):
    print "binman:- there's ya bins mate."
    return np.arange(np.floor(np.min(x)),np.ceil(np.max(x)),1.0)

#find the morphologies and how many are in each bin
def w(x):
    return len(np.where(cat[mor]==x)[0])
    
def dw(m,f):
    return len(np.where((cat[mor]==m) & (f>0.0))[0])
    
def l10(x):
    return np.log10(x)
############ control ########################

#create bins
optical = cat[opt] #optical data 
#bins = binman(optical.copy())
bins = np.arange(10,22)

#create fig to hold plots
fig = plt.figure(figsize = (8.,8.),facecolor='w',edgecolor='w')

a = np.arange(0.5,1,0.5/len(bands))
#loop through all bands
for i in range(1,len(bands)+1):
    flux = cat[bands[i-1]] #fir data for that band
    opt_det = rtn_det(flux,optical)
    
    ################ optical ###########################
    #hist(optical,bins)
    f1 = plt.subplot(len(bands),1,i)
    
    
    #f1.hist(optical,bins,normed=False, color="black", histtype='bar', alpha=1.0,log=False)
    #f1.hist(opt_det,bins,normed=False,histtype='bar',color="aqua",log=False,alpha=a[i-1])
    
    o_y,o_x = np.histogram(optical,bins=bins)
    od_y,od_x = np.histogram(opt_det,bins=bins) 
    
    f1.bar(bins[:-1], l10(o_y)+1, color="black",  alpha=1.0, width=1)
    f1.bar(bins[:-1], l10(od_y)+1,color="aqua",alpha=a[i-1], width=1)
    f1.tick_params(axis='x', labelbottom='off')
    #remove x lables
    f1.set_ylim(0,3)
    f1.yaxis.set_major_locator(MaxNLocator(3))
    f1.set_yticklabels(['0','1','10',''])     
    f1.set_ylabel("log(N)")
    #f1.grid()
    #f1.set_yticks([0,10,100])    
    
    #add text
    f1.text(10.1,30, bands[i-1][-3:]+"$\mu m$", fontsize=15)
    print bands[i-1],opt_det
    ############### Morphologica detection rates ########
    #f2 = plt.subplot(len(bands),2,2*i)
    #totals = [w("1"),w("2"),w("3"),w("4")]
    #det = [dw("1",flux),dw("2",flux),dw("3",flux),dw("4",flux)]
    
    #f2.bar(np.arange(1,5),totals,align = 'center',color ="black",width=1.0,log=True, alpha = 1.0) #totals
    #f2.bar(np.arange(1,5),det,   align = 'center',color ="purple",width=1.0,log=True, alpha = 1.0) #det
    #f2.tick_params(axis='x', labelbottom='off')
    #f2.set_ylim(1,380)   
    #f2.set_ylabel("N")   
    
f1.set_xlabel("$m_{bt}$")
f1.tick_params(axis='both',labelleft='on', labelbottom='on')

#f2.set_xticklabels(mLables) 
#f2.set_xlabel("Morphological Type")
#f2.tick_params(axis='both',labelleft='on', labelbottom='on')
#plt.subplots_adjust(left=0.14, bottom=0.06, right=0.98, top=0.992, wspace=0.0, hspace=0.0)


#fig.savefig(pj(folder,"detection_rates_optical.pdf"))
plt.subplots_adjust(bottom=0.07, left=0.08, right=0.98, top=0.98, wspace=0.0, hspace=0.0)

fig.savefig(pj('/Users/chrisfuller/Dropbox/','detection_rates_optical.pdf'))

plt.show()


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
cat = at.Table("/Users/chrisfuller/Documents/phd/herchel/coma/final_outputs/coma_supercluster_cal12_pacscorrected.fits")

opt ="r" #name of optical magnitude column
#bands = ["F500","F350","F250","F160","F100"]
bands =["F100","F160","F250","F350","F500"]

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
    x = np.array(x, dtype = np.float)
    #temp = np.log10(x)
    temp = x
    #temp += 1.
    temp[np.where(temp== float('-inf'))[0]] = 0.0
    return temp
############ control ########################

#create bins
optical = cat[opt] #optical data 
#bins = binman(optical.copy())
wid = 0.5
bins = np.arange(11.5,19.,wid)

#create fig to hold plots
fig = plt.figure(figsize = (4.,6),facecolor='w',edgecolor='w')

a = np.arange(0.5,1,0.5/len(bands))
#loop through all bands
for i in range(1,len(bands)+1):
    flux = cat[bands[i-1]] #fir data for that band
    opt_det = rtn_det(flux,optical)
    
    ################ optical ###########################
    #hist(optical,bins)
    f1 = plt.subplot(len(bands),1,i)
    
    
    #f1.hist(optical,bins,normed=False, color="black", histtype='bar', alpha=1.0,log=True)
    #f1.hist(opt_det,bins,normed=False,histtype='bar',color="aqua",log=True,alpha=a[i-1])
    
    o_y,o_x = np.histogram(optical,bins=bins)
    od_y,od_x = np.histogram(opt_det,bins=bins)
    
    f1.bar(bins[:-1], l10(o_y ), width = wid, color="black",  alpha=1.0)
    f1.bar(bins[:-1], l10(od_y), width = wid, color="aqua", alpha=a[i-1])


    f1.tick_params(axis='x', labelbottom='off')
    #remove x lables
    f1.set_ylim(0,380)   
    f1.set_ylabel("N")
    #f1.grid()

    f1.yaxis.set_major_locator(MaxNLocator(5))
    #f1.set_yticklabels(['0','1','10','100','1000'])    
    
    #add text
    #f1.text(10.1,np.log10(1000), , fontsize=15)
    f1.text(0.02, 0.9, bands[i-1][-3:]+"$\mu m$", transform=f1.transAxes, fontsize=10, verticalalignment='top')
    print bands[i-1],' Detected: ', od_y
    print bands[i-1],' Total   : ', o_y

f1.set_xlabel("$m_{r}$")
f1.tick_params(axis='both',labelleft='on', labelbottom='on')

plt.subplots_adjust(left=0.16, bottom=0.09, right=0.97, top=0.98, wspace=0.0, hspace=0.0)
fig.savefig(pj("/Users/chrisfuller/Desktop/","detection_rates_optical_coma.pdf"))
plt.show()



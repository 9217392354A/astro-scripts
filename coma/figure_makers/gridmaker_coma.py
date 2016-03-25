#morph converter
#Chris Fuller Feb 2014
#This program takes surface density and then looks at the fraction detected for a given morphological type
#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
import matplotlib.pylab  as pl
from matplotlib.colors import Colormap

#Inputs
cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12_pacscorrected.fits")


###################### functions ######################
def w(cat,x):
    return where(cat['goldmine']==x)[0]
    
def dw(f):
    f = nan_to_num(f)
    return where((f>0.0))
    

#p(total_cluster[i],detected_cluster[i])
def p(a,b):
    per = np.round(float(b)*100.0/a, decimals = 0)
    error = np.sqrt(b) / a
    return str(np.int(per))+'\\,$\\pm$\\,'+str(np.int(np.round(error*100.0, decimals=0)))
    
def er(a,b):
    return np.sqrt(a) *100.0 / b



def lstr(x):
	return str(np.around(np.log10(x)))

def double_plot(subs, bin_limit, A, col, xlims, ylims, xlabs, ylabs,cn, lab):
	#create bins
	start = bin_limit[0]
	stop = bin_limit[1]
	step = bin_limit[2]
	print bin_limit 
	b = np.arange(start,stop, step) #bin edeges
	print b
	#galaxyies detected
	A_detected  = A[col][where(np.nan_to_num(A['F250']) > 0.0)[0]]
	
	#totals
	A_totals    = A[col]
	
	#caculate fraction detected
	detected_A_hist, bin_edges=  histogram(A_detected, bins = b)
	total_A_hist,bin_edges    =  histogram(A_totals, bins = b)
	per_A = array(detected_A_hist, dtype=float)/ total_A_hist
	
	
	print total_A_hist, b

	################## plotting ##########################

	sub1 = subs[0]
	sub2 = subs[1]


	#plot fraction detected
	sub1.plot(bin_edges[:-1]+step*0.5, per_A*100.0 , color='b', alpha=cn)
	print per_A
	#plot fraction of sample
	sub2.plot(bin_edges[:-1]+step*0.5, total_A_hist *100.0 / len(A_totals) , color='b', alpha=cn)

	#set unfied limits
	sub1.set_ylim(ylims[0],ylims[1])
	sub1.set_xlim(xlims[0],xlims[1])
	#set unfied limits
	sub2.set_ylim(ylims[0],ylims[1])
	sub2.set_xlim(xlims[0],xlims[1])


	#hide x?
	if xlabs == False:
		sub1.tick_params(axis='x', labelbottom='off')

	#hide y?
	if ylabs == False:
		sub1.tick_params(axis='y', labelleft='off')


	# xlabes?
	if xlabs == True:
		sub2.set_xlabel('p(S)')

	# ylables?
	if ylabs == True:
		sub1.set_ylabel('Detected (%)')
		sub2.set_ylabel('Total (%)')




############################ end functions ################################

#####loop over lots of different densityies and overplot########
sigma5 = cat.SIGMA5 

#shaped of gray
no_bins = 10
colours = np.linspace(0.3, 1., no_bins) 

#create density bins
density_bins = np.logspace(np.min(np.log10(sigma5)),np.max(np.log10(sigma5)),no_bins)

#create figure to plot over 
fig= plt.figure(figsize = (4.0,6.5),facecolor='w',edgecolor='w')

#create subplots
p1 = plt.subplot(211)
p2 = plt.subplot(212)

#loop through all density bins
for i in range(0,len(density_bins)-1):

	#density limits
	den_min = density_bins[i]#*0.2 
	den_max = density_bins[i+1]#*1.8

	#select cat
	selection = cat.where((sigma5 > den_min) & (sigma5 < den_max))

	#send cat to plotta
	double_plot([p1,p2], [0.,1.1,0.25], selection, 'pS0', [-0.09, 1.09], [0,100], False, True, colours[i], lstr(den_min) + '<log10sigma5<' + lstr(den_max) )




plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
#fig.savefig(pj('/Users/chrisfuller/Desktop/','per_det_accros_hubble_seq_zoo.pdf'))
plt.show()



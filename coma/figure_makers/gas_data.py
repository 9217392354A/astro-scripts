#creates plot of atomic gas from literature
# Chris Fuller, June 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input nametest.fits

cat = Table(pj(folder,fname))

alfaCat = cat.where(np.nan_to_num(cat.GMASS) > 0.0)
litgCat = cat.where((np.nan_to_num(cat.HI_ALL2) > 0.0) & (np.nan_to_num(cat.GMASS) == 0.0))

gas1 = alfaCat.HI_ALL2
gas2 = litgCat.HI_ALL2
gas3 = cat.HI_ALL2



#clean
gas1 = gas1[gas1>0.0]
gas2 = gas2[gas2>0.0]
gas3 = gas3[gas3>0.0]

bins = np.linspace(min(gas3), max(gas3), 10)

#create histograms
hist1, binEdges = np.histogram(gas1, bins=bins)
hist2, binEdges = np.histogram(gas2, bins=bins)
hist3, binEdges = np.histogram(gas3, bins=bins)


#plot
fig = plt.figure(figsize = (8.0,4.5),facecolor='w',edgecolor='w')
sub = fig.add_subplot(1,1,1)

w = binEdges[1] - binEdges[0] 


sub.bar(binEdges[:-1], hist1, width=w, fc ='None', ec = 'b', hatch='/')
sub.bar(binEdges[:-1], hist2, width=w, fc ='None', ec = 'r', hatch='\\')
sub.bar(binEdges[:-1], hist3, width=w, fc ='None', ec = 'k', lw=2)

sub.set_ylabel('N')
sub.set_xlabel('$\log_{10} (M_{gas} / $M$_{\odot}$)')

plt.subplots_adjust(left=0.11, bottom=0.12, right=0.95, top=0.97, wspace=0.0, hspace=0.0)
fig.savefig('/Users/chrisfuller/Dropbox/phd/thesis/main/chapter-coma/gas-data.pdf')

plt.show()



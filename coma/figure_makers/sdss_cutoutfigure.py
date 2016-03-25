#figure maker
# Chris Fuller, April 2013

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import atpy as at
from os.path import join as pj
import os

print 'reading in cats'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = at.Table(pj(folder,cat_name))

cat = cat.where((cat.goldmine < 2) & (cat.goldmine > -10) & (cat.late == 1))

cat.sort('goldmine')
print cat.goldmine
#folder
outFolder = '/Users/chrisfuller/Desktop/'

os.chdir(outFolder)


fig, subs = plt.subplots(nrows=2, ncols=4, sharex=False, sharey=False, squeeze=False, figsize = (8.5,4.5), facecolor='w',edgecolor='w')

def clear_frame(ax=None):
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.itervalues():
        spine.set_visible(False)

cat1 = cat.where(cat.goldmine==0)
cat1.sort('pEdg')

for i in range(4):
	#left
    name = cat1.OBJECT[i]
    coord = cat1.GRA2000[i], cat1.GDEC2000[i]

    outFile = name  + "-SDSS.jpg"
    img=mpimg.imread(outFile)

    
    subl=subs[0,i]
    subl.imshow(img, interpolation='nearest')
    clear_frame(subl)
    subl.text(0.05, 0.95,  name, transform=subl.transAxes, fontsize=12, verticalalignment='top', color='white')

    #right
    #subr=subs[i,1]

ii = range(4,8)

cat2 = cat.where(cat.goldmine==1)
cat2.sort('pEdg')
for j in range(4):
	i = j
	#left
	name = cat2.OBJECT[j]
	coord = cat2.GRA2000[j], cat2.GDEC2000[j]

	outFile = name  + "-SDSS.jpg"
	img=mpimg.imread(outFile)


	subr=subs[1,i]
	subr.imshow(img, interpolation='nearest')
	clear_frame(subr)
	subr.text(0.05, 0.95,  name, transform=subr.transAxes, fontsize=12, verticalalignment='top', color='white')

	#right
	#subr=subs[i,1]

subs[0,0].text(0.05, 0.15,  'goldmine:E', transform=subs[0,0].transAxes, fontsize=13, verticalalignment='top', color='white')
subs[1,0].text(0.05, 0.15,  'goldmine:S0', transform=subs[1,0].transAxes, fontsize=13, verticalalignment='top', color='white')

plt.subplots_adjust(left=0.0, bottom=0.0, right=1., top=1., wspace=0.01, hspace=0.0)
fig.savefig(pj('/Users/chrisfuller/Desktop/','gold_early_zoo_late.pdf'))
plt.show()

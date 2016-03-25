#plot goldmine hist of early late and inter
# Chris Fuller, March 2014

#import
print 'importing modules...'
from atpy import Table
from numpy import histogram, arange, sqrt
from os.path import join as pj
from matplotlib.pyplot import subplots, show, subplots_adjust, legend
from matplotlib.ticker import MaxNLocator

#Inputs
print 'reading in cats'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(pj(folder,cat_name))

early = cat.where(cat.early == 1)
late  = cat.where(cat.late  == 1)
inter = cat.where(cat.inter == 1)


cats = [early, late, inter]
cols = ['r', 'b', 'g']

morph = ['early', 'late', 'inter']

fig, subs = subplots(nrows=1, ncols=1, sharex=True, sharey=False, squeeze=True, figsize = (8.0,4.5), facecolor='w',edgecolor='w')

sub1 = subs

bins = arange(0,10,1)

total, _ = histogram(cat.goldmine, bins=bins)

for i in range(3):
	hist_N, _ = histogram(cats[i].goldmine, bins=bins)
	sub1.errorbar(bins[:-1], hist_N, yerr=sqrt(hist_N), c=cols[i], ls='-', label=morph[i])
	print morph[i], hist_N

sub1.set_xlabel('Morphological Type')
sub1.set_ylabel('N')
sub1.set_xticklabels( ('E', 'S0', 'S0/Sa', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd') )
subplots_adjust(left=0.1, bottom=0.12, right=0.96, top=0.98, wspace=0.0, hspace=0.0)
#legend()




fig.savefig(pj('/Users/chrisfuller/Dropbox/','goldmine_vs_zoo.pdf'))

show()






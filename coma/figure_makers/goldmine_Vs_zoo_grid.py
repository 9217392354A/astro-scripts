#plot goldmine hist of early late and inter
# Chris Fuller, March 2014

#import
print 'importing modules...'
from atpy import Table
from numpy import histogram, arange, sqrt, linspace, min, max, around, nan_to_num
from os.path import join as pj
from matplotlib.pyplot import subplots, show, subplots_adjust, legend
from matplotlib.ticker import MaxNLocator

#Inputs
print 'reading in cats'
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat_name = 'coma_supercluster_cal12.fits' #input name
cat_master = Table(pj(folder,cat_name))

cols = ['r', 'b', 'g']
N=2

fig, subs = subplots(nrows=N, ncols=1, sharex=True, sharey=False, squeeze=True, figsize = (4.5,8.5), facecolor='w',edgecolor='w')

sub1 = subs[N-1]

bins = arange(0,10,1)

smass_bins = linspace(min(cat_master.SIGMA10), max(cat_master.SIGMA10), N+1)


for j in range(N):
	lower = smass_bins[j]
	upper = smass_bins[j+1]

	cat  = cat_master.where((cat_master.SIGMA10 >= lower) & (cat_master.SIGMA10 <= upper)  & ( nan_to_num(cat_master.F250) > 0.0))

	sub = subs[j]

	sub.plot(cat.pS0, cat.D2S, 'kx') #yerr=sqrt(hist_N), c=cols[i], ls='-')#, label=morph[i])

	sub.text(0.02, 0.9,  '$'+ str(around(lower)) + '< $log_{10} (\SIGMA_{10})$ < ' + str(around(upper)) + '$', transform=sub.transAxes, fontsize=12, verticalalignment='top')
sub1.set_xlabel('Morphological Type')
sub1.set_ylabel('$log(M_{dust}/M_{stars})$')
#sub1.set_xticklabels( ('E', 'S0', 'S0/Sa', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd') )
subplots_adjust(left=0.11, bottom=0.12, right=0.96, top=0.98, wspace=0.0, hspace=0.0)
legend()




#fig.savefig(pj('/Users/chrisfuller/Desktop/','goldmine_vs_zoo_total.pdf'))

show()






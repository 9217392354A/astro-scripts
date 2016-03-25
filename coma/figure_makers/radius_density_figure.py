# Radius vs Density plot
# Chris Fuller, April 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MaxNLocator


#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))
cat.VELOCITY_1 = cat.VELOCITY_1/1000.0

if False:
	cat = cat.where((cat.RADIUS_VIR > 1.0))


types = ['late', 'inter', 'early']
colours = ['b', 'g', 'r']
sigmas = ['SIGMA1', 'SIGMA5', 'SIGMA10']
sig_labs = ['$\log_{10}$($\Sigma_{1}$) (Mpc$^{-2}$)', '$\log_{10}$($\Sigma_{5}$) (Mpc$^{-2}$)', '$\log_{10}$($\Sigma_{10}$) (Mpc$^{-2}$)']

fig, subs = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True, squeeze=True, figsize = (8,9), facecolor='w',edgecolor='w')


for j in range(3):
	sub = subs[j]
	for i in range(3):
		selection = cat.where(cat[types[i]]==1)

		x = np.log10(selection.RADIUS_VIR)
		y = selection[sigmas[j]]

		sub.scatter(x, y, s=10, marker = 'o', color = colours[i], alpha=0.5)
	sub.set_ylim(-1.5, 4.3)
	sub.set_xlim(-2.7, 1.1)
	sub.axvline(x=np.log10(1.0), ls = '--', c='k')

	sub.set_ylabel(sig_labs[j])

subs[2].set_xlabel('$\log_{10}$(R/R$_{Virial}$)')
plt.subplots_adjust(left=0.08, bottom=0.06, right=0.99, top=0.99, wspace=0.0, hspace=0.0)
#fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/coma/radius_vs_sigmas.pdf')
plt.savefig(pj('/Users/chrisfuller/Desktop/','radius_vs_sigmas.pdf'))
plt.show()
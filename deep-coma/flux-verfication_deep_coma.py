# flux verfication plotting tool
# Chris Fuller, 2013

import numpy as np
import matplotlib.pyplot as plt
import atpy as at
from os.path import join as pj

#import stuff
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/flux_verication/"
cat = at.Table(pj(folder,"flux_verfication-all-v2.fits"),type='fits')






fa = list(reversed(["F500","F350","F250","F160","F100"]))

fb = list(reversed([["F500_NGP"],["F350_NGP"],["F250_NGP"], ['F160_DEEP', 'F160_DEEP_MY-METHOD'],['F100_DEEP', 'F100_DEEP_MY-METHOD', 'F100_IRAS-FAINT', 'F100_IRAS-POINT']]))
colour = list(reversed([["black"],["black"],["black"], ['red','purple'],['red', 'purple', 'green', 'blue']]))

names = list(reversed(["500$\mu m$","350$\mu m$","250$\mu m$","160$\mu m$","100$\mu m$"]))


#function for finding errors in a straight line of two arras x and y
# assuming least squares and equal error on each point
def errorInLine(x_data,y_data):
    x_data = np.array(x_data, dtype=np.float)
    y_data = np.array(y_data, dtype=np.float)
    
    n = len(x_data)
    D = np.sum(x_data**2) - 1./n * np.sum(x_data)**2
    x_bar = np.mean(x_data)
    p_coeff, residuals, _, _, _ = np.polyfit(x_data, y_data, 1, full=True)
    dm_squared = 1./(n-2)*residuals/D
    dc_squared = 1./(n-2)*(D/n + x_bar**2)*residuals/D
    
    print fa[j],' vs.', fb[j][i], ' m=', np.round(p_coeff[0], decimals=3), ' dm=', np.round(np.sqrt(dm_squared[0]), decimals=3),' c=', np.round(p_coeff[1], decimals=3), ' dc=', np.round(np.sqrt(dc_squared[0]), decimals=3)

    return p_coeff
#make empty lists to hold fluxes
fluxA = []
fluxB = []
fits = []

fig = plt.figure(figsize = (4.5,4.5),facecolor='w',edgecolor='w')
ax = fig.add_subplot(111)

# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

#loop through all bands
for j in range(0,2):
	sub = fig.add_subplot(1,2,j+1)
	#loop through comparsion sample
	sub.plot(np.logspace(-3,3), np.logspace(-3,3), '--k', alpha=0.5)
	for i in range(0,len(fb[j])):
		

		#find where both are detected
		if fb[j][i].split('_')[1] != 'NGP': w = np.where((cat[fa[j]] > 0.0) & (np.nan_to_num(cat[fb[j][i]]) > 0.0))[0]
		if fb[j][i].split('_')[1] == 'NGP': 
			#create empty list
			w_point = []

			#determine what extendedness parameter to use
			if fa[j][-3:] == '250': ext_index = 2
			if fa[j][-3:] == '350': ext_index = 1
			if fa[j][-3:] == '500': ext_index = 0

			#loop through all souces and determine if point or exteneded
			for k in range(0, len(cat[fa[j]])):
				if cat.EXTENDEDNESS[k][ext_index] == 'P': w_point.append(1)
				if cat.EXTENDEDNESS[k][ext_index] == 'E': w_point.append(0)

			#turn w_point into numpy array
			w_point = np.array(w_point, dtype=np.int8) 


			#find where both are greater than 0 and sn>3 and point as ngp exteneded not relaiable in current release
			my_fluxes = cat[fa[j]]
			lit_fluxes = np.nan_to_num(cat[fb[j][i]])
			lit_sn = lit_fluxes / cat['E'+fb[j][i][1:]]

			w = np.where((my_fluxes > 0.0) & (lit_fluxes > 0.0) & (lit_sn > 5.0) & (w_point == 1))[0]

		


		x,y, e_x, e_y = cat[fa[j]][w], cat[fb[j][i]][w], cat['E'+fa[j][1:]][w], cat['E'+fb[j][i][1:]][w]
		c = colour[j][i] 
		fit__ =  errorInLine(x,y)

		sub.errorbar(x, y, yerr=e_y, xerr=e_x, c=c,marker='x', alpha=0.3, ls='none')
		#plt.loglog(x,y,'x', c=c)


		#xd,yd = log10(xdata),log10(ydata)
		#polycoef = polyfit(xd, yd, 1)
		#yfit = 10**( polycoef[0]*xd+polycoef[1] )
		
		fit = np.polyfit(np.log10(x),np.log10(y),deg=1)
		p = np.poly1d(fit)
		#x = myrange(x,y)
		x = np.linspace(np.log10(np.min(x)*0.9),np.log10(np.max(x)*1.1))
		sub.loglog(10**x,10**p(x), c=c, label=names[j]+", m = "+str(np.round(fit[0],decimals=3))+fb[j][i])
		
		#if j == 0 or j==2:
		#	sub.tick_params(axis='xaxis', labelbottom='off')
		if j ==0:
			sub.set_xlabel("Fuller Flux (Jy)")
			sub.set_ylabel("Hickinbottom Flux (Jy)")

		if j == 1:
			sub.tick_params(axis='y', labelleft='off')
			sub.set_xlabel("Fuller Flux(Jy)")



		


	sub.set_xlim(xmin=10.0**-2.25, xmax=10**1.2 )
	sub.set_ylim(ymin=10.0**-2.25, ymax=10**1.2 )

	sub.text(0.05, 0.95, names[j], transform=sub.transAxes, fontsize=14, verticalalignment='top')

plt.subplots_adjust(left=0.135, bottom=0.08, right=0.98, top=0.98, wspace=0.0, hspace=0.0)

plt.savefig(pj("/Users/chrisfuller/Desktop/","flux_verification.pdf"),dpi=600.0)
plt.show()

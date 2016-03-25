# programs to invistigate trends with density
# Chris Fuller, March 2014

#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'

c1 = 'coma_supercluster_cal12.fits' #input name
c3 = 'fornax_at_100mpc.fits'
c2 = 'virgo_at_100mpc.fits'

no_density = 5 #number of sigma and radial density bins between logmin and logmax
sfactor = 0.0
bin_type = 'fixed'  
N = 40 #number in each bin

cats = [Table(pj(folder,c1)), Table(pj(folder2,c2)), Table(pj(folder2,c3))] 
cat_names = ['Coma', 'Virgo', 'Fornax']
morphsys = ['zoo', 'goldmine', 'goldmine']


#xcols  = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']
xcols  = ['SIGMA1_1000', 'SIGMA5_1000', 'SIGMA10_1000', 'RADIUS_VIR']
xlabs = ['$\log_{10}(\Sigma 1)($Mpc$^{-2})$','$\log_{10}(\Sigma 5)($Mpc$^{-2})$','$\log_{10}(\Sigma 10) ($Mpc$^{-2})$', 'Projected Cluster Radius (Mpc)']
xlims = [[-1.2, 4.1], [-1.2,3.3], [-1.2,3.3], [-2.7, 1.5]] 

#quanties on the yaxis
#quanties on the yaxis
if True: 
	ycols = [['DMASS'], ['DMASS', 'SMASS']]
	ylabs = ['$Log_{10} (M_{dust} / $M$_{\odot}$)', '$Log_{10} (M_{dust} / M_{star}$)']
	ylims = [[5.0, 9.1], [-4.8, -1.4 ]]
	yclean = ylims


if False:#full columns
	ycols = [['DMASS', 'SMASS'], ['SRF'], ['sSFR'], ['GMASS', 'SMASS'], ['METAL'],['g','r'], ['SMASS']]
	ylabs = ['$Log_{10} (M_{dust} / M{star}$)', '$Log_{10} (SRF)$', '$Log_{10} (sSFR)$', '$ Log_{10} (M_{HI} / M_{star}) $', '$Log_{10} (O/H)$', 'g-r', '$Log_{10}(M_{star}/$M$_{\odot})$']
	ylims = [[-4.8, -1.4  ],[-3.1, 1.1 ],[-13.2,-9.1 ], [-2.,2. ], [7.9, 9.5], [0.1,0.95], [7.6,11.5]]
	clean_y = [[-4.8, -1.4],[-3.1, 1.1],[-13.2,-9.1], [-20.,20.], [7.9, 9.5], [-200.0,200.0], [-200.0,200.0]]



#create figure and subplots
fig, subs = plt.subplots(nrows=len(ycols), ncols=len(xcols), sharex=False, sharey=False, squeeze=False, figsize = (8.,4.), facecolor='w',edgecolor='w')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  Functions  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def adaptive_histogram(x, y, N):
	#intialise list to hold y values
	y_mean, x_mid, y_error, y_n, x_mean = [], [], [], [], []

	#sort arrays
	ii = np.argsort(x)

	#create new arrays sorted
	xx = x[ii]
	yy = y[ii]

	#empty list of bins
	bins = []

	#caculate how many bins needed
	number = np.int(np.ceil(len(xx)/ np.float(N)))

	print 'Number of bins: ', number
	#create intial lower and upper limits
	lower_i = 0
	upper_i = N

	#loop through and add N to each bin
	for i in range(number):
		#if not final bin
		if i != number -1:
			#create index upper and lower limits
			lower_i, upper_i = N*i, N*(i+1) - 1

		#if final bin
		else:
			lower_i, upper_i = N*i, len(ii) - 1 

		#physical bin widths
		lower, upper = xx[lower_i], xx[upper_i]

		#find width of bin
		width = upper-lower

		#select y values that fall in bins
		y_vals = yy[lower_i:upper_i+1]

		#mean, error in the mean, and midbin val
		y_n.append(np.float(len(y_vals)))
		y_mean.append(np.mean(y_vals))
		y_error.append(np.std(y_vals)/ np.sqrt(len(y_vals)))
		x_mid.append((lower+upper)/2.0)
		x_mean.append(np.mean(xx[lower_i:upper_i+1]))
	
	return y_n, y_mean, y_error, x_mid, x_mean

def bin_plotter(x,y, bin_lims, subplot, colour):
	#intialise list to hold y values
	y_mean, x_mean, y_error, y_n = [], [], [], []

	# caculate for fixed bin width
	if bin_type == 'fixed':
		#create bins
		bins = np.linspace(bin_lims[0], bin_lims[1], no_density+1)
		#loop through bins and bin up
		for i in range(len(bins)-1):
			#limits
			lower, upper = bins[i], bins[i+1]
			width = upper-lower
			#select y values that fall in bins
			y_vals = y[(x > lower-width*sfactor ) & (x < upper+width*sfactor)]

			#mean, error in the mean, and midbin val
			y_n.append(np.float(len(y_vals)))
			y_mean.append(np.mean(y_vals))
			y_error.append(np.std(y_vals)/ np.sqrt(len(y_vals)+0.001))
			x_mean.append((lower+upper)/2.0)


	if bin_type == 'adaptive':
		print 'adaptive method...'
		y_n, y_mean, y_error, x_mid, x_mean = adaptive_histogram(x, y, N)		

	print 'N', y_n
	print 'Err:', y_error

	#now plot up the values
	subplot.errorbar(x_mean, y_mean, yerr = y_error, color = colour, ls=' ')
	subplot.scatter(x_mean, y_mean, color=colour, marker = 'o')
	#subplot.errorbar(x_mid, y_mean, yerr = y_error, color = colour)
	#subplot.scatter(x , y , s = 20, c = colour, marker = '+', alpha = 1.)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#loop through all the clusters 
for cluster_i in range(1):
	cat = cats[cluster_i] #pull out cat from list of cats


	# create two catalogues early and late type galxies, and morphsys determines the morphological sep method
	if morphsys[cluster_i] == 'zoo':
		#create early and late galaxies
		early = cat.where(cat['pE0'] > 0.8)
		late  = cat.where(cat['pS0'] > 0.8)
		inter = cat.where((cat['pS0'] < 0.8) & (cat['pE0'] < 0.8) )

	if morphsys[cluster_i] == 'goldmine':
		#create early and late galaxies
		early = cat.where(cat['goldmine'] >= 2)
		late  = cat.where(cat['goldmine'] <  2) 




	#loop through all y values
	for j in range(len(ycols)):
		ycol = ycols[j]

		#two columns
		if len(ycol) == 2:
			y_early = early[ycol[0]] - early[ycol[1]] #this will only work for log10 otherwise divide is needed
			y_late  =  late[ycol[0]] -  late[ycol[1]]
			y_inter = inter[ycol[0]] - inter[ycol[1]]

		#single col
		if len(ycol) == 1:
			y_early = early[ycol[0]]
			y_late =   late[ycol[0]]
			y_inter = inter[ycol[0]]

		if ycol[0] == 'DMASS': #little hack that bassicly only plots galaxies detected at 250
			w_early = np.where(early[ycol[0]]>0.0)[0]
			w_late  = np.where( late[ycol[0]]>0.0)[0]
			w_inter = np.where(inter[ycol[0]]>0.0)[0]

		else: #if not equal to dmass then check if Ben is hiding in any of the arragys
			w_early = np.where( (clean_y[j][1] > early[ycol[0]]) & ( clean_y[j][0] < early[ycol[0]]) & ( early[ycol[0]] != 999999.0))[0]
			w_late  = np.where( (clean_y[j][1] >  late[ycol[0]]) & ( clean_y[j][0] <  late[ycol[0]]) & (  late[ycol[0]] != 999999.0))[0]
			w_inter = np.where( (clean_y[j][1] > inter[ycol[0]]) & ( clean_y[j][0] < inter[ycol[0]]) & ( inter[ycol[0]] != 999999.0))[0]


			#if not then it will return a full list of indicies
			#w_early = np.where(early[ycol[0]] != 999999.0)[0]
			#w_late  = np.where( late[ycol[0]] != 999999.0)[0]
			#w_inter = np.where(inter[ycol[0]] != 999999.0)[0]

		#loop through x cols and plot
		for i in range(len(xcols)):
			xcol = xcols[i]

			#extract xcols
			x_early = early[xcol]
			x_late =   late[xcol]
			x_inter = inter[xcol]



			sub = subs[j,i] #select subplot

			#plot x and y early and late
			bin_plotter(x=x_late[w_late],y=y_late[w_late], bin_lims=xlims[i], subplot = sub, colour = 'b')
			bin_plotter(x=x_inter[w_inter],y=y_inter[w_inter], bin_lims=xlims[i], subplot = sub, colour = 'g')
			bin_plotter(x=x_early[w_early],y=y_early[w_early], bin_lims=xlims[i], subplot = sub, colour = 'r')
			
			#plotting options
			sub.xaxis.set_major_locator(	MaxNLocator(4)	)
			sub.yaxis.set_major_locator(	MaxNLocator(4)	)

			#lable x on final y
			if j == len(ycols)-1:
				sub.set_xlabel(xlabs[i])

			#label y on final x
			if i == 0:
				sub.set_ylabel(ylabs[j])

			#if not on bottom switch off x axis
			if j != len(ycols)-1:
				sub.tick_params(axis = 'x', labelbottom = 'off')

			#if not on the left switch off y axis
			if i != 0:
				sub.tick_params(axis = 'y', labelleft = 'off')

			#set xlims
			#sub.set_xlim(xmin = xlims[i][0], xmax = xlims[i][1])

			#set ylims
			if ylims[j][0] != '': sub.set_ylim(ymin = ylims[j][0], ymax = ylims[j][1])




##### plotting options
plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)



fig.savefig('/Users/chrisfuller/Desktop/trends_density.pdf')
plt.show()
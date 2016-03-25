# programs to invistigate trends with density
# Chris Fuller, March 2014

#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from lmfit import minimize, Parameters,report_fit
from scipy.stats import pearsonr 
#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'

c1 = 'coma_supercluster_cal12.fits' #input name
c3 = 'fornax_at_100mpc.fits'
c2 = 'virgo_at_100mpc.fits'

cats = [Table(pj(folder,c1)), Table(pj(folder2,c2)), Table(pj(folder2,c3))] 
cat_names = ['Coma', 'Virgo', 'Fornax']
morphsys = ['zoo', 'goldmine', 'goldmine']

#quanties on x axcis
xcols  = ['SIGMA1', 'SIGMA5', 'SIGMA10', 'RADIUS_VIR']
xlabs = ['$Log_{10}(\Sigma 1)($Mpc$^{-2})$','$Log_{10}(\Sigma 5)($Mpc$^{-2})$','$Log_{10}(\Sigma 10) ($Mpc$^{-2})$', 'Projected Cluster Radius (Mpc)']
xlims = [[-1.2, 4.1], [-1.2,3.3], [-1.2,3.3], [-2.7, 1.5]] 

if False:
	#quanties on the yaxis
	ycols = [['DMASS_250', 'SMASS'], ['SRF'], ['sSFR'], ['GMASS', 'SMASS'], ['METAL']]
	ylabs = ['$Log_{10} (M_{dust} / M{star}$)', '$Log_{10} (SRF)$', '$Log_{10} (sSFR)$', '$ Log_{10} (M_{HI} / M_{star}) $', '$Log_{10} (O/H)$']
	ylims = [[-4.8, -1.4  + 2.0],[-3.1, 1.1 + 2.0],[-13.2,-9.1 + 2.0], [-2.,2. + 2.0], [7.9, 9.5 + 2.0]]
	clean_y = [[-4.8, -1.4],[-3.1, 1.1],[-13.2,-9.1], [-20.,20.], [7.9, 9.5]]


if True:
	#quanties on the yaxis
	ycols = [['DMASS', 'SMASS']]
	ylabs = ['$Log_{10} (M_{dust} / M{star}$)']
	ylims = [[-4.8, -1.4  + 2.0]]
	clean_y = [[-4.8, -1.4]]

	
#create figure and subplots
fig, subs = plt.subplots(nrows=len(ycols), ncols=len(xcols), sharex=False, sharey=False, squeeze=False, figsize = (10.5,10.5), facecolor='w',edgecolor='w')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  Functions  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def residual(params, x, y_data):
    m = params['m'].value
    c = params['c'].value
    y_model = m*x + c
    return (y_data-y_model)**2

def lmfitter(x , y):
    params = Parameters()
    params.add('m', value=-1., vary=True)
    params.add('c', value=1., vary=True)
    out = minimize(residual, params, args=(x, y))
    #report_fit(params)
    return out

def model(params, x):
    m = params['m'].value
    c = params['c'].value
    return m*x + c

def fit_line(subplot, x_data , y_data, colour, x_loc, y_loc, text):
	#fit the data
	out = lmfitter(x_data , y_data)

	#plot straight line to data
	subplot.plot(x_data , model(out.params, x_data), color = colour)

	m, dm  = out.params['m'].value, out.params['m'].stderr

	#do ppc test
	pcc = pearsonr(x_data, y_data)

	s = text + 'pcc = ' + str(np.round(pcc[0],decimals=3)) + ' p=' + str(np.round(pcc[1],decimals=3)) + 'm=' + str(np.round(m,decimals=3)) + '+-' + str(np.round(dm,decimals=3))
	print s
	#show text over the top
	subplot.text(x_loc, y_loc, s, transform=subplot.transAxes, fontsize=7, color = colour, verticalalignment='top')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#loop through all the clusters 
for cluster_i in range(1):
	cat = cats[cluster_i] #pull out cat from list of cats
	cat.RADIUS_VIR = np.log10(cat.RADIUS_VIR)

	# create two catalogues early and late type galxies, and morphsys determines the morphological sep method
	if morphsys[cluster_i] == 'zoo':
		#create early and late galaxies
		early = cat.where(cat['early'] == 1.0)
		late  = cat.where(cat['late'] == 1.0)
		inter = cat.where(cat['inter'] == 1.0)

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

		if ycol[0] == 'DMASS_250': #little hack that bassicly only plots galaxies detected at 250
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
			sub.scatter(x_early[w_early], y_early[w_early], s = 5, c = 'r', marker = '+', label = 'early' + cat_names[cluster_i], alpha = 0.3)
			sub.scatter(x_inter[w_inter], y_inter[w_inter], s = 5, c = 'g', marker = '+', label = 'inter' + cat_names[cluster_i], alpha = 0.3)
			sub.scatter(x_late[w_late] ,   y_late[w_late] , s = 5, c = 'b', marker = '+', label = 'late ' + cat_names[cluster_i], alpha = 0.3)



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
			sub.set_xlim(xmin = xlims[i][0], xmax = xlims[i][1])

			#set ylims
			if ylims[j][0] != '': sub.set_ylim(ymin = ylims[j][0], ymax = ylims[j][1])

			if i != len(xcols)-1: sub.invert_xaxis()

						#plot straight line on figures, for some it wont be possible to plot
			if True:
				fit_line(subplot = sub, x_data = x_early[w_early] , y_data = y_early[w_early], colour='r',x_loc= 0.02, y_loc= 0.95, text = 'early:')

			#except:
			#	print 'error', ycols[j]

			try:
				fit_line(subplot = sub, x_data = x_inter[w_inter] , y_data = y_inter[w_inter], colour='g',x_loc= 0.02, y_loc=0.88, text = 'inter:')
				
			except:
				print 'error', ycols[j]

			try:
				fit_line(subplot = sub, x_data = x_late[w_late] , y_data = y_late[w_late], colour='b',x_loc= 0.02, y_loc=0.81, text = 'late: ')
				
			except:
				print 'error', ycols[j]


##### plotting options
plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)



fig.savefig('/Users/chrisfuller/Desktop/trends_density.pdf')
plt.show()
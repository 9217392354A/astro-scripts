# programs to invistigate trends with density
# Chris Fuller, March 2014

#import moduals
print 'reading in python modules...'
print 'atpy...'
#from atpy import Table
import atpy
print 'numpy...'
import numpy as np
print 'os.path...'
from os.path import join as pj
print 'matplotlib.pyplot...'
import matplotlib.pyplot as plt
print 'MaxNLocator...'
from matplotlib.ticker import MaxNLocator
#remove numpy runtime warings
np.seterr(invalid='ignore')

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'

c1 = 'coma_supercluster_cal12.fits' #input name
c3 = 'fornax_input.fits'#_at_100mpc.fits'
c2 = 'virgo_input.fits'#_at_100mpc.fits'

N_params = 6 #number of bins

print 'reading in cats'
cats = [atpy.Table(pj(folder,c1)), atpy.Table(pj(folder2,c2)), atpy.Table(pj(folder2,c3))] 
cat_names = ['Coma', 'Virgo', 'Fornax']
morphsys = ['zoo', 'goldmine', 'goldmine']
cluster_ls = ['-', '--', '--']
cluster_col = ['purple', 'black', 'orange']

#quanties on x axcis
density  = 'RADIUS_VIR'	# ['SIGMA10', 'RADIUS_VIR']
density_lab = '$R/$R$_{virial}$'

#quanties on the yaxis
if True: 
	paras = [['DMASS'], ['DMASS', 'SMASS']]
	ylabs = ['$Log_{10} (M_{dust} / $M$_{\odot}$)', '$Log_{10} (M_{dust} / M_{star}$)']
	ylims = [[5.5, 9.1], [-4.8, -1.4 ]]

if False:
	paras = [ ['SRF'], ['sSFR'], ['METAL'],['g','r'], ['SMASS']]
	ylabs = ['$Log_{10} (SRF)$', '$Log_{10} (sSFR)$', '$Log_{10} (O/H)$', 'g-r', '$Log_{10}(M_{star}/$M$_{\odot})$']
	ylims = [[-3.1, 1.1 ],[-13.2,-9.1], [7.9, 9.5], [0.1,0.95], [7.6,11.5]]	

if False:#full columns
	paras = [['DMASS', 'SMASS'], ['SRF'], ['sSFR'], ['GMASS', 'SMASS'], ['METAL'],['g','r'], ['SMASS']]
	ylabs = ['$Log_{10} (M_{dust} / M_{star}$)', '$Log_{10} (SRF)$', '$Log_{10} (sSFR)$', '$ Log_{10} (M_{HI} / M_{star}) $', '$Log_{10} (O/H)$', 'g-r', '$Log_{10}(M_{star}/$M$_{\odot})$']
	ylims = [[-4.8, -1.4 ],[-3.1, 1.1 ],[-13.2,-9.1], [-2.,2.], [7.9, 9.5], [0.1,0.95], [7.6,11.5]]

#create bins that are universal for each cluster
bins_density = [0, 1.4, 1.45] #, 7.0]

#create figure and subplots
fig, subs = plt.subplots(nrows=len(bins_density)-1, ncols=len(paras), sharex=False, sharey=True, squeeze=True, figsize = (3.0*len(paras),3.5*len(bins_density)-1), facecolor='w',edgecolor='w')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  Functions  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def histo(y, bin_lims, subplot, colour, style, a1, a2):
	#create histogram
	hist, bin_edges = np.histogram(y, bins=bin_lims) #bin_lims)

	#cacualte other parameters for plotting
	width = 0.5 * (bin_edges[1] - bin_edges[0])

	#plot
	subplot.axvspan(np.mean(y) - (np.std(y)/np.sqrt(len(y)*1.0)), np.mean(y) + (np.std(y)/np.sqrt(len(y)*1.0)), color = colour, alpha=0.1)
	subplot.axvline(x=np.mean(y), ls='--', color = colour, alpha = a2)
	subplot.bar(bin_edges[:-1], hist, width*2.0, edgecolor=colour, fc ='None',hatch=style, alpha = a1 ) #, ls = style)
	
	
	
	#subplot.errorbar(x_mid, y_mean, yerr = y_error, color = colour)
	#subplot.scatter(x , y , s = 20, c = colour, marker = '+', alpha = 0.3)

def bin_creator(a, N):

	if len(a) == 1: new = a
	else:
		new = a[0]
 		for i in range(1,len(a)): new = np.append(new, a[i])

 	#new is now a list of all the values that need to be binned
 	return np.linspace(np.min(new), np.max(new), N)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#loop through all the clusters 
for cluster_i in range(1):
	print 'starting cluster...', cat_names[cluster_i]
	cat = cats[cluster_i] #pull out cat from list of cats


	# create two catalogues early and late type galxies, and morphsys determines the morphological sep method
	if morphsys[cluster_i] == 'zoo':
		#create early and late galaxies
		early = cat.where(cat['pE0'] >= 0.8)
		late  = cat.where(cat['pS0'] >= 0.8)
		inter = cat.where((cat['pS0'] <= 0.8) & (cat['pE0'] <= 0.8) )

	if morphsys[cluster_i] == 'goldmine':
		#create early and late galaxies
		early = cat.where(cat['goldmine'] >= 2)
		late  = cat.where(cat['goldmine'] <  2) 




	#loop through all y values
	for j in range(len(paras)):
		hdr = paras[j]
		num_p  = 0
		#two columns
		if len(hdr) == 2:
			print 'two parameters detected'
			y_early = early[hdr[0]] - early[hdr[1]] #this will only work for log10 otherwise divide is needed
			y_late  =  late[hdr[0]] -  late[hdr[1]]
			y_inter = inter[hdr[0]] - inter[hdr[1]]
			num_p = 2
			lablab = hdr[0] + ' - ' + hdr[1] 

		#single col
		if len(hdr) == 1:
			y_early = early[hdr[0]]
			y_late =   late[hdr[0]]
			y_inter = inter[hdr[0]]
			num_p = 1
			lablab = hdr[0]

		print  'starting ' + lablab 
		print  'cleaning data...'

		#find where values are not null or not detected
		clean_early = np.where( (ylims[j][0] <   y_early) & ( ylims[j][1] >  y_early))[0]
		clean_late  = np.where( (ylims[j][0] <    y_late) & ( ylims[j][1] >   y_late))[0]
		clean_inter = np.where( (ylims[j][0] <   y_inter) & ( ylims[j][1] >  y_inter))[0]

		#extract densitys
		x_early = early[density][clean_early]
		x_late =   late[density][clean_late]
		x_inter = inter[density][clean_inter]

		#extract paramater (y)
		y_early = y_early[clean_early]
		y_late =   y_late[clean_late]
		y_inter = y_inter[clean_inter]



		# make parameter bins
		bins_params = bin_creator([y_late, y_inter, y_early], N_params)

		print bins_params
		#loop through x cols and plot
		for i in range(len(bins_density)-1):
			#sep into bins acording to density
			lower = bins_density[i]
			upper = bins_density[i+1]

			print  lablab + ' ' + str(np.round(lower, decimals=2)) + ' $\leq$ ' +  density_lab + ' $\leq$ ' + str(np.round(upper, decimals=2))

			try:
				sub = subs[i,j] #select subplot

			except:
				sub = subs[i]
			

			#plot x and y early and late
			histo(y=y_late [(x_late  > lower) & (x_late  < upper)], bin_lims=bins_params, subplot = sub, colour = cluster_col[cluster_i], style = '/', a1 = 1. , a2= 1. )
			#histo(y=y_early[(x_early > lower) & (x_early < upper)], bin_lims=bins_params, subplot = sub, colour = 'r', style = '\\', a1 = 1. , a2= 1.   )
			#histo(y=y_inter[(x_inter > lower) & (x_inter < upper)], bin_lims=bins_params, subplot = sub, colour = 'g', style = '', a1= 0.3, a2= 0.3)
			
			#plotting options
			sub.xaxis.set_major_locator(	MaxNLocator(4)	)
			sub.yaxis.set_major_locator(	MaxNLocator(4)	)

			if j == 0: sub.text(0.02, 0.9, str(np.round(lower, decimals=2)) + ' $\leq$ ' +  density_lab + ' $\leq$ ' + str(np.round(upper, decimals=2)) , transform=sub.transAxes, fontsize=12, verticalalignment='top')

			if i == len(bins_density) - 2: sub.set_xlabel(ylabs[j])

			if j == 0: sub.set_ylabel('N')

			if i != len(bins_density) - 2: sub.tick_params(axis='x', labelbottom='off')

			
##### plotting options
plt.subplots_adjust(left=0.11, bottom=0.08, right=0.98, top=0.99, wspace=0., hspace=0.0)



fig.savefig('/Users/chrisfuller/Desktop/trends_density_DUST-threebin.pdf')
plt.show()
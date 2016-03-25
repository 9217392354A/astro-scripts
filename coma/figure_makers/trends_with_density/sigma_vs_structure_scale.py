# sigma vs structure size
# program that uses nearest neigbour stats to tests how corelatted various components are with different
# density scales

#this program works with an extended optical catlogue that is crossmatched with our FIR cat


# programs to invistigate trends with density
# Chris Fuller, March 2014

#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import spatial
from scipy.stats import spearmanr, spearmanr
#import pdb

#pdb.set_trace()

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'

c1 = 'coma_supercluster_cal12.fits' #input name
c3 = 'fornax_at_100mpc.fits'
c2 = 'virgo_at_100mpc.fits'


cats = [Table(pj(folder,c1)), '', '']#Table(pj(folder2,c2)), Table(pj(folder2,c3))] 
cat_names = ['Coma', 'Virgo', 'Fornax']
morphsys = ['zoo', 'goldmine', 'goldmine']
x_mpc =  [1.77, 0.25, 0.30]

#switchs
N_start = 1
N_end = 40

#quanties on x axcis
xcols  = ['DENSITY']
xlabs = ['$Log_{10}(\Sigma N)($Mpc$^{-2})$']
xlims = [[-1.2, 4.1], [-1.2,3.3], [-1.2,3.3], [0.0, 7.1]]

#quanties on the yaxis
if False: 
	ycols = [['RADIUS']]
	ylabs = ['$Log_{10}(M_{star}/$M$_{\odot})$']
	ylims = [[0.,7.]]
	clean_y = [[-200.0,200.0]]

else:#full columns
	ycols = [['DMASS_250', 'SMASS'], ['SRF'], ['sSFR'], ['GMASS', 'SMASS'], ['METAL'],['g','r'], ['SMASS'], ['RADIUS']]
	ylabs = ['$Log_{10} (M_{dust} / M{star}$)', '$Log_{10} (SRF)$', '$Log_{10} (sSFR)$', '$ Log_{10} (M_{HI} / M_{star}) $', '$Log_{10} (O/H)$', '$m_{g} - m_{r}$', '$Log_{10}(M_{star}/$M$_{\odot})$', 'RADIUS']
	ylims = [[-4.8, -1.4  ],[-3.1, 1.1 ],[-13.2,-9.1 ], [-2.,2. ], [7.9, 9.5], [0.1,0.95], [7.6,11.5], [0.0, 7.1]]
	clean_y = [[-4.8, -1.4],[-3.1, 1.1],[-13.2,-9.1], [-20.,20.], [7.9, 9.5], [-200.0,200.0], [-200.0,200.0], [-200.0,200.0]]



#create figure and subplots
fig, subs = plt.subplots(nrows=2, ncols=len(ycols), sharex=True, sharey=False, squeeze=True, figsize = (17.5,8.5), facecolor='w',edgecolor='w')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  Functions  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def plotter_pcc(w, y, t, subplot1, subplot2, colour):
	#initlise values to hold x and y 
	xx, y_pcc, y_pstat, nn = [], [], [], []

	#create array of Ns
	Ns = np.arange(N_start, N_end+1, 1)

	for N in Ns:
		temp_x = t['DENSITY' + str(N)][w]
		
		#do pcc test
		pcc, pstat = spearmanr(temp_x, y)

		#append lists
		y_pcc.append(pcc)
		y_pstat.append(pstat) 

	subplot1.scatter(Ns, y_pcc, marker = 'x', color = colour, alpha = 0.5)
	subplot2.scatter(Ns, y_pstat, marker = 'x', color = colour, alpha = 0.5)

	if True:
		#fit a polynomial
		p1 = np.poly1d(np.polyfit(Ns, y_pcc, 1))
		p2 = np.poly1d(np.polyfit(Ns, y_pstat, 1))

		xxx = np.arange(N_start, N_end+1, 0.01)
		subplot1.plot(xxx, p1(xxx), color = colour)
		subplot2.plot(xxx, p2(xxx), color = colour)


		subplot2.axhline(y=0.05, ls='--')

		subplot1.set_ylim(-1,1)
		subplot2.semilogy()



def plotter_m(w, y, t, subplot1, subplot2, colour):
	#initlise values to hold x and y 
	xx, y_m, y_dm, nn = [], [], [], []

	#create array of Ns
	Ns = np.arange(N_start, N_end+1, 1)

	for N in Ns:
		temp_x = t['DENSITY' + str(N)][w]
		
		#do pcc test
		p1 = np.poly1d(np.polyfit(temp_x, y, 1))

		#append lists
		y_m.append(p1[0])
		y_dm.append(p1[1]) 

	subplot1.scatter(Ns, y_m, marker = 'x', color = colour, alpha = 0.5)
	subplot2.scatter(Ns, y_dm, marker = 'x', color = colour, alpha = 0.5)



def distance(ra1, dec1, ra2, dec2):
	delta_ra = (ra1 - ra2) * np.cos(np.radians((dec1+dec2)/2.0))
	delta_dec = (dec1 - dec2)

	return np.sqrt(delta_ra**2.0 + delta_dec**2.0)

def sigma_n(d, N):
	return np.log10(N / (np.pi*(d)**2.0) )

def nearest_neigbours_n(t, N_start, N_end, x):
	#create array on Ns
	Ns = np.arange(N_start, N_end+1, 1)

	#create KD tree
	tree= spatial.cKDTree(zip(t['GRA2000'], t['GDEC2000']))

	#find nearest neigbours of random cat to kd tree
	seps_list, seps_indexs = tree.query(zip(t['GRA2000'], t['GDEC2000']), k=N_end+1)



	#add new cols to table to hold sigmans
	for N in Ns: t.add_empty_column('DENSITY' + str(N), dtype = np.float)

	#now loop through each galaxy and find with SigmaN is
	for i in range(len(t)):
		#indexes of nearest neibours
		index = seps_indexs[i][1:]

		#ra of galaxy in cat
		ra, dec = t['GRA2000'][i], t['GDEC2000'][i]

		#create list the same lenth to cacualte distance
		ra_long, dec_long = np.array([ra]*len(index)), np.array([dec]*len(index))


		#caculate distance in mpc
		seps = distance(ra_long, dec_long, t['GRA2000'][index],t['GDEC2000'][index]) * x

		#caculte sigman
		sigmas = sigma_n(seps, Ns)

		#now loop through sigma and place them in table
		for j in range(len(Ns)):
			N = Ns[j]
			sigma = sigmas[j]
			col_head = 'DENSITY'+str(N)

			#add value to header
			t[col_head][i] = sigma

	return t





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#loop through all the clusters 
for cluster_i in range(1):
	raw_cat = cats[cluster_i] #pull out cat from list of cats

	#caculate nearest neigbour stats for all galaxies
	cat = nearest_neigbours_n(raw_cat, N_start, N_end, x_mpc[cluster_i])

	#cat = raw_cat 

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

		if ycol[0] == 'DMASS_250': #little hack that bassicly only plots galaxies detected at 250
			w_early = np.where(early[ycol[0]]>0.0)[0]
			w_late  = np.where( late[ycol[0]]>0.0)[0]
			w_inter = np.where(inter[ycol[0]]>0.0)[0]

		else: #if not equal to dmass then check if Ben is hiding in any of the arragys
			w_early = np.where( (clean_y[j][1] > early[ycol[0]]) & ( clean_y[j][0] < early[ycol[0]]) & ( early[ycol[0]] != 999999.0))[0]
			w_late  = np.where( (clean_y[j][1] >  late[ycol[0]]) & ( clean_y[j][0] <  late[ycol[0]]) & (  late[ycol[0]] != 999999.0))[0]
			w_inter = np.where( (clean_y[j][1] > inter[ycol[0]]) & ( clean_y[j][0] < inter[ycol[0]]) & ( inter[ycol[0]] != 999999.0))[0]


		#select subplot
		sub1 = subs[0,j]
		sub2 = subs[1,j]


		#labeling subs	
		sub1.xaxis.set_major_locator(	MaxNLocator(4)	)
		sub1.yaxis.set_major_locator(	MaxNLocator(4)	)

		sub2.xaxis.set_major_locator(	MaxNLocator(4)	)
		sub2.yaxis.set_major_locator(	MaxNLocator(4)	)

		sub1.text(0.05, 0.95, ylabs[j], transform=sub1.transAxes, fontsize=14, verticalalignment='top')
		#sub2.text(0.05, 0.95, ylabs[j], transform=sub1.transAxes, fontsize=14, verticalalignment='top')

		


		#plot pcc vs sigma n
		plotter_pcc(w_late, y_late[w_late], cat, sub1, sub2, 'b')
		plotter_pcc(w_inter, y_inter[w_inter], cat, sub1, sub2, 'g')
		plotter_pcc(w_early, y_early[w_early], cat, sub1, sub2, 'r')




##### plotting options




subs[0][0].set_ylabel('SRC')
subs[1][0].set_ylabel('p-value')
#subs[0][0].set_ylim(-1.1,1.1)
#subs[1][0].set_ylim(0,1.1)
subs[0][0].set_xlim(0.0, float(N_end))

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.99 , wspace=0.0, hspace=0.0)


cat.write('/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/full_density.fits',  overwrite=True)
#fig.savefig('/Users/chrisfuller/Desktop/trends_density.pdf')
#plt.show()
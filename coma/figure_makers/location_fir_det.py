# programs to invistigate location of dected vs undeteced galaxies in the FIR 
# Chris Fuller, April 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MaxNLocator


#Inputs
folder = "/Users/chrisfuller/Documents/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(pj(folder,fname))
cat.VELOCITY_1 = cat.VELOCITY_1/1000.0

if True:
	cat = cat.where((cat.RADIUS_VIR <= 1.0))


types = ['early', 'inter', 'late']
colours = ['r', 'g', 'b']
labels = ['Early', 'Intermediate', 'Late']
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # # # # # # # # # # # # # # #   Functions     # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def residual(params, x, y_data, y_error):
	y_model = model(params, x)
	return (y_data-y_model)/y_error

def lmfitter(x, y, y_error):
    params = Parameters()
    params.add('sigma', value=905.2/1000.0, vary=True)
    params.add('mean', value=6984.5/1000.0, vary=True)
    params.add('hmax', value=1000., vary=True)



    # remove inf values in errors
    out = minimize(residual, params, args=(x, y, y_error))
    #report_fit(out.params)
    return out

def model(params, x):
    sigma = params['sigma'].value
    mean = params['mean'].value
    hmax = params['hmax'].value
	
    return hmax*np.exp(-(x-mean)**2/(2*sigma**2))

def s(a,b):
	a = str(np.int(a*1000.0))
	b = str(np.int(b*1000.0))
	line = a + '$\pm$' + b 
	return line

def s1(a,b):
	a = str(np.round(a, decimals=2))
	b = str(np.round(b, decimals=2))
	line = a + '$\pm$' + b 
	return line

def ch(a):
	return "(" + str(np.round(a, decimals=2)) + ")"

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#create a figure
fig, subs = plt.subplots(nrows=3, ncols=2, sharex=False, sharey=False, squeeze=False, figsize = (8.,8.5), facecolor='w',edgecolor='w')

#create bins for hist smass
bins_smass = np.linspace(8.0, 12.0, 10)
w = (bins_smass[1]-bins_smass[0])

#create bins for hist vel
bins_VELOCITY_1 = np.linspace(4., 10., 15)
w_vel = (bins_VELOCITY_1[1]-bins_VELOCITY_1[0])


#create stuff for circles
b = np.arange(0, 18, 1.*(2.0))

#create list of angles for cirplot
an = np.linspace(0,2*np.pi,100)

master_x = []
master_y = []

for bi in b:
	master_x.append((bi)*np.cos(an) + 194.9531)
	master_y.append((bi)*np.sin(an) + 27.9807)

#loop through types
for j in range(3):
	print '\n'*3
	print 'type: ', types[j]
	#left = subs[j,0]
	mid = subs[j,0]
	right = subs[j,1]

	#select galaxies of type
	selection = cat.where(cat[types[j]] == 1)

	#select detected and undetected
	detected = selection.where(selection.D250 == 1)
	undetected = selection.where(selection.D250 == 0)


	#### LEFT ####
	"""
	#histogram of mass's of detected and un detected galalxies	
	h_de, _ = np.histogram(detected.SMASS, bins=bins_smass)
	h_nd, _ = np.histogram(undetected.SMASS, bins=bins_smass)

	left.bar(bins_smass[:-1], h_nd, width=w, color='k') 
	left.bar(bins_smass[:-1], h_de, width=w, edgecolor=colours[j], facecolor='None', hatch ='\\', lw=4)

	left.axvline(x=9.2, ls='--', c='grey', lw=3)

	left.set_xlim(8.0, 12.0)
	left.text(0.05, 0.95, labels[j] , transform=left.transAxes, fontsize=20, verticalalignment='top', color='grey')

	left.set_ylim(0,max(h_nd)+5)

	"""
	#### MID ####

	#create list of angles for cirplot
	an = np.linspace(0,2*np.pi,100)

	r_det, r_deterr = np.mean(detected.RADIUS_DEG), np.std(detected.RADIUS_DEG)/np.sqrt(len(detected.RADIUS_DEG)*1.0)
	r_un, r_unerr = np.mean(undetected.RADIUS_DEG),  np.std(undetected.RADIUS_DEG)/np.sqrt(len(undetected.RADIUS_DEG)*1.0)

	x_det = (r_det)*np.cos(an) + 194.9531
	y_det = (r_det)*np.sin(an) + 27.9807

	x_undet = (r_un)*np.cos(an) + 194.9531
	y_undet = (r_un)*np.sin(an) + 27.9807

	#mid.plot(master_x, master_y, color='pink', linestyle='dashed', lw=1, alpha =0.2)

	for val in range(len(master_x)):
		mid.plot(master_x[val], master_y[val], color='grey', linestyle='dashed', lw=1)

	mid.scatter(selection.GRA2000, selection.GDEC2000, s=20, edgecolor = 'k', facecolor='None')
	mid.scatter(  detected.GRA2000,   detected.GDEC2000, s=20, color = colours[j])

	#mid.text(0.05, 0.95, s1(r_un/1.8,r_unerr/1.8), transform=mid.transAxes, fontsize=12, verticalalignment='top', color = 'grey')
	#mid.text(0.05, 0.85, s1(r_det/1.8, r_deterr/1.8), transform=mid.transAxes, fontsize=12, verticalalignment='top', color = 'purple')


	mid.set_ylim(min(cat.GDEC2000), max(cat.GDEC2000))
	mid.set_xlim(min(cat.GRA2000),  max(cat.GRA2000))
	mid.invert_xaxis()

	print 'radius undetected', s1(r_un/1.8, r_unerr/1.8)
	print 'radius detected', s1(r_det/1.8, r_deterr/1.8)
	#### RIGHT ####

	#histogram of mass's of detected and un detected galalxies	
	h_de, _ = np.histogram(detected.VELOCITY_1, bins=bins_VELOCITY_1)
	h_nd, _ = np.histogram(undetected.VELOCITY_1, bins=bins_VELOCITY_1)

	right.bar(bins_VELOCITY_1[:-1], h_nd, width=w_vel, color='k') 
	right.bar(bins_VELOCITY_1[:-1], h_de, width=w_vel, edgecolor=colours[j], facecolor='None', hatch ='\\', lw=4)

	fit_undetected = lmfitter(bins_VELOCITY_1[:-1]+w_vel*0.5 , h_nd , np.sqrt(h_nd+1.0))
	fit_detected   = lmfitter(bins_VELOCITY_1[:-1]+w_vel*0.5 , h_de , np.sqrt(h_de+1.0))

	#caculate gaussian 
	xx = np.linspace(min(cat.VELOCITY_1), max(cat.VELOCITY_1), 10000)
	yy_det = model(fit_detected.params, xx)
	yy_udet = model(fit_undetected.params, xx)

	right.plot(xx, yy_udet, color='grey', ls='-', lw=3)
	right.plot(xx, yy_det, color='purple', ls='-', lw=3)

	sigma_det, sigma_det_err = fit_detected.params['sigma'].value, fit_detected.params['sigma'].stderr 
	sigma_undet, sigma_undet_err = fit_undetected.params['sigma'].value, fit_undetected.params['sigma'].stderr 

	#right.text(0.05, 0.85, s(sigma_det, sigma_det_err), transform=right.transAxes, fontsize=12, verticalalignment='top', color = 'purple')
	#right.text(0.05, 0.95, s(sigma_undet, sigma_undet_err), transform=right.transAxes, fontsize=12, verticalalignment='top', color = 'grey')

	right.set_xlim(bins_VELOCITY_1[0], bins_VELOCITY_1[-1])
	print 'detected vel dis:', s(sigma_det, sigma_det_err), ch(fit_detected.redchi)
	print 'undetected vel dis:', s(sigma_undet, sigma_undet_err), ch(fit_undetected.redchi)

	right.set_ylim(0,max(h_nd)+3)

	if j != 2:
		#left.tick_params(axis='both', labelbottom='off')
		mid.tick_params(axis='both', labelbottom='off')
		right.tick_params(axis='both', labelbottom='off')

	mid.set_xlabel('RA (J2000)')
	mid.set_ylabel('Dec(J2000)')

	if j ==2:
		#left.set_xlabel('Stellar Mass ($\log_{10}(M_{star}/$M$_{\odot})$)')
		right.set_xlabel('Velocity (x10$^{3}$km s$^{-1}$)')
	
	#left.set_ylabel('N')
	right.set_ylabel('N')

	#left.xaxis.set_major_locator(	MaxNLocator(6)	)
	right.xaxis.set_major_locator(	MaxNLocator(6)	)
	mid.xaxis.set_major_locator(	MaxNLocator(4)	)

	#left.yaxis.set_major_locator(	MaxNLocator(6)	)
	right.yaxis.set_major_locator(	MaxNLocator(6)	)
	mid.yaxis.set_major_locator(	MaxNLocator(4)	)


plt.subplots_adjust(left=0.1, bottom=0.07, right=0.97, top=0.99, wspace=0.3, hspace=0.0)
fig.savefig('/Users/chrisfuller/Desktop/location_cluster.pdf')
plt.show()



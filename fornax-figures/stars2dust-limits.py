#program to analise detected and undetected galaxies

# Chris Fuller August  2013

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#user variables
dust_lim = 4.6


#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/plots/stellar_mass/"
cat = Table(pj(folder,"stellar_mass.fits"))

#input cols
log_mass = cat['SMASS'] 
morphology = cat['MTYPE_PGC_2'] 
radius = cat['RADIUS']/60.0
ra = cat['GRA2000']
dec = cat['GDEC2000']


#this line changes dS0 into the later type by asigning them a morphology of 21
"""
w_s0 = np.where(morphology == -3)
morphology[w_s0] = 21

print 'where s0', w_s0
print 'new morph:' ,morphology[w_s0]
"""


#incides of each morphology type
dwarf = np.where((morphology >= -2) & (morphology <= -1))[0]
early = np.where((morphology >= 0) & (morphology <= 1))[0]
late  = np.where((morphology >= 2) & (morphology <= 9))[0]
later = np.where((morphology >= 10) & (morphology <= 21))[0]

#detected at 250
dwarf_detected = np.where((morphology >= -3) & (morphology <= -1) & (np.nan_to_num(cat['F250'])>0))[0]
early_detected = np.where((morphology >= 0) & (morphology <= 1) & (np.nan_to_num(cat['F250'])>0))[0]
late_detected = np.where((morphology >= 2) & (morphology <= 9) & (np.nan_to_num(cat['F250'])>0))[0]
later_detected = np.where((morphology >= 10) & (morphology <= 20) & (np.nan_to_num(cat['F250'])>0))[0]

#undetected at 250
dwarf_undetected = np.where((morphology >= -3) & (morphology <= -1) & (np.nan_to_num(cat['F250'])==0))[0]
early_undetected = np.where((morphology >= 0) & (morphology <= 1) & (np.nan_to_num(cat['F250'])==0))[0]
late_undetected = np.where((morphology >= 2) & (morphology <= 9) & (np.nan_to_num(cat['F250'])==0))[0]
later_undetected = np.where((morphology >=10) & (morphology <= 20) & (np.nan_to_num(cat['F250'])==0))[0]

#put them all into lists for loopage
types = [dwarf, early, late, later]
detected = [dwarf_detected, early_detected, late_detected, later_detected]
undetected = [dwarf_undetected, early_undetected, late_undetected, later_undetected]
types_names = ['Dwarf', 'Early', 'Late', 'Irregulars']
colours = ['yellow', 'red', 'blue', 'green']

###################################################### Functions ##############################################################
#bin maker so that they all have the same area
def binmaker(maxval, initial_step):
	#area of initial bin
	r0 = initial_step**2

	r_old = 0
	r_new = initial_step
	count = 0

	rlist = [0,initial_step]
	#loop through untill maxval
	while r_old <= maxval:
		count += 1

		#create new picture
		r_new = np.sqrt(r0 + r_old**2)

		#append new value
		rlist.append(r_new)
		#print r_old, np.pi*(r_new**2 - r_old**2)
		r_old = r_new

		if count == 20:
			r_old = 1000000000.0
			break

	return np.array(rlist,dtype=np.float)

# find the centers of the bins
def bin_param(bins):
	cen   = []
	width = [] 
	area  = []
	for i in range(0,len(bins)-1):
		#find bin width
		bin_width = (bins[i+1]-bins[i])
		width.append(bin_width)
		#add to centre of bin
		cen.append(bins[i] + bin_width*0.5)
		#find bin area
		area.append(np.pi*(bins[i+1]**2 - bins[i]**2))

	return cen, width, area


def print_lay(a,b):

	try:
	    val1 = np.mean(a)
	    error1 = np.std(a)/len(a)
	    val2 = np.mean(b)
	    error2 = np.std(b)/len(b)

	    val3= np.mean(np.append(a,b))
	    error3 = np.std(np.append(a,b)) / len(np.append(a,b))

	    print str(np.round(val1, decimals=3))+'\\,$\\pm$\\,'+str(np.round(error1, decimals=3)) + '&' +str(np.round(val2, decimals=3))+'\\,$\\pm$\\,'+str(np.round(error2, decimals=3))+ '&' + str(np.round(val3, decimals=3))+'\\,$\\pm$\\,'+str(np.round(error3, decimals=3)) +'\\\\'

	except: 
		val1 = np.mean(a)
		error1 = np.std(a)/len(a)
	    #print str(np.round(val1, decimals=2))+'\\,$\\pm$\\,'+str(np.round(error1, decimals=2)) 

##############################################################################################################################

											# part 1 #
					#caculated limits for dust to stars for all of these galaxies

##############################################################################################################################

#loop through each type
for i in range(0, len(types)):
	#if we assume that they all have a dust mass of the limit then what is the lower limit on stars to dust ratio 
	s2d = log_mass[undetected[i]] - 4.6 

	#now select galaxies with a stars to dust less than our milkway
	try:
		print 'undetected', types_names[i], " mean: ", np.mean(s2d[np.where(s2d>3.)[0]]), " sd: ", np.std(s2d[np.where(s2d>3.)[0]]), " range: ", np.min(s2d[np.where(s2d>3.)[0]]), " - ", np.max(s2d[np.where(s2d>3.)[0]])
	except:
		print 'none undetected for ', types_names[i]

#now log at where the galaxies greater than stars to dust of 3 are detected
for i in range(0,len(types)):
	s2d = log_mass[undetected[i]] - 4.6 
	rad = radius[undetected[i]]

	s2d2 = log_mass[detected[i]] - 4.6 
	rad2= radius[detected[i]]

	s2d2 = log_mass - 4.6 
	rad2= radius
	#now selecte all galaxies that have a s2d of >3

	r_vir = 2.33 #virial radius in degrees 0.7Mpc

	#print types_names[i]
	#print 'mean-undetected stars:dust >3: ',np.mean(rad[np.where(s2d>3)[0]]), np.std(rad[np.where(s2d>3)[0]])/len(np.where(s2d>3)[0])
	#print 'mean-detected stars:dust >3: ',np.mean(rad2[np.where(s2d2>3)[0]]), np.std(rad2[np.where(s2d2>3)[0]])/len(np.where(s2d2>3)[0])


	print_lay(rad2[np.where(s2d2>3)[0]] / r_vir, rad[np.where(s2d>3)[0]]/ r_vir)
##############################################################################################################################

											# part 2 #
								# fraction of galaxies detected in each bin	

##############################################################################################################################








##############################################################################################################################

											# part 3 #
								#look at where they are in the cluster 	

##############################################################################################################################

delta = 5.
#make radius bins!!
#b = binmaker(max(radius), 30.0) #for all with equal area
b = np.arange(0, delta, 1.)
cen, widths, area = bin_param(b)

#create list of angles for cirplot
an = np.linspace(0,2*np.pi,100)

master_x = []
master_y = []

for bi in b:
	master_x.append((bi)*np.cos(an) + 54.6289)
	master_y.append((bi)*np.sin(an) - 35.4545)



figM = plt.figure(figsize=(7.0, 9.5), facecolor='w', edgecolor='k')

for i in range(0, len(types)):
	j = i+1
	#create two submarines
	sub1 = plt.subplot(4,2,j*2)
    #sub2 = plt.subplot(4,3,j*3-1)
	#sub3 = plt.subplot(4,3,j*3)

	#plot up ra and dec
	sub1.scatter(ra[undetected[i]],dec[undetected[i]], c='black', marker = 'o', s=80,facecolors='none')
	sub1.scatter( ra[detected[i]], dec[detected[i]], c=colours[i], marker = 'o', s=80)

	#plot up circles
	for k in range(0,len(b)):
		sub1.plot(master_x[k], master_y[k], color='grey', linestyle='dashed')

	#now caculated the histogtram
	hist_undetected, bin_edges = np.histogram(radius[undetected[i]], bins=b)
	hist_detected, bin_edges = np.histogram(radius[detected[i]], bins=b)

	#plot histograms
#sub2.bar(cen, hist_undetected/area, width=widths, color = 'k', align='center')
#sub2.bar(cen, hist_detected/area, width=widths, color = colours[i], align='center')

	#sub2.plot(cen,hist_undetected/area, color='k')	
	#sub2.plot(cen,hist_detected/area, color=colours[i])
	#plot fraction detected
	#fd = (hist_detected*1. / hist_undetected*1.)*100.0 +0.01
	#sub3.plot(cen, fd)

	#set plot limits
	x_cen =  54.6289
	y_cen =  -35.4545
	sub1.set_ylim( y_cen - delta, y_cen + delta )
	sub1.set_xlim( x_cen + delta, x_cen - delta )

#sub2.set_xlim(0,delta)
	#sub3.set_xlim(0,delta)
	#sub3.set_ylim(0,101)
	
	if i != 0:
		print i
		#sub2.set_ylim(0,12)

	# play around with lables
	if i != 3:
		print types_names[i]
		sub1.set_xticklabels([])
	#sub2.set_xticklabels([])
		#sub3.set_xticklabels([])
		sub1.set_yticklabels([])
		#sub2.set_yticklabels([])
		#sub3.set_yticklabels([])

	if i == 3:
		sub1.set_xlabel('RA (J2000)')
		sub1.set_ylabel('Dec (J2000)')


width = 0.5

binsa = np.arange(6.,12.0,width)
for i in range(0,len(types)):
    sub2 = plt.subplot(4,2,(i*2+1))
    hist, bin_edges = np.histogram(log_mass[types[i]], bins=binsa) 
    hist_d, bin_edges = np.histogram(log_mass[detected[i]], bins=binsa) 
    bin_centers = bin_edges[:-1]+ ((bin_edges[1] - bin_edges[0]) / 2.0)
    #plot undetected
    sub2.bar(bin_centers,hist, width=width,align ='center', color='black')
    sub2.bar(bin_centers,hist_d, width=width,align ='center', color=colours[i])
    sub2.set_xlim(5.,12.)
    #sub2.set_ylim(0,3)
        
    #plot lines of limit

    #sub2.axvline(x=6.2, color='k', ls='--')
    #sub2.axvline(x=7.2, color='k', ls='--')
    sub2.axvline(x=8.1, color='k', ls='--')
    #sub2.axvline(x=9.2, color='k', ls='--')
    sub2.axvline(x=10.1, color='k', ls='--')
    #sub2.axvline(x=11.2, color='k', ls='--')

    if i == 2:
        #sub2.text(6.2,8,'10$^1$')
        #sub2.text(7.2,8,'10$^2$')
        sub2.text(8.2,8,'10$^{-3}$')
        #sub2.text(9.2,8,'10$^4$')
        sub2.text(10.2,8,'10$^{-5}$')
        #sub2.text(11.2,8,'10$^6$')

    if i != 0: 
        sub2.set_ylim(0.,12.9)
        hight = 12*0.95
    if i == 0:
    	print types_names[i],  
    	hight = 55.0#*0.95
    if i != 3: 
        sub2.set_xticklabels([])
        
    if i == 3: sub2.set_xlabel('$\log_{10}$$(M_{star} / M_{\odot})$')
    sub2.set_ylabel('N')

    sub2.text(0.023, 0.975, types_names[i], transform=sub2.transAxes, fontsize=14, weight = 'semibold', verticalalignment='top')
    #sub2.text(5.3,names[i] )

#plt.subplots_adjust(left=0.125,bottom=0.075, top = 0.985, hspace=0.0,wspace=0.0)



plt.subplots_adjust(left=0.12, right=0.93 ,bottom=0.055, top = 0.985, hspace=0.0,wspace=0.3)
#figM.savefig(pj('/Users/chrisfuller/Desktop/','smass-detected-location.pdf'))
plt.show()



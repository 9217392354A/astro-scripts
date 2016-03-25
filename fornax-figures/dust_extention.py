# program to examine dust extention in galaxies
# Chris Fuller, 2013

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/plots/upper_limits/"
virgo  = Table(pj(folder,"virgo.fits" ))
fornax = Table(pj(folder,"fornax.fits"))

#print virgo.columns
#print fornax.columns


###################################### Functions ########################################################

# caculates dust extention
def dust_extention_precursors(cat):
	# select galaxies that are detected
	w = np.where(np.nan_to_num(cat.F250) > 0.0)[0]


	# morphology
	temp = cat.goldmine
	morph = temp[w].copy()
	R_fir = cat.R250[w]
	R_opt = cat.FULLMAJAX[w]
	mass = cat.MSTAR[w]
	dmass = cat.DMASS[w]
	ext = cat.EXTENDEDNESS[w]

	#print np.mean(np.nan_to_num(mass)), np.nan_to_num(np.mean(dmass))
	#print mass, dmass
	return R_fir, R_opt, morph, mass, dmass, ext

#dust extention
def dust_extention_main(D_fir, D_opt, morph, name, ext):
	#create an arrray that is the same lenght as our optical array and then seed it with the beam
	D_beam = D_fir.copy()
	
	D_beam[:] = beam

	#only select galaxies 5x bigger than beam

	w_big = np.where(D_opt>3*beam)[0]
	print 'number of galaxies bigger than the beam:',len(w_big), '  no of gal: ',len(D_opt)
	#split into early and late
	early = np.where((morph>=0) & (morph<= 1)& (D_opt>3*beam))[0]
	late = np.where((morph>=3) & (morph<= 20)& (D_opt>3*beam))[0]



	print name, ' early', ' min_ave = ', de(D_beam[early], D_opt[early]), ' ave = ', de(D_fir[early], D_opt[early]), ' point/extended/total (percentage point) : ', point_ext(ext[early], D_fir[early])
	print name, ' late', ' min_ave = ', de(D_beam[late], D_opt[late]), ' ave = ', de(D_fir[late], D_opt[late]), ' point/extended/total (percentage point) : ', point_ext(ext[late], D_fir[late])

	return early, late


# caculated dust extention
def de(fir,opt):
	return np.round(np.mean(fir/opt), decimals=3)

def point_ext(x, y):
	point = 0
	extended = 0


	for i in range(0,len(x)):
		val =  x[i].strip()[2]
		d_fir = float(y[i])

		if val == 'P' or d_fir <= beam : point +=1
		elif val== 'E': extended += 1

	per_point = point*100. / float(point + extended)
	return str(point) + '/' + str(extended) + '/' + str(point + extended) + '  (' + str(per_point) +')'

###################################### Main Program #####################################################

#need to caculated dust extention in each galaxy cluster
D_fornax_fir, D_fornax_opt, morph_fornax, mass_fornax, dmass_fornax, ext_fornax = dust_extention_precursors(fornax)
#D_virgo_fir, D_virgo_opt, morph_virgo, mass_virgo, dmass_virgo, ext_virgo = dust_extention_precursors(virgo)


# 250 beam width
beam = 17.8 / 60.0 # arcmin

# caculate what dust extention would be if all sources for early and late galaxies were point sources
early_fornax, late_fornax = dust_extention_main(D_fornax_fir, D_fornax_opt, morph_fornax, 'Fornax', ext_fornax)

"""
early_virgo, late_virgo = dust_extention_main(D_virgo_fir, D_virgo_opt, morph_virgo, 'Virgo', ext_virgo)

# create figures and sub figures
figM = plt.figure(figsize=(8.5,11.5), facecolor='w', edgecolor='k')
sub1 = plt.subplot(3,2,1)
sub2 = plt.subplot(3,2,2)
sub3 = plt.subplot(3,2,3)
sub4 = plt.subplot(3,2,4)
sub5 = plt.subplot(3,2,5)
sub6 = plt.subplot(3,2,6)

# dust mass vs dust extension

sub1.plot(D_fornax_fir[early_fornax]/D_fornax_opt[early_fornax], mass_fornax[early_fornax],'+', color='red')
sub1.plot(D_fornax_fir[late_fornax]/D_fornax_opt[late_fornax], mass_fornax[late_fornax], '+', color='blue')
sub1.set_ylim(6,12)
sub1.set_xlim(0, 1.5)

sub2.plot(D_virgo_fir[early_virgo]/D_virgo_opt[early_virgo], mass_virgo[early_virgo],'+', color='red')
sub2.plot(D_virgo_fir[late_virgo]/D_virgo_opt[late_virgo], mass_virgo[late_virgo], '+', color='blue')
sub2.set_ylim(8,12)
sub2.set_xlim(0, 1.5)


sub3.plot(mass_fornax[early_fornax], dmass_fornax[early_fornax]/D_fornax_fir[early_fornax]**3,'+', color='red')
sub3.plot(mass_fornax[late_fornax], dmass_fornax[late_fornax]/D_fornax_fir[late_fornax]**3,'+', color='blue')
#sub1.set_ylim(6,12)
#sub1.set_xlim(0, 1.5)
sub3.semilogy()


sub4.plot(mass_virgo[early_virgo], dmass_virgo[early_virgo]/D_virgo_fir[early_virgo]**3,'+', color='red')
sub4.plot(mass_virgo[late_virgo], dmass_virgo[late_virgo]/D_virgo_fir[late_virgo]**3,'+', color='blue')
#sub2.set_ylim(8,12)
#sub2.set_xlim(0, 1.5)
sub4.semilogy()

sub5.plot(mass_fornax[early_fornax], D_fornax_fir[early_fornax],'+', color='red')
sub5.plot(mass_fornax[late_fornax], D_fornax_fir[late_fornax],'+', color='blue')
#sub1.set_ylim(6,12)
#sub1.set_xlim(0, 1.5)
sub5.semilogy()


sub6.plot(mass_virgo[early_virgo], D_virgo_fir[early_virgo],'+', color='red')
sub6.plot(mass_virgo[late_virgo], D_virgo_fir[late_virgo],'+', color='blue')
#sub2.set_ylim(8,12)
#sub2.set_xlim(0, 1.5)
sub6.semilogy()
#plt.show()

"""



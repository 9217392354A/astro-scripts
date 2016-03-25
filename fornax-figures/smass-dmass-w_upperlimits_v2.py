# program to convert virgo and fornax catalouges to upper limits or 250um fluxes and stellars masses
# Chris Fuller, 2013

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
from scipy.stats import pearsonr 

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/plots/upper_limits/"
virgo  = Table(pj(folder,"virgo.fits" ))
fornax = Table(pj(folder,"fornax.fits"))

###################### Functions ############################
# fornax mag to mass BT mag
def fornax_mag2mass(mag):
	return -0.51*mag + 16.6

# virgo mag to mass B band
def virgo_mag2mass(mag):
	return -0.46772188*mag + 15.819291

# convert from 250um flux
def flux2mass(f):
	return 0.789*np.log10(f)+6.51


#fill in the blanks 
def fill_in_the_blanks(orginal_mass, fit_mass):
	w = np.where(np.nan_to_num(orginal_mass)==0.0)[0]

	orginal_mass[w] = fit_mass[w]
	return orginal_mass

#fill in the blanks dust mass
def fill_in_the_blanks_dust(x,flux):
	dmass = x.copy()
	fluxes = np.nan_to_num(flux.copy())
	what_dust = []
	#loop through each

	for i in range(0,len(dmass)):
		if dmass[i] != 0.: # if it has a dust mass asign a 1 to the what dust and move on
			what_dust.append(1)
			continue 

		elif (dmass[i] == 0.) and (fluxes[i] !=0.): # if it has no dust sed fit but has a 250um detection then use 250um to gen dust
			what_dust.append(2)
			dmass[i] = flux2mass(fluxes[i])

		elif (dmass[i] == 0.) and (fluxes[i] ==0.): # no dust mass or 250 just us min dust value
			what_dust.append(3)
			dmass[i] = 4.6
		else: print 'errorr!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

	#convert what dust into numpy
	return dmass, np.array(what_dust)


#def myplo(sub1,dwarf,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 1,'x','yellow')
def myplot(can,morph_select, morph, dmass, what_dust, stars, dust_type,marker, colour):
	w = np.where((morph >= morph_select[0]) & (morph <= morph_select[1]) & (what_dust == dust_type))[0]

	#X and Y def
	X = stars[w]
	Y = dmass[w] - stars[w]
	can.scatter(X,Y,marker = 'o',facecolors='none', edgecolors=colour, s=20, alpha=1)
#facecolors='none', edgecolors='r'
#def myplo(sub1,dwarf,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 1,'x','yellow')
def myplot2(can,morph_select, morph, dmass, what_dust, stars, dust_type,marker, colour):
	w = np.where((morph >= morph_select[0]) & (morph <= morph_select[1]) & (what_dust == dust_type))[0]

	#X and Y def
	X = stars[w]
	Y = dmass[w] 
	can.scatter(X,Y,marker = 'o',facecolors='none', edgecolors=colour, s=20, alpha=1)
#draw 1 2 1 line that has the c cept
def draw_line(can,c):
	#make equation for line
	y = np.arange(-10,100,1)
	x = y/1. - c 

	#now plot line
	can.plot(x,y,'k--', alpha=1.)

	#add label
	yy = 4.3
	xx = yy - c +0.1
	can.text(xx, yy, '10$^'+str(line*-1)+'$')



def draw_curve_1(can,x,y,order,col):
	y = y-x
	#fit line
	fit = np.polyfit(x,y,order)
	err = errorInLine(x,y)

	print 'fit:',fit,'err', err
	#gen line
	p = np.poly1d(fit)

	x = np.arange(6,12,0.01)
	#plot
	can.plot(x, p(x), '-', color=col)

def draw_curve(can,x,y,order,break_val,col):
	"""
	min_vals =np.arange(6.5,14.5,2.)
	for mi in min_vals:
		ma = mi + 2.

		w = np.where((x>mi) & (x<ma))
		try:
			draw_curve_1(can,x[w],y[w],order)
		except: continue
	"""
	#select upper
	w1 = np.where(x>break_val)

	#select lower
	w2 = np.where(x<break_val)

	draw_curve_1(can,x[w1],y[w1],order,col)
	#draw_curve_1(can,x[w2],y[w2],order)
	


#selected_data_excluding   (early,3,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass)
def selected_data_excluding(morph_select,dust_type,morph, cat_all_dmass, what_dust_cat, cat_smass):
	# first for cat
	w =  np.where(((morph < morph_select[0]) | (morph > morph_select[1])) & (what_dust_cat != dust_type))[0]
	return cat_smass[w], cat_all_dmass[w]

def selected_data_including(morph_select,dust_type,morph, cat_all_dmass, what_dust_cat, cat_smass):
	# first for cat
	w =  np.where(((morph == morph_select[0]) | (morph == morph_select[1])) & (what_dust_cat != dust_type))[0]
	return cat_smass[w], cat_all_dmass[w]


def errorInLine(x_data,y_data):
    x_data = np.array(x_data, dtype=np.float)
    y_data = np.array(y_data, dtype=np.float)
    
    n = len(x_data)
    D = np.sum(x_data**2) - 1./n * np.sum(x_data)**2
    x_bar = np.mean(x_data)
    p_coeff, residuals, _, _, _ = np.polyfit(x_data, y_data, 1, full=True)
    dm_squared = 1./(n-2)*residuals/D
    dc_squared = 1./(n-2)*(D/n + x_bar**2)*residuals/D
    
    return np.sqrt(dm_squared[0]),np.sqrt(dc_squared[0])

def draw_path(xx,yy):
	y = []
	sd = []
	d=0.6
	d2 = d
	x = np.arange(7.4,11.6,0.55)
	#x= np.array([ 6.81,  7.81000004,   8.76000004,   9.71000004,  10.66000004,  11.61000004], dtype=np.float)

	for i in range(0,len(x)-1):
		w= np.where((xx>=x[i]) & (xx<=x[i+1]))[0]
		y.append(np.median(yy[w]))
		sd.append(np.std(yy[w])/np.sqrt(len(w)))
		
	x = x[:-1]+(x[1]-x[0])/2.

	x = np.array(x, dtype = np.float)
	y = np.array(y, dtype = np.float)
	sd = np.array(sd, dtype = np.float)

	#sub1.scatter(xx,yy, c='k', s=100, alpha=0.1)
	sub1.plot(x,y,c='#00EEEE',linewidth = 2, ls='-')
	sub1.plot(x,y+sd*3.0 , c='#00EEEE',linewidth = 2, ls='--')
	sub1.plot(x,y-sd*3.0 , c='#00EEEE',linewidth = 2, ls='--')

def print_stats(name,x):
	print name,' mean:', np.median(x), ' +- ' ,np.std(x)/len(x)

def compare(x_virgo_ , y_virgo_ ,x_fornax_, y_fornax_ ):

	s2d_virgo_ = -y_virgo_ + x_virgo_
	s2d_fornax_ = -y_fornax_ + x_fornax_

	print_stats('virgo_',s2d_virgo_)
	print_stats('fornax_',s2d_fornax_)


def stats(d_mass, s_mass):
	s2d = s_mass - d_mass
	return ' ' +  str(np.mean(s2d)) + '+-' + str(np.std(s2d)/len(s2d))

def small_comparsion(M, vals):
	# remove upper limits from dust mass and steller mass
	w = np.where(vals[2] != 3)[0]

	morph = vals[0][w]
	dmass = vals[1][w]
	what = vals[2][w]
	smass = vals[3][w]

	#print np.max(what), 'should be 2, not 3'

	split = 8.5
	names = ['dwarf', 'early', 'late', 'later']
	for i in range(0,len(M)):
		types = M[i]
		wl = np.where((smass <= split) & (morph>=types[0]) & (morph<=types[1]))[0]
		wu = np.where((smass >= split) & (morph>=types[0]) & (morph<=types[1]))[0]

		#print out the comparsion of each 
		print types,' upper: ',stats(dmass[wu], smass[wu])
		print types,' lower: ',stats(dmass[wl], smass[wl])



def big_comparsion(MORPH, v, f):
#[virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass]
#[fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass]
	print 'Virgo'
	small_comparsion(MORPH, v)
	print 'Fornax'
	small_comparsion(MORPH, f)
###################### Main ############################

# distance ratio for virgo
D = virgo.DIST__MPC_**2 / 17.0**2

# stellar mass
fornax_smass = fill_in_the_blanks(fornax['MSTAR'], fornax_mag2mass(fornax['BTmag']))
virgo_smass = fill_in_the_blanks(virgo['MSTAR'], virgo_mag2mass(virgo['B']))

# dust mass
fornax_dmass = np.nan_to_num(fornax['DMASS']).copy()
virgo_dmass = np.nan_to_num(virgo['DMASS']).copy()

#fill in the blanks dust mass
fornax_all_dmass, what_dust_fornax = fill_in_the_blanks_dust(fornax_dmass,fornax['F250'])
virgo_all_dmass, what_dust_virgo = fill_in_the_blanks_dust(virgo_dmass,virgo['F250']*D)

# create figures and sub figures
figM = plt.figure(figsize=(8., 9.5), facecolor='w', edgecolor='k')
sub1 = plt.subplot(2,1,1)
sub2 = plt.subplot(2,1,2)

# select morphological catagries

early = [-3,2]
late = [3,20]



################### plot 1 ############################
#draw on lines for stars to dust lines
#lines = np.arange(-1, -8, -1)
#for line in lines: draw_line(sub1,line)
w__fornax = np.where(what_dust_fornax!=3)[0]
sub1.scatter(fornax_smass[w__fornax], fornax_all_dmass[w__fornax] - fornax_smass[w__fornax], marker='s', c='k', s=50, alpha=.5)


#virgo early
myplot(sub1,early,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 1,'o','red')
#myplot(sub1,early,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 2,'o','red')
#myplot(sub1,early,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 3,'v','red')

#virgo late
myplot(sub1,late,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 1,'o','blue')
#myplot(sub1,late,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 2,'o','blue')
#myplot(sub1,late,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 3,'v','blue')

#plot a squar round all fornax galaxies
#fornax early
myplot(sub1,early,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 1,'o','red')
myplot(sub1,early,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 2,'+','red')
#myplot(sub1,early,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 3,'v','red')

#fornax late
myplot(sub1,late,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 1,'o','blue')
myplot(sub1,late,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 2,'+','blue')
#myplot(sub1,late,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 3,'v','blue')


#draw line to represent min dust mass that could be detected
XXX = np.arange(1,18.,0.1)
sub1.plot(XXX,5.05-XXX,'k--')
sub1.tick_params(axis='bottom', labelbottom='off')

###############################plot 2#######################################################


#draw on lines for stars to dust lines
#lines = np.arange(-1, -8, -1)
#for line in lines: draw_line(sub1,line)
w__fornax = np.where(what_dust_fornax!=3)[0]
sub2.scatter(fornax_smass[w__fornax], fornax_all_dmass[w__fornax], marker='s', c='k', s=50, alpha=.5)

#virgo early
myplot2(sub2,early,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 1,'o','red')
#myplot2(sub2,early,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 2,'o','red')
#myplot2(sub2,early,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 3,'v','red')

#virgo late
myplot2(sub2,late,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 1,'o','blue')
#myplot2(sub2,late,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 2,'o','blue')
#myplot2(sub2,late,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass, 3,'v','blue')


#plot a squar round all fornax galaxies
#fornax early
myplot2(sub2,early,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 1,'o','red')
myplot2(sub2,early,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 2,'+','red')
#myplot2(sub2,early,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 3,'v','red')

#fornax late
myplot2(sub2,late,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 1,'o','blue')
myplot2(sub2,late,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 2,'+','blue')
#myplot2(sub2,late,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 3,'v','blue')

#myplot2(sub2,later,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass, 3,'v','green')#plot vertical line at stellar mass of 8
#sub1.axvline(x=8.5, color='k', ls='--')
sub2.axhline(y=5.05, color='k', ls='--')


#fit line to "main sequence"
x_virgo_late,y_virgo_late = selected_data_excluding([-3,1],3,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass)
#draw_curve(sub1, x_virgo_late , y_virgo_late , 1, 8.5,'blue')
#draw_curve_1(sub1,x_virgo_late,y_virgo_late,1,'blue')
#draw_path(x,y)

#fit line to "main sequence"
x_fornax_late,y_fornax_late = selected_data_excluding([-3,1],3,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass)
#compare(x_virgo_late , y_virgo_late ,x_fornax_late, y_fornax_late )
#append virgo data to fornax data for path fit
x_late = np.append(x_virgo_late, x_fornax_late)
y_late = np.append(y_virgo_late, y_fornax_late)
#draw_curve_1(sub1,x_fornax,y_fornax,1,'#00FF00')
#sub1.scatter(x,y,marker='+',s=200, alpha=0.2)
#draw_path(x,y-x)
#big_comparsion([dwarf, early, late, later], [virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass] , [fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass])
x_virgo,y_virgo = selected_data_including([0,1],3,virgo['goldmine'], virgo_all_dmass, what_dust_virgo, virgo_smass)
x_fornax,y_fornax = selected_data_including([0,1],3,fornax['goldmine'], fornax_all_dmass, what_dust_fornax, fornax_smass)
x_early = np.append(x_virgo, x_fornax)
y_early = np.append(y_virgo, y_fornax)


#min dust mass
print flux2mass(15.E-3)






print ''
print 'late pearson:' , pearsonr(x_late, y_late)
print 'early pearson:' , pearsonr(x_early, y_early)
sub1.set_xlim(6.4,11.5)
sub1.set_ylim(-5.7,-1)
sub2.set_xlim(6.4,11.5)
sub2.set_xlabel('$\log_{10}$$(M_{star} / \mathrm{M_{\odot})}$')
sub1.set_ylabel('$\log_{10}$$(M_{dust} / M_{stars})$')
sub2.set_ylabel('$\log_{10}$$(M_{dust}/  \mathrm{M_{\odot})}$')

plt.subplots_adjust(left=0.09, bottom=0.065, right=0.98, top=0.98, wspace=0.0, hspace=0.0)
figM.savefig(pj('/Users/chrisfuller/Desktop/','smass_vs_dmass_2.pdf'))
plt.show()


# 250um to dust mass
# Chris Fuller, Sep 2013

#import stuff
import numpy as np
from atpy import Table
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/plots/upper_limits/"
virgo  = Table(pj(folder,"virgo.fits" ))
fornax = Table(pj(folder,"fornax.fits"))


# # # # # # # # # # # # # # # # # # # # Functions # # # # # # # # # # # # # # # # # # # # # # # # # 
#fit 250 to dust mass
def produce_fit(sub,x_un,y_un,colour):
	x = np.nan_to_num(x_un)
	y = np.nan_to_num(y_un)

	w = np.where((x>0.01) & (y>0.01))[0]
	#print 'len where ', len(w), 'len x ', len(x)
	#print 'min ', np.min([w]), ' after log min ', np.min(np.log10(x[w]))
	x = np.log10(x)
	fit = np.polyfit(x[w],y[w],1)

	p = np.poly1d(fit)
	#now create plot
	xx = np.arange(-2.0, 3.0)
	sub.plot(xx,p(xx), ls='--', c=colour)
	sub.plot(x[w],y[w],'+', c=colour)

	#print 'fit ', fit


def line_only(sub, x1,y1,x2,y2):
	x1 = np.nan_to_num(x1)
	y1 = np.nan_to_num(y1)
	x2 = np.nan_to_num(x2)
	y2 = np.nan_to_num(y2)

	x = np.append(x1,x2)
	y = np.append(y1,y2)

	w = np.where((x>0.01) & (y>0.01))[0]
	
	x = np.log10(x)
	fit = np.polyfit(x[w],y[w],1)

	p = np.poly1d(fit)
	#now create plot
	xx = np.arange(-2., 2.,0.0001)
	sub.plot(xx,p(xx), ls='-', c='k')
	#sub.plot(x,y,'x', c='k')


	print 'fit ', fit
# # # # # # # # # # # # # # # # # # # # Main Program # # # # # # # # # # # # # # # # # # # # # # # 
#create fig
fig= plt.figure(figsize = (8.,8.),facecolor='w',edgecolor='w')

#sub fig
p1 = plt.subplot(111)
#p2 = plt.subplot(212)

#filter virgo galaxies so we only include galaxies that are a 17MPC

D = virgo.DIST__MPC_**2 / 17.0**2
print D


produce_fit(p1,virgo.F250*D,virgo.DMASS,'blue')
produce_fit(p1,fornax.F250,fornax.DMASS,'red')
line_only(p1, virgo.F250*D,virgo.DMASS,fornax.F250,fornax.DMASS)

p1.set_xlabel('$Log_{10}(S_{250})$ (Jy)')
p1.set_ylabel('$Log_{10}(M_{dust}/M_{\odot})$')
p1.set_ylim(4,9)
p1.set_xlim(-2.5,2.5)	
# old fit  [ 0.75986014  6.60244077]
# new fit  [ 0.78038038  6.47259655]
# new fit  [ 0.78937274  6.48619262]
#produce_fit(p2,virgo.F250,virgo.TEMP,'blue')
#produce_fit(p2,fornax.F250,fornax.T,'red')
#line_only(p2, virgo.F250,virgo.TEMP,fornax.F250,fornax.T)




#plt.subplots_adjust(left=0.12, bottom=0.12, right=0.95, top=0.95, wspace=0.0, hspace=0.0)
plt.subplots_adjust(bottom=0.07, left=0.08, right=0.98, top=0.98)

fig.savefig(pj('/Users/chrisfuller/Dropbox/','250_to_dmass.pdf'))

plt.show()



# reads in sevral samples and compares whether there is a stastical difference 
# in dust mass between different samples and morphologies
# Chris Fuller, May 2013

import numpy as np
from scipy.stats import ks_2samp, mannwhitneyu, f_oneway
from os import chdir
from atpy import Table
import matplotlib.pyplot as plt

#root folder
folder = '/Users/chrisfuller/Dropbox/phd/plots/statdiff/'
chdir(folder)

fornax = Table('fornax_stellar_updated.fits', type='fits')
virgo  = Table('virgo_stellar.fits' , type='fits')


#get sep into early and late
#
vType = virgo.TYPE
vE = np.where((vType<3) & (vType >= -1))[0]
vL = np.where(virgo.TYPE>=3)[0]

#fornax
fE = np.where(fornax.Morphology == "E")[0]
fL = np.where(fornax.Morphology == "L")[0]


########### functions ##############

#print line peformes various statisial tests and then prints them out to a line in 
#the terminal
def pl(n1,n2,x,y):
    u1 , u2 = np.mean(x) , np.mean(y) #caculate means
    s1 , s2 = np.std(x) / np.sqrt(len(x))  , np.std(y) / np.sqrt(len(y))  #caculate error
    
    #ks test
    ks = ks_2samp(x,y)
    
    #mwu test
    mwu = mannwhitneyu(x,y)
    
    #FTEST?
    f = f_oneway(x,y)

    #N
    N1 = len(x)
    N2 = len(y)
    
    #result = n1 +r(u1,2),s1,u2,s2,ks[0],ks[1],mwu[0],mwu[1],f[0],f[1]
    a = n1 + '&'+ n2 + '&' + r(u1,2) +'('+r(s1,2)+')'+ '&' + r(u2,2) +'('+r(s2,2)+')' + '&'
    b = r(ks[0],3) +'&' + r(ks[1],3) + '&' + str(N1) + '&'  + str(N2) + '\\\\'#+ '&' + r(mwu[0],5) +'&' + r(mwu[1],5) + '&' + r(f[0],5) +'&' + r(f[1],5) +'\\\\'
    print a + b
def r(x,d):
    return str(np.round(x,decimals=d))
    
def linefit(x,y,l):
    fit = np.polyfit(x,y,1.0)
    err_m, err_c = errorInLine(x,y)
    dd = 2
    print l +"&" + r(fit[0],dd) + "\,$\pm$\," + r(err_m,dd) + "&" + r(fit[1],dd) + "$\,\pm\,$" + r(err_c,dd)  + "\\\\"    
    p = np.poly1d(fit)
    return x , p(x)

def myplot(x,y,l,mark, mark2):
        xl,yl = linefit(x,y,l)
        f1.plot(x,y, mark ,label=l)
        f1.plot(xl,yl,mark2, label=l)
        
        
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
    
    return np.sqrt(dm_squared[0]),np.sqrt(dc_squared[0])
############ control ###############

#mass
vMass = virgo.mass
fMass = fornax.MASS

#dust mass / stellar mass
vd2s =   virgo.M_STAR - vMass
fd2s =   fornax.SMASS - fMass

#temp
vTemp = virgo.temp
fTemp = fornax.T

#mag
#vMag = virgo.BT
#fMag = fornax.BTmag


print '\\begin{table*}'
print '\\centering'
print '\\begin{tabular}{cccccc}'#cccc}'
print '\\hline'
print '\\\\'

a = 'Sample 1 & Sample 2 & $\mu_{1}$($\sigma_{1}$) & $\mu_{2}$($\sigma_{2}$) &' 
b = '\multicolumn{2}{c}{K-S test} & \multicolumn{2}{c}{Sample Size} \\\\' #& \multicolumn{2}{c}{M-W U-test} & \multicolumn{2}{c}{F-test}\\\\'
print a+b
print '~&~&~&~&Value& $P_{value}$& $N_{1}$& $N_{2}$' #&Value& $P_{value}$&Value& $P_{value}$ \\\\ \\hline'
print '\\\\ \\hline'
print'\multicolumn{6}{l}{\\textbf{Dust Mass (log($M_{Dust}$/$M_{\odot}$))}} \\\\'


pl("Virgo Early","Virgo Late",vMass[vE],vMass[vL])
pl("Fornax Early","Fornax Late",fMass[fE],fMass[fL])
pl("Virgo Early","Fornax Early",vMass[vE],fMass[fE])
pl("Virgo Late","Fornax Late",vMass[vL],fMass[fL])

print '\\\\'
print'\multicolumn{6}{l}{\\textbf{ Stellar Mass / Dust Mass (log($M_{Stellar}$/$M_{\odot}$) - log($M_{Dust}$/$M_{\odot}$))}} \\\\'
pl("Virgo Early","Virgo Late",vd2s[vE],vd2s[vL])
pl("Fornax Early","Fornax Late",fd2s[fE],fd2s[fL])
pl("Virgo Early","Fornax Early",vd2s[vE],fd2s[fE])
pl("Virgo Late","Fornax Late",vd2s[vL],fd2s[fL])


print '\\\\'
print'\multicolumn{6}{l}{\\textbf{Dust Temp. (K)}} \\\\'
pl("Virgo Early","Virgo Late",vTemp[vE],vTemp[vL])
pl("Fornax Early","Fornax Late",fTemp[fE],fTemp[fL])
pl("Virgo Early","Fornax Early",vTemp[vE],fTemp[fE])
pl("Virgo Late","Fornax Late",vTemp[vL],fTemp[fL])


print '\\hline'
print '\\end{tabular}'
print '\\end{table*}'







print ''
print '####################################'
print ''

print 'fornax early range mass:', str(np.min(fMass[fE])) ," - ", str(np.max(fMass[fE]))
print 'fornax late  range mass:', str(np.min(fMass[fL])) ," - ", str(np.max(fMass[fL]))

print ''
print 'fornax early range temp:', str(np.min(fTemp[fE])) ," - ", str(np.max(fTemp[fE]))
print 'fornax late  range temp:', str(np.min(fTemp[fL])) ," - ", str(np.max(fTemp[fL]))
print ''
print '####################################'









print '\\begin{table}'
print '\\centering'
print '\\begin{tabular}{ccc}'#cccc}'
print '\\hline'
print '~&Gradient, M& Intercept, C \\\\ \\hline'


fig = plt.figure(figsize = (4.5,4.5),facecolor='w',edgecolor='w')
f1 = plt.subplot(1,1,1)

c= +3.0
YY = np.arange(4.5,9.0)
XX = np.arange(4.5+c,9.0+c)
f1.plot(XX,YY, 'k-', linewidth=1, alpha=0.7)

myplot(fornax.SMASS[fE],fMass[fE],'Fornax Early','or','r-')
myplot(fornax.SMASS[fL],fMass[fL],'Fornax Late', 'ob','b--')

myplot(virgo.M_STAR[vE],vMass[vE],'Virgo Early','xr','r--.')
myplot(virgo.M_STAR[vL],vMass[vL],'Virgo Late','xb','b-.')


f1.set_ylabel('Dust Mass, log($M_{Dust}$/$M_{\odot}$)')
f1.set_xlabel('Stellar Mass, log($M_{Star}$/$M_{\odot}$)')
plt.subplots_adjust(left=0.15, bottom=None, right=None, top=None, wspace=None, hspace=None)
#fig.savefig('/Users/chrisfuller/Dropbox/phd/papers/fornax/stellardustmass.pdf')
#plt.show()


'''


print '\\hline'
print '\\end{tabular}'
print '\\end{table}'

print ''
print '####################################'
print ''
print ''







#mass and temp
bm=np.arange(3,9,.2)
print bm
fig = plt.figure(figsize = (4,8),facecolor='w',edgecolor='w')

f1 = plt.subplot(2,1,1)
f2 = plt.subplot(2,1,2)

#mass 
print fMass[fE]
print fMass[fL]
f1.hist(fMass[fE],bins=bm,histtype='step',color='red')#, normed=True)
f1.hist(fMass[fL],bins=bm,histtype='step',color='blue')#, normed=True)
#f1.hist(vd2s[vE],bins=bm,histtype='step',color='purple')#, normed=True)
#f1.hist(vd2s[vL],bins=bm,histtype='step',color='green')#, normed=True)
f1.axvline(np.mean(fMass[fE]),color='red',ls='--')
f1.axvline(np.mean(fMass[fL]),color='blue',ls='--')
#f1.axvline(np.mean(vMass[vE]),color='red',ls='--')
#f1.axvline(np.mean(vMass[vL]),color='blue',ls='--')
f1.set_xlim(4,10)
#f1.set_ylim(0,5)
f1.set_xlabel('log($M_{Dust}$/$M_{\odot}$)' )
f1.set_ylabel('N')


bt = np.arange(6,26,3.0)
#temp
f2.hist(fTemp[fE],bins=bt,histtype='step',color='red')#, normed=True)
f2.hist(fTemp[fL],bins=bt,histtype='step',color='blue')#, normed=True)
#f2.hist(vTemp[vE],bins=bt,histtype='step',color='purple')#, normed=True)
#f2.hist(vTemp[vL],bins=bt,histtype='step',color='green')#, normed=True)
f2.axvline(np.mean(fTemp[fE]),color='red',ls='--')
f2.axvline(np.mean(fTemp[fL]),color='blue',ls='--')
#f2.axvline(np.mean(vTemp[fE]),color='red',ls='--')
#f2.axvline(np.mean(vTemp[fL]),color='blue',ls='--')
f2.set_xlim(6,26)
f2.set_ylim(0,4)
f2.set_xlabel('Dust Temp. (K))' )
f2.set_ylabel('N')
ya = f2.get_yaxis()
ya.set_major_locator(plt.MaxNLocator(integer=True))

plt.subplots_adjust(hspace=0.3,wspace=0.0)
plt.show()
'''


# reads in sevral samples and compares whether there is a stastical difference 
# in dust mass between different samples and morphologies
# Chris Fuller, May 2013

import numpy as np
from scipy.stats import ks_2samp, mannwhitneyu, f_oneway
from os import chdir
from atpy import Table
import matplotlib.pyplot as plt
from os.path import join as pj

#Filament = Table('Filament_stellar_updated.fits', type='fits')
#Cluster  = Table('Cluster_stellar.fits' , type='fits')

#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

subCat = cat.where(cat.DMASS_TYPE == 2)

cluster = subCat.where(subCat.RADIUS_VIR <= 1.0)
filament = subCat.where(subCat.RADIUS_VIR > 1.0)


#get sep into early and late
#cluster
cE = np.where(cluster.early == 1)[0]
cI = np.where(cluster.inter == 1)[0]
cL = np.where(cluster.late == 1)[0]

#filament
fE = np.where(filament.early == 1)[0]
fI = np.where(filament.inter == 1)[0]
fL = np.where(filament.late == 1)[0]


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
    
    #result = n1 +r(u1,2),s1,u2,s2,ks[0],ks[1],mwu[0],mwu[1],f[0],f[1]
    a = n1 + '&'+ n2 + '&' + r(u1,2) +'('+r(s1,2)+')'+ '&' + r(u2,2) +'('+r(s2,2)+')' + '&'
    b = r(ks[0],3) +'&' + r(ks[1],3) + '\\\\'#+ '&' + r(mwu[0],5) +'&' + r(mwu[1],5) + '&' + r(f[0],5) +'&' + r(f[1],5) +'\\\\'
    print a + b

def pl1(n1,n2,x,y):
    u1 , u2 = np.mean(x) , np.mean(y) #caculate means
    s1 , s2 = np.std(x) / np.sqrt(len(x))  , np.std(y) / np.sqrt(len(y))  #caculate error

    diff = abs(u1-u2)/ np.sqrt(s1**2 +s2**2)
    
    #ks test
    ks = ks_2samp(x,y)
    
    #mwu test
    mwu = mannwhitneyu(x,y)
    
    #FTEST?
    f = f_oneway(x,y)
    
    #result = n1 +r(u1,2),s1,u2,s2,ks[0],ks[1],mwu[0],mwu[1],f[0],f[1]
    a = n1 + '&' + r(u1,2) +'('+r(s1,2)+')'+ '&' + r(u2,2) +'('+r(s2,2)+')' + '&' + r(diff,1)+ '\\sigma' + '&'
    b = r(ks[0],3) +'&' + r(ks[1],3) + '\\\\'#+ '&' + r(mwu[0],5) +'&' + r(mwu[1],5) + '&' + r(f[0],5) +'&' + r(f[1],5) +'\\\\'
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
cMass = cluster.SMASS
fMass = filament.SMASS

#dust mass / stellar mass
cd2s =   cluster.DUST_STARS
fd2s =   filament.DUST_STARS

#temp
cTemp = cluster.DTEMP
fTemp = filament.DTEMP

#mag
#vMag = Cluster.BT
#fMag = Filament.BTmag


print '\\begin{table*}'
print '\\centering'
print '\\begin{tabular}{cccccc}'#cccc}'
print '\\hline'
print '\\\\'

a = 'Sample 1 & Sample 2 & $\mu_{1}$($\sigma_{1}$) & $\mu_{2}$($\sigma_{2}$) &' 
b = '\multicolumn{2}{c}{K-S test} \\\\' #& \multicolumn{2}{c}{M-W U-test} & \multicolumn{2}{c}{F-test}\\\\'
print a+b
print '~&~&~&~&Value& $P_{value}$' #&Value& $P_{value}$&Value& $P_{value}$ \\\\ \\hline'
print '\\\\ \\hline'
print'\multicolumn{6}{l}{\\textbf{Dust Mass (log($M_{Dust}$/$M_{\odot}$))}} \\\\'


#pl("Cluster Early","Cluster Late",cMass[cE],cMass[cL])
#pl("Filament Early","Filament Late",fMass[fE],fMass[fL])
pl1("Cluster Early","Filament Early",cMass[cE],fMass[fE])
pl1("Cluster Intermediate","Filament Intermediate",cMass[cI],fMass[fI])
pl1("Cluster Late","Filament Late",cMass[cL],fMass[fL])

print '\\\\'
print'\multicolumn{6}{l}{\\textbf{ Stellar Mass / Dust Mass (log($M_{Stellar}$/$M_{\odot}$) - log($M_{Dust}$/$M_{\odot}$))}} \\\\'
#pl1("Cluster Early","Cluster Late",cd2s[cE],cd2s[cL])
#pl1("Filament Early","Filament Late",fd2s[fE],fd2s[fL])
pl1("Cluster Early","Filament Early",cd2s[cE],fd2s[fE])
pl1("Cluster Intermediate","Filament Intermediate",cd2s[cI],fd2s[fI])
pl1("Cluster Late","Filament Late",cd2s[cL],fd2s[fL])


print '\\\\'
print'\multicolumn{6}{l}{\\textbf{Dust Temp. (K)}} \\\\'
#pl1("Cluster Early","Cluster Late",cTemp[cE],cTemp[cL])
#pl1("Filament Early","Filament Late",fTemp[fE],fTemp[fL])
pl1("Cluster Early","Filament Early",cTemp[cE],fTemp[fE])
pl1("Cluster Intermediate","Filament Intermediate",cTemp[cI],fTemp[fI])
pl1("Cluster Late","Filament Late",cTemp[cL],fTemp[fL])


print '\\hline'
print '\\end{tabular}'
print '\\end{table*}'







print ''
print '####################################'
print ''

print 'Filament early range mass:', str(np.min(fMass[fE])) ," - ", str(np.max(fMass[fE]))
print 'Filament late  range mass:', str(np.min(fMass[fL])) ," - ", str(np.max(fMass[fL]))

print ''
print 'Filament early range temp:', str(np.min(fTemp[fE])) ," - ", str(np.max(fTemp[fE]))
print 'Filament late  range temp:', str(np.min(fTemp[fL])) ," - ", str(np.max(fTemp[fL]))
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

myplot(Filament.SMASS[fE],fMass[fE],'Filament Early','or','r-')
myplot(Filament.SMASS[fL],fMass[fL],'Filament Late', 'ob','b--')

myplot(Cluster.M_STAR[cE],cMass[cE],'Cluster Early','xr','r--.')
myplot(Cluster.M_STAR[cL],cMass[cL],'Cluster Late','xb','b-.')


f1.set_ylabel('Dust Mass, log($M_{Dust}$/$M_{\odot}$)')
f1.set_xlabel('Stellar Mass, log($M_{Star}$/$M_{\odot}$)')
plt.subplots_adjust(left=0.15, bottom=None, right=None, top=None, wspace=None, hspace=None)
fig.sacefig('/Users/chrisfuller/Dropbox/phd/papers/Filament/stellardustmass.pdf')
plt.show()


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
#f1.hist(cd2s[cE],bins=bm,histtype='step',color='purple')#, normed=True)
#f1.hist(cd2s[cL],bins=bm,histtype='step',color='green')#, normed=True)
f1.axcline(np.mean(fMass[fE]),color='red',ls='--')
f1.axcline(np.mean(fMass[fL]),color='blue',ls='--')
#f1.axcline(np.mean(cMass[cE]),color='red',ls='--')
#f1.axcline(np.mean(cMass[cL]),color='blue',ls='--')
f1.set_xlim(4,10)
#f1.set_ylim(0,5)
f1.set_xlabel('log($M_{Dust}$/$M_{\odot}$)' )
f1.set_ylabel('N')


bt = np.arange(6,26,3.0)
#temp
f2.hist(fTemp[fE],bins=bt,histtype='step',color='red')#, normed=True)
f2.hist(fTemp[fL],bins=bt,histtype='step',color='blue')#, normed=True)
#f2.hist(cTemp[cE],bins=bt,histtype='step',color='purple')#, normed=True)
#f2.hist(cTemp[cL],bins=bt,histtype='step',color='green')#, normed=True)
f2.axcline(np.mean(fTemp[fE]),color='red',ls='--')
f2.axcline(np.mean(fTemp[fL]),color='blue',ls='--')
#f2.axcline(np.mean(cTemp[fE]),color='red',ls='--')
#f2.axcline(np.mean(cTemp[fL]),color='blue',ls='--')
f2.set_xlim(6,26)
f2.set_ylim(0,4)
f2.set_xlabel('Dust Temp. (K))' )
f2.set_ylabel('N')
ya = f2.get_yaxis()
ya.set_major_locator(plt.MaxNLocator(integer=True))

plt.subplots_adjust(hspace=0.3,wspace=0.0)
plt.show()
'''


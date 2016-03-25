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
foldercoma = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
cat = Table(pj(foldercoma,'coma_supercluster_cal12_pacscorrected.fits'))
subCat = cat.where(cat.DMASS_TYPE == 2)
cluster = subCat.where(subCat.RADIUS_VIR <= 1.0)
filament = subCat.where(subCat.RADIUS_VIR > 1.0)

#root folder
folder = '/Users/chrisfuller/Dropbox/phd/plots/statdiff/'
chdir(folder)

fornax = Table('fornax_stellar_updated.fits', type='fits')
virgo  = Table('virgo_stellar.fits' , type='fits')
hrsField = Table('/Users/chrisfuller/Dropbox/phd/herchel/HRS/HRS_field_inc_outerVirgo.fits')

#get sep into early and late
hrsType = hrsField.Type
hrsE = np.where((hrsType<3) & (hrsType >= -1))[0]
hrsL = np.where(hrsType >= 3)[0]


#virgo
vType = virgo.TYPE
virgoE = np.where((vType<3) & (vType >= -1))[0]
virgoL = np.where(virgo.TYPE>=3)[0]

#fornax
fornaxE = np.where(fornax.Morphology == "E")[0]
fornaxL = np.where(fornax.Morphology == "L")[0]


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
    b = r(ks[0],3) +'&' + r(ks[1],3) + '&' + str(N1) + '&'  + str(N2) + '\\\\'#+ '&' + r(mwu[0],5) +'&' + r(mwu[1],5) + '&' + r(f[0],5) +'&' + r(f[1],5) +'\\\\'
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

    #N
    N1 = len(x)
    N2 = len(y)
    
    #result = n1 +r(u1,2),s1,u2,s2,ks[0],ks[1],mwu[0],mwu[1],f[0],f[1]
    a = n1 + '&' + n2 + '&' + r(u1,2) +'('+r(s1,2)+')'+ '&' + r(u2,2) +'('+r(s2,2)+')' + '&' #+ r(diff,1)+ '\\sigma' + '&'
    b = r(ks[0],3) +'&' + r(ks[1],3) + '&' + str(N1) + '&'  + str(N2) + '\\\\'#+ '&' + r(mwu[0],5) +'&' + r(mwu[1],5) + '&' + r(f[0],5) +'&' + r(f[1],5) +'\\\\'
    print a + b




def r(x,d):
    return str(np.round(x,decimals=d))


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
cMass = cluster.DMASS_250 #cluster.DMASS
fMass = filament.DMASS_250 #filament.DMASS
vMass = virgo.mass
forMass = fornax.MASS
hrsMass = hrsField.DMASS

#dust mass / stellar mass
cd2s =   cMass - cluster.SMASS#cluster.DUST_STARS
fd2s =   fMass - filament.SMASS#filament.DUST_STARS

vd2s =    vMass - virgo.M_STAR
ford2s =   forMass- fornax.SMASS
hrsd2s = hrsMass - hrsField.SMASS.astype(np.float) + 0.36


#temp
cTemp = cluster.DTEMP
fTemp = filament.DTEMP
vTemp = virgo.temp
forTemp = fornax.T
hrsTemp = hrsField.TEMP

#mag
#vMag = Cluster.BT
#fMag = Filament.BTmag


print '\\begin{table}'
print '\\centering'
print '\\begin{tabular}{cccccc}'#cccc}'
print '\\hline'
print '\\\\'

a = 'Sample 1 & Sample 2 & $\mu_{1}$($\sigma_{1}$) & $\mu_{2}$($\sigma_{2}$) & $\\frac{\mu_{2} - \mu_{2}}{(\sigma^{2}_{1} + \sigma^{2}_{2})^{1/2}}$&' 
b = '\multicolumn{2}{c}{K-S test} \\\\' #& \multicolumn{2}{c}{M-W U-test} & \multicolumn{2}{c}{F-test}\\\\'
print a+b
print '~&~&~&~&Value& $P_{value}$' #&Value& $P_{value}$&Value& $P_{value}$ \\\\ \\hline'
print '\\\\ \\hline'
print'\multicolumn{6}{l}{\\textbf{Dust Mass (log($M_{Dust}$/$M_{\odot}$))}} \\\\'

"""
pl1("Coma Early","Filament Early",cMass[cE],fMass[fE])
pl1("Coma Uncertain","Filament Uncertain",cMass[cI],fMass[fI])
pl1("Coma Late","Filament Late",cMass[cL],fMass[fL])
print '\\\\'
pl1("Coma Early","Virgo Early",cMass[cE],vMass[virgoE])
pl1("Coma Early","Fornax Early",cMass[cE],vMass[fornaxE])
pl1("Coma Late","Virgo Late",cMass[cL],vMass[virgoL])
pl1("Coma Late","Fornax Late",cMass[cL],vMass[fornaxL])
"""

print '\\\\'
print'\multicolumn{6}{l}{\\textbf{ Stellar Mass / Dust Mass (log($M_{Stellar}$/$M_{\odot}$) - log($M_{Dust}$/$M_{\odot}$))}} \\\\'
#pl1("Coma Early","Coma Late",cd2s[cE],cd2s[cL])
#pl1("Filament Early","Filament Late",fd2s[fE],fd2s[fL])
#pl1("Coma Early","Filament Early",cd2s[cE],fd2s[fE])
#pl1("Coma Uncertain","Filament Uncertain",cd2s[cI],fd2s[fI])
#pl1("Coma Late","Filament Late",cd2s[cL],fd2s[fL])
print '\\\\'
pl1("Coma Early","Virgo Early",cd2s[cE],vd2s[virgoE])
pl1("Coma Early","Fornax Early",cd2s[cE],ford2s[fornaxE])
print '\\\\'
pl1("Coma Late","Virgo Late",cd2s[cL],vd2s[virgoL])
pl1("Coma Late","Fornax Late",cd2s[cL],ford2s[fornaxL])
#pl1("Coma Late","HRS Late",cd2s[cL],hrsd2s[hrsL])
#pl1("Coma Early","HRS Early",cd2s[cE],hrsd2s[hrsE])
print '\\\\'
pl1("Filament Late","HRS Late",fd2s[fL],hrsd2s[hrsL])
pl1("Filament Early","HRS Early",fd2s[fE],hrsd2s[hrsE])

"""
print '\\\\'
print'\multicolumn{6}{l}{\\textbf{Dust Temp. (K)}} \\\\'
pl1("Coma Early","Coma Late",cTemp[cE],cTemp[cL])
pl1("Filament Early","Filament Late",fTemp[fE],fTemp[fL])
pl1("Coma Early","Filament Early",cTemp[cE],fTemp[fE])
pl1("Coma Uncertain","Filament Uncertain",cTemp[cI],fTemp[fI])
pl1("Coma Late","Filament Late",cTemp[cL],fTemp[fL])
pl1("Coma Early","Virgo Early",cTemp[cE],vTemp[virgoE])
pl1("Coma Early","Fornax Early",cTemp[cE],forTemp[fornaxE])
pl1("Coma Late","Virgo Late",cTemp[cL],vTemp[virgoL])
pl1("Coma Late","Fornax Late",cTemp[cL],forTemp[fornaxL])

pl1("Coma Late","HRS Late",cTemp[cL],hrsTemp[hrsL])
pl1("Coma Early","HRS Early",cTemp[cE],hrsTemp[hrsE])
pl1("Filament Late","HRS Late",fTemp[cL],hrsTemp[hrsL])
pl1("Filament Early","HRS Early",fTemp[cE],hrsTemp[hrsE])
"""
print '\\hline'
print '\\end{tabular}'
print '\\end{table}'







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









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
fname = 'coma_supercluster_cal12_pacscorrected.fits' #input name
cat = Table(pj(folder,fname))

subCat = cat.where(cat.DMASS_TYPE == 2)


cluster = subCat.where(subCat.RADIUS_VIR <= 1.0)
filament = subCat.where(subCat.RADIUS_VIR > 1.0)

#gross sep
EE = np.where(subCat.early == 1)[0]
II = np.where(subCat.inter == 1)[0]
LL = np.where(subCat.late == 1)[0]


#get sep into early and late
#cluster
cE = np.where(cluster.early == 1)[0]
cI = np.where(cluster.inter == 1)[0]
cL = np.where(cluster.late == 1)[0]

#filament
fE = np.where(filament.early == 1)[0]
fI = np.where(filament.inter == 1)[0]
fL = np.where(filament.late == 1)[0]


#hi def def
hi_nor = np.where(subCat.DefHI >= 0.5)[0]
hi_def = np.where(subCat.DefHI < 0.5)[0]

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
cMass = cluster.DMASS
fMass = filament.DMASS
Mass = subCat.DMASS

#dust mass / stellar mass
cd2s =   cluster.DUST_STARS
fd2s =   filament.DUST_STARS
d2s = subCat.DUST_STARS

#temp
cTemp = cluster.DTEMP
fTemp = filament.DTEMP
Temp = subCat.DTEMP

#mag
#vMag = Cluster.BT
#fMag = Filament.BTmag


print '\\begin{table}'
print '\\centering'
print '\\begin{tabular}{cccccc}'#cccc}'
print '\\hline'
print '\\\\'

a = 'Sample 1 & Sample 2 & $\mu_{1}$($\sigma_{1}$) & $\mu_{2}$($\sigma_{2}$) & $\\frac{\mu_{2} - \mu_{2}}{(\sigma^{2}_{1} + \sigma^{2}_{2})^{1/2}}$&' 
b = '\multicolumn{2}{c}{K-S test} & \multicolumn{2}{c}{Sample Size} \\\\' #& \multicolumn{2}{c}{M-W U-test} & \multicolumn{2}{c}{F-test}\\\\'
print a+b
print '~&~&~&~&Value& $P_{value}$& $N_{1}$& $N_{2}$' #&Value& $P_{value}$&Value& $P_{value}$ \\\\ \\hline'
print '\\\\ \\hline'
print'\multicolumn{6}{l}{\\textbf{Dust Mass (log($M_{Dust}$/$M_{\odot}$))}} \\\\'



#pl("Cluster Early","Cluster Late",cMass[cE],cMass[cL])
#pl("Filament Early","Filament Late",fMass[fE],fMass[fL])
pl1("Cluster Early","Filament Early",cMass[cE],fMass[fE])
pl1("Cluster Uncertain","Filament Uncertain",cMass[cI],fMass[fI])
pl1("Cluster Late","Filament Late",cMass[cL],fMass[fL])
pl1("HI-Normal","HI-Deficent",Mass[hi_nor],Mass[hi_def])

print '\\\\'
print'\multicolumn{6}{l}{\\textbf{ Stellar Mass / Dust Mass (log($M_{Stellar}$/$M_{\odot}$) - log($M_{Dust}$/$M_{\odot}$))}} \\\\'
#pl1("Cluster Early","Cluster Late",cd2s[cE],cd2s[cL])
#pl1("Filament Early","Filament Late",fd2s[fE],fd2s[fL])
pl1("Cluster Early","Filament Early",cd2s[cE],fd2s[fE])
pl1("Cluster Uncertain","Filament Uncertain",cd2s[cI],fd2s[fI])
pl1("Cluster Late","Filament Late",cd2s[cL],fd2s[fL])
pl1("HI-Normal","HI-Deficent",d2s[hi_nor],d2s[hi_def])

print '\\\\'
print'\multicolumn{6}{l}{\\textbf{Dust Temp. (K)}} \\\\'
#pl1("Cluster Early","Cluster Late",cTemp[cE],cTemp[cL])
#pl1("Filament Early","Filament Late",fTemp[fE],fTemp[fL])
pl1("Cluster Early","Filament Early",cTemp[cE],fTemp[fE])
pl1("Cluster Uncertain","Filament Uncertain",cTemp[cI],fTemp[fI])
pl1("Cluster Late","Filament Late",cTemp[cL],fTemp[fL])
pl1("HI-Normal","HI-Deficent",Temp[hi_nor],Temp[hi_def])

print '\\hline'
print '\\end{tabular}'
print '\\end{table}'














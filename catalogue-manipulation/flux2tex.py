# Taking a table of fluxes and printing them into a latex table
# Chris Fuller, May 2013

import numpy as np
from atpy import Table

folder = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/'
cat = Table(folder+"HeFoCS-fluxes-260313-mybgsub-v2-handmeasure-bad-detections-removed-morphology.fits",type='fits')
start = 0
n = 30
n2 = 4
bands = ['500','350','250','160','100']
d=3
######### function ############
def pl(l):
    s = ''
    i = 0 
    for val in l:
        if i != len(l)-1: x = val + ' &'
        else: x = val + ' \\\\' 
        i +=1
        s = s+x
    print s
    
def dtest(flux):
    flux = np.nan_to_num(flux)
    if float(flux) > 0.0: return 1
    else: return 0
    
def r(x,d):
    return str(np.round(x,decimals=d))
######## control ##############




h1 = ['OBJECT','RA','Dec.','$S_{500} (\sigma S_{500})$','$S_{350} (\sigma S_{350})$','$S_{250} (\sigma S_{250})$','$S_{160} (\sigma S_{160})$','$S_{100} (\sigma S_{100})$']
h2 = ['~','h:m:s','d:m:s','(Jy)','(Jy)','(Jy)','(Jy)','(Jy)']
h3 = ['~','(J2000)','(J2000)','~','~','~','~','~']

space = ['~','.','.','.','.','.','.','.']

headers = [h1,h2,h3]

print '\\begin{table*}'
print '\\centering'
print '\\begin{tabular}{'+'c'*len(h1)+'}'
print '\\hline'

for header in headers: pl(header) #print header

print '\\hline'

for i in range(start,n):
    line = ['FCC '+str(cat.OBJECT[i]),str(cat.RA[i]),str(cat.DEC[i])]
    for band in bands:
        if dtest(cat["F"+band][i])==1: y = r(cat["F"+band][i],d) +" (" + r(cat["E"+band][i],d) + ")"
        
        elif dtest(cat["F"+band][i])== 0 and (band[0] =='1') and cat.HEVICS_PLW[i] ==3 : y = '-'
        
        else: y = '\\textless(' + r(cat["E"+band][i],d) + ')'
        
        if int(cat['remove-'+band][i]) == 1: y = y+"*" 
        
        line.append(y)
    pl(line)
  
for j in range(0,n2): pl(space)   

print '\\hline'
print '\\end{tabular}'
print '\\end{table*}'        
        
    
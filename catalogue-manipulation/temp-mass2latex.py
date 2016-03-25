# Taking a table of fluxes and printing them into a latex table
# Chris Fuller, May 2013

import numpy as np
from atpy import Table

folder = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/'
cat = Table(folder+"sed.fits",type='fits')
d=2
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
    
def val(x,e,d):    
    return r(x,d) + ' ('+r(e,d)+')' 
######## control ##############




h1 = ['OBJECT','RA','Dec.','Type','Dust Temperature','Duss Mass', 'Stellar Mass']
h2 = ['~','h:m:s','d:m:s','~','K','log($M_{Dust}$/M$_{\odot}$)', 'log($M_{Stars}$/M$_{\odot}$)']
h3 = ['~','(J2000)','(J2000)','~','~','~','~']

space = ['~','.','.','.','.']

headers = [h1,h2,h3]

print '\\begin{table*}'
print '\\centering'
print '\\begin{tabular}{'+'c'*len(h1)+'}'
print '\\hline'

for header in headers: pl(header) #print header

print '\\hline'

for i in range(0,len(cat)):
    line = ['FCC '+str(cat.OBJECT[i]),str(cat.RA[i]),str(cat.DEC[i]),str(cat.TYPE[i]).rstrip()]
    m , me = cat.MASS[i], cat.MERROR[i]
    t , te = cat.T[i] , cat.TERROR[i]
    s = cat.SMASS[i]
    line.append(val(t,te,d))
    line.append(val(m,me,d))
    line.append(r(s,d))
    pl(line)
   

print '\\hline'
print '\\end{tabular}'
print '\\end{table*}'        
        
    
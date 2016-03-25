#Program to create cataloge tables in laytex
"""
\begin{center}
\begin{longtable}{|l|l|l|}
\caption[Feasible triples for a highly variable Grid]{Feasible triples for 
highly variable Grid, MLMMH.} \label{grid_mlmmh} \\

\hline \multicolumn{1}{|c|}{\textbf{Time (s)}} & \multicolumn{1}{c|}{\textbf{Triple chosen}} & \multicolumn{1}{c|}{\textbf{Other feasible triples}} \\ \hline 
\endfirsthead

\multicolumn{3}{c}%
{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
\hline \multicolumn{1}{|c}{\textbf{Time (s)}} &
\multicolumn{1}{c|}{\textbf{Triple chosen}} &
\multicolumn{1}{c|}{\textbf{Other feasible triples}} \\ \hline 
\endhead

\hline \multicolumn{3}{c}{{Continued on next page}} \\ \hline
\endfoot

\hline \hline
\endlastfoot
"""

# Taking a table of fluxes and printing them into a latex table
# Chris Fuller, May 2013

import numpy as np
from atpy import Table

#folder = '/Users/chrisfuller/Desktop/'
folder = '/Users/chrisfuller/Dropbox/phd/herchel/fornax/final_outputs/'
cat = Table(folder+"fornax_paper.fits",type='fits')
start = 0
n = len(cat)
n2 = 6
bands = ['500','350','250','160','100']
d=0
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
    if d == 0:
        return str(np.int(x))

    else:
        return str(np.round(x,decimals=d))
######## control ##############




h1 = ['OBJECT','RA','Dec.','$S_{500}$', '$E_{500}$','$S_{350}$', '$E_{350}$','$S_{250}$', '$E_{250}$','$S_{160}$', '$E_{160}$','$S_{100}$', '$E_{100}$']
#h2 = ['~','h:m:s','d:m:s','(Jy)','(Jy)','(Jy)','(Jy)','(Jy)','(Jy)','(Jy)','(Jy)','(Jy)','(Jy)']
h2 = ['~','h:m:s','d:m:s','(mJy)','(mJy)','(mJy)','(mJy)','(mJy)','(mJy)','(mJy)','(mJy)','(mJy)','(mJy)']
h3 = ['~','(J2000)','(J2000)','~','~','~','~','~','~','~','~','~','~']

space = ['.']*len(h1)

headers = [h1,h2,h3]
print '\\begin{longtable}{'+'c'*len(h1)+'}'
print '\\hline'

for header in headers: pl(header) #print header

print '\\hline'
print '\\endfirsthead'
print '\\hline'

for header in headers: pl(header) #print header

print '\\hline'
print '\\endhead'

print "\\hline \\multicolumn{3}{r}{{Continued on next page}} \\\\ \\hline"
print "\\endfoot"

print "\\hline \\hline"
print "\\endlastfoot"

for i in range(start,n):
    line = ['FCC '+str(cat.OBJECT[i]),str(cat.RA[i]).split('.')[0],str(cat.DEC[i]).split('.')[0]]
    for band in bands:
        y = r(cat["F"+band][i]*1000.0,d) +" & " + r(cat["E"+band][i]*1000.0,d) 
        
        if dtest(cat["F"+band][i])== 0 and (band[0] =='1') and cat.HEVICS_PLW[i] ==3 : y = '-&-'
        

        #else: y = '\\textless(' + r(cat["E"+band][i],d) + ')'
        
        #if int(cat['remove-'+band][i]) == 1: y = y+"*" 
        
        line.append(y)
    pl(line)
  
#for j in range(0,n2): pl(space)   

print '\\hline'
print '\\end{tabular}'
print '\\end{longtable}'        
        
    
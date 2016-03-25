#morph converter
#Chris Fuller 11th July 2013

#This program takes a catalogue and converts it into a goldmine type morphological classification

#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt
import matplotlib.pylab  as pl

#Inputs

cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12.fits")

fornax = cat.where((cat['RADIUS'] <  6048.0)& (cat.goldmine != -32768))
virgo  = cat.where((cat['RADIUS'] >= 6048.0)& (cat.goldmine != -32768))




###################### functions ######################
def w(cat,x):
    return where(cat['goldmine']==x)[0]
    
def dw(f):
    f = nan_to_num(f)
    return where((f>0.0))
    

#p(total_fornax[i],detected_fornax[i])
def p(a,b):
    per = np.round(float(b)*100.0/a, decimals = 0)
    error = np.sqrt(b) / a
    return str(np.int(per))+'\\,$\\pm$\\,'+str(np.int(np.round(error*100.0, decimals=0)))
    
def er(a,b):
    return np.sqrt(a) *100.0 / b


    
####################### main ##########################
    
#find detected
fornax_m_d = fornax['goldmine'][dw(fornax['F250'])]
virgo_m_d  =  virgo['goldmine'][dw( virgo['F250'])]

#find totals
fornax_m_total = fornax['goldmine']
virgo_m_total  =  virgo['goldmine']



################### plotting ##########################

bins_all = arange(-3,21)

bins_early_late_later_dwalf = [0,2,10,21]
bina = bins_early_late_later_dwalf


mList =  ['','Early','', 'Late','', 'Irregulars']
mList_nos =  ['E/S', 'Sa/Sb/Sc/Sd', 'Sm/Im/dS']

total_fornax,_    =  histogram(fornax_m_total, bins = bina)
detected_fornax,_ =  histogram(fornax_m_d, bins = bina)

total_virgo,_    =  histogram( virgo_m_total, bins = bina)
detected_virgo,_ =  histogram( virgo_m_d, bins = bina)

per_fornax = array(detected_fornax, dtype=float)*100. / total_fornax
per_virgo = array(detected_virgo, dtype=float)*100. / total_virgo

all_fornax = len(fornax['goldmine'])
all_virgo = len(virgo['goldmine'])

fig= plt.figure(figsize = (4.5,7.5),facecolor='w',edgecolor='w')


p1 = plt.subplot(211)
p2 = plt.subplot(212)
#p3 = plt.subplot(313)

#plot fraction detected
p1.errorbar(arange(1,len(bina)), per_fornax, yerr=er(detected_fornax,total_fornax), color='red', label='Fornax Cluster')
p1.errorbar(arange(1,len(bina)), per_virgo, yerr=er(detected_virgo,total_virgo), color='blue', label='Virgo Cluster')

#plot fraction compostion for each cluster
p2.errorbar(arange(1,len(bina)), total_fornax*100./all_fornax, yerr=er(total_fornax,all_fornax), color='red', label='Fornax Cluster', ls='--')
p2.errorbar(arange(1,len(bina)), total_virgo*100./all_virgo, yerr=er(total_virgo,all_virgo), color='blue', label='Virgo Cluster', ls='--')


p1.set_xlim(0.5,3.5)
p1.set_ylim(0,110)
#p1.set_ylabel('Percentage detected (%)')
#p1.set_xticklabels(mList)
p1.tick_params(axis='x', labelbottom='off')
p1.set_ylabel('Detected (%)')

p2.set_xlim(0.5,3.5)
p2.set_ylim(0,110)
p2.set_xlabel('Morphological Type')
#p2.tick_params(axis='x', labelbottom='off')
p2.set_ylabel('Total (%)')
p2.set_xticklabels(mList)

#p3.bar(arange(1,len(bina))-0.2, total_fornax, width=0.2, color='r')
#p3.bar(arange(1,len(bina)), total_virgo, width=0.2, color='b')
#p3.set_xlim(0.5,4.5)
#p3.set_ylim(0,340)
#p3.set_xlabel('Morphological Type')
#p3.set_xticklabels(mList)
#p1.grid()
#subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.subplots_adjust(left=0.150, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
fig.savefig(pj('/Users/chrisfuller/Desktop/','per_det_accros_hubble_seq_goldmine.pdf'))
plt.show()

"""
print 'Morphological  &  \\multicolumn{3}{c}{Virgo} & \\multicolumn{3}{c}{Fornax} \\\\'
print 'Type & Total & Detected & \\% & Total & Detected & \\%\\\\'
for i in range(0,len(mList)-1):

    print mList_nos[i],'&',total_virgo[i],'&',detected_virgo[i],'&',p(total_virgo[i],detected_virgo[i]),'&',total_fornax[i],'&',detected_fornax[i],'&',p(total_fornax[i],detected_fornax[i]),'\\\\'
"""
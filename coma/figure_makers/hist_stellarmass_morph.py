#morph converter
#Chris Fuller Sept 2014


#import moduals
from atpy import Table
from numpy import nan_to_num, where, arange, histogram, log10, array
import numpy as np
from os.path import join as pj
import matplotlib.pyplot as plt

#Inputs
cat = Table("/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/coma_supercluster_cal12.fits")


###################### functions ######################
def w(cat,x):
    return where(cat['goldmine']==x)[0]
    
def dw(f):
    f = nan_to_num(f)
    return where((f>0.0))
    

#p(total_cluster[i],detected_cluster[i])
def p(a,b):
    per = np.round(float(b)*100.0/a, decimals = 0)
    error = np.sqrt(b) / a
    return str(np.int(per))+'\\,$\\pm$\\,'+str(np.int(np.round(error*100.0, decimals=0)))
    
def er(a,b):
    return np.sqrt(a) #*100.0 / b


    
####################### main ##########################
    

early = cat.where(cat.early==1)
inter = cat.where(cat.inter==1)
late = cat.where(cat.late==1)


################### plotting ##########################




early_hist,bin_edges    =  histogram(early.SMASS, bins = 10)
inter_hist,bin_edges    =  histogram(inter.SMASS, bins = 10)
late_hist,bin_edges    =  histogram(late.SMASS, bins = 10)


#per_early = np.cumsum(array(early_hist, dtype=float)*1.0) / len(cat)
#per_late = np.cumsum(array(late_hist, dtype=float)*1.0) / len(cat)
#per_inter = np.cumsum(array(inter_hist, dtype=float)*1.0) / len(cat)

#per_early = np.log10(array(early_hist, dtype=float)*1.0 )
#per_late = np.log10(array(late_hist, dtype=float)*1.0 )
#per_inter = np.log10(array(inter_hist, dtype=float)*1.0 )


fig= plt.figure(figsize = (8.0,4.5),facecolor='w',edgecolor='w')



p1 = plt.subplot(111)

#p1.text(0.02, 0.9, texts, transform=p1.transAxes, fontsize=20, verticalalignment='top')
#plot fraction detected
p1.errorbar(bin_edges[:-1], early_hist, yerr=er(early_hist,len(cat)), color='red')
p1.errorbar(bin_edges[:-1], inter_hist, yerr=er(inter_hist,len(cat)), color='green')
p1.errorbar(bin_edges[:-1], late_hist, yerr=er(late_hist,len(cat)), color='blue')







#p1.set_xlim(-0.09,1.09)
#p1.set_ylim(0,110)
#p1.set_ylabel('Percentage detected (%)')
#p1.set_xticklabels(mList)
#p1.tick_params(axis='x', labelbottom='off')
p1.set_ylabel('Fraction')
p1.set_xlabel('$\log_{10}$($M_{star}$/M$_{\odot}$)')
#p2.set_xlim(-0.09,1.09)
#p2.set_ylim(0,110)
#p2.set_xlabel('p(S)')
#p2.tick_params(axis='x', labelbottom='off')
#p2.set_ylabel('Total (%)')



plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
#fig.savefig(pj('/Users/chrisfuller/Desktop/','optical_test_' + str(m_lim) + sel_type +'.pdf'))
plt.show()
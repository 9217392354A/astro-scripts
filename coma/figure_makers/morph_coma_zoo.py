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

m_lim = 17.77
sel_type = 'lt' 
texts = '$m_{r} <' + str(m_lim) +'$'



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
    return np.sqrt(a) *100.0 / b


    
####################### main ##########################
    

early = cat.where(cat.early==1)
inter = cat.where(cat.inter==1)
late = cat.where(cat.late==1)


################### plotting ##########################




early_hist,bin_edges    =  histogram(early.SMASS, bins = 10)
inter_hist,bin_edges    =  histogram(inter.SMASS, bins = 10)
late_hist,bin_edges    =  histogram(late.SMASS, bins = 10)


per_early = array(early_hist, dtype=float)*1.0 / len(cat)
per_late = array(late_hist, dtype=float)*1.0 / len(cat)
per_inter = array(inter_hist, dtype=float)*1.0 / len(cat)

fig= plt.figure(figsize = (8.0,8.5),facecolor='w',edgecolor='w')



p1 = plt.subplot(111)

#p1.text(0.02, 0.9, texts, transform=p1.transAxes, fontsize=20, verticalalignment='top')
#plot fraction detected
p1.errorbar(bin_edges[:-1], per_early, yerr=er(per_early,len(cat)), color='red')
p1.errorbar(bin_edges[:-1], per_inter, yerr=er(per_inter,len(cat)), color='green')
p1.errorbar(bin_edges[:-1], per_late, yerr=er(per_late,len(cat)), color='blue')







#p1.set_xlim(-0.09,1.09)
#p1.set_ylim(0,110)
#p1.set_ylabel('Percentage detected (%)')
#p1.set_xticklabels(mList)
#p1.tick_params(axis='x', labelbottom='off')
p1.set_ylabel('Fraction')
p1.set_xlabel('$\log_{10}$($M_{star}$/M$_{\odot}$)')
p2.set_xlim(-0.09,1.09)
p2.set_ylim(0,110)
p2.set_xlabel('p(S)')
#p2.tick_params(axis='x', labelbottom='off')
p2.set_ylabel('Total (%)')

###################################################################################################################
################################################  pS0 - end #######################################################
###################################################################################################################
###################################################################################################################
################################################  pE0 - start #####################################################
###################################################################################################################

#find detected

cluster_m_d  = cat['pE0'][where((cat['RADIUS_VIR'] <  1.0) & (np.nan_to_num(cat['F250']) > 0.0))[0]]
filament_m_d = cat['pE0'][where((cat['RADIUS_VIR'] >= 1.0) & (np.nan_to_num(cat['F250']) > 0.0))[0]]

#find totals
cluster_m_total  = cat['pE0'][where((cat['RADIUS_VIR'] <  1.0))[0]]
filament_m_total = cat['pE0'][where((cat['RADIUS_VIR'] >= 1.0))[0]]



################### plotting ##########################
total_cluster,bin_edges    =  histogram(cluster_m_total, bins = bins_pe)
detected_cluster, bin_edges=  histogram(cluster_m_d, bins = bins_pe)

total_filament, bin_edges   =  histogram( filament_m_total, bins = bins_pe)
detected_filament, bin_edges=  histogram( filament_m_d, bins = bins_pe)

per_cluster = array(detected_cluster, dtype=float)*100. / total_cluster
per_filament = array(detected_filament, dtype=float)*100. / total_filament

p3 = plt.subplot(232)
p4 = plt.subplot(235)

#plot fraction detected
p3.errorbar(bin_edges[:-1]+delta*0.5, per_cluster, yerr=er(detected_cluster,total_cluster), color='red', label='Cluster')
p3.errorbar(bin_edges[:-1]+delta*0.5, per_filament, yerr=er(detected_filament,total_filament), color='blue', label='Filament')

#plot fraction compostion for each cluster
p4.errorbar(bin_edges[:-1]+delta*0.5, total_cluster*100./all_cluster, yerr=er(total_cluster,all_cluster), color='red', label='Cluster', ls='--')
p4.errorbar(bin_edges[:-1]+delta*0.5, total_filament*100./all_filament, yerr=er(total_filament,all_filament), color='blue', label='Filament', ls='--')


p3.set_xlim(-0.09,1.09)
p3.set_ylim(0,110)
#p3.set_ylabel('Percentage detected (%)')
#p3.set_xticklabels(mList)
p3.tick_params(axis='both', labelbottom='off', labelleft='off')
#p3.set_ylabel('Detected (%)')

p4.set_xlim(-0.09,1.09)
p4.set_ylim(0,110)
p4.set_xlabel('p(E)')
p4.tick_params(axis='y', labelleft='off')
#p4.set_ylabel('Total (%)')

p3.invert_xaxis()
p4.invert_xaxis()

###################################################################################################################
################################################  pE0 - end #######################################################
###################################################################################################################



###################################################################################################################
################################################  goldmine - start ################################################
###################################################################################################################


cluster = cat.where((cat['RADIUS_VIR'] <  1.0)& (cat.goldmine != -32768))
filament  = cat.where((cat['RADIUS_VIR'] >= 1.0)& (cat.goldmine != -32768))
 
#find detected
cluster_m_d = cluster['goldmine'][dw(cluster['F250'])]
filament_m_d  =  filament['goldmine'][dw( filament['F250'])]

#find totals
cluster_m_total = cluster['goldmine']
filament_m_total  =  filament['goldmine']



################### plotting ##########################

bins_early_late_later_dwalf = [0,2,10,21]
bina = bins_early_late_later_dwalf


mList =  ['','Early','', 'Late','', 'Irregulars']
mList_nos =  ['E/S', 'Sa/Sb/Sc/Sd', 'Sm/Im/dS']

total_cluster,_    =  histogram(cluster_m_total, bins = bina)
detected_cluster,_ =  histogram(cluster_m_d, bins = bina)

total_filament,_    =  histogram( filament_m_total, bins = bina)
detected_filament,_ =  histogram( filament_m_d, bins = bina)

per_cluster = array(detected_cluster, dtype=float)*100. / total_cluster
per_filament = array(detected_filament, dtype=float)*100. / total_filament

all_cluster = len(cluster['goldmine'])
all_filament = len(filament['goldmine'])

p5 = plt.subplot(233)
p6 = plt.subplot(236)
#p3 = plt.subplot(313)

#plot fraction detected
p5.errorbar(arange(1,len(bina)), per_cluster, yerr=er(detected_cluster,total_cluster), color='red', label='Cluster')
p5.errorbar(arange(1,len(bina)), per_filament, yerr=er(detected_filament,total_filament), color='blue', label='Filament')

#plot fraction compostion for each cluster
p6.errorbar(arange(1,len(bina)), total_cluster*100./all_cluster, yerr=er(total_cluster,all_cluster), color='red', label='Cluster', ls='--')
p6.errorbar(arange(1,len(bina)), total_filament*100./all_filament, yerr=er(total_filament,all_filament), color='blue', label='Filament', ls='--')


p5.set_xlim(0.5,3.5)
p5.set_ylim(0,110)
#p5.set_ylabel('Percentage detected (%)')
#p5.set_xticklabels(mList)
p5.tick_params(axis='both', labelleft = 'off', labelbottom='off')

p6.set_xlim(0.5,3.5)
p6.set_ylim(0,110)
p6.set_xlabel('Goldmine Morphology')
p6.tick_params(axis='y', labelleft='off')
p6.set_xticklabels(mList)
###################################################################################################################
################################################  goldmine - end ##################################################
###################################################################################################################





"""

###################################################################################################################
################################################  Virgo and Fornax - start ################################################
###################################################################################################################


#Inputs
folder2 = '/Users/chrisfuller/Dropbox/phd/herchel/coma/aux_data/'

c3 = 'fornax_at_100mpc_goldmine.fits'
c2 = 'virgo_at_100mpc_goldmine.fits'

virgo  = Table(pj(folder2,c2))
fornax = Table(pj(folder2,c3))
 
#find detected
fornax_m_d = fornax['goldmine'][dw(fornax['F250'])]
virgo_m_d  =  virgo['goldmine'][dw( virgo['F250'])]

#find totals
fornax_m_total = fornax['goldmine']
virgo_m_total  =  virgo['goldmine']



################### plotting ##########################

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

p5 = plt.subplot(233)
p6 = plt.subplot(236)
#p3 = plt.subplot(313)

#plot fraction detected
p5.errorbar(arange(1,len(bina)), per_fornax, yerr=er(detected_fornax,total_fornax), color='purple', label='fornax', ls='-.',lw=2)
p5.errorbar(arange(1,len(bina)), per_virgo, yerr=er(detected_virgo,total_virgo), color='green', label='virgo', ls='-.',lw=2)

#plot fraction compostion for each fornax
p6.errorbar(arange(1,len(bina)), total_fornax*100./all_fornax, yerr=er(total_fornax,all_fornax), color='purple', label='fornax', ls='-.',lw=2)
p6.errorbar(arange(1,len(bina)), total_virgo*100./all_virgo, yerr=er(total_virgo,all_virgo), color='green', label='virgo', ls='-.')


p5.set_xlim(0.5,3.5)
p5.set_ylim(0,110)
#p5.set_ylabel('Percentage detected (%)')
#p5.set_xticklabels(mList)
p5.tick_params(axis='both', labelleft = 'off', labelbottom='off')

p6.set_xlim(0.5,3.5)
p6.set_ylim(0,110)
p6.set_xlabel('Goldmine Morphology')
p6.tick_params(axis='y', labelleft='off')
p6.set_xticklabels(mList)
###################################################################################################################
################################################  fornax and Virgo - end ##################################################
###################################################################################################################
"""
p6.legend()

plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.99, wspace=0.0, hspace=0.0)
#fig.savefig(pj('/Users/chrisfuller/Desktop/','optical_test_' + str(m_lim) + sel_type +'.pdf'))
plt.show()










"""
print 'Morphological  &  \\multicolumn{3}{c}{filament} & \\multicolumn{3}{c}{cluster} \\\\'
print 'Type & Total & Detected & \\% & Total & Detected & \\%\\\\'
for i in range(0,len(bins_pe)-1):

    print bins_pe[i],'&',total_filament[i],'&',detected_filament[i],'&',p(total_filament[i],detected_filament[i]),'&',total_cluster[i],'&',detected_cluster[i],'&',p(total_cluster[i],detected_cluster[i]),'\\\\'
"""
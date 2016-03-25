#program to create hi limits from alfacat
# Chris Fuller, July 2014

#import mods
from atpy import Table
import numpy as np
from os.path import join as pj



#Inputs
folder = "/Users/chrisfuller/Dropbox/phd/herchel/coma/final_outputs/" # input/output folder
fname = 'coma_supercluster_cal12.fits' #input name
cat = Table(pj(folder,fname))

HiComp = []

for i in range(len(cat)):
	w50 = cat.W50_1[i]
	s21 = cat['Si(HI)'][i]
		
	if 0.0 == np.nan_to_num(w50):
		result = -1 

	if w50 > 2.5:
		test90 = 0.5*np.log10(w50) - 1.11
		test50 = test90 - 0.130
		test25 = test90 - 0.202

	elif w50 <= 2.5:
	 	test90 = np.log10(w50) - 2.36
	 	test50 = test90 - 0.130
	 	test25 = test90 - 0.202


	if np.log10(s21) >= test90:
		result = 90

	elif np.log10(s21) >= test50:
		result = 50

	elif np.log10(s21) >= test25:
		result = 25

	else:
		result = 0

	HiComp.append(result)

hiComp = np.array(HiComp)

w90 = np.where(hiComp == 90)[0]
w50 = np.where(hiComp == 50)[0]
w25 = np.where(hiComp == 25)[0]
wn90= np.where(hiComp == 0)[0]

import matplotlib.pyplot as plt

xx1 = np.linspace(1.2,2.5,40)
yy1a = 0.5*xx1 - 1.11

xx2 = np.linspace(2.5,3.0,40)
yy2a = xx2 - 2.36


yy1b = yy1a - 0.130
yy2b = yy2a - 0.130


yy1c = yy1a - 0.202
yy2c = yy2a - 0.202

plt.plot(np.log10(cat.W50_1[w90]), np.log10(cat['Si(HI)'][w90]), 'ok')
plt.plot(np.log10(cat.W50_1[w50]), np.log10(cat['Si(HI)'][w50]), 'og')
plt.plot(np.log10(cat.W50_1[w25]), np.log10(cat['Si(HI)'][w25]), 'ob')
plt.plot(np.log10(cat.W50_1[wn90]), np.log10(cat['Si(HI)'][wn90]), 'or')
plt.plot(xx1,yy1a, '-k', xx2, yy2a, '-k', label='90$\%$ completeness')
plt.plot(xx1,yy1b, '--g', xx2, yy2b, '--g', label='50$\%$ completeness')
plt.plot(xx1,yy1c, '--b', xx2, yy2c, '--b', label='25$\%$ completeness')
#plt.legend()
plt.show()





cat.add_column('hiComp', hiComp)
cat.write(pj(folder,'test-v2.fits'))
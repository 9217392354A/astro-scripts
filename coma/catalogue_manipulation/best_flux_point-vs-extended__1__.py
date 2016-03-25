# program for finding the best flux between point and exteneded source measurement methods.
# program written by Chris Fuller, Jan 2014

#import modules
import numpy as np
from os.path import join as pj
import atpy as at
from copy import copy, deepcopy
import matplotlib.pyplot as plt

#switches
five_band = True #if True it only selects galaxies that are in all five bands ie inpacs and hevics_plw == 1



#folder stuff
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/initial_outputs/'
outfolder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/catalogue_creation_phase/'
outname = 'ngp+20140331__best-flux__.fits'
#load catalogues one where all are forced as point to serve as the base catalogue
extended = at.Table(pj(folder,"ngp+atlas+fluxes-20140328-EXTENDED-mybgsub-cal12-rerun.fits"),type='fits')
point    = at.Table(pj(folder,"ngp+atlas+fluxes-20140328-POINT-mybgsub-cal12-rerun.fits")   ,type='fits')

#bands
bands = ['500','350','250','160','100']
beams =  [36.0, 24.5, 18.2, 13.4, 9.4] #beams in arcsec
colour = ['b', 'r', 'k', 'g', 'c']

# # # # # # # # # # # # # # # # # # # # # # # # Functions # # # # # # # # # # # # # # # # # # # # # # # #
def remove_badgal(cat):
	#find which rows have pacs and spire
	return cat.where((cat['IN_PACS'] == 1) & (cat['HEVICS_PLW'] == 1))

def mod(x):
	return np.sqrt(x**2)

def remove_single_from_cat(index,t,band):
	cols = ['F', 'SN', 'R']
	#now times the error by 3 for all galaxies removed
	if t['F'+band][index] != 0.0:
		#print 'OBJECT: ', t['OBJECT'][index], 'BAND: ', band, 'FLUX: ', t['F'+band][index], 'ERROR: ', t['E'+band][index]
		t['E'+band][index] = t['E'+band][index]*3.0

	#loop through each col
	for col in cols:
		#set col to zero
		t[col+band][index] = 0.0 

def set_ext(jj,type_val,ext):
	y=[] #create empty list
	for kk in range(0,len(ext)):
		if kk == jj: y.append(type_val)
		else: y.append(ext[kk])

	return ''.join(y)

# # # # # # # # # # # # # # # # # # # # # # # # Functions # # # # # # # # # # # # # # # # # # # # # # # #

#remove galaxies that arn't in all 5 bands
if five_band == True:
	extended = remove_badgal(extended)
	point    = remove_badgal(point)

#test that both are the same length
if len(extended) != len(point):
	raise 'Tables different lenghts.... #sigh#'


#create new cat
new = deepcopy(point)
counts = []

#repeat analsis for each band
for j in range(0, len(bands)):
	band = bands[j]
	beam = beams[j]
	print band
	print ''
	count = 0
	#work rowwise down each band
	for i in range(0,len(point)):
		#fluxes
		F_point = point['F'+band][i]*10**3
		F_extended = extended['F'+band][i]*10**3

		#fluxes
		E_point = point['E'+band][i]*10**3
		E_extended = extended['E'+band][i]*10**3

		#extendedness
		T_point = point['EXTENDEDNESS'][i][j]
		T_extended = extended['EXTENDEDNESS'][i][j]

		#Radii
		R_extended = extended['R'+band][i]*60.0
		a = extended['FULLMAJAX'][i]*60.0
		#if undetected as a pure point then skip to next source
		if (F_point == 0.0) or (F_extended == 0.0): 
			remove_single_from_cat(i,new,band)
			#count +=1
			


		#if both are point sources skip to next source
		elif ((T_point == 'P') and (T_extended == 'P')): continue
 
		#if sources are different inside 1sigma skip to the next source and use point method
		elif ((F_extended - F_point) > np.sqrt(E_point**2.0 + E_extended**2.0)) and (R_extended < 5.0*a) and (R_extended > beam): #use extended flux
			#print np.around(R_extended), np.around(a), beam
			count +=1
			new['F'+band][i] = extended['F'+band][i]
			new['E'+band][i] = extended['E'+band][i]
			new['R'+band][i] = extended['R'+band][i]
			new['SN'+band][i] = extended['SN'+band][i]
			new['EXTENDEDNESS'][i] = set_ext(j,'E',new['EXTENDEDNESS'][i])
	counts.append(count)

print bands
print counts
new.write(pj(outfolder, outname), overwrite=True)

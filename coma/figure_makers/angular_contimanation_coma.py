# program written to find the angular contamination by using a simulated optical catalogue

# Chris Fuller, Jan 2013

#import modules
import numpy as np
from os.path import join as pj
import atpy as at
from copy import copy, deepcopy
import matplotlib.pyplot as plt
from random import random
from pylab import bar
from time import time
import datetime


#i/o
print 'reading in cats . . .'
folder = '/Users/chrisfuller/Dropbox/phd/herchel/coma/source_measurement/positions/'
cat_mine = at.Table(pj(folder, 'ngp+20140109__best-flux__.fits'),type='fits')
catA = at.Table(pj(folder, 'ngp_psf.fits'),type='fits')



#key parameters
r_search = 10.0 / 3600.0 # search radius in degrees
coor_names_A = ['GRA2000', 'GDEC2000'] # these are the colum names that containe ra and dec
n=10 #number of interations of monty carlo
n_bins = 35 #number of radial bins

# # # # # # # # # # # # # # # Function # # # # # # # # # # # # # # # # # # # # # # # # # #

def finished():
	from pygame import mixer#
	mixer.init() #you must initialize the mixer
	alert=mixer.Sound('/Users/chrisfuller/Dropbox/personal/sounds/R2D2e.wav')

	import smtplib

	fromaddr = 'farpoint.fuller@gmail.com'
	toaddrs  = 'farpoint.fuller@gmail.com'
	msg = 'So long and thanks for all the fish!'


	# Credentials (if needed)
	username = 'farpoint.fuller'
	password = '1L=&8_a$>7/$*5#'

	# The actual mail send
	server = smtplib.SMTP('smtp.gmail.com:587')
	server.starttls()
	server.login(username,password)
	server.sendmail(fromaddr, toaddrs, msg)
	server.quit()

	alert.play()

def email(msg):
	import smtplib

	fromaddr = 'farpoint.fuller@gmail.com'
	toaddrs  = 'farpoint.fuller@gmail.com'


	# Credentials (if needed)
	username = 'farpoint.fuller'
	password = '1L=&8_a$>7/$*5#'

	# The actual mail send
	server = smtplib.SMTP('smtp.gmail.com:587')
	server.starttls()
	server.login(username,password)
	server.sendmail(fromaddr, toaddrs, msg)
	server.quit()

#function to produce new cat with matched columns
def x_matcher(catA, coor_names_A, catB, coor_names_B, r_search):
	print 'starting x matcher . . .'
	###### part 1 ######
	# first create new cat structure to hold old cat and new matches
	t = at.Table() #empty cat object

	colsA =  catA.columns #get columns of the first table
	colsB =  catB.columns #get columns of the second table

	cols_names_b = colsB.keys #col names for cat b

	#loop through each column adding it to the table 
	for j in range(0, len(colsA)): t.add_column(colsA.keys[j],catA[colsA.keys[j]],dtype = colsA[colsA.keys[j]])

	#loop through each column adding an empty one to the new empty table
	for k in range(0, len(colsB)): t.add_empty_column(colsB.keys[k],dtype = colsB[colsB.keys[k]])

	#add column for seporation
	t.add_empty_column('radial_sep', dtype = np.float)

	#if verbous: print t.columns #print list of new cols

	###### part 2 #######
	# find nearest neighbours

	#create arrays of ra and dec for cat B, save making new ones each loop
	ra_2_long = catB[coor_names_B[0]]
	dec_2_long = catB[coor_names_B[1]]

	#loop through all members of catA
	for i in range(0, len(t)):

		####  wiget #####
		#if i%100 == 0: print np.round(i*100.0/len(t), decimals=2), '%'

		ra = t[coor_names_A[0]][i]
		dec = t[coor_names_A[1]][i]
		
		#create a search box around the source
		temp_cat = search_box(ra, dec, ra_2_long , dec_2_long, r_search, catB)

		#create new list of ra and dec
		ra_2 = temp_cat[coor_names_B[0]]
		dec_2 = temp_cat[coor_names_B[1]]

		#ra1 and dec1
		ra_1 = np.array([t[coor_names_A[0]][i]]*len(ra_2), dtype=np.float)
		dec_1 = np.array([t[coor_names_A[1]][i]]*len(ra_2), dtype=np.float)

		#if verbous: print 'caculating distance for all soruces'
		#caculate distance to all sources from ra1 and dec1
		radius = distance(ra_1, dec_1, ra_2, dec_2)

		#find where radius is closer than search radius
		w = np.where(radius < r_search)[0]	
	
		#now we need to find the nearest neigbourgh
		if len(w) == 0: #no source found that is whithin search radius
			#if verbous: print 'No source found within ', np.around(r_search*3600.0), ' moving on to next source'
			t['radial_sep'][i] = -1.0

		elif len(w) >= 1: # a single or multipul source found
			#add all the seps to a temp cat
			distances = np.array(radius[w], dtype=np.float)
			keys = np.array(w, dtype=np.int)

			d = np.sort(distances)[0]
			 
			#reshape into an array
			temp = np.append(distances.reshape(len(keys),1) , keys.reshape(len(keys),1),1)

			#sort to find the nearest source
			temp = deepcopy(temp[np.argsort(temp[:,0])])

			#now the top row second element will hold the key for the galaxy
			# from catB
			key = temp[0,1]
			rad = temp[0,0]
		
			#now add source to new cat
			t['radial_sep'][i] = d*3600.0

			#loop through each column adding to new cat
			if True: 
				for j in range(0, len(cols_names_b)): t[cols_names_b[j]][i] = temp_cat[cols_names_b[j]][key]
	return deepcopy(t.where(t['radial_sep']>0.0)) #.write(pj(folder, 'test.fits'), overwrite=True)

#distance equation designed to do arraywise caculations
def distance(ra1, dec1, ra2, dec2):
	delta_ra = (ra1 - ra2) * np.cos(np.radians((dec1+dec2)/2.0))
	delta_dec = (dec1 - dec2)

	return np.sqrt(delta_ra**2.0 + delta_dec**2.0)

#returns a cat with galaxies inside a search box
def search_box(box_ra, box_dec, box_ra_2 , box_dec_2, r, box_cat):
	#factor of box greater than search radius
	factor = 2.0
	
	max_ra = box_ra + r*factor
	min_ra = box_ra - r*factor

	max_dec = box_dec + r*factor
	min_dec = box_dec - r*factor

	new_cat = box_cat.where((box_ra_2 < max_ra) & (box_ra_2 > (min_ra)) & (box_dec_2 < max_dec) & (box_dec_2 > min_dec))

	return new_cat

def random_cat(cen_x, cen_y, size, length):
	print 'creating random cat . . .'
	rand1, rand2 =[], [] #create empty lists

	#create list of random numbers
	for i in range(0,length):
		rand1.append(random())
		rand2.append(random())

	# turn them into numpy arrays
	rand1 = np.array(rand1, dtype= np.float)
	rand2 = np.array(rand2, dtype= np.float)

	#create lists of ra and dec
	ra_rand =  cen_x + (rand1 - 0.5)*size
	dec_rand = cen_y + (rand2 - 0.5)*size

	# create new cat structure
	t = at.Table() #empty cat object

	#add cols for ra dec
	t.add_column('ra_rand',ra_rand, dtype=np.float)
	t.add_column('dec_rand', dec_rand, dtype=np.float)

	return t

def let_chris_know(time):
	if time < 60.*5.0: 
		print 'I should finish in under 5 minutes you may as well hold fast'
	if time > 60.*5.0 and time < 60.*30.0: 
		msg = 'this is going to take' + str(datetime.timedelta(seconds=loop_time*n)) +'/n'
		msg += 'go and have lunch'
		email( msg)

	if time > 60.*30.0 and time < 60.*60.0*2.0: 
		msg = 'this is going to take' + str(datetime.timedelta(seconds=loop_time*n)) +'/n'
		msg += 'this better be worth my time....'
		email( msg)

	if time > 60.*60.0*2.0 and time < 60.*60.0*10.0: 
		msg = 'this is going to take' + str(datetime.timedelta(seconds=loop_time*n)) +'/n'
		msg += 'go to bed meatbag... leave this to the bigboys'
		email( msg)

	if time < 60.*60.0*10.0: 
		msg = 'this is going to take' + str(datetime.timedelta(seconds=loop_time*n)) +'/n'
		msg += 'go to bed meatbag... leave this to the bigboys'
		email( msg)

	if str(time) == '7.5 million years': raise '42'

# for saving an array to a text file .csv
def savetext(folder,filename,array):
    outfile = open(pj(folder, filename), 'w')
    for row in array:
        count = 0
        for element in row:
            count += 1
            if count < len(row):
                x = element + ","
            else:
                x = element
            outfile.write(x)
        outfile.write("\n")
    outfile.close()

 # # # # # # # # # # # # # # # Main Program # # # # # # # # # # # # # # # # # # # # # # # #

####### part 1 ###### create plot of random radial distribution ####### part 1 ###### 
start = time()

print 'simulation random distribution'
for i in range(0,n):
	print "starting cycle ",i+1, ' of ',n
	rand_cat = random_cat(200.0, 29.0, 11.0, len(cat_mine))
	new_cat = x_matcher(rand_cat, ['ra_rand', 'dec_rand'], catA, ['ra', 'dec'], r_search )

	#list of seps
	seps_list = new_cat.radial_sep


	#find sources with matches
	hist, binedges = np.histogram(seps_list, bins=np.arange(0.0,r_search*3600.0, r_search*3600.0/n_bins) , range=None, normed=False, weights=None, density=None)

	#now do mean of sample
	if i == 0:
		rand_hist = hist
		loop_time = time() - start
		print 'program will complete in ', str(datetime.timedelta(seconds=loop_time*n))
		try:let_chris_know(loop_time*n)
		except: continue
	if i > 0: 
		rand_hist = rand_hist + hist
#find average random hist
rand_hist = rand_hist.astype(np.float)/np.float(n)

####### part 2 ###### perform x_match on my optical cat ####### part 2 ###### 
print 'caculating my distribution'
new_cat_mine = x_matcher(cat_mine, ['GRA2000', 'GDEC2000'], catA, ['ra', 'dec'], r_search )
#list of seps
#find sources with matches
mine_hist, binedges = np.histogram(new_cat_mine.radial_sep, bins=np.arange(0.0,r_search*3600.0, r_search*3600.0/n_bins) , range=None, normed=False, weights=None, density=None)

#new_cat.write(pj(folder, 'test.fits'), overwrite=True)


####### part 3 ###### caculate the % contamination per bin ####### part 3 ###### 
conta_hist_fraction = rand_hist.astype(np.float)/ mine_hist.astype(np.float)

#fit polynomial to contamination plot
p_con = np.poly1d(np.polyfit(binedges[:-1]+r_search*3600.0/n_bins,conta_hist_fraction,2.))

X_con = np.arange(0.0,r_search*3600.0, 0.001)
Y_con = p_con(X_con)

#save all relevent details
a = rand_hist.reshape(len(conta_hist_fraction),1)
b = mine_hist.reshape(len(conta_hist_fraction),1)
c = conta_hist_fraction.reshape(len(conta_hist_fraction),1)
d = binedges[:-1].reshape(len(conta_hist_fraction),1) + r_search*3600.0/n_bins

a = np.append(a, b, 1)
a = np.append(a, c, 1)
a = np.append(a, d, 1)

header = np.array(['random_histogram', 'my_optical_hist', 'contamination_fraction', 'angular_distance_arcsec'], dtype=np.str).reshape(1,4)

a = np.append(header,a,0)
savetext(folder,'values_'+str(n)+'_.csv',a)


####### part 4 ####### plot results ###### part 4 #########

fig = plt.figure(figsize = (4.6,4.6),facecolor='w',edgecolor='w')
sub1 = fig.add_subplot(2,1,1)
sub2 = fig.add_subplot(2,1,2)

sub1.bar(binedges[:-1], rand_hist, width=r_search*3600.0/n_bins)#, fc ='None',hatch='/')
sub1.bar(binedges[:-1], mine_hist, width=r_search*3600.0/n_bins, fc ='None',hatch='/')
sub1.tick_params(axis='x', labelbottom='off')
sub1.set_ylabel('Number of galaxies')

sub2.bar(binedges[:-1], conta_hist_fraction, width=r_search*3600.0/n_bins, fc ='b')
sub2.plot(X_con, Y_con, 'k--')
sub2.set_xlabel('Angular Separation, (arcsec)')
sub2.set_ylabel('Fractional Contamination')
sub2.set_ylim(ymin=0)
plt.subplots_adjust(left=0.12, bottom=0.06, right=0.97, top=0.98,hspace=0.0)

fig.savefig('/Users/chrisfuller/Desktop/angular_contamination.pdf', dpi=600.0)
finished()
plt.show()

#program to measure noise varaition between two maps

#import stuff
from pyfits import fits
from numpy import where, isnan, std, median
import os




folder = '/home/scuba2/spx6cff/coma/maps/'

imps = [s for s in os.listdir(folder) if ".fits" in s]
print 'list of images'
print '--------------'*3
for rebs in imps: print rebs
print '--------------'*3
os.chdir(folder)

for im in imps:

	#import fits file
	print '--------------'*3
	print 'reading in image....', im
	hdulist = fits.open(im)

	#get image
	try:image =  hdulist[1].data
	except: image = hdulist[0].data
	#find pixels in image that arn't nan
	w = where(isnan(image)== False)
	pix_val = image[w].flatten()

	hdulist.close()

	#now sigma clip uptill noise is found

	old_sigma = std(pix_val)
	new_sigma = 1000.0


	safe = 0

	print 'sigma clipping starting....'
	while (old_sigma - new_sigma)**2  > 0.0000000000001:
		old_sigma = new_sigma

		# on the 10th iteration break as something has gone wrong
		if safe >= 100: break
		safe +=1

		#caculate sigma
		new_sigma = std(pix_val)

		#median 
		med = 0.0 #median(pix_val)

		#remove everything outside 3sigma
		w_clip = where((pix_val < (med + 3.*new_sigma)) & (pix_val > (med - 3.*new_sigma)))

		pix_val = pix_val[w_clip[0]]
		
	print im, ' ---- ', new_sigma, safe
	print '--------------'*3

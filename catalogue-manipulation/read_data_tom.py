# read in data 
#import mods

#from atpy import Table
from os.path import join as pj
import numpy as np




fname= 'plot_galaxyspectrum_template.cold.dat'
folder = "/Users/chrisfuller/Dropbox/" # input/output folder



#function to make a table that is easy to play about with
def makeMeATable(path):

	#create data array
	data = np.loadtxt(path, dtype=np.float,  skiprows=12, delimiter="   ", unpack=False)

	#create empty dic
	templates = {}


	#find how many templates there are
	Ntemps = len(data[0,:]) - 1
	Nwaves = len(data[:,0])

	wavelengths = data[:,0]


	#loop through templates
	for i in range(Ntemps):

		templates[i] =  data[:,i]


	return wavelengths, templates



waves, temps = makeMeATable(pj(folder, fname))

print 'finished'
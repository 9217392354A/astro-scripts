#program to fit a gaussian to data
# Chris Fuller, March 2014

#import modules
import numpy as np
from scipy.optimize import fmin
import matplotlib.pyplot as plt 

#### functions ######
#gaussian function with a cisquared likehood
def gauss(x,*args): 
	#X, Y, and Y error data extracted from args   
    xDat = args[0]
    yDat = args[1]
    yError=args[2]

    #gaussian parameters   
    hmax = x[0]
    hmin = x[1]
    mean =x[2]
    sigma = x[3]

    #loop throguh data and sum ciaquared
    A = (1.0/np.sqrt(2*np.pi*sigma**2))
    B = -((xDat-mean)**2)/(2.0*sigma**2)
    fx = (hmax*(A*np.exp(B)))+hmin
    
    return np.sum(((yDat - fx)**2)/ (yError**2)) 

#function to genorate a gaussian
def genGauss(hmax, hmin, mean, sigma, x):
    y = [] 
    for i in range(0,len(x)):
        A = (1.0/np.sqrt(2*np.pi*sigma**2.0))
        B = -((x[i]-mean)**2.0)/(2.0*sigma**2.0)
        temp = (hmax*A*np.exp(B))
        temp += hmin
        y.append(temp)
    return y

##### end functions ####

#first genorate gaussian distribution of data
x_data = np.arange(-100.0,100.0,0.1) 
y_data = genGauss(10.0, 1.0, 19.0, 10.0, x_data) + np.random.rand(1,len(x_data)).flatten()*0.03


#intertive clipping data loop thingy
sigma_new = 0.0
sigma_old = 10000000.0

#create initial guess
sigma = 1.0
mean = 200.0
hmin = 1.0
hmax = 1.0
guess = [hmax,hmin,mean,sigma]

fit = fmin(gauss,guess,args=(x_data, y_data, np.ones(len(x_data)).flatten()))

# work through 
while np.abs(sigma_old - sigma_new)*100.0 / sigma_old > 70.0:
	fit =  fmin(gauss,fit,args=(x_data, y_data, np.ones(len(x_data)).flatten()))
	sigma_new = fit[3]

	mu = fit[2] #find mean

	w = np.where((x_data > mu - sigma_new*2.0) & (x_data < mu + sigma_new*2.0))[0] #select data inside 2sigma

	x_data = x_data[w]
	y_data = y_data[w] #select where try in ydata

	sigma_old = sigma_new




#create initial guess
sigma = 1.0
mean = 1.0
hmin = 1.0
hmax = 1.0
guess = [hmax,hmin,mean,sigma]

#fit gaussian 
fit = fmin(gauss,guess,args=(x_data, y_data, np.ones(len(x_data)).flatten()))

#genorate result from fit
y_fit = genGauss(fit[0], fit[1], fit[2], fit[3],  x_data)

#plot result
plt.plot(x_data, y_data, 'xk')
plt.plot(x_data, y_fit, '-r')
plt.show()





#import mods
import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters



def model(params, x):
	phiStar = params['phiStar'].value
	LStar = params['LStar'].value
	alpha = params['alpha'].value + 1

	L = (x/LStar)
	phi = phiStar*(L**alpha)*np.exp(-1.0*L)
	return phi


def residual(params, x, y_data, y_error):
	return (y_data-model(params, x))/y_error



colours = ['g','r','b']


#guess = np.array([1.0E8, 1.0E9, 1.0E10])
#guess = np.array([1., 0.1, 0.01])
guess  = np.array([-2., -1., 0.])
for i in range(3):


	#model
	params = Parameters()
	params.add('phiStar', value=1.0, vary=True)
	params.add('LStar', value=1.0E10, vary=True)
	params.add('alpha', value=guess[i], vary=True)



	x_data = np.linspace(6,14)
	#out = minimize(residual, params, args=(10**x_data, 10**y_data, 10**y_error))
	y_line = np.log10(model(params, 10**x_data))



	#plt.plot(x_data, y_data, 'kx')
	#plt.axhline(y=guess[i], ls='--', c=colours[i])
	plt.plot(x_data, y_line, ls='--', c=colours[i], label=guess[i])

plt.ylim(-4, 1)
plt.xlim(6, 12)

plt.legend()
plt.show()

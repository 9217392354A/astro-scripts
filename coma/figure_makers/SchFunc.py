#import mods
import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit

def residual(params, x, y_data, y_error):
	return y_data-model(params, x)

def lmfitter(x , y, y_error):
	params = Parameters()
	params.add('phiStar', value=0.1, vary=True)
	params.add('LStar', value=1.0, vary=True)
	params.add('alpha', value=-1, vary=True)


	# remove inf values in errors
	out = minimize(residual, params, args=(x, y, y_error))
	report_fit(params)
	return out

def model(params, x):
	phiStar = params['phiStar'].value
	LStar = params['LStar'].value
	alpha = params['alpha'].value

	L = (x/LStar)
	return phiStar*(L**alpha)*np.exp(-L)
	




x_data = np.array(		[5.47673083e-04,   1.49203807e-03,   4.06479280e-03,   1.10738062e-02, 3.01686188e-02,   8.21890453e-02,   2.23909460e-01,   6.10001568e-01, 1.66184096e+00,   4.52739062e+00], dtype=np.float64)

y_data = np.array(		[ 0.01751313,  0.09632224,  0.74430823,  2.05779335,  2.21541156,  0.83187391, 0.43782837,  0.22767075,  0.14010508,  0.00875657], dtype=np.float64)

y_error = np.array(		[ 0.01238366,  0.02904225,  0.08073156,  0.13423564,  0.13928173,  0.08534846, 0.06191828,  0.04464991,  0.03502627,  0.00875657], dtype=np.float64)

x_mod = np.linspace(x_data[0], x_data[-1], 200)

if True:
	#fit the data
	#out = lmfitter(x_data , y_data, y_error)
	#y_mod =  model(out.params, x_mod)

	plt.errorbar(x_data, y_data, yerr=y_error, c='r', ls='', label='Raw Data')
	plt.plot(x_data, y_data, 'xr', label='Raw Data' )
	#plt.plot(x_mod, y_mod, 'k--', label='Model Fit')


if True:
	x_mod = np.logspace(-2, 1, 200)
	Lstar = 0.1
	phiStar = 0.3
	x = x_mod/Lstar

	y_modMine = phiStar*(x**0.8)*np.exp(-x)
	plt.plot(x_mod, y_modMine, '-g', label='Model Mine')


plt.loglog()
plt.legend(loc=5)
plt.show()

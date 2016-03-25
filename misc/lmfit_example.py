
#import modules
import numpy as np
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, report_fit



def residual(params, x, y_data, y_error):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    e = params['e'].value
   
    y_model = a*x**4 + b*x**3 + c*x**2 + d*x + e

    return (y_data-y_model)/y_error

def lmfitter(x , y, y_error):
    params = Parameters()
    params.add('a', value=-1., vary=True)
    params.add('b', value=1., vary=True)
    params.add('c', value=1., vary=True)
    params.add('d', value=1., vary=True)
    params.add('e', value=0.0, vary=False)


    # remove inf values in errors
    y_error[y_error ==  np.float('inf')] = 1.0  
    out = minimize(residual, params, args=(x, y, y_error))
    report_fit(params)
    return out

def model(params, x):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    e = params['e'].value
   
    return a*x**4 + b*x**3 + c*x**2 + d*x + e

def rmodel(params, x):
    x_long = np.round(np.arange(0.,10.0, 0.001), decimals = 2)
    y_long = np.round(model(params, x_long), decimals = 2)

    #round all to the same decima  place
    sol = x_long[y_long == np.round(x, decimals=2)][0]
    print x_long[y_long == np.round(x, decimals=2)]
    return sol



xx = np.array([ 0.33333333,  0.66666667,  1.        ,  1.33333333,  1.66666667,
        2.        ,  2.33333333,  2.66666667,  3.        ,  3.33333333,
        3.66666667,  4.        ,  4.33333333,  4.66666667])

yy = np.array([  0.23809524,   0.40243902,   0.36792453,   0.70175439,
         0.7037037 ,   1.29487179,   1.67647059,   3.16      ,
         4.4       ,   4.97727273,   8.24137931,   8.92857143,
        10.22222222,  12.86956522])

#yy_error = np.np.array([1.0]*len(xx), dtype= np.float)

yy_error = np.array([ 0.03673889,  0.04444196,  0.03573599,  0.06572532,  0.06771392,
        0.14661536,  0.20330192,  0.44689149,  0.69570109,  0.7503521 ,
        1.53038572,  1.6873414 ,  1.96726758,  2.68348985])

#yy_error = np.array([1.0]*len(yy), dtype= np.float)
print 'error'

fit = lmfitter(xx , yy, yy_error)

contam = rmodel(fit.params, 5.0)
plt.errorbar(xx, yy, yerr=yy_error)
plt.plot(xx,model(fit.params, xx), '-g')
plt.plot(xx,yy, 'ro')
plt.show()

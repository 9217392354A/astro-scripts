# Quick program to solve an equation numerically
# Chris Fuller Sept 2013

from numpy import log, log10, logspace
from scipy.optimize import fmin


def eq(M):
	f = M/(dust+stars+M)
	return f*log(1/f)*p*n

def gas_frac(x):
	val = eq(x)
	return (dust- val)

def flux2gal():
	A = 2.36*10**5
	D = 17.5
	dv = 5
	flux = 15.*10.**-3.
	return A*D**2 * flux * dv





p = 0.004
n = 0.5
dust = 10**5.2
stars = 10**6.76

print log10(flux2gal())
#for stars in logspace(1,10):
#max_dust = fmin(gas_frac, 10.0**7)
#print log10(stars),log10(max_dust)


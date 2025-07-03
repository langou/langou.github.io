# Program 2.3 Broyden's Method II
# 2e update to Python 2/21/2012
# Input: x0 initial vector, k = max steps
# Output: solution x
# Example usage: broyden2(f,[1,2],20)

from numpy import *

def broyden2(f,x0in,k):
	# Allows input of x0 as list or 1D array
	x0 = array(matrix(x0in,dtype=float).T)
	n = x0.shape[0]
	b = eye(n)					# initial b
	for i in range(k):
		x = x0 - dot( b, f(x0) )
		deltax = x-x0; deltaf = f(x)-f(x0)
		bdeltaf = dot( b, deltaf )
		tmp1 = dot( deltax-bdeltaf , deltax.T )
		tmp2 = dot( tmp1, b)
		tmp3 = dot( deltax.T, dot( b, deltaf ) )
		b += tmp2/tmp3
		x0 = x
		print i, x[:,0].tolist()
	return x[:,0].tolist()

def myF(x):
	[u,v] = x
	F = array( [u*u - 4.0*v*v - 4.0, (u-1.0)**2.0 + v*v - 4.0 ] )
	return F

x0 = [3,1]
print 0,x0
sol = broyden2( myF, x0, 12 )
print sol


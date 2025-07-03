# Translation to Python/Numpy of Sauer's
# Program 3.1 Newton Divided Difference Interpolation Method
# Computes coefficients of interpolating polynomial
# Input: x and y are vectors containing the x and y coordinates
#        of the n data points
# Output: coefficients c of interpolating polynomial in nested form
# Use with nest.m to evaluate interpolating polynomial

from numpy import *

def newtdd(x,y):
	n = len(x)
	v = zeros((n,n))
	for j in range(n):
		v[j,0] = y[j]			# Fill in y column of Newton triangle
	for i in range(1,n):		# For column i,
		for j in range(n-i):	# 1:n+1-i		# fill in column from top to bottom
#			print j,i," ",j+1,i-1," ", j,i-1," ",j+i," ",j
			v[j,i] = (v[j+1,i-1]-v[j,i-1])/(x[j+i]-x[j])
	c = v[0,:].copy()			# Read along top of triangle for output coefficients
	return c

'''
# test_newtdd.py
from newtdd import newtdd
from nest   import nest
from numpy  import linspace
from pylab  import plot, show
xdata = [0.0,1,2,3,3.4]
ydata = [0.0,3,12,6,2]
plot(xdata,ydata,'ro')
interpolant = newtdd(xdata,ydata)
x = linspace(-1,4,100,endpoint=True)
y = nest(interpolant,x,xdata)
plot(x,y,'b')
show()
'''

# Program 6.2 Euler's Method for Solving Initial Value Problems
# Python translation JR 1/17/12
# Update for 2e 2/24/12
# ydot evaluates rhs of differential equation
# Inputs: interval inter, initial vector y0, number of steps n
# Output: time steps t, solution y
# Example usage: y=euler2([0,1],[0,1],10);

from numpy import *
from pylab import plot, show

def euler2(interval,y0,n):
	h = float(interval[1]-interval[0])/n
	t = [i*h for i in range(n+1)]
	y = empty((n+1,2)); y[0,:] = y0
	for i in range(n):
		y[i+1,:] = eulerstep(t[i],y[i,:],h) 
	plot(t,y[:,0])
	plot(t,y[:,1])
	show()
	
def eulerstep(t,y,h):
	return y + h*ydot(t,y)

def ydot(t,y):
	return array( [ y[1]**2-2*y[0], y[0]-y[1]-t*y[1]**2 ] )


# try it
euler2( [0,1], [0.,1.], 10 )



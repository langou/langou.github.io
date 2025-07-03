# Translation to Python of Sauer's
# Program 1.2 Fixed-Point Iteration
# Computes approximate solution of g(x)=x
# Input: inline function g, starting guess x0, 
#        number of steps k
# Output: Approximate solution
def fpi(g,x0,k):
	x=x0; print x
	for i in range(k):
		x = g(x); print x
	return x

'''
# test_fpi.py
def g1(x): 
	return 1.0 - x**3.0
def g2(x): 
	return (1.0 - x)**(1.0/3.0)
def g3(x):
	return (2.0*x**3.0 + 1.0)/(3.0*x**2.0 + 1.0);

from fpi import fpi
xc = fpi(g1,0.5,15); print 'For g1, get', xc
xc = fpi(g2,0.5,15); print 'For g2, get', xc
xc = fpi(g3,0.5,15); print 'For g3, get', xc
'''

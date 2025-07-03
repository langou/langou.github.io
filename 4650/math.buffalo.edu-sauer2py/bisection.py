# bisection.py (note: name "bisect.py" conflicts with some other package)
# Python/NumPy translation of Sauer's
# Program 1.1 Bisection Method
# Computes approximate solution of f(x)=0
# Input: function f; a,b such that f(a)*f(b)<0, 
#      and tolerance tol
# Output: Approximate solution xc

def sign(x):
	if x < 0:
		return -1
	elif x > 0:
		return 1
	else:
		return 0

def bisect(f,a,b,tol):
	fa = f(a)
	fb = f(b)
	if sign(fa)*sign(fb) >= 0:
		print 'f(a)f(b)<0 not satisfied!'
		quit()
	while (b-a)/2. > tol:
		c = (a+b)/2.
		fc = f(c)
		if fc == 0:						# c is a solution, done
			return c
		if sign(fc)*sign(fa) < 0:	# a and c make the new interval
			b = c
			fb = fc
		else:								# c and b make the new interval
			a = c
			fa = fc
	return (a+b)/2.					# new midpoint is best estimate


'''
# test_bisect.py
from bisection import bisect
from math import cos
def f(x):
	return cos(x)-x
tol = 1.0e-5
print bisect(f,0,1,tol)
'''


# Program 7.1 Nonlinear Finite Difference Method for BVP
# Translated to Python by JR 2/18/2012
# Uses Multivariate Newton's Method to solve nonlinear equation
# Inputs: interval inter, boundary values bv, number of steps n
# Output: solution w
# Example usage: w = nlbvpfd([0,1],[1,4],40)

from numpy import *
from pylab import plot, show

def nlbvpfd(inter,bv,n):
	global h,ya,yb				# needed in f and jac functions
	[a,b] = inter; [ya,yb] = bv
	h = (b-a)/float(n+1)		# h is step size
	w = zeros(n)				# initialize solution array w
	for i in range(20):			# loop of Newton step
		w -= linalg.solve( jac(w), f(w) )
								# plot w with boundary data
	plot( linspace(a,b,n+2), hstack((ya,w,yb)), linewidth=2, alpha=0.55 )
	return w

def f(w):
	global h,ya,yb
	n = len(w); h2=h*h
	y = zeros(n)
	y[ 0] =     ya     - (2+h2)*w[ 0] + h2*w[ 0]**2 + w[1]
	y[-1] =     w[ -2] - (2+h2)*w[-1] + h2*w[-1]**2 + yb
	for i in range(1,n-1):
		y[i] =  w[i-1] - (2+h2)*w[ i] + h2*w[ i]**2 + w[i+1]
	return y

def jac(w):
	global h,ya,yb
	n = len(w); h2=h*h
	a = zeros((n,n))
	for i in range(n):
		a[i,i] = 2*h2*w[i] - 2 - h2
	for i in range(n-1):
		a[i,i+1] = 1
		a[i+1,i] = 1
	return a

w = nlbvpfd([0,1],[1,4],  3)
set_printoptions(precision=14); print w
w = nlbvpfd([0,1],[1,4],500)
show()


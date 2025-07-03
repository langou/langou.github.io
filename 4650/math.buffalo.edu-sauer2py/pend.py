# Program 6.3 Animation program for pendulum using IVP solver
# Translated to Python by JR 2/17/2012.
# Inputs: 	ab = [a,b] time interval,
# 			initial values ic = [angle,velocity], number of steps n
# Calls a one-step method such as trapstep
# Example usage: pend([0,10],[pi/2,0],.05)

from numpy import *
from pylab import *

def trapstep(t,x,h):
	#one step of the Trapezoid Method
	z1 = ydot(t,x)
	g = x + h*z1
	z2 = ydot(t+h,g)
	return x + h*(z1+z2)/2.

def ydot(t,x):
	g = 9.81; length = 1
	return array( [x[1], -(g/length)*sin(x[0])-0.25*x[1]] )

def pend(ab,ic,n):
	subplot(111,aspect='equal')
	h = float(ab[1]-ab[0])/n	# plot n points in total
	y = zeros((n+1,2)); t = empty(n+1)					
	y[0,:] = ic; t[0] = ab[0]	# enter initial conds in y
	ion(); figure(figsize=(5,5))
	Rp = 1.5;
	xbob = sin(y[0,0]); ybob = -cos(y[0,0]-pi/2)
	xrod = [0,xbob]; yrod = [0,ybob]
	plot([-Rp,Rp,Rp],[-Rp,-Rp,Rp],'w')
	rod, = plot(xrod,yrod,'b',linewidth=3)
	bob, = plot(xbob,ybob,'ro',markersize=15)
	for k in range(n):
		t[k+1] = t[k] + h
		y[k+1,:] = trapstep( t[k],y[k,:], h )
		xbob = cos(y[k+1,0]-pi/2); ybob = sin(y[k+1,0]-pi/2)
		xrod = [0,xbob]; yrod = [0,ybob]
		rod.set_xdata( xrod ); rod.set_ydata( yrod )
		bob.set_xdata( xbob ); bob.set_ydata( ybob )
		draw()
	tmp = raw_input('Enter to leave! ')

pend([0,10],[pi/2,-8],200)

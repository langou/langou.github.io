# Program 6.4 Plotting program for one-body problem 
# Translated to Python by JR 2/17/2012. 2e update 2/21/2012.
#  Inputs: inter=[a b] time interval, initial conditions
#  ic = [x0 vx0 y0 vy0], x position, x velocity, y pos, y vel,
#  number of steps n, p = steps per point plotted
#  Calls a one-step method such as trapstep.m
#  Example usage: orbit([0, 100],[0, 1, 2, 0],10000,5)

from numpy import *
from pylab import *

def eulerstep(t,x,h):
	# one step of the Euler method
	return x + h*ydot(t,x)

def ydot(t,x):
	m2 = 3.; g = 1; mg2 = m2*g; px2 = 0; py2 = 0
	px1 = x[0]; py1 = x[2]; vx1 = x[1]; vy1 = x[3]
	dist = sqrt((px2-px1)**2+(py2-py1)**2)
	z = zeros(4)
	z[0] = vx1; z[1] = mg2*(px2-px1)/dist**3
	z[2] = vy1; z[3] = mg2*(py2-py1)/dist**3
	return z

def orbit(inter,ic,n,p):
	h = (inter[1]-inter[0])/float(n)				# plot n points 
	[x0,vx0,y0,vy0] = ic							# grab initial conds
	y = zeros((p+1,4)); t = empty(p+1)
	y[0,:] = [x0, vx0, y0, vy0]; t[0] = inter[0]	# build y vector
	ion(); figure(figsize=(5,5))
	Rp = 5; 
	plot([-Rp,Rp,Rp],[-Rp,-Rp,Rp],'w')
	sun, = plot(0,0,'ro')
	tail,= plot(y[1:,0],y[1:,2],'k')
	sat, = plot(y[0,0],y[0,2],'bo')
	for k in range(n/p):
		for i in range(p):
			t[i+1] = t[i] + h
			y[i+1,:] = eulerstep( t[i],y[i,:], h )
		y[0,:] = y[p,:]; t[0] = t[p]
		tail.set_xdata( y[1:,0] )
		tail.set_ydata( y[1:,2] )
		sat.set_xdata( y[0,0] )
		sat.set_ydata( y[0,2] )
		draw()
	show()
	tmp=raw_input('Enter to continue ')

orbit([0, 18.5],[0, 1, 2, 1],1000,250)
"""
orbit([0, 100],[0, 1, 2, 0],10000,5)
"""


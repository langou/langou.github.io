# Program 6.6 Animation program for bridge using IVP solver
# Translated to Python by JR 2/17/2012. Update for 2e 2/24/12.
# Inputs: inter = time interval inter, 
#   ic = [y[0,0] y[0,1] y[0,2] y[0,3]],
#   number of steps n, p = steps per point plotted
# Calls a one-step method such as trapstep.m
# Example usage: tacoma([0,1000],[0,0,0.001,0],25000,3)

from numpy import *
from pylab import *
from time  import sleep

def tacoma(inter,ic,n,p):
	[a,b] = inter; h = float(b-a)/n 	# plot n points
	y = zeros((p+1,4)); t = empty(p+1)
	y[0,:] = ic							# enter initial conds in y
	t[0] = a; len = 6
	ion(); figure(figsize=(5,5)); Rp = 9; 
	plot([-Rp,Rp,Rp],[-Rp,-Rp,Rp],'w')
	c = len*cos(y[0,2]); s = len*sin(y[0,2])
	road,   = plot([-c, c],[-s-y[0,0], s-y[0,0]],linewidth=5) 
	lcable, = plot([-c,-c],[-s-y[0,0],8]) 
	rcable, = plot([ c, c],[ s-y[0,0],8]) 
	for k in range(n):
		for i in range(p):
			t[i+1] =  t[i]+h; 
			y[i+1,:] = trapstep(t[i],y[i,:],h); 
		y[0,:] = y[p,:]; t[0] = t[p]
		c = len*cos(y[0,2]); s = len*sin(y[0,2])
		road.set_xdata(  [-c, c]);   road.set_ydata([-s-y[0,0], s-y[0,0]])
		lcable.set_xdata([-c,-c]); lcable.set_ydata([-s-y[0,0],8])
		rcable.set_xdata([ c, c]); rcable.set_ydata([ s-y[0,0],8])
		draw(); #sleep(h)

def trapstep(t,x,h): #one step of the Trapezoid Method
	z1 = ydot(t,x)
	g = x + h*z1
	z2 = ydot(t+h,g)
	return x + h*(z1+z2)/2.

def ydot(t,y):
	leng = 6;  a = 0.2;  W = 80;  omega = 2*pi*38/60.
	a1 = exp(a*(y[0]-leng*sin(y[2])))
	a2 = exp(a*(y[0]+leng*sin(y[2])))
	ydot = empty(4)
	ydot[0] = y[1]
	ydot[1] = -0.01*y[1]-0.4*(a1+a2-2)/a+0.2*W*sin(omega*t)
	ydot[2] = y[3]
	ydot[3] = -0.01*y[3]+1.2*cos(y[2])*(a1-a2)/(leng*a)
	return ydot

tacoma([0,1000],[0,0,0.001,0],25000,3)


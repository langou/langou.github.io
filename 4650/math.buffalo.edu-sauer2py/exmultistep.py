# Program 6.7 Multistep method plotting program
# Translated to Python by JR 2/17/2012. 2e update 2/24/12.
# Inputs: [inter[0],inter[1]] time interval,
#  ic = [y0]  initial condition,
#  h = stepsize, s = number of (multi)steps, e.g. 2 for 2-step method
# Output: solution y
# Calls a multistep method such as ab2step.m
# Example usage: y = exmultistep([0,1],[1],20,2)

from numpy import *
from pylab import plot, show

def trapstep(t,x,h):				# one step of the Trapezoid Method
	z1 = ydot(t,x)
	g = x + h*z1
	z2 = ydot(t+h,g)
	return x + h*(z1+z2)/2.

def ab2step(t,i,y,f,h):				# one step of the Adams-Bashforth 2-step method
	return y[i,:]+h*(3*f[i,:]/2-f[i-1,:]/2)

def unstable2step(t,i,y,f,h): 		# one step of an unstable 2-step method
	return -y[i,:]+2*y[i-1,:]+h*(5*f[i,:]/2+f[i-1,:]/2)

def weaklystable2step(t,i,y,f,h): 	# one step of a  weakly-stable 2-step method
	return y[i-1,:]+h*2*f[i,:]

def exmultistep(inter,ic,n,s,multistepper=ab2step):
	h = float(inter[1]-inter[0])/n
	y = empty((n+1,len(ic))); t = empty(n+1)
	f = empty((n+1,len(ic)))
	y[0,:] = ic; t[0] = inter[0]
	for i in range(s-1):	# start-up phase, using one-step method
		t[i+1] = t[i]+h
		y[i+1,:] = trapstep(t[i],y[i,:],h)
		f[i,:] = ydot(t[i],y[i,:])
	for i in range(s-1,n):	# multistep method loop
		t[i+1] = t[i]+h
		f[i,:] = ydot(t[i],y[i,:])
		y[i+1,:] = multistepper(t[i],i,y,f,h)
	plot(t,y)
	return t,y

def ydot(t,y):  # IVP from section 6.1
	return t*y+t**3

set_printoptions(precision=4)#,linewidth=100)
t,y = exmultistep([0,1],[1],20,2)
print 't = ', t; 
print 'ab2step';           print 'y =', y.T
t,y = exmultistep([0,1],[1],20,2,unstable2step)
print 'unstable2step';     print 'y =', y.T
t,y = exmultistep([0,1],[1],20,2,weaklystable2step)
print 'weaklystable2step'; print 'y =', y.T
show()

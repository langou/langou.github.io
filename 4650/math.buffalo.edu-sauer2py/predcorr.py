# Program 6.8 Adams-Bashforth-Moulton second-order p-c
# Python update for 2e JR 2/24/2012
# Inputs: 
#  ydot(t,y) function that evaluates RHS of DE
#  [inter[0],inter[1]] time interval, 
#  ic = [y0]  initial condition
#  number of steps n, number of (multi)steps s for explicit method
# Output: time steps t, solution y
# Calls multistep methods such as ab2step.m and am1step.m
# Example usage: predcorr(ydot,[0,1],1,20,2)

from numpy import *

def predcorr(ydot,inter,ic,n,s):
	h = float(inter[1]-inter[0])/n
	t = empty(  n+1 )
	y = empty( (n+1,len(ic)) )
	f = empty( (n+1,len(ic)) )
	# Start-up phase
	y[0,:] = ic; t[0] = inter[0]
	for i in range(s-1):	# start-up phase, using one-step method
		t[i+1]   = t[i]+h
		y[i+1,:] = trapstep(ydot,t[i],y[i,:],h)
		f[i,:]   = ydot(t[i],y[i,:])
	for i in range(s-1,n):				# multistep method loop
		t[i+1]   = t[i]+h
		f[i,:]   = ydot(t[i],y[i,:])
		y[i+1,:] = ab2step(t[i],i,y,f,h)  # predict
		f[i+1,:] = ydot(t[i+1],y[i+1,:])
		y[i+1,:] = am1step(t[i],i,y,f,h)  # correct
	return t,y

def trapstep(ydot,t,x,h):
# one step of the Trapezoid Method from section 6.2
	z1 = ydot(t,x) 
	g = x+h*z1 
	z2 = ydot(t+h,g) 
	y = x+h*(z1+z2)/2
	return y

def ab2step(t,i,y,f,h):
# one step of the Adams-Bashforth 2-step method
	return y[i,:]+h*(3*f[i,:]-f[i-1,:])/2

def am1step(t,i,y,f,h):
# one step of the Adams-Moulton 1-step method
	return y[i,:]+h*(f[i+1,:]+f[i,:])/2

"""
def ydot(t,y):  # IVP
	return t*y + t**3
"""

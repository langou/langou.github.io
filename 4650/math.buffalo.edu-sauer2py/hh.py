# Program 6.5 Hodgkin-Huxley equations
# Inputs: time interval inter,
#  ic  =  initial voltage v and 3 gating variables, steps n
# Output: solution y
# Calls a one-step method such as rk4step.m
# Example usage: y = hh([0,100],[-65,0,0.3,0.6],2000)

from numpy import *
from pylab import *

def hh(inter,ic,n):
	global pa,pb,pulse
	inp = input('pulse start, end, muamps in [ ], e.g. [50,51,7]: ')
	[pa,pb,pulse] = inp
	[a,b] = inter; h = float(b-a)/n		# plot n points in total
	y = empty((n+1,4)); t = empty(n+1)
	y[0,:] = ic							# enter initial conditions in y
	t[0] = a
	for i in range(1,n):
		t[i+1] = t[i]+h
		y[i+1,:] = rk4step(t[i],y[i,:],h)
	subplot(3,1,1)
	plot([a,pa,pa,pb,pb,b],[0,0,pulse,pulse,0,0])
	grid(True); axis([0,100,0,2*pulse])
	ylabel('input pulse')
	subplot(3,1,2)
	plot(t,y[:,0]); grid(True); axis([0,100,-100,100])
	ylabel('voltage (mV)')
	subplot(3,1,3)
	plot(t,y[:,1],t,y[:,2],t,y[:,3]); grid(True); axis([0,100,0,1])
	xlabel('time (msec)'); ylabel('gating variables')
	legend(['m','n','h'])
	show()
	return y

def rk4step(t,w,h):
	#one step of the Runge-Kutta order 4 method
	s1 = ydot(t,w)
	s2 = ydot(t+h/2,w+h*s1/2.)
	s3 = ydot(t+h/2,w+h*s2/2.)
	s4 = ydot(t+h,  w+h*s3)
	return w + h*(s1+2*s2+2*s3+s4)/6.

def ydot(t,w):
	global pa,pb,pulse
	c = 1; g1 = 120; g2 = 36; g3 = 0.3; T = (pa+pb)/2.; leng = pb-pa
	e0 = -65; e1 = 50; e2 = -77; e3 = -54.4
	inp = pulse*(1-sign(abs(t-T)-leng/2.))/2.
	# square pulse input on interval [pa,pb] of pulse muamps
	[v,m,n,h] = w
	z = zeros(4)
	z[0] = (inp-g1*m*m*m*h*(v-e1)-g2*n*n*n*n*(v-e2)-g3*(v-e3))/c
	v  =  v-e0 
	z[1] = (1-m)*(2.5-0.1*v)/(exp(2.5-0.1*v)-1)-m*4*exp(-v/18)
	z[2] = (1-n)*(0.1-0.01*v)/(exp(1-0.1*v)-1)-n*0.125*exp(-v/80)
	z[3] = (1-h)*0.07*exp(-v/20) - h/(exp(3-0.1*v)+1)
	return z

y = hh([0,100],[-65,0,0.3,0.6],2000)



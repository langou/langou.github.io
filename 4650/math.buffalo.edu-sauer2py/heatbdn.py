# Program 8.3 Backward difference method for heat equation 
# with Neumann boundary conditions 
# Translation of heatbdn.m to Python/NumPy/Matplotlib JR 2/25/2012.
# Input: space interval [xl,xr], time interval [yb,yt],
#        number of space steps M, number of time steps N
#        function for initial conditions f
# Output: solution w
# Example usage: w=heatbdn(0,1,0,1,20,20,f)

from numpy import *

def heatbdn(xl,xr,yb,yt,M,N,f):
	D = 1.									# diffusion coefficient
	h = float(xr-xl)/M; k=float(yt-yb)/N; m=M+1; n=N
	sigma = D*k/(h*h)
	a  = diag(1+2*sigma*ones(m)) + diag(-sigma*ones(m-1),1) 
	a += diag(-sigma*ones(m-1),-1)			# define matrix a
	a[  0,:] = hstack([-3, 4, -1, zeros(m-3)]) 	# Neumann conditions
	a[m-1,:] = hstack([zeros(m-3), -1, 4, -3]) 
	xvals = linspace(xl,xr,M+1)
	tvals = linspace(yb,yt,N+1)	
	w = zeros((m,N+1))   					# 2nd index is time index
	w[:,0] = f(xvals)						# initial conditions
	for j in range(n):
		b = w[:,j].copy(); b[0]=0; b[m-1]=0
		w[:,j+1] = linalg.solve(a,b) 
	print w
	[x,t] = meshgrid(xvals,tvals)			# 3-D plot of solution w
	from mpl_toolkits.mplot3d import axes3d, Axes3D
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_wireframe(x, t, w.T, rstride=1, cstride=1)
	ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('w')
	plt.show()
	return w


'''
# test_heatbdn.py
from numpy import *
from heatbdn import heatbdn

def f(x):
	return sin(2*pi*x)**2

set_printoptions(precision=14,linewidth=250)

w = heatbdn(0,1,0,0.2,4,5,f)
print w
'''

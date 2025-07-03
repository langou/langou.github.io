# Program 8.1 Forward difference method for heat equation
# Translation of fdheat.m to Python/NumPy/Matplotlib
# Python update for 2e JR 2/24/2012.
# Input: space interval [xl,xr], time interval [yb,yt],
#        number of space steps M, number of time steps N,
#        functions for initial and L & R boundary conditions
# Output: solution w
# Example usage: w=heatfd(0,1,0,1,10,250,f,l,r)

from numpy import *

def heatfd(xl,xr,yb,yt,M,N,f,l,r):
	D = 1.							# diffusion coefficient
	h = float(xr-xl)/M; k=float(yt-yb)/N; m=M-1; n=N
	sigma = D*k/(h*h)
	a  = diag(1-2*sigma*ones(m)) + diag(sigma*ones(m-1),1) 
	a += diag(sigma*ones(m-1),-1)	# define matrix a
	x = linspace(xl,xr,M+1)
	t = linspace(yb,yt,N+1)	
	w = zeros((m,N+1))   			# 2nd index is time index
	w[:,0] = f(x[1:-1])				# initial conditions
	lside = l(t); rside = r(t)		# boundary conditions
	for j in range(n):
		w[:,j+1] = dot(a,w[:,j]) \
               + sigma*hstack((lside[j],zeros(m-2),rside[j]))
	w = vstack((lside,w,rside))
	print w.T
	[X,T] = meshgrid(x,t)			# 3-D plot of solution w
	from mpl_toolkits.mplot3d import axes3d, Axes3D
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_wireframe(X, T, w.T)
	# ax.plot_surface(X, T, w.T, cstride=1,rstride=50 )
	ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('w')
	plt.show()
	return w


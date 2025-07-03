# Program 13.3 Nelder-Mead Search
# Input: inline function f, best guess xbar (column vector),
#    initial search radius rad and number of steps k
# Output: matrix x whose columns are vertices of simplex,
#    function values y of those vertices
# Translation to Python of neldermead.m
from numpy import *
def neldermead(f,xbar,rad,k):
	xbarf = array(xbar,dtype=float) # in case input array is list or has int data type
	n = xbarf.shape[0]
	x = empty((n,n+1))
	x[:,0] = xbarf				# each column of x is a simplex vertex
	x[:,1:n+1] = xbar*ones((1,n))+rad*eye(n,n)
	y = empty(n+1)
	for j in range(n+1):
		y[j] = f( x[:,j] )   # evaluate obj function f at each vertex
	oy = argsort(y)			# sort the function values in ascending order
	y = y[oy]
	x = x[:,oy]					# and rank the vertices the same way
	for i in range(k):
		xbar = mean( x[:,0:n], axis=1)	# xbar is the centroid of the face
		xh = x[:,n].copy()					# omitting the worst vertex xh
		xr = 2*xbar - xh; yr = f(xr)
		if yr < y[n-1]:
			if yr < y[0]:		# try expansion xe
				xe = 3*xbar - 2*xh; ye = f(xe)
				if ye < yr:		# accept expansion
					x[:,n] = xe; y[n] = f(xe)
				else:				# accept reflection
					x[:,n] = xr; y[n] = f(xr)
			else:					# xr is middle of pack, accept reflection
				x[:,n] = xr; y[n] = f(xr)
		else:						# xr is still the worst vertex, contract
			if yr < y[n]:		# try outside contraction xoc
				xoc = 1.5*xbar - 0.5*xh; yoc = f(xoc)
				if yoc < yr:	# accept outside contraction
					x[:,n] = xoc; y[n] = f(xoc)
				else:				# shrink simplex toward best point
					for j in range(1,n+1):
						x[:,j] = 0.5*x[:,0] + 0.5*x[:,j]; y[j] = f(x[:,j])
			else:					# xr is even worse than the previous worst
				xic = 0.5*xbar + 0.5*xh; yic = f(xic)
				if yic < y[n]:	# accept inside contraction
					x[:,n] = xic; y[n] = yic
				else:				# shrink simplex toward best point
					for j in range(1,n+1):
						x[:,j] = 0.5*x[:,0]+0.5*x[:,j]; y[j] = f(x[:,j])
		oy = argsort(y)		# re-sort the function values in ascending order
		y = y[oy]
		x = x[:,oy]				# and rank the vertices the same way
		print i,y
	return x

'''
# test_neldermead.py

from neldermead import neldermead

def f(x):
	return sum((x-[2,3,4,5,6])**2)

print neldermead(f,[7,7,7,7,7],0.1,150)
'''


# Program 3.5 Calculation of spline coefficients
# Calculates coefficients of cubic spline
# Input: x,y vectors of data points
#   plus two optional extra data v1, vn
# Output: matrix of coefficients b1,c1,d1;b2,c2,d2;...

from numpy import *

def splinecoeff(xin,yin,option=1,v1=0,vn=0):
	x = array(xin,dtype=float)	# in case inputs are integer lists
	y = array(yin,dtype=float) 
	n = len(x); 
	A = zeros((n,n))			# matrix A is nxn
	r = zeros(n)
	dx = x[1:] - x[:-1]			# define the deltas
	dy = y[1:] - y[:-1]
	for i in range(1,n-1):			# load the A matrix
		A[i,i-1:i+2] = hstack( (dx[i-1], 2*(dx[i-1]+dx[i]), dx[i]) )
		r[i] = 3*(dy[i]/dx[i] - dy[i-1]/dx[i-1]) # right-hand side
# Set endpoint conditions
	if   option==1:					# natural spline conditions
		A[ 0, 0]  =  1.
		A[-1,-1]  =  1.
	elif option==2:					# curvature-adj conditions
		A[ 0, 0] = 2; r[ 0] = v1
		A[-1,-1] = 2; r[-1] = vn
	elif option==3:					# clamped
		A[ 0, :2] = [2*dx[  0],  dx[  0]]; r[ 0] = 3*(dy[  0]/dx[  0]-v1) 
		A[-1,-2:] = [  dx[n-2],2*dx[n-2]]; r[-1] = 3*(vn-dy[n-2]/dx[n-2])
	elif option==4:					# parabol-term conditions, for n> = 3
		A[ 0, :2] = [1,-1]
		A[-1,-2:] = [1,-1]
	elif option==5: 				# not-a-knot for n> = 4
		A[ 0, :3] = [dx[  1], -(dx[  0]+dx[  1]), dx[  0]]
		A[-1,-3:] = [dx[n-2], -(dx[n-3]+dx[n-2]), dx[n-3]]
	coeff = zeros((n,3))
	coeff[:,1] = linalg.solve(A,r)	# solve for c coefficients
	for i in range(n-1):			# solve for b and d
		coeff[i,2] = (coeff[i+1,1]-coeff[i,1])/(3.*dx[i])
		coeff[i,0] = dy[i]/dx[i]-dx[i]*(2.*coeff[i,1]+coeff[i+1,1])/3.
	coeff = coeff[0:n-1,:]
	print  coeff
	return coeff
"""
splinecoeff([1,2,4,7],[3,4,1,2])
splinecoeff([1,2,4,7],[3,4,1,2],2)
splinecoeff([1,2,4,7],[3,4,1,2],3,0,0)
splinecoeff([1,2,4,7],[3,4,1,2],4)
splinecoeff([1,2,4,7],[3,4,1,2],5)
"""

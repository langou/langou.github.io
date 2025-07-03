# Program 8.5 Finite difference solver for 2D Poisson equation 
# with Dirichlet boundary conditions on a rectangle 
# Python translation JR 2/25/2012.
# Input: rectangle domain [xl,xr]x[yb,yt] with MxN space steps 
# Output: matrix w holding solution values 
# Example usage: w = poisson(0,1,1,2,4,4) 
from numpy import *
from mymesh import mymesh
def poisson(xl,xr,yb,yt,M,N):
	def f(x,y): return 0 			# define input function data 
	def g1(x) : return log(x**2+1) 	# define boundary values 
	def g2(x) : return log(x**2+4)	# Example 8.8 is shown 
	def g3(y) : return 2*log(y)
	def g4(y) : return log(y**2+1) 
	m = M+1; n = N+1;  mn = m*n 
	h = float(xr-xl)/M; h2 = h**2; k = float(yt-yb)/N; k2 = k**2 
	x = linspace(xl,xr,m) # set mesh values 
	y = linspace(yb,yt,n)
	A = zeros((mn,mn)); b = zeros(mn) 
	for i in range(1,m-1):  # interior points 
		for j in range(1,n-1):
			A[i+j*m,i-1+j*m] = 1./h2
			A[i+j*m,i+1+j*m] = 1./h2
			A[i+j*m,i  +j*m] =-2./h2 - 2./k2
			A[i+j*m,i  +(j-1)*m] = 1./k2
			A[i+j*m,i  +(j+1)*m] = 1./k2
			b[i+j*m] = f(x[i],y[j]) 
	for i in range(m): 		# bottom and top boundary points 
		j = 0  ; A[i+j*m,i+j*m] = 1; b[i+j*m] = g1(x[i]) 
		j = n-1; A[i+j*m,i+j*m] = 1; b[i+j*m] = g2(x[i]) 
	for j in range(1,n-1):	# left and right boundary points 
		i = 0  ; A[i+j*m,i+j*m] = 1; b[i+j*m] = g3(y[j]) 
		i = m-1; A[i+j*m,i+j*m] = 1; b[i+j*m] = g4(y[j]) 
	v = linalg.solve(A,b)	# solve for solution in v labeling 
	w = reshape(v,(m,n),order='F') #translate from v to w
	print w
	mymesh(x,y,w.T,'x','y','w') 
	return w

set_printoptions(precision=15,linewidth=140)
poisson(0,1,1,2,3,4) 

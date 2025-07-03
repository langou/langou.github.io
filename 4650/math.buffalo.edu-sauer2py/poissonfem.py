# Program 8.6 Finite element solver for 2D PDE 
#   with Dirichlet boundary conditions on a rectangle 
# Python/NumPy translation JR 2/25/2012
# Input: rectangle domain [xl,xr]x[yb,yt] with MxN space steps 
# Output: matrix w holding solution values 
# Example usage: w = poissonfem(0,1,1,2,4,4) 

from numpy import *
from mymesh import mymesh
def poissonfem(xl,xr,yb,yt,M,N):
	def f(x,y): return 0 # define input function data 
	def r(x,y): return 0 
	def g1(x) : return log(x**2 + 1)	# define boundary values on bottom 
	def g2(x) : return log(x**2 + 4)	# top 
	def g3(y) : return 2*log(y) 		# left side 
	def g4(y) : return log(y**2 + 1)	# right side 
	m = M + 1; n = N + 1;  mn = m*n 
	h = float(xr - xl)/M; h2 = h**2; ho3 = h/3.
	k = float(yt - yb)/N; k2 = k**2; ko3 = k/3.; hk = h*k 
	x = linspace(xl,xr,m) # set mesh values 
	y = linspace(yb,yt,n)
	A = zeros((mn,mn)); b = zeros(mn) 
	for i in range(1,m-1):  # interior points 
		for j in range(1,n-1):
			rsum  = r(x[i]-2*ho3,y[j]-ko3) + r(x[i]-ho3,y[j]-2*ko3) + r(x[i]+ho3,y[j]-ko3) 
			rsum += r(x[i]+2*ho3,y[j]+ko3) + r(x[i]+ho3,y[j]+2*ko3) + r(x[i]-ho3,y[j]+ko3) 
			A[i+j*m,i  +    j*m] = 2*(h2 + k2)/hk - hk*rsum/18. 
			A[i+j*m,i-1+    j*m] = -k/h - hk*(r(x[i]-  ho3,y[j]+  ko3) + r(x[i]-2*ho3,y[j]-  ko3))/18. 
			A[i+j*m,i-1+(j-1)*m] =      - hk*(r(x[i]-2*ho3,y[j]-  ko3) + r(x[i]-  ho3,y[j]-2*ko3))/18. 
			A[i+j*m,i  +(j-1)*m] = -h/k - hk*(r(x[i]-  ho3,y[j]-2*ko3) + r(x[i]+  ho3,y[j]-  ko3))/18. 
			A[i+j*m,i+1+    j*m] = -k/h - hk*(r(x[i]+  ho3,y[j]-  ko3) + r(x[i]+2*ho3,y[j]+  ko3))/18. 
			A[i+j*m,i+1+(j+1)*m] =      - hk*(r(x[i]+2*ho3,y[j]+  ko3) + r(x[i]+  ho3,y[j]+2*ko3))/18. 
			A[i+j*m,i  +(j+1)*m] = -h/k - hk*(r(x[i]+  ho3,y[j]+2*ko3) + r(x[i]-  ho3,y[j]+  ko3))/18. 
			fsum  = f(x[i]-2*ho3,y[j]-ko3) + f(x[i]-ho3,y[j]-2*ko3) + f(x[i]+ho3,y[j]-ko3) 
			fsum += f(x[i]+2*ho3,y[j]+ko3) + f(x[i]+ho3,y[j]+2*ko3) + f(x[i]-ho3,y[j]+ko3) 
			b[i+j*m] = -h*k*fsum/6 
	for i in range(m):		# boundary points 
		j =   0; A[i+j*m,i+j*m] = 1; b[i+j*m] = g1(x[i]) 
		j = n-1; A[i+j*m,i+j*m] = 1; b[i+j*m] = g2(x[i]) 
	for j in range(1,n-1): 
		i =   0; A[i+j*m,i+j*m] = 1; b[i+j*m] = g3(y[j]) 
		i = m-1; A[i+j*m,i+j*m] = 1; b[i+j*m] = g4(y[j]) 
	v = linalg.solve(A,b)	# solve for solution in v numbering 
	w = reshape(v,(m,n),order='F')
	print w 
	mymesh(x,y,w.T,'x','y','w') 
	return w

set_printoptions(precision=15,linewidth=150)
w = poissonfem(0,1,1,2,5,4) 


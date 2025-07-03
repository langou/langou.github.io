# rqi.py
# Python translation of Sauer Matlab, JR 1/19/12. Updated for 2e 2/24/12.
# (Note this can readily fail if A has an integer eigenvalue.)
# Program 12.3 Rayleigh Quotient Iteration
# Input: square numpy array A, (nonzero) 1d numpy array x, shift s, steps k
# Output: eigenvalue lam and eigenvector x

from numpy import *

def rqi( A, x, k ):
	I = eye(A.shape[0])
	for j in range(k):
		u = x/linalg.norm(x)		# normalize
		lam = dot(u,dot(A,u))		# Rayleigh quotient
#		print j,u, lam
#		print A-lam*I
		x = linalg.solve(A-lam*I,u)	# inverse power iteration
	u = x/linalg.norm(x)
	lam = dot(u,dot(A,u))
	return lam,x/linalg.norm(x,inf)

# try it out
A = array([[11.,-12,-6],[5,-5,-4],[-1,0,3]])
for i in range(10):
	x = random.rand(3)-0.5; lam,v = rqi( A, x, 4 ); print lam, v


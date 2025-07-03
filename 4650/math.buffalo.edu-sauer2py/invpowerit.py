# Program 12.2 Inverse Power Iteration
# Python translation by JR 1/19/12
# Computes eigenvalue of square matrix nearest to input s
# Input: square numpy array A, (nonzero) 1d numpy array x, shift s, steps k
# Output: eigenvalue lam nearest to s, eigenvector 

from numpy import *

def invpowerit( A, x, s, k ):
	As = A - s*eye(A.shape[0])
	for i in range(k):
		u = x/linalg.norm(x)
		x = linalg.solve(As,u)
		lam = dot(u,x)
	return 1./lam + s, x/linalg.norm(x,inf)

# try it
A = diag([1,2,3,4]) + 0.01*random.rand(4,4)
x = ones(4)
lam,u = invpowerit( A, x, 1., 25 ); print lam, u
lam,u = invpowerit( A, x, 2., 25 ); print lam, u
lam,u = invpowerit( A, x, 3., 25 ); print lam, u
lam,u = invpowerit( A, x, 4., 25 ); print lam, u




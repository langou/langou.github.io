# Program 12.4 Normalized Simultaneous Iteration
# Translated to Python/NumPy by JR 2/18/2012
# Computes eigenvalues/vectors of symmetric matrix
# Input: matrix A, number of steps k
# Output: eigenvalues lam and eigenvector matrix Q

from numpy import *

def nsi(A,k):
	m = A.shape[0]
	Q = eye(m)
	for j in range(k):
		Q,R = linalg.qr( dot(A,Q) )	# QR factorization
	lam = diag( dot(Q.T,dot(A,Q)) )	# Rayleigh quotient
	return lam,Q

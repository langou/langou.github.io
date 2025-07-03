# Program 12.5 Unshifted QR Algorithm
# Translated to Python/NumPy by JR 2/18/2012
# Computes eigenvalues/vectors of symmetric matrix
# Input: matrix A, number of steps k
# Output: eigenvalues lam and eigenvector matrix Qbar

from numpy import *

def unshiftedqr(A,k):
	m = A.shape[0]
	Q = eye(m)
	Qbar = Q.copy(); R = A.copy()
	for j in range(k):
		Q,R = linalg.qr( dot(R,Q) )	# QR factorization
		Qbar = dot(Qbar,Q)
	lam = diag( dot(R,Q) )	# Rayleigh quotient
	return lam,Qbar

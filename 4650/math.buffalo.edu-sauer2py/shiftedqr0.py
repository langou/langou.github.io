# Program 12.6 Shifted QR Algorithm, preliminary version
# Translated to Python/NumPy by JR 2/18/2012
# Input: matrix a
# Output: eigenvalues lam
# Computes eigenvalues of matrices without equal magnitude eigenvalues

from numpy import *

def shiftedqr0(a):
	tol = 1.e-14
	m = a.shape[0]; lam = zeros(m)
	n = m-1
	while n>0:
		while max(abs(a[n,0:n])) > tol:
			mu = a[n,n];				# define shift mu
			q,r = linalg.qr( a - mu*eye(n+1) )
			a = dot(r,q) + mu*eye(n+1)
		lam[n] = a[n,n]					# declare eigenvalue
		a = a[:n,:n]					# deflate
		n -= 1							# decrement n
	lam[0] = a[0,0]						# 1x1 matrix remains
	return lam

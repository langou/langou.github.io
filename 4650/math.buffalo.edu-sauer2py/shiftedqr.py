# Program 12.7 Shifted QR Algorithm, general version
# Computes real and complex eigenvalues of square matrix
# Translated to Python/NumPy by JR 2/18/2012
# Input: matrix a
# Output: eigenvalues lam

from numpy import *

def shiftedqr(a):
	tol = 1.e-14; kounttol = 500
	m = a.shape[0]; lam = zeros(m,dtype=complex)
	n = m-1
	while n>0:
		kount = 0
		while max(abs(a[n,0:n])) > tol and kount < kounttol:
			kount += 1     	# keep track of number of qr's
			mu = a[n,n];	# define shift mu
			q,r = linalg.qr( a - mu*eye(n+1) )
			a = dot(r,q) + mu*eye(n+1)
		if kount < kounttol:	# have isolated 1x1 block
			lam[n]=a[n,n]		# declare eigenvalue
			n -= 1				# decrement n
			a = a[:n+1,:n+1]	# deflate by 1
		else:					# have isolated 2x2 block
			disc = complex(( a[-2,-2] - a[-1,-1] )**2 + 4*a[-1,-2]*a[-2,-1])
			lam[n  ]=( a[-2,-2] + a[-1,-1] + sqrt(disc) )/2.
			lam[n-1]=( a[-2,-2] + a[-1,-1] - sqrt(disc) )/2.
			n -= 2
			a = a[:n+1,:n+1]	# deflate by 2
	if n==0: lam[0] = a[0,0]	# only a 1x1 block remains
	return lam

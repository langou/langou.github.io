# Program 12.8 Upper Hessenberg form
# Translated to Python/NumPy by JR 2/18/2012
# Input: matrix a
# Output: hessenberg form matrix a and reflectors v
# Usage: [a,v] = hessen(a) yields similar matrix a of
#    Hessenberg form and a matrix v whose columns hold
#    the v's defining the Householder reflectors.

from numpy import *

def hessen(a):
	[m,n]  =  a.shape
	print m,n
	v = zeros((m,m))
	eps = spacing(1)
	for k in range(m-2):
		x = a[k+1:m,k]
		vk = v[:m-k-1,k] 
		vk[:] = -sign(x[0]+eps)*linalg.norm(x)*eye(m-k-1,1)[:,0] - x
		vk /= linalg.norm( vk )
		vkvk = outer(vk,vk)
		a[k+1:m,k:m  ] -= 2*dot( vkvk, a[k+1:m,k:m] )
		a[:    ,k+1:m] -= 2*dot( a[:,k+1:m], vkvk   )
	return a,v

# Test it:
a  =  ones((4,6)) 
a[0,0]=2; a[1,1]=3; a[2,2]=7; a[3,3]=13; a[3,5]=23;
set_printoptions(precision=5,linewidth=120)
print a
[a,v] = hessen(a)
print 'ah =\n', a
print 'v =\n', v


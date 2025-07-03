# Program 7.2 Finite element solution of linear BVP
# Inputs: interval inter, boundary values bv, number of steps n
# Output: solution values c
# Example usage c = bvpfem( [0, 1], [1, 3], 9 )

from numpy import ones, zeros
import scipy.sparse as sp
import scipy.sparse.linalg as spl

def bvpfem(inter,bv,n):
	a = inter[0]; b = inter[1]; ya = bv[0]; yb = bv[1];
	h = (b-a)/float(n+1);
	alpha = (8./3.)*h+2./h; beta = (2./3.)*h-1./h;
	e = ones(n);
	M = sp.dia_matrix( ([beta*e, alpha*e, beta*e],[-1,0,1]), shape=(n,n) ).tocsr()
	d = zeros(n);
	d[ 0] =  -ya*beta;
	d[-1] =  -yb*beta;
	c = spl.spsolve(M,d);
	return c

# test it:
# print bvpfem( [0, 1], [1, 3], 49 )

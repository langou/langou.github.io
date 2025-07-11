# Sauer Program 2.1 Sparse matrix setup
# Translated to Python. JR 1/17/12

from numpy import *
from scipy.sparse import *

def sparsesetup(n):
	e = ones(n)
	n2 = n/2
	a = dia_matrix( ([-e,3*e,-e],[-1,0,1]), shape=(n,n) ).tocsr()
	c = coo_matrix( (e/2, (range(n),range(n-1,-1,-1))), shape=(n,n)).tocsr()
	a = a + c
	a[n2  ,n2-1] = -1; 
	a[n2-1,n2  ] = -1; 
	b = zeros(n)
	b[0]=2.5; b[-1]=2.5; b[1:n-1]=1.5; b[n2-1:n2+1]=1.0 
	return a,b

a,b = sparsesetup(8)
print a
print a.todense()
print b






















'''
Original Matlab:
% Program 2.1 Sparse matrix setup
% Input: n = size of system
% Outputs: sparse matrix a, r.h.s. b
function [a,b] = sparsesetup(n)
e = ones(n,1);
n2=n/2;
a = spdiags([-e 3*e -e],-1:1,n,n);    % Entries of a
c=spdiags([e/2],0,n,n);c=fliplr(c);a=a+c;
a(n2+1,n2) = -1; a(n2,n2+1) = -1;     % Fix up 2 entries
b=zeros(n,1);                         % Entries of r.h.s. b
b(1)=2.5;b(n)=2.5;b(2:n-1)=1.5;b(n2:n2+1)=1;
'''


# Sauer Program 2.2  Jacobi Method
# Translated to Python. JR 1/17/12
# Inputs: full or sparse matrix a, r.h.s. b,
#         number of Jacobi iterations k
# Output: solution x
from numpy import *
from sparsesetup import sparsesetup

def jacobi(a,b,k):
	n = len(b)			
	d = [a[i,i] for i in range(n)]	# extract diagonal of a (diag won't work for sparse matrix)
	r = a - diag(d)			# r is the remainder
	x = zeros((n,1))			# initialize vector x
	bb = reshape(b,(n,1))
	dd = reshape(d,(n,1))
	for j in range(k):		# loop for Jacobi iteration
		x=(bb-dot(r,x))/dd
		# print x.T
	return x

a,b = sparsesetup(6)		# test with sparse matrix
print jacobi(a,b,6).T

da = a.todense()
print jacobi(da,b,6).T	# test with regular dense matrix


'''
Original Matlab:
% Program 2.2  Jacobi Method
% Inputs: full or sparse matrix a, r.h.s. b,
%         number of Jacobi iterations k
% Output: solution x
function x = jacobi(a,b,k)
n=length(b);      % find n
d=diag(a);        % extract diagonal of a
r=a-diag(d);      % r is the remainder
x=zeros(n,1);     % initialize vector x
for j=1:k         % loop for Jacobi iteration
  x=(b-r*x)./d;
end               % End of Jacobi iteration loop
'''

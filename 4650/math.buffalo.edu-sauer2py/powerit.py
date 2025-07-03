# Program 12.1 Power Iteration
# Python translation JR 1/18/12. 2e update 2/21/12.
# Computes dominant eigenvector of square matrix
# Input: matrix A, initial (nonzero) vector x, number of steps k
# Output: dominant eigenvalue, eigenvector 

from numpy import *

def powerit(A,x,k):
	for j in range(k):
		u = x/linalg.norm(x)	# normalize vector
		x = dot(A,u)			# power step
		lam = dot(u,x)			# Rayleigh quotient
		#print lam,u
	return lam,x/linalg.norm(x,inf)

# example
A = array([[1.,3.],[2.,2.]])
x = random.rand(2)
lam,v = powerit(A,x,20)
print "Dominant eigenvalue:",lam
print "Corresponding eigenvector:",v


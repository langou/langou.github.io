# Program 9.1  Quasi-random number generator 
# Halton sequence in base p
# Input: prime number p, random numbers required n
# Output: array u of quasi-random numbers in [0,1]
# Example usage: halton(2,100)
from numpy import *
def halton(p,n):
	eps = finfo(double).eps
	b = zeros( ceil(log(n)/log(p)) )   # largest number of digits
	u = empty( n )
	for j in range(n):
		i = 0
		b[0] += 1                 		  # add one to current integer
		while b[i] > p-1+eps:           # this loop does carrying
			b[i] = 0                     #    in  base p
			i = i+1
			b[i] += 1
		u[j] = 0
		for k in range(len(b)):            # add up reversed digits
			u[j] += b[k]*p**-(k+1)
	return u

"""
# test it

u = halton(2,100)
print u
from pylab import plot,show
plot( u[:-1],u[1:],'o' )
show()
"""

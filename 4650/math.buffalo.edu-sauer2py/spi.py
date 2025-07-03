# Program 13.2 Successive Parabolic Interpolation
# Input: inline function f, initial guesses r,s,t, number of steps k
# Output: approximate minimum x
# Translation to Python of spi.m
#  Storing and returning all the iterates.
#  Will crash if iterates enough times to get denominator zero.
def spi(f,r,s,t,k):
	x = [r,s,t]
	fr=f(r); fs=f(s); ft=f(t)
	for i in range(3,k+3):
		x.append((r+s)/2.-(fs-fr)*(t-r)*(t-s)/(2.*((s-r)*(ft-fs)-(fs-fr)*(t-s))))
		t=s; s=r; r=x[-1]
		ft=fs; fs=fr; fr=f(r)		# single function evaluation
	return x








































'''
# test_spi.py
from spi import spi

def f(x):
	return (x-3.7)**2 + (x-3.7)**3

print spi(f,3,4,5,20)
'''

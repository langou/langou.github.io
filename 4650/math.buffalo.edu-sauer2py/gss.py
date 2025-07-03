# Program 13.1 Golden Section Search for minimum of f(x)
# Start with unimodal f(x) and minimum in [a,b]
# Input: inline function f, interval [a,b], number of steps k
# Output: approximate minimum y
# Translation to Python of gss.m
from math import sqrt
def gss(f,a,b,k):
	g=(sqrt(5.)-1.)/2.
	x1 = a+(1-g)*(b-a)
	x2 = a+g*(b-a)
	f1=f(x1); f2=f(x2)
	for i in range(k):
	  if f1 < f2:				# if f(x1) < f(x2), replace b with x2
		 b=x2; x2=x1; x1=a+(1-g)*(b-a)
		 f2=f1; f1=f(x1)		# single function evaluation
	  else:						# otherwise, replace a with x1
		 a=x1; x1=x2; x2=a+g*(b-a)
		 f1=f2; f2=f(x2)		# single function evaluation
	return (a+b)/2.


'''
# test_gss.py
from gss import gss

def f(x):
	return (x-3.7)**2

print gss(f,0,5,30)
'''


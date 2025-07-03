# Program 5.2 Adaptive Quadrature
# Computes approximation to definite integral
# Inputs: inline function f, interval [a0,b0], 
#         error tolerance tol0
# Output: approximate definite integral
# Translated to Python 2/15/12 by JR
# This is a less literal translation than most of the other codes,
# in order to work with lists whose maximum length is unknown.
# But the arithmetic performed should be identical to that of the Matlab version.

def trap(f,a,b):
	s = (f(a)+f(b))*(b-a)/2.
	return s

def adapquad(f,a0,b0,tol0):
	sum = 0.; n = 0; alist = [a0]; blist = [b0]; tols = [tol0]; apps = [trap(f,a0,b0)]
	while len(alist) > 0:		# n is index of lalistt element of lists
		if n>15:
			quit()
		a = alist.pop(); b = blist.pop(); oldtol = tols.pop(); oldapp = apps.pop()
		c = (a+b)/2.
		leftapp = trap(f,a,c); rightapp = trap(f,c,b)
		errest = oldapp - (leftapp + rightapp)
		if abs(errest) < 3.*oldtol:
			sum  += leftapp+rightapp		# success
		else:
			alist.append(a); blist.append(c);		# divide into two intervals
			alist.append(c); blist.append(b);		
			tols.append(oldtol/2.); tols.append(oldtol/2.); 
			apps.append(leftapp); apps.append(rightapp); 
	return sum



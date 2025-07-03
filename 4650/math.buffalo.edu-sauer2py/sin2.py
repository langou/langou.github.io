# sin2.py
# Translation to Python/Numpy JR 1/22/12
# Program 3.4 Building a sin calculator key, attempt #2
# Approximates sin curve with degree 9 polynomial
# Input:  x
# Output: approximation for sin(x), correct to 10 decimal places

from numpy  import *
from nest   import nest
from newtdd import newtdd

def sin2(xin):
	x = array(xin)  # in case input is plain list
	# First calculate the interpolating polynomial and 
	# store coefficients
	n = 10
	b = pi/4.*(1. + cos(array(range(1,2*n,2))*pi/2/n)) 
	# b holds Chebyshev base points
	yb = sin(b)
	c = newtdd( b, yb )
	# For each input x, move x to the fundamental domain 
	# and evaluate the interpolating polynomial
	x1 = x % (2*pi)
	s = ones_like(x)       # Correct the sign of sin
	s [x1>pi] = -1
	x1[x1>pi] = 2*pi - x1[x1>pi]
	x1[x1>pi/2] = pi - x1[x1>pi/2]
	y = s*nest(c,x1,b)
	return y

# try it out
x =[17.0]
print sin2(x), sin(x), sin2(x)-sin(x)
x = linspace(0,10,500)
yinterp = sin2(x)
ytrue   = sin(x)
from pylab import plot, show
plot(x,yinterp,'r')
plot(x,ytrue,'b')
show()









'''
Original Matlab
n=10;
b=pi/4+(pi/4)*cos((1:2:2*n-1)*pi/(2*n));
yb=sin(b);                      % b holds Chebyshev base points
c=newtdd(b,yb,n);
%For each input x, move x to the fundamental domain and evaluate
%   the interpolating polynomial
s=1;                            % Correct the sign of sin
x1=mod(x,2*pi);
if x1>pi
  x1 = 2*pi-x1;
  s = -1;
end
if x1 > pi/2
  x1 = pi-x1;
end
y = s*nest(n-1,c,x1,b);
'''


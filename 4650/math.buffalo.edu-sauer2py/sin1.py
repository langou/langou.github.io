# sin1.py
# Translation to Python/Numpy JR 1/22/12
# Program 3.3 Building a sin calculator key, attempt #1
# Approximates sin curve with degree 3 polynomial
#     (Caution: do not use to build bridges,
#     at least until we have discussed accuracy.)
# Input:  x - a list or 1D NumPy array of evaluation points
# Output: approximation for sin(x)

from numpy  import *
from nest   import nest
from newtdd import newtdd

def sin1(xin):
	x = array(xin)  # in case input is plain list
	# First calculate the interpolating polynomial and 
	# store coefficients
	b = linspace(0,pi/2,4,endpoint=False)   # b holds base points
	yb = sin(b)
	c = newtdd( b, yb )
	# For each input x, move x to the fundamental domain and evaluate
	# the interpolating polynomial
	s = ones_like(x)       # Correct the sign of sin
	x1 = x % (2*pi)
	s [x1>pi] = -1
	x1[x1>pi] = 2*pi - x1[x1>pi]
	x1[x1>pi/2] = pi - x1[x1>pi/2]
	y = s*nest(c,x1,b)
	return y

# try it out
x =[17.0]
print sin1(x), sin(x), sin1(x)-sin(x)
x = linspace(0,10,500)
yinterp = sin1(x)
ytrue   = sin(x)
from pylab import plot, show
plot(x,yinterp,'r')
plot(x,ytrue,'b')
show()


'''
Original Matlab
% Program 3.3 Building a sin calculator key, attempt #1
% Approximates sin curve with degree 3 polynomial
%     (Caution: do not use to build bridges,
%     at least until we have discussed accuracy.)
% Input:  x
% Output: approximation for sin(x)
function y=sin1(x)
% First calculate the interpolating polynomial and 
%   store coefficients
b=pi*(0:3)/6;yb=sin(b);    % b holds base points
c=newtdd(b,yb,4);
%For each input x, move x to the fundamental domain and evaluate
%      the interpolating polynomial
s=1;                       % Correct the sign of sin
x1=mod(x,2*pi);
if x1>pi
  x1 = 2*pi-x1;
  s = -1;
end
if x1 > pi/2
  x1 = pi-x1;
end
y = s*nest(3,c,x1,b);
'''


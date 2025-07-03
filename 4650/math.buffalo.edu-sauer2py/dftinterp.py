'''
Original Matlab code:
%Program 10.1 Fourier interpolation
%Interpolate n data points on [c,d] with trig function P(t)
%  and plot interpolant at p (>= n) evenly spaced points.
%Input: interval [c,d], data points x,
%  even number of data points n, even number p>=n
%Output: data points of interpolant xp
function xp=dftinterp(inter,x,n,p)
c=inter(1);d=inter(2);
t=c+(d-c)*(0:n-1)/n;      % n evenly-spaced time points
tp=c+(d-c)*(0:p-1)/p;     % p evenly-spaced time points 
y=fft(x);                 % apply DFT
yp=zeros(p,1);            % yp will hold coefficients for ifft
yp(1:n/2+1)=y(1:n/2+1);   % move n frequencies from n to p
yp(p-n/2+2:p)=y(n/2+2:n); % same for upper tier
xp=real(ifft(yp))*(p/n);  % invert fft to recover data 
plot(t,x,'o',tp,xp)       % plot data points and interpolant
'''

# dftinterp.py
# Python/Numpy translation of Tim Sauer's Matlab, John Ringland 9/2/11.
# Program 10.1 Fourier interpolation
# Interpolate n data points on [c,d] with trig function P(t)
#   and plot interpolant at p (>= n) evenly spaced points.
# Input: interval [c,d], data points x,
#   even number of data points n, even number p>=n
# Output: data points of interpolant xp

from numpy import *
from scipy import fft, ifft
from pylab import plot, show

def dftinterp( inter, x, n, p ):
	t  = linspace( inter[0], inter[1], n, endpoint=False )  # n evenly-spaced time points
	tp = linspace( inter[0], inter[1], p, endpoint=False )  # p evenly-spaced time points
	y = fft(x)                                              # apply DFT
	yp = zeros(p,dtype=complex)                             # yp will hold coefficients for ifft
	yp[:n/2+1]   = y[:n/2+1]                                # move n frequencies from n to p
	yp[p-n/2+1:] = y[n/2+1:]                                # same for upper tier
	xp = real( ifft(yp) )*float(p)/n                        # invert fft to recover data 	
	plot( t, x, 'o', tp, xp )                               # plot data points and interpolant
	show()                                                   

dftinterp( [0,1], [-2.2, -2.8, -6.1, -3.9, 0.0, 1.1, -0.6, -1.1], 8, 100 ) 

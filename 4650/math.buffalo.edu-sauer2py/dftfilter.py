# Python translation of Sauer's
# Program 10.2 Least squares trigonometric fit
# Least squares fit of n data points on [c,d] with trig function
#    where 2 <= m <= n. Plot best fit at p (>= n) points
# Input: interval [c,d], data points x, even number m,
#  even number of data points n, even number p>=n
# Output: filtered points xp
from numpy import *
from scipy import fft, ifft
from pylab import plot, show

def dftfilter(inter,x,m,n,p):
	c = inter[0]; d = inter[1]
	t  = linspace(c,d,n,endpoint=False)	# time points for data (n)
	tp = linspace(c,d,p,endpoint=False)	# time points for interpolant (p)
	y = fft(x)									# compute interpolation coefficients
	yp = zeros( p, dtype=complex )		# yp will hold coefficients for ifft
	yp[:m/2] = y[:m/2]						# keep only first m frequencies 
	yp[ m/2] = real(y[m/2])					# since m is even, keep cos term only
	if m < n:                  			# unless at the maximum frequency,
	  yp[p-m/2] = yp[m/2]					# add complex conjugate to corresponding place in upper tier
	yp[p-m/2+1:p] = y[n-m/2+1:n]			# more conjugates for upper tier
	xp = real(ifft(yp))*(float(p)/n)		# invert fft to recover data
	plot(t,x,'o',tp,xp)						# plot data and least square approx
	show()
	return xp

'''
# test_dftfilter.py
from dftfilter import dftfilter
print dftfilter([0,1],[-2.2,-2.8,-6.1,-3.9,0.0,1.1,-0.6,-1.1],4,8,200)
'''

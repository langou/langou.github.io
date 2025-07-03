# Program 3.6 Cubic spline plot
# Python update for 2e JR 2/24/2012
# Plots spline from data points
# Input: x,y vectors of data points, k number of plotted points per segment
#        option, v1, vn, endpoint condition info
# Output: x1,y1 spline values at plotted points

from numpy import *
from pylab import plot,show, clf
from splinecoeff import splinecoeff

def splineplot(x,y,k,option=1,v1=0,vn=0):
	n = len(x)
	coeff = splinecoeff(x,y,option,v1,vn)
	x1 = empty((n-1)*k+1)
	y1 = empty((n-1)*k+1)
	for i in range(n-1):
		xs = linspace(x[i],x[i+1],k+1)
		dx = xs - x[i]
		ys = coeff[i,2]*dx	# evaluate using nested multiplication
		ys = (ys+coeff[i,1])*dx
		ys = (ys+coeff[i,0])*dx + y[i]
		#plot([x[i],x[i+1]],[y[i],y[i+1]],'o',xs,ys)
		plot(xs,ys)
		x1[i*k:(i+1)*k] = xs[:-1]
		y1[i*k:(i+1)*k] = ys[:-1]
	x1[-1] = x[-1]; y1[-1] = y[-1]
	plot(x,y,'o')#,x1,y1)
	return x1,y1

#splineplot([1,2,4,7],[3,4,1,2],100)		# natural spline conditions
#splineplot([1,2,4,7],[3,4,1,2],100,2)		# curvature-adj conditions
splineplot( [1,2,4,7],[3,4,1,2],100,3,0,0)	# clamped
#splineplot([1,2,4,7],[3,4,1,2],100,4)		# parabol-term conditions, for n> = 3
#splineplot([1,2,4,7],[3,4,1,2],100,5)		# not-a-knot for n> = 4
show()

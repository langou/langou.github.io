# Program 3.2.  Polynomial Interpolation Program
# Translated to Python by JR, 2/17/2012.
# Left click in figure window to locate data point.
#     Continue, to add more points.
#     Kill plot window to terminate program.

from numpy  import *
from pylab  import plot, ginput, grid, show, ylim
from newtdd import newtdd
from nest   import nest

def clickinterp():
	xl=-3; xr=3; yb=-3; yt=3
	plot([xl,xr],[0,0],'k',[0,0],[yb,yt],'k'); grid(True); ylim(yb,yt)
	x = linspace(xl,xr,300)			# define x coordinates of curve
	xlist = []; ylist = []
	pt = ginput(1)[0]; 				# get one mouse click
	xlist.append(pt[0]); ylist.append(pt[1])
	while(True):
		pt = ginput(1)[0]; 			# get one mouse click
		xlist.append(pt[0]); ylist.append(pt[1]); k=len(xlist)
		c = newtdd(xlist,ylist)		# get interpolation coeffs
		y = nest(c,x,xlist)		# get y coordinates of curve
		plot(xlist,ylist,'o',x,y,[xl,xr],[0,0],'k',[0,0],[yb,yt],'k')

clickinterp()


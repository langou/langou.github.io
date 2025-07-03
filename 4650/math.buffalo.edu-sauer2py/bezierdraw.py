# Program 3.7 to Python by JR 1/16/12 
# Translation of Sauer to Python by JR 1/16/2012
# Renamed from draw to bezierdraw 2/24/2012
# Freehand draw program using B\'ezier curves
#	Left click in figure window to locate endpoint, and click three
#	more times to specify 2 control points and the next spline endpoint.
#	Continue with groups of 3 points to add more to
#	the curve. Kill the window to terminate the program.

from numpy import *
from pylab import plot, ginput, show, subplot

def bezierdraw():
	subplot(111,aspect='equal')
	plot([ 0,0],[-1,1],'k')
	plot([-1,1],[ 0,0],'k')
	points = ginput(1)
	t=linspace(0.,1.,50,endpoint=True)
	while( True ):
		points += ginput(3) 	# get 3 new points
		x = [pt[0] for pt in points[-4:]]
		y = [pt[1] for pt in points[-4:]]
		plot(x[0:2],y[0:2],'c' ,x[1],y[1],'cs')
		plot(x[2:4],y[2:4],'c' ,x[2],y[2],'cs')
		plot(x[0]  ,y[0]  ,'bo',x[3],y[3],'bo')
		bx=3*(x[1]-x[0]);    by=3*(y[1]-y[0])	# spline equations ...
		cx=3*(x[2]-x[1])-bx; cy=3*(y[2]-y[1])-by
		dx=x[3]-x[0]-bx-cx;  dy=y[3]-y[0]-by-cy
		xp=x[0]+t*(bx+t*(cx+t*dx))				# Horner's method
		yp=y[0]+t*(by+t*(cy+t*dy))
		plot(xp,yp)                				# Plot spline curve and 
	show()

bezierdraw()



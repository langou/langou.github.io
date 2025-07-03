# bounce.py
# Translated to Python by JR 2/17/2012.
# Illustrates Matlab animation using the drawnow command
# Usage: Type "python bounce.py"

from pylab import *
from time import sleep

hx0=.005; hy0=.0039; hx=hx0; hy=hy0; 
xl=.02; xr=.98; yb=xl; yt=xr; x=.1; y=.1; 
ion(); figure(figsize=(5,5))
Rp = 1.0; 
plot([0,Rp,Rp,0],[0,0,Rp,Rp],'k')
ball, = plot( x, y, 'ro', markersize=20 )

while True:
	if x < xl:
		hx = hx0 
	if x > xr:
		hx = -hx0 
	if y < yb:
		hy = hy0
	if y > yt:
		hy = -hy0 
	x += hx; y += hy
	ball.set_xdata(x); ball.set_ydata(y)
	draw(); #sleep(0.01)


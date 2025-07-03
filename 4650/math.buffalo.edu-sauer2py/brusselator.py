# Program 8.9 Backward difference method with Newton iteration
#             for the Brusselator 
# Translation to Python/NumPy, JR 2/26/2012.
# input: space region [xl,xr]x[yb,yt], time interval [t0,te], 
#    M,N space steps in x and y directions, tsteps time steps 
# output: solution mesh [x,y,p,q]
# Example usage: [x,y,p,q] = brusselator(0,40,0,40,0,20,40,40,20)

from numpy import *
from pylab import ion, figure, plot, contour, draw
from matplotlib.pyplot import cla

def brusselator(xl,xr,yb,yt,tb,te,M,N,tsteps):
	Dp = 1.; Dq = 8.; C = 4.5; K = 9.
#	def fp(X,Y): return (C   + 0.1)*ones_like(X)
#	def fq(X,Y): return (K/C + 0.2)*ones_like(X)
	def fp(X,Y): return (C   + 0.1)*(1+sin(X))
	def fq(X,Y): return (K/C + 0.2)*(1+cos(Y))
	delt = float(te-tb)/tsteps
	m = M+1; n = N+1; mn = m*n; mn2 = 2*mn
	h = float(xr-xl)/M; k = float(yt-yb)/N; h2 = h**2; k2 = k**2
	x = linspace(xl,xr,m); y = linspace(yb,yt,n)
	Y,X = meshgrid(y,x)		# Careful: this way to get 1st index for x
	p = fp(X,Y)          	# Define initial conditions
	q = fq(X,Y)
	ion()
	figure(figsize=(7,7))
	plot([xl,xr,xr],[yb,yb,yt],'r')
	for tstep in range(tsteps):
		v = hstack([reshape(p,mn,order='F'),reshape(q,mn,order='F')])
		pold = p.copy(); qold = q.copy()
		for it in range(3):
			DF1 = zeros((mn2,mn2)); DF3 = zeros((mn2,mn2))
			# Note: for efficiency, take the DF1 code outside the it and tstep loops
			b = zeros(mn2)
			for i in range(1,m-1):
				for j in range(1,n-1):
					ipjm = i+j*m; mnpipjm = mn + ipjm
					DF1[    ipjm,    ipjm-1 ] = -Dp/h2
					DF1[    ipjm,    ipjm   ] =  Dp*(2./h2+2./k2) + K + 1. + 1./(1.*delt)
					DF1[    ipjm,    ipjm+1 ] = -Dp/h2
					DF1[    ipjm,    ipjm-m ] = -Dp/k2
					DF1[    ipjm,    ipjm+m ] = -Dp/k2
					DF1[ mnpipjm, mnpipjm-1 ] = -Dq/h2
					DF1[ mnpipjm, mnpipjm   ] =  Dq*(2./h2+2./k2)          + 1./(1.*delt)
					DF1[ mnpipjm, mnpipjm+1 ] = -Dq/h2
					DF1[ mnpipjm, mnpipjm-m ] = -Dq/k2
					DF1[ mnpipjm, mnpipjm+m ] = -Dq/k2
					DF1[ mnpipjm,    ipjm   ] = -K
					DF3[    ipjm,    ipjm   ] = -2.*p[i,j]*q[i,j]
					DF3[    ipjm, mnpipjm   ] = -p[i,j]**2
					DF3[ mnpipjm,    ipjm   ] =  2.*p[i,j]*q[i,j]
					DF3[ mnpipjm, mnpipjm   ] =  p[i,j]**2
					b[   ipjm] = -pold[i,j]/(1.*delt) - C      
					b[mnpipjm] = -qold[i,j]/(1.*delt) 

			for i in range(m):   	# bottom and top Neumann conditions
				j = 0;  ipjm = i + j*m
				DF1[    ipjm,    ipjm     ] =  3.
				DF1[    ipjm,    ipjm+  m ] = -4.
				DF1[    ipjm,    ipjm+2*m ] =  1.
				j = n-1; ipjm = i + j*m; 
				DF1[    ipjm,    ipjm     ] =  3.
				DF1[    ipjm,    ipjm-  m ] = -4.
				DF1[    ipjm,    ipjm-2*m ] =  1.
				j = 0;   mnpipjm = mn + i + j*m
				DF1[ mnpipjm, mnpipjm     ] =  3.
				DF1[ mnpipjm, mnpipjm+  m ] = -4.
				DF1[ mnpipjm, mnpipjm+2*m ] =  1.
				j = n-1; mnpipjm = mn + i + j*m
				DF1[ mnpipjm, mnpipjm     ] =  3.
				DF1[ mnpipjm, mnpipjm-  m ] = -4.
				DF1[ mnpipjm, mnpipjm-2*m ] =  1.

			for j in range(1,n-1):   #left and right
				i = 0;  ipjm = i + j*m
				DF1[    ipjm,    ipjm   ] =  3.
				DF1[    ipjm,    ipjm+1 ] = -4.
				DF1[    ipjm,    ipjm+2 ] =  1.
				i = m-1;ipjm = i + j*m; 
				DF1[    ipjm,    ipjm   ] =  3.
				DF1[    ipjm,    ipjm-1 ] = -4.
				DF1[    ipjm,    ipjm-2 ] =  1.
				i = 0;   mnpipjm = mn + i + j*m
				DF1[ mnpipjm, mnpipjm   ] =  3.
				DF1[ mnpipjm, mnpipjm+1 ] = -4.
				DF1[ mnpipjm, mnpipjm+2 ] =  1.
				i = m-1; mnpipjm = mn + i + j*m
				DF1[ mnpipjm, mnpipjm   ] =  3.
				DF1[ mnpipjm, mnpipjm-1 ] = -4.
				DF1[ mnpipjm, mnpipjm-2 ] =  1.
			DF = DF1 + DF3

			F = dot(DF1+DF3/3.,v) + b
			v -= linalg.solve(DF,F)
			p = reshape(v[:mn],(m,n),order='F')
			q = reshape(v[mn:],(m,n),order='F')
		cla(); contour(x,y,p.T); draw()
	return x,y,p,q

set_printoptions(precision=15,linewidth=150)
x,y,p,q = brusselator(0,40,0,30,0,20,3,4,1)
print 'x\n',x; print 'y\n',y; print 'p\n',p; print 'q\n',q; raw_input(" ")



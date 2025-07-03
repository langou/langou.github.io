# Program 8.7 Implicit Newton solver for Burgers equation 
# Python translation by JR 2/24/2012.
# input: space interval [xl,xr], time interval [tb,te], 
#        number of space steps M, number of time steps N 
# output: solution w 
# Example usage: w = burgers(0,1,0,1,10,10) 
from numpy import *

def burgers(xl,xr,tb,te,M,N):
	alf = 5; bet = 4; D = .05 
	def f(x): return 2*D*bet*pi*sin(pi*x)/(alf+bet*cos(pi*x)) 
	def l(t): return 0*t
	def r(t): return 0*t
	t = linspace(tb,te,N+1); L = l(t); R = r(t)
	h = (xr-xl)/float(M);  k = (te-tb)/float(N);  m = M+1;  n = N 
	sigma = D*k/(h*h) 
	w = empty((m,n+1)) 		# row index space, column index time
	x = linspace(xl,xr,m)
	w[:,0] = f(x) 			# initial conditions
	w1 = w[:,0].copy(); print w1 
	for j in range(n):
		for it in range(3): # Newton iteration 
			DF1 =  diag(1+2*sigma*ones(m)) + diag(-sigma*ones(m-1),1) 
			DF1 += diag(-sigma*ones(m-1),-1) 
			DF2  = diag( hstack( [0, k*w1[2:m  ]/(2*h), 0] )   )	\
			      -diag( hstack( [0, k*w1[0:m-2]/(2*h), 0] )   ) 
			DF2 += diag( hstack( [0, k*w1[1:m-1]/(2*h)]    ), 1)	\
				  -diag( hstack(    [k*w1[1:m-1]/(2*h), 0] ),-1) 
			DF = DF1+DF2 
			F = -w[:,j] + dot(DF1+DF2/2.,w1) 		# Using Lemma 8.11 
			DF[0 ,:] = hstack([1., zeros(m-1)]) 	# Dirichlet conditions for DF 
			DF[-1,:] = hstack([zeros(m-1), 1.]) 
			F[0] = w1[0]-L[j+1]; F[-1] = w1[-1]-R[j+1] # Dirichlet conditions for F. 
												   # Original should be l(j*k), r(j*k)?
			w1 -= linalg.solve(DF,F) 
		w[:,j+1] = w1 
		print w1
	mesh(x,t,w.T)  
	return w 

def mesh(xvals,yvals,w):
	[x,y] = meshgrid(xvals,yvals)			# 3-D plot of solution w
	from mpl_toolkits.mplot3d import axes3d, Axes3D
	import matplotlib.pyplot as plt
	from matplotlib import cm
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	#ax.plot_wireframe(x, y, w, rstride=1, cstride=1)
	ax.plot_surface(x, y, w, rstride=4, cstride=1, cmap=cm.jet, antialiased=True,linewidth=0.1)
	ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('w')
	#ax.set_zlim3d(0,0.5)   					
	plt.show()

set_printoptions(precision=5,linewidth=150,suppress=True); 
w = burgers(10,11,0,1,10,11)
#w = burgers(10,11,0,1,50,100)


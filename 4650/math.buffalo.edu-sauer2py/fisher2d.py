# Program 8.8 Backward difference method with Newton iteration
#             for Fisher's equation with two-dim domain
# input: space region [xl,xr]x[yb,yt], time interval [t0,te], 
# M,N space steps in x and y directions, tsteps time steps 
# output: solution mesh [x,y,w]
# Example usage: [x,y,w] = fisher2d(0,1,0,1,0,2,20,20,40) 

from numpy import *
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

def fisher2d(xl,xr,yb,yt,t0,te,M,N,tsteps):
	def f(x,y): return 1 + cos(pi*x)*cos(pi*y)
	delt = float(te-t0)/tsteps
	D = 1.
	m = M+1; n = N+1; mn = m*n
	h = float(xr-xl)/M; k = float(yt-yb)/N
	x = linspace(xl,xr,m); y = linspace(yb,yt,n)
	Y,X = meshgrid(y,x) # careful: this way to get 1st index for x
	w = f(X,Y)			# define initial u
	plt.ion()
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for tstep in range(tsteps):
		v = reshape(w,mn,order='F')
		wold = w.copy()
		for it in range(3):
			b = zeros(mn); DF1 = zeros((mn,mn)); DF2 = zeros((mn,mn))
			for i in range(1,m-1):
				for j in range(1,n-1):
					DF1[ i+j*m , i-1+    j*m] = -1./h**2
					DF1[ i+j*m , i+1+    j*m] = -1./h**2
					DF1[ i+j*m , i  +    j*m] =  2./h**2 + 2./k**2 - 1. + 1./(1.*delt)
					DF1[ i+j*m , i  +(j-1)*m] = -1./k**2
					DF1[ i+j*m , i  +(j+1)*m] = -1./k**2
					DF2[ i+j*m , i  +    j*m] =  2.*w[i,j] 
					b[   i+j*m ] = -wold[i,j]/(1.*delt)  
			for i in range(m):		# bottom and top
				j = 0 
				DF1[ i+j*m , i+(j  )*m] = -3./(2.*k)
				DF1[ i+j*m , i+(j+1)*m] =  4./(2.*k)
				DF1[ i+j*m , i+(j+2)*m] = -1./(2.*k)
				j = n-1
				DF1[ i+j*m , i+(j  )*m] =  3./(2.*k)
				DF1[ i+j*m , i+(j-1)*m] = -4./(2.*k)
				DF1[ i+j*m , i+(j-2)*m] =  1./(2.*k)
			for j in range(1,n-1):	# left and right  
				i = 0 
				DF1[ i+j*m , i+  j*m] = -3./(2.*h)
				DF1[ i+j*m , i+1+j*m] =  4./(2.*h)
				DF1[ i+j*m , i+2+j*m] = -1./(2.*h)
				i = m-1 
				DF1[ i+j*m , i  +j*m] =  3./(2.*h)
				DF1[ i+j*m , i-1+j*m] = -4./(2.*h)
				DF1[ i+j*m , i-2+j*m] =  1./(2.*h)
			DF = DF1 + DF2
			F = dot(DF1 + DF2/2.,v) + b
			v -= linalg.solve(DF,F)
			w = reshape(v,(m,n),order='F')
			plt.cla()
			ax.plot_surface(X,Y,w,rstride=1,cstride=1,cmap=cm.jet,antialiased=True,linewidth=0.1)
			ax.set_zlim3d(0,2) 
			ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('w'); plt.draw()
	return x,y,w
set_printoptions(precision=14,linewidth=200)
x,y,w = fisher2d(-1,1,-1,1,0,0.02,40,40,2)  
#print w
raw_input("Enter to leave: ")


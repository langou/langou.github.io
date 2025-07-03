# Program 6.1 Euler's Method for Solving Initial Value Problems
# Python translation JR 1/17/12. Update for 2e 2/24/12.
# ydot evaluates rhs of differential equation
# Input: interval [a,b], initial value y0, number of steps n
# Output: time steps t, solution y
# Example usage: y=euler([0 1],1,10)

from pylab import plot, show

def euler(interval,y0,n):
	h = float(interval[1]-interval[0])/n
	t = [i*h for i in range(n+1)]
	y = [y0]
	for i in range(n):
		y.append( eulerstep(t[i],y[i],h) )
	plot(t,y)
	show()
	
def eulerstep(t,y,h):
	return y + h*ydot(t,y)

def ydot(t,y):
	# right-hand side of differential equation
	return t*y + t**3

euler( [0,1], 1., 10 )


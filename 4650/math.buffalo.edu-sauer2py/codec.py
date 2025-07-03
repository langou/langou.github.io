# Program 11.1 Audio codec
# Translated to Python by JR, 2/17/2012
from numpy import *
from scipy.io import wavfile
import os
from time import sleep
n  = 32								# length of window
nb = 127							# number of windows; must be > 1
b = 4; L = 5						# quantization information
q = 2*L/float(2**b-1)				# b bits on interval [-L, L]
J,I = meshgrid(range(2*n),range(n))	# form the MDCT matrix
M = cos((I+0.5)*(J+0.5+n/2)*pi/n)
M *= sqrt(2./n)
N = M.T								# inverse MDCT
Fs = 8192; f = 1					# Fs = sampling rate, f is multiple of base freq.
Fs2 = 4096
x = cos( array( range(1,(Fs2+1)) )*pi*64.*f/4096.)		# test signal
#from scikits.audiolab import play	
#play(x, Fs)						# Could not import alsa backend
soundamp = 2**14					# Largest allowed value is 2**15-1
wavfile.write('codec.py.wav',Fs, array( soundamp*x/max(x), dtype=int16) )
os.system('play codec.py.wav&')		# In Linux, package sox provides play

out = empty(n*(nb-1))
w   = empty((2*n,nb))
for k in range(nb):					# loop over windows
	x0 = x[k*n:2*n+k*n]
	y0 = dot(M,x0)
	y1 = around(y0/q)				# transform components quantized
	y2 = y1*q						#                and dequantized
	w[:,k] = dot(N,y2)				# invert the MDCT
	if k>0:
		w2 = w[n:2*n,k-1]; w3 = w[0:n,k]
		out[(k-1)*n:k*n] = (w2+w3)/2.		# collect the reconstructed signal
sleep(2.)
wavfile.write('codec.py.out.wav',Fs, array( soundamp*out/max(out), dtype=int16) )
os.system('play codec.py.out.wav&')		

from pylab import plot, show
plot(x)
plot(out)
show()

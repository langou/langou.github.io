# Program 10.3 Wiener filter
# Translated to Python by JR, 2/17/2012
from numpy import *
from scipy.io import wavfile
import os
def weiner(sourcefile):
	os.system('play '+sourcefile)
	rate,data = wavfile.read(sourcefile)	
	y = array(data,dtype=float)				# y is clean signal
	n = y.shape[0]						
	c = y[:min(n,2000000),0]*0.5 # work with first 2M samples of channel 0
	p = 1.3									# parameter for cutoff
	noise = std(c)*.50						# 50 percent noise
	n = len(c)								# n is length of signal
	r = noise*random.randn(n)				# pure noise
	x = c+r									# noisy signal
	wavfile.write('noisy.wav',rate,array(x,dtype=int16))
	os.system('play noisy.wav')
	fx = fft.fft(x); sfx = real(conj(fx)*fx)# take fft of signal, and
	sfcapprox = maximum(sfx-n*(p*noise)**2,0)# apply cutoff
	phi = sfcapprox/sfx						# define phi as derived
	xout = real(fft.ifft(phi*fx))			# invert the fft
	wavfile.write('noisycutoff.wav',rate,array(xout,dtype=int16))
	os.system('play noisycutoff.wav')		# compare x and xout

weiner('daft.wav')


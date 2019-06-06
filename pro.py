#! /usr/bin/python2.7
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import wfdb
import math

#Kernel function for Segmentation
def kernel(X,t):
	if t==0:
		print("Error: t can't be zero")
	else:
		g=[]

		for x in X:
			g.append(math.exp((x*x)/(2*t))/(math.sqrt(2*np.pi*t)))

		return g

#convolution function for finite signals
def convolution(f,g,M,x):
	p=0
	for n in range(-int(M),int(M)+1):
		if x-n>0 and x-n<4000:
			p+= f[x-n]*g[n]
		else:
			p+=0

	return p

def lin_interpolate(x1,y1,x2,y2,d,t):
	m=[]
	if t==True:
		t_x = np.arange(x1,x2,(x2-x1)/d)
		for x in t_x:
			m.append((((y2-y1)*(x-x1))/(x2-x1)) + y1)
		return m , t_x
	else:
		t_x = np.arange(x1+(x2-x1)/(d-1),x2,(x2-x1)/(d-1))
		for x in t_x:
			m.append((((y2-y1)*(x-x1))/(x2-x1)) + y1)
		return m , t_x


#importing data from a patients record
signals , fields = wfdb.rdsamp("s0010_re",sampto=4000)


freq = np.arange(0,np.pi,np.pi/4000)

print(freq)

#separated v1 for testing
v1 = []

for x in range(len(signals)):
	v1.append(signals[x][6])


#plt.plot(v1)
#plt.show()
print("Signal appended")

fc = 0.5  # Cut-off frequency of the filter
fs = 1000
w = fc / (fs / 2) # Normalizing the frequency
b, a = signal.butter(5, w, 'high')
output = signal.filtfilt(b, a, v1) #High pass filtering with the given cut-off frequency

print("filtered")

#FBSE application
fbse = []

N=len(output)
for x in range(N):
	k=0
	for i in range(N):
		k = k+output[i]*np.sin(x*np.pi*i/N)
	fbse.append(2*x*np.pi*k/N)
fbse1 = abs(np.array(fbse))

print("Transformed")

t = []

nstep = 6

for x in range(1,nstep+1):
	t.append(x*0.5)


L = []

for x in range(len(freq)):
	c=[]
	for p in t:
		c.append(convolution(fbse1,kernel(freq,p),6*p+1,x))
	L.append(c)



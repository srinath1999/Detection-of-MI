#! /usr/bin/python2.7
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import wfdb

#importing data from a patients record
signals , fields = wfdb.rdsamp("s0010_re",sampto=3999)


#separated v1 for testing
v1 = []

for x in range(len(signals)):
	v1.append(signals[x][6])


plt.plot(v1)
plt.show()


fc = 0.5  # Cut-off frequency of the filter
fs = 1000
w = fc / (fs / 2) # Normalizing the frequency
b, a = signal.butter(5, w, 'high')
output = signal.filtfilt(b, a, v1) #High pass filtering with the given cut-off frequency

plt.plot(output)
plt.show()

#FBSE application
fbse = []

N=len(output)
for x in range(N):
	k=0
	for i in range(N):
		k = k+output[i]*np.sin(x*np.pi*i/N)
	fbse.append(2*x*np.pi*k/N)
fbse1 = abs(np.array(fbse))


#inversing the FBSE
infbse = np.zeros(N)
for x in range(1,N):
	for i in range(1,N):
		infbse[x] = infbse[x] + fbse[i-1]*np.sin(i*np.pi*x/N)*N/(i*np.pi*N)
plt.plot(infbse)

plt.show()

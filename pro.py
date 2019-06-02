#! /usr/bin/python2.7
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

f = []
with open('s0054lre.dat', 'r') as fp:
    hex_list = ["{:02x}".format(ord(c)) for c in fp.read()]
def s16(value):
    return -(value & 0x8000) | (value & 0x7fff)
for x in xrange(len(hex_list)//2):
	f.append(s16(int(hex_list[2*x]+hex_list[2*x+1],16)))
signala = []
for x in xrange(4000):
	signala.append(f[12*x+9])
#plt.plot(signala)
signala=np.array(signala)
#plt.show()
fc = 0.5  # Cut-off frequency of the filter
fs = 1000
w = fc / (fs / 2) # Normalize the frequency
b, a = signal.butter(5, w, 'high')
output = signal.filtfilt(b, a, signala)
#plt.plot( output, label='filtered')
#plt.show()

fbse = []

N=len(output)
for x in xrange(1,N):
	k=0
	for i in xrange(N):
		k = k+output[i]*np.sin(x*np.pi*i/N)
	fbse.append(2*x*np.pi*k/N)
fbse1 = abs(np.array(fbse))
#plt.plot(fbse1)
#plt.show()
#1107,1350
fbse = np.array(fbse)
"""
t = np.linspace(0.0,1.0,50)
fbse[0:1116] = 0
for x in xrange(50):
	fbse[1116+x]=t[x]*fbse[1116+x]
	fbse[1350-x]=t[x]*fbse[1350-x]
fbse[1350:] = 0
#plt.plot(fbse)
#plt.show()"""
infbse = np.zeros(N)
for x in xrange(1,N):
	for i in xrange(1,N):
		infbse[x] = infbse[x] + fbse[i-1]*np.sin(i*np.pi*x/N)*N/(i*np.pi*x)
#plt.plot(infbse)
plt.plot(output)
plt.show()
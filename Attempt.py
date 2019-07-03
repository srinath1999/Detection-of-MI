import ewtpy as ew
import wfdb
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import timeit
import os


def entropy(arr,bins=10):
    pk =plt.hist(arr,bins=bins)
    return st.entropy(pk[0])

start = timeit.default_timer()

data_dir = "Data/"

folders = []

for x in os.walk(data_dir):
	if(x[1]!=[]):
		folders = x[1]



for currentfolder in folders:
	curr_dir = data_dir+currentfolder+"/"+currentfolder
	signals , fields = wfdb.rdsamp(curr_dir,sampto=4000)

	f_v = open(curr_dir+"_kurtosis_vector.txt","w+")
	p_v = open(curr_dir+ "_entropy_vector.txt","w+")
	y_v = open(curr_dir+"_skewness_vector.txt","w+")

	for i in range(12):
		sig = []
		for p in range(len(signals)):
			sig.append(signals[p][i])

		sig = np.array(sig)

		ewt, mfb,boundaries = ew.EWT1D(sig,N=9,type = "fbse")

		for p in range(9):
			sub_band = []
			for m in range(len(ewt)):
				sub_band.append(ewt[m][p])
			np.array(sub_band)
			f_v.write(str(st.kurtosis(sub_band))+'\n')
			p_v.write(str(entropy(sub_band))+"\n")
			y_v.write(str(st.skew(sub_band))+'\n')


stop = timeit.default_timer()

print('Time: ', stop - start)  



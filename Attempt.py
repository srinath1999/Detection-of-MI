import ewtpy as ew
import wfdb
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import timeit

start = timeit.default_timer()

file= "s0010_re"

signals , fields = wfdb.rdsamp(file,sampto=4000)

f_v = open(file +"_kurtosis_vector.txt","w+")
p_v = open(file + "_entropy_vector.txt","w+")
y_v = open(file+"_skewness_vector.txt","w+")

kurtosis_vector = []


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
		f_v.write(str(st.kurtosis(sub_band))+"\n")
		#p_v.write("Entropy Value here")
		y_v.write(str(st.skew(sub_band))+"\n")


stop = timeit.default_timer()

print('Time: ', stop - start)  




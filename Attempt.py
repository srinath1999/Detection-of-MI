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



data_dir = "Data/PATHS"

names = [line.rstrip('\n') for line in open(data_dir)]

heal_path = "Data/HEALTHY"

heal_names = [line.rstrip('\n') for line in open(heal_path)]

for name in names:
	name = name.split("/")
	signals , fields = wfdb.rdsamp(name[1],sampto=4000,pb_dir="ptbdb/"+name[0])

	outfile = "Data/"+name[1]+".npy"

	features = np.zeros([3,108])
	n=0
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

			features[0][n] = st.kurtosis(sub_band)
			features[1][n] = entropy(sub_band)
			features[2][n] = st.skew(sub_band)

			n +=1

	np.save(outfile,features)
	
stop = timeit.default_timer()

print('Time: ', stop - start)  



import ewtpy as ew
import wfdb
import numpy as np
import matplotlib.pyplot as plt

signals , fields = wfdb.rdsamp("s0010_re",sampto=4000)

v1=[]

for p in range(len(signals)):
	v1.append(signals[p][6])

v1 = np.array(v1)

plt.plot(v1)
plt.show()

ewt,mfb,boundaries = ew.EWT1D(v1,N=9,type = "fbse")




print(boundaries)

plt.plot(ewt)
plt.show()


plt.plot(mfb)
plt.show()
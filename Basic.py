import ewtpy
import wfdb
import numpy as np


#----------------------------------Reading Multiple Records----------------------------------#


def ReadSamp( Names ,Path = None, hop = None, channels = None ) :

	"""
	=============================

	Read the data samples from multiple files and return the data samples in the form 
	of a list of 3-D numpy arrays.

	----------

	Inputs:
		-Names	: The names of the different patient files which are to be read 
		without any extention.

	Optional Inputs:
		-Path	: A variable which contains the path of the folder in which all the 
		wfdb files(.dat, .hea and .atr) are there.
		-hop	:  The length of the hop for each file.

	Returns:
		-Data	: A list of numpy arrays.

	----------

	Notes:
	-The output is a list of 3 dimentional numpy arrays, each one of the 3 Dimentional 
	array represent data from a single patient.
	-The first level of abstraction of the numpy array is the number of channels we are
	taking the data of.
	-The next level of abstraction of the numpy array is the number of hops, the number 
	of windows we are having, if we consider a certain window length. If its not 
	provided then the length of the window will be the whole length of the file.
	-The third level is the hop length of the window we considered to segment the file.
	
	----------

	***I M P***

	Not every file's length can be integral multiple of hop length, so to have uniformity 
	in the data extracted, last few observations are omitted to make the number of readings 
	integral multiple of the hop size.

	=============================

	Example:
	
	>>> Data = ReadSamp( ['r1', 'r2', 'r3'], Path = "Desktop/DataFolder", hop = 4000, channels = [0, 4, 14] )

	"""
	
	if Path is not None:
		if Path[-1] != '/':
			Path = Path + "/"

	Data = []

	for name in Names:

		Signal, Fields = wfdb.rdsamp( Path + name, channels = channels )


		if hop is not None:
			#Number of hop
			nhops = int( len(Signal)/hop )
		else:
			nhops = 1
			hop = Signal.shape[0]

		if nhops is 0:
			raise ValueError('The hop size is greater than the length of the signal')

		arr = np.zeros(shape = (Signal.shape[1], nhops , hop))

		ini = 0
		for i in range(Signal.shape[1]):
			for j in range(nhops):
				for k in range(hop):
					arr[i,j,k] = Signal[(j*hop + k), i]

		Data.append(arr)

	return Data





if __name__ == "__main__" :

	Data = ReadSamp( ["s0014lre"], Path = "Patient1")

	print(Data)
	print()
	print(Data[0].shape)
	print()
	print(Data[0])
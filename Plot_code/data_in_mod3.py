# data_in_mod2.py
# Module that contains functions for reading in data
# Focus is on reading in difference between different run types

# to do: currently casting to floats, check if this is sufficient precision
# to do: check if data_dif_read should delete temporary array, or if it is done when they go out of scope

#### Include packages #####
import math
import mpmath as mp
import numpy as np

##### Path to data files #####
path_n = '/home/morrison/Lattice/Plots/'
##### File name variables #####
ener = 'energy_spec'
spec = 'spectrum'
zeta = 'zeta'
##### Run type variables #####
# wdv:		with delta V
# wodv:		without delta V
# hwdv:		homogeneous with delta V
# hwodv:	homogeneous without delta V
# lon:		longitudinal
# hlon:		homogeneous longitudinal
wdv = '_wdv1.out'
wodv = '_wodv.out'
hwdv = '_hwdv1.out'
hwodv = '_hwodv.out'
lon = '_lon1.out'
hlon = '_hlon1.out'
# nl:			number of lattice sites
# nl_h:		number of lattice sites for homogeneous run
nl = 64**3
nl_h = 2**3
# db:			number of lines in output data block
# db_h:		number of lines in output data block of homogeneous run
db = 58
db_h = 4

##### Run type variable list #####
rt_list = [hlon,hwodv,hwdv,lon,wodv,wdv]
nl_list = [nl_h,nl_h,nl_h,nl,nl,nl]
db_list = [db_h,db_h,db_h,db,db,db]

def data_read(f_in, rt, f_col, db_size=0, db_row=1, scl=1, root=False):
	# ar_out: 	np array that will hold the data read in
	# f_in:			taken from file name variables where data is
	# rt:				run type
	# f_col:		column number of interest in input file, counting from 0
	# db_size:	size of data block in file for rt
	# db_row:		row of interest in data block for rt
	# scl:			scaling factor that might be usefule eg for dividing by N_lat
	# root:			option to take sqrt of input
	
	##### Open data file #####
	data_f_n = path_n + f_in + rt
	data_f = open(data_f_n, 'r')
	##### Read in data #####
	data = []
	i = 0
	if db_size != 0:
		for line in data_f:
			i = i+1
			if i % (db_size) == db_row:
				line = line.strip()
				row = line.split()
				data.append(float(row[f_col]))
	elif db_size == 0:
		for line in data_f:
			line = line.strip()
			row = line.split()
			data.append(float(row[f_col]))
	##### Close data file #####
	data_f.close()
	##### Cast data to array #####
	ar_out = np.array(data)
	if root == True:
		ar_out = np.sqrt(ar_out)
	ar_out = ar_out/scl
#	print ar_out
	return ar_out

def data_dif_read(f_in, rt1, rt2, f_col, db_size1=1, db_size2=1, db_row1=0, db_row2=0, scl1=1, scl2=1, root=False):
	# ar_out: 	np array that will hold the data read in
	# f_in:			taken from file name variables where data is
	# rt1:			run type 1
	# rt2:			run type 2, rt1-rt2 computed
	# f_col:		column number of interest in input file
	# db_size1:	size of data block in file for rt1
	# db_size2:	size of data block in file for rt2
	# db_row1:	row of interest in data block for rt1
	# db_row2:	row of interest in data block for rt2

	##### Create data arrays #####
	#ar_data1 = np.array([])
	#ar_data2 = np.array([])
	##### Populate data arrays #####
	ar_data1 = data_read(f_in, rt1, f_col, db_size1, db_row1)
	ar_data2 = data_read(f_in, rt2, f_col, db_size2, db_row2)
	##### Difference arrays with possible options #####
	if root==True:
		ar_out = np.sqrt(ar_data1) - np.sqrt(ar_data2)
	if root==False:
		ar_out = ar_data1 - ar_data2
	return ar_out

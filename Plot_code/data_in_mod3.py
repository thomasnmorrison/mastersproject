# data_in_mod3.py
# Module that contains functions for reading in data
# Focus is on reading in difference between different run types

# to do: currently casting to floats, check if this is sufficient precision
# to do: check if data_dif_read should delete temporary array, or if it is done when they go out of scope
# to do: add an epsilon funtion
# to do: add a function that locates an index for a give phi value

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
lat = 'lat_sample'
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
l_lattice = 0.5 #different for some earlier runs
nl = 32**3 #64**3
nl_h = 2**3
# db:			number of lines in output data block
# db_h:		number of lines in output data block of homogeneous run
db = 58
db_h = 4

##### Run type variable list #####
rt_list = [hlon,hwodv,hwdv,lon,wodv,wdv]
#nl_list = [nl_h,nl_h,nl_h,nl,nl,nl]
#db_list = [db_h,db_h,db_h,db,db,db]

# Sampled runtype variable list
run_number = 11
wdv_s = '_wdv'#'_wdv8_sampled.out'#
wodv_s = '_wodv'#'_wodv8_sampled.out'#
hwdv_s = '_hwdv'#'_hwdv8_sampled.out'#
hwodv_s = '_hwodv'#'_hwodv8_sampled.out'#
lon_s = '_lon1_sampled.out'
hlon_s = '_hlon1_sampled.out'

n_sample = 32
n_sample_h = 1

#rt_s_list = [hwodv_s,wodv_s,hwdv_s,wdv_s,hlon_s,lon_s]
rt_s_list = [hwodv_s,wodv_s,hwdv_s,wdv_s]
for i in range(0,len(rt_s_list)):
	rt_s_list[i] = rt_s_list[i]+str(run_number)+'_sampled.out'
rt_s_list.append(hlon_s); rt_s_list.append(lon_s)
nl_list = [nl_h,nl,nl_h,nl,nl_h,nl]
db_list = [db_h,db,db_h,db,db_h,db]
db_s_list = [n_sample_h+1, n_sample+1, n_sample_h+1, n_sample+1, n_sample_h+1, n_sample+1]

##### Delta V parameters #####
# each entry corresponds to a run number
a2_list = [100,100,100,100,100,100,5.0e4,2.0e4,1.0e4,2.0e4,2.0e4]
b_1_list = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.3,0.2,0.2]
phi_p_list = [6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7]

a2 = a2_list[run_number-1]
b_1 = b_1_list[run_number-1]
phi_p = phi_p_list[run_number-1]

##### Function to read in data #####
def data_read(f_in, rt, f_col, db_size=0, db_row=1, scl=1, root=False, cut_off=0):
	# ar_out: 	np array that will hold the data read in
	# f_in:			taken from file name variables where data is
	# rt:				run type
	# f_col:		column number of interest in input file, counting from 0
	# db_size:	size of data block in file for rt
	# db_row:		row of interest in data block for rt
	# scl:			scaling factor that might be usefule eg for dividing by N_lat
	# root:			option to take sqrt of input
	# cut_off:	if not equal to zero will read in a number of data points equal to the value of cut_off
	
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
			if cut_off != 0 and i>=cut_off:
				break
	elif db_size == 0:
		for line in data_f:
			i = i+1
			line = line.strip()
			row = line.split()
			data.append(float(row[f_col]))
			if cut_off != 0 and i>=cut_off:
				break
	##### Close data file #####
	data_f.close()
	##### Cast data to array #####
	ar_out = np.array(data)
	if root == True:
		ar_out = np.sqrt(ar_out)
	ar_out = ar_out/scl
#	print ar_out
	return ar_out

##### Function to return an array with H #####
def get_hub(f_in, rt):
	# f_in:			taken from file name variables where data is
	# rt:				run type
	hub = data_read(f_in,rt,8)
	hub = np.sqrt(-1./3.*hub)
	return hub

def get_fld_dot(f_in, rt_index, col, a_scl, db_size=1, db_row=0):
	# f_in:			taken from file name variables where data is
	# rt:				run type
	# col: 			data column
	fld_dot = data_read(f_in, rt_s_list[rt_index], col, db_size, db_row)
	fld_dot = fld_dot/a_scl[rt_index]/a_scl[rt_index]/a_scl[rt_index]
	return fld_dot

##### Function to return the index location of the end of inflation #####
# to be more accurate I could fit a_scl*hub and interpolate the maximum
def get_end_infl(a_scl, hub):
	# a_scl:	array of a_scl for a particular run
	# hub:		array of H for a particular run
	end_infl = np.argmax(a_scl*hub)
	return end_infl

def index_finder(phi, phi_c):
	i = np.argmin(np.absolute(phi-phi_c))
	return i

##### Function to determine at what index location phi=phi_p	#####
# this could also be improved using some interpolation
# this function only returns the first instance, if it is used post inflation it will need to be modified
def get_p(phi,phi_s):
	# phi:		array with the values of phi for a particular run
	# phi_s:	value of phi_s to be found
	# p:			index location of phi=phi_p (round down)
	p = 0
	for i in range(0,len(phi)):
		if phi[i] >= phi_s:
			p = i
		else:
			break
	return p

##### Function to determine k_p asociated with phi_p #####
# k_p is not properly calculated
def get_k_p(a_scl,hub,end_infl,p):
	# a_scl:		array of a_scl for a particular run
	# hub:			array of H for a particular run
	# end_infl:	index location for end of inflation
	# p:				index location for phi = phi_p (or whatever other index location)
	# k_p:				output for k_p asociated with phi_p
	k_p = hub[p]*a_scl[p]/a_scl[end_infl]
	return k_p

# This may require taking derivatives
def get_epsilon(a_scl,hub):
	# a_scl:		array of a_scl for a particular run
	# hub:			array of H for a particular run
	epsilon = 0 
	return epsilon

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

# data_in_mod4.py
# Module that contains functions for reading in data from trapped/trapped-instability runs

# to do: clean out antiquated functions and variables
# to do: initialize db based on nl
# to do: initialize n_sample based on run number
# to do: fix longitudinal only run referencing
# to do: add the capability to use .out files named after the job-run id
# to do: add the capability to read in the run parameters from a file

#### Include packages #####
import math
import mpmath as mp
import numpy as np

##### Path to data files #####
path_n = '/home/morrison/Lattice/Plots/'#'/home/morrison/Lattice/Plots/Trapped_plots/Trapped_runs'#
##### File name variables #####
ener = 'energy_spec'
spec = 'spectrum'
zeta_n = 'zeta'
lat = 'lat_sample'

# nl_h:		number of lattice sites for homogeneous run
nl_h = 2**3
# db:			number of lines in output data block
# db_h:		number of lines in output data block of homogeneous run
db = 30#58 # 
db_h = 4

n_sample = 32#512#
n_sample_h = 1

##########################
##### Run parameters #####
##########################
l_lattice = [] # the parameter len used in each run
n_lattice = [1] # the number number of lattice sites used in each run
mpl = [1] # the value of the mpl parameter used in each run
m2 = []
lambda_chi = []
g2 = []
beta2 = []
phi_p = []

run_offset = 60 # this number is for file naming, 60 corresponds to the start of trapped inflation runs. eg g2_list[0] corresponds to run 0+run_offset.

l_lattice_list = [0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.053,0.053,0.053,0.053,0.053,0.64,0.64,0.64,0.64,0.64,0.64,0.64,0.64] # up to run 80
n_lattice_list = [32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3,32**3] # up to run 80
mpl_l = [1.e3/(math.sqrt(8.0*math.pi)),1.e6/(math.sqrt(8.0*math.pi)),1.e3,1.e6]
mpl_list = [mpl_l[0],mpl_l[0],mpl_l[0],mpl_l[0],mpl_l[0],mpl_l[1],mpl_l[2],mpl_l[3],mpl_l[2],mpl_l[2],mpl_l[2],mpl_l[2],mpl_l[2],mpl_l[2],mpl_l[2],mpl_l[2],mpl_l[2],
mpl_l[2],mpl_l[2],mpl_l[2],mpl_l[2]] # up to run 80
##### Potential parameters #####
m2_list = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0] # up to run 80
lambda_chi_list = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] # up to run 80
g2_list = [0.0,1.0e6/(8.0*math.pi),1.0e6/(8.0*math.pi),1.0e6/(8.0*math.pi),1.0e6/(8.0*math.pi),0.0,0.0,0.0,0.0,1.e6,1.e6,1.e6,1.e6,0.0,1.e6,1.e6,1.e6,1.e6,
1.e6,1.e5,1.e6] # up to run 80
beta2_list = [0.0,0.0,1.e-8*g2_list[2],1.e-7*g2_list[3],1.e-6*g2_list[4],0.0,0.0,0.0,0.0,0.0,1.e-6*g2_list[10],1.e-4*g2_list[11],1.e-2*g2_list[12],0.0,0.0,
8.e-4*g2_list[15],9.e-4*g2_list[16],1.e-4*g2_list[17],7.5e-4*g2_list[18],2.5e-3*g2_list[19],9.2e-4*g2_list[20]] # up to run 80
phi_p_l = [3.2*math.sqrt(8.0*math.pi),3.2]
phi_p_list = [phi_p_l[0],phi_p_l[0],phi_p_l[0],phi_p_l[0],phi_p_l[0],phi_p_l[0],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],
phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1],phi_p_l[1]] # up to run 80

##############################
##### Run initialization #####
############################## 
wdv_s = '_wdv'# wdv:		with delta V
wodv_s = '_wodv'# wodv:		without delta V
hwdv_s = '_hwdv'# hwdv:		homogeneous with delta V
hwodv_s = '_hwodv'# hwodv:	homogeneous without delta V
lon_s = '_lon1_sampled.out'# lon:		longitudinal
hlon_s = '_hlon1_sampled.out'# hlon:		homogeneous longitudinal
rt_s_list_base = [hwodv_s,wodv_s,hwdv_s,wdv_s]
rt_s_list = []

def clear_rt():
	rt_s_list[:]=[]
	return

def init_rt(rn):
	for i in range(0,len(rt_s_list_base)):
		rt_s_list.append(rt_s_list_base[i]+str(rn)+'_sampled.out')
	rt_s_list.append(hlon_s); rt_s_list.append(lon_s)
	#print rt_s_list
	return

def init_param(rn):
	l_lattice[:] = []; n_lattice[:] = []; mpl[:] = []; m2[:] = []; lambda_chi[:] = []; g2[:] = []; beta2[:] = []; phi_p[:] = []
	l_lattice.append(l_lattice_list[rn-run_offset])
	n_lattice.append(n_lattice_list[rn-run_offset])
	mpl.append(mpl_list[rn-run_offset])
	m2.append(m2_list[rn-run_offset])
	lambda_chi.append(lambda_chi_list[rn-run_offset])
	g2.append(g2_list[rn-run_offset])
	beta2.append(beta2_list[rn-run_offset])
	phi_p.append(phi_p_list[rn-run_offset])
	return

#nl_list = [nl_h,nl,nl_h,nl,nl_h,nl] # nl is not defined
db_list = [db_h,db,db_h,db,db_h,db]
db_s_list = [n_sample_h+1, n_sample+1, n_sample_h+1, n_sample+1, n_sample_h+1, n_sample+1]

#####################
##### Data read #####
#####################

##### Function to read in data #####
# remove root function
def data_read(f_in, rt, f_col, db_size=0, db_row=1, scl=1, cut_off=0, index_list=[]):
	# ar_out: 	np array that will hold the data read in
	# f_in:			taken from file name variables where data is
	# rt:				run type
	# f_col:		column number of interest in input file, counting from 0
	# db_size:	size of data block in file for rt
	# db_row:		row of interest in data block for rt
	# scl:			scaling factor that might be usefule eg for dividing by N_lat
	# cut_off:	if not equal to zero will read in a number of data points equal to the value of cut_off
	
	##### Open data file #####
	data_f_n = path_n + f_in + rt
	data_f = open(data_f_n, 'r')
	##### Read in data #####
	data = []
	i = 0
	n = 0
	m=0
	if (len(index_list)==0):
		if db_size != 0:
			for line in data_f:
				i = i+1
				if i % (db_size) == db_row:
					line = line.strip()
					row = line.split()
					data.append(float(row[f_col]))
					n = n+1
				if cut_off != 0 and n>=cut_off:
					break
		elif db_size == 0:
			for line in data_f:
				i = i+1
				line = line.strip()
				row = line.split()
				data.append(float(row[f_col]))
				n = n+1
				if cut_off != 0 and n>=cut_off:
					break
	else:
		if db_size != 0:
			for line in data_f:
				i = i+1
				if i % (db_size) == db_row:
					if n == index_list[m]:
						line = line.strip()
						row = line.split()
						data.append(float(row[f_col]))
						m = m+1
					n = n+1
				if cut_off != 0 and n>=cut_off:
					break
				if m>=len(index_list):
					break
		elif db_size == 0:
			for line in data_f:
				i = i+1
				if n == index_list[m]:
					line = line.strip()
					row = line.split()
					data.append(float(row[f_col]))
					m = m+1
				n = n+1
				if cut_off != 0 and n>=cut_off:
					break
				if m>=len(index_list):
					break
	##### Close data file #####
	data_f.close()
	##### Cast data to array #####
	ar_out = np.array(data)
	ar_out = ar_out/scl
#	print ar_out
	return ar_out

##### Function to return an array with H #####
def get_hub(f_in, rt, cut_off1=0, units=True):
	# f_in:			taken from file name variables where data is
	# rt:				run type
	hub = data_read(f_in,rt,8,cut_off=cut_off1)
	hub = np.sqrt(-1./3.*hub)
	if units == True:
		print('Hubble units=True')
		hub = hub/mpl[0]
	return hub

def get_fld_dot(f_in, rt_index, col, a_scl, db_size=1, db_row=0, cut_off1=0, units=True):
	# f_in:			taken from file name variables where data is
	# rt:				run type
	# col: 			data column
	fld_dot = data_read(f_in, rt_s_list[rt_index], col, db_size, db_row,cut_off=cut_off1)
	fld_dot = fld_dot/a_scl[rt_index]/a_scl[rt_index]/a_scl[rt_index]
	if units == True:
		fld_dot = fld_dot/mpl[0]
	return fld_dot

def get_spec(f_in, rt_index, col, db_in, cut_off1=0, units=True):
	spec = []
	print 'in spec'
	for i in range (1,db_in):
		print i
		spec.append(data_read(f_in, rt_s_list[rt_index], col, db_size=db_in, db_row=i, scl=nl, cut_off=cut_off1))
	return spec

##### Function to return the index location of the end of inflation #####
# to be more accurate I could fit a_scl*hub and interpolate the maximum
# if inflation ends after the run this will return the last index of the run
# if the run is not in inflation this will return the first index of the run
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

##### Function to initialize all data #####
# I'm leaving this function here for now, but never use this, it is just too slow.
def data_init(a_scl, hub, tau, phi, phi_dot, chi, chi_dot, zeta, zeta_mean, rho_in=[], hub_machine=[],cut=0):
	
	for i in range(0,len(rt_s_list)):
		a_scl.append(data_read(ener, rt_s_list[i], 1, cut_off=cut))
		hub.append(get_hub(ener, rt_s_list[i], cut_off1=cut))
		hub_machine.append(get_hub(ener, rt_s_list[i], cut_off1=cut, units=False))
		tau.append(data_read(ener, rt_s_list[i], 0, cut_off=cut))
		zeta_mean.append(data_read(zeta_n, rt_s_list[i], 3, cut_off=cut))
		rho_in.append(data_read(ener, rt_s_list[i], 2, cut_off=cut))

	for i in range(0,len(rt_s_list)):
		phi_temp2 = []	
		chi_temp2 = []
		phi_dot_temp2 = []
		chi_dot_temp2 = []
		zeta_temp2 = []
		#for j in range(1,n_sample+1):
		for j in range(1,db_s_list[i]):
			print(j)
			if i in range(0,4):
				phi_temp1 = []
				chi_temp1 = []
				phi_dot_temp1 = []
				chi_dot_temp1 = []
				zeta_temp1 = []
				#print rt_s_list[i]
				#print 'i=',i
				#print 'j=',j
				phi_temp1.append(data_read(lat, rt_s_list[i], 1, db_size=db_s_list[i], db_row=j,cut_off=cut))
				chi_temp1.append(data_read(lat, rt_s_list[i], 2, db_size=db_s_list[i], db_row=j,cut_off=cut))
				phi_dot_temp1.append(get_fld_dot(lat, i, 3, a_scl, db_size=db_s_list[i], db_row=j,cut_off1=cut))
				chi_dot_temp1.append(get_fld_dot(lat, i , 4, a_scl, db_size=db_s_list[i], db_row=j,cut_off1=cut))
				zeta_temp1.append(data_read(lat, rt_s_list[i], 5, db_size=db_s_list[i], db_row=j,cut_off=cut))
				phi_temp2.append(phi_temp1)
				chi_temp2.append(chi_temp1)
				phi_dot_temp2.append(phi_dot_temp1)
				chi_dot_temp2.append(chi_dot_temp1)
				zeta_temp2.append(zeta_temp1)
			if i in range(4,6):
				phi_temp1 = []
				phi_dot_temp1 = []
				zeta_temp1 = []
				phi_temp1.append(data_read(lat, rt_s_list[i], 1, db_size=db_s_list[i], db_row=j,cut_off=cut))
				phi_dot_temp1.append(get_fld_dot(lat, i, 2, a_scl, db_size=db_s_list[i], db_row=j,cut_off1=cut))
				zeta_temp1.append(data_read(lat, rt_s_list[i], 3, db_size=db_s_list[i], db_row=j,cut_off=cut))
				phi_temp2.append(phi_temp1)
				zeta_temp2.append(zeta_temp1)
		phi.append(phi_temp2)
		chi.append(chi_temp2)
		phi_dot.append(phi_dot_temp2)
		chi_dot.append(chi_dot_temp2)
		zeta.append(zeta_temp2)

	return

# function to read in a_scl and hub for a run
def data_init_scl_short(a_scl, hub, rt, cut=0):
	#a_scl.append(data_read(ener, rt_s_list[rt], 1, cut_off=cut))
	#hub.append(get_hub(ener, rt_s_list[rt], cut_off1=cut))
	a_scl[:] = data_read(ener, rt_s_list[rt], 1, cut_off=cut)
	hub[:] = get_hub(ener, rt_s_list[rt], cut_off1=cut)
	return

# function to read the homogeneous components of the fields
def data_init_fld_hom_short(f1,f2,rt,cut=0):
	#f1.append(data_read(ener, rt_s_list[rt],9,cut_off=cut))
	#f2.append(data_read(ener, rt_s_list[rt],10,cut_off=cut))
	f1[:] = data_read(ener, rt_s_list[rt],9,cut_off=cut)
	f2[:] = data_read(ener, rt_s_list[rt],10,cut_off=cut)
	return

# function to read in zeta moments
# old version
#def data_init_zeta_short(a_scl, hub, zeta_mean, rt, zeta_moment_in,cut=0,moment_max=0):
#	i=rt
#	zeta_mean.append(data_read(zeta_n, rt_s_list[i], 3, cut_off=cut))
#
#	for j in range(0,moment_max):
#		if i in range(0,4):
#			zeta_moment_temp1 = []
#			#zeta_moment_temp2 = []
#			zeta_moment_temp1.append(data_read(zeta_n, rt_s_list[i], 4+j, cut_off=cut))
#			zeta_moment_temp2.append(zeta_moment_temp1)
#	zeta_moment_in.append(zeta_moment_temp2)
#	return

# function to read in zeta moments
# new version
# Organized as zeta_moment_in[moment number-1][time step]
# Moments 2 and up are reduced, not raw
def data_init_zeta_moment(zeta_moment_in, rt,cut=0, moment_max=0):
	#zeta_mean[:] = data_read(zeta_n, rt_s_list[rt], 3, cut_off=cut)
	zeta_moment_in.append(data_read(zeta_n, rt_s_list[rt], 3, cut_off=cut))
	for j in range(0,moment_max):
		if rt in range(0,4):
			zeta_moment_in.append(data_read(zeta_n, rt_s_list[rt], 4+j, cut_off=cut))
	return

# function to read in zeta at samples lattice sites
def data_init_zeta_sample(zeta, rt, cut=0, sample_max=n_sample):
	i=rt
	zeta_temp1 = []
	zeta_temp2 = []
	for j in range(1,n_sample+1):
		if i in range(0,4):
			zeta_temp1 = []
			zeta_temp1.append(data_read(lat, rt_s_list[i], 5, db_size=db_s_list[i], db_row=j,cut_off=cut))
			zeta_temp2.append(zeta_temp1)
	zeta.append(zeta_temp2)
	return

# function to read in spectrum info
# need to pass in a_scl after it is read in
# Output is of the form P_ff[k_number][time step]
# Previously output was of the form P_ff[k number][0][time step]
def data_init_spec(k_spec, P_f1f1, P_df1df1, P_df1f1_re, P_df1f1_im, P_f2f2, P_df2df2, P_df2f2_re, P_df2f2_im, P_zz, rt, a_scl, cut=0):
	##### Read k_spec #####
	k_spec.append(data_read(spec, rt_s_list[rt], 0, cut_off=db-1))
	##### Read Spectral #####
	for i in range(1,db):
		P_phiphi_temp = []
		P_dphidphi_temp = []
		P_dphiphi_re_temp = []
		P_dphiphi_im_temp = []
		P_chichi_temp = []
		P_dchidchi_temp = []
		P_dchichi_re_temp = []
		P_dchichi_im_temp = []
		P_zz_temp = []
		#P_dzdz_temp = []
		#P_dzz_re_temp = []
		#P_dzz_im_temp = []
		# Read data (out of date)
		#P_phiphi_temp.append(data_read(spec, rt_s_list[rt], 1, db_size=db, db_row=i, cut_off=cut))
		#P_dphidphi_temp.append(data_read(spec, rt_s_list[rt], 2, db_size=db, db_row=i, cut_off=cut))
		#P_dphiphi_re_temp.append(data_read(spec, rt_s_list[rt], 3, db_size=db, db_row=i, cut_off=cut))
		#P_dphiphi_im_temp.append(data_read(spec, rt_s_list[rt], 4, db_size=db, db_row=i, cut_off=cut))
		#P_chichi_temp.append(data_read(spec, rt_s_list[rt], 5, db_size=db, db_row=i, cut_off=cut))
		#P_dchidchi_temp.append(data_read(spec, rt_s_list[rt], 6, db_size=db, db_row=i, cut_off=cut))
		#P_dchichi_re_temp.append(data_read(spec, rt_s_list[rt], 7, db_size=db, db_row=i, cut_off=cut))
		#P_dchichi_im_temp.append(data_read(spec, rt_s_list[rt], 8, db_size=db, db_row=i, cut_off=cut))
		#P_zz_temp.append(data_read(spec, rt_s_list[rt], 9, db_size=db, db_row=i, cut_off=cut))
		### Read data ###
		P_phiphi_temp[:] = data_read(spec, rt_s_list[rt], 1, db_size=db, db_row=i, cut_off=cut)
		P_dphidphi_temp[:] = data_read(spec, rt_s_list[rt], 2, db_size=db, db_row=i, cut_off=cut)
		P_dphiphi_re_temp[:] = data_read(spec, rt_s_list[rt], 3, db_size=db, db_row=i, cut_off=cut)
		P_dphiphi_im_temp[:] = data_read(spec, rt_s_list[rt], 4, db_size=db, db_row=i, cut_off=cut)
		P_chichi_temp[:] = data_read(spec, rt_s_list[rt], 5, db_size=db, db_row=i, cut_off=cut)
		P_dchidchi_temp[:] = data_read(spec, rt_s_list[rt], 6, db_size=db, db_row=i, cut_off=cut)
		P_dchichi_re_temp[:] = data_read(spec, rt_s_list[rt], 7, db_size=db, db_row=i, cut_off=cut)
		P_dchichi_im_temp[:] = data_read(spec, rt_s_list[rt], 8, db_size=db, db_row=i, cut_off=cut)
		P_zz_temp[:] = data_read(spec, rt_s_list[rt], 9, db_size=db, db_row=i, cut_off=cut)
		# Scale by factors of nlat and a_scl (out of date)
		#print(len(P_phiphi_temp))
		#print(len(P_phiphi_temp[0]))
		#P_phiphi_temp[0] = P_phiphi_temp[0]/(n_lattice[0]**2)
		#P_dphidphi_temp[0] = P_dphidphi_temp[0]/(n_lattice[0]**2*np.power(a_scl,6))
		#P_dphiphi_re_temp[0] = P_dphiphi_re_temp[0]/(n_lattice[0]**2*np.power(a_scl,3))
		#P_dphiphi_im_temp[0] = P_dphiphi_im_temp[0]/(n_lattice[0]**2*np.power(a_scl,3))
		#P_chichi_temp[0] = P_chichi_temp[0]/(n_lattice[0]**2)
		#P_dchidchi_temp[0] = P_dchidchi_temp[0]/(n_lattice[0]**2*np.power(a_scl,6))
		#P_dchichi_re_temp[0] = P_dchichi_re_temp[0]/(n_lattice[0]**2*np.power(a_scl,3))
		#P_dchichi_im_temp[0] = P_dchichi_im_temp[0]/(n_lattice[0]**2*np.power(a_scl,3))
		#P_zz_temp[0] = P_zz_temp[0]/(n_lattice[0]**2)
		### Scaling ###
		P_phiphi_temp_ar = np.array(P_phiphi_temp)/(n_lattice[0]**2)
		P_dphidphi_temp_ar = np.array(P_dphidphi_temp)/(n_lattice[0]**2*np.power(a_scl,6))
		P_dphiphi_re_temp_ar = np.array(P_dphiphi_re_temp)/(n_lattice[0]**2*np.power(a_scl,3))
		P_dphiphi_im_temp_ar = np.array(P_dphiphi_im_temp)/(n_lattice[0]**2*np.power(a_scl,3))
		P_chichi_temp_ar = np.array(P_chichi_temp)/(n_lattice[0]**2)
		P_dchidchi_temp_ar = np.array(P_dchidchi_temp)/(n_lattice[0]**2*np.power(a_scl,6))
		P_dchichi_re_temp_ar = np.array(P_dchichi_re_temp)/(n_lattice[0]**2*np.power(a_scl,3))
		P_dchichi_im_temp_ar = np.array(P_dchichi_im_temp)/(n_lattice[0]**2*np.power(a_scl,3))
		P_zz_temp_ar = np.array(P_zz_temp)/(n_lattice[0]**2)
		# Append data (out of date)
		#P_f1f1.append(np.array(P_phiphi_temp))
		#P_df1df1.append(np.array(P_dphidphi_temp))
		#P_df1f1_re.append(np.array(P_dphiphi_re_temp))
		#P_df1f1_im.append(np.array(P_dphiphi_im_temp))
		#P_f2f2.append(np.array(P_chichi_temp))
		#P_df2df2.append(np.array(P_dchidchi_temp))
		#P_df2f2_re.append(np.array(P_dchichi_re_temp))
		#P_df2f2_im.append(np.array(P_dchichi_im_temp))
		#P_zz.append(np.array(P_zz_temp))
		### Append data ###
		P_f1f1.append(P_phiphi_temp_ar)
		P_df1df1.append(P_dphidphi_temp_ar)
		P_df1f1_re.append(P_dphiphi_re_temp_ar)
		P_df1f1_im.append(P_dphiphi_im_temp_ar)
		P_f2f2.append(P_chichi_temp_ar)
		P_df2df2.append(P_dchidchi_temp_ar)
		P_df2f2_re.append(P_dchichi_re_temp_ar)
		P_df2f2_im.append(P_dchichi_im_temp_ar)
		P_zz.append(P_zz_temp_ar)
		#print(len(P_phiphi))
		#print(len(P_phiphi[0]))
		#print(len(P_phiphi[0][0]))
		# Calc Det P and n_k_phi
		#Det_P.append(np.absolute(P_f1f1[i-1][0]*P_df1df1[i-1][0] - (np.power(P_df1f1_re[i-1][0],2) + np.power(P_dphiphi_im[i-1][0],2)))) # Assuming power spectrum is Hermitian
		#Det_P_chi.append(np.absolute(P_chichi[i-1][0]*P_dchidchi[i-1][0] - (np.power(P_dchichi_re[i-1][0],2) + np.power(P_dchichi_im[i-1][0],2)))) # Assuming power spectrum is Hermitian
		#print(len(Det_P))
		#print(len(Det_P[0]))
		#print(len(Det_P[0][0]))
		#n_k_phi.append(0.5*np.log(Det_P[i-1][0]))
		#n_k_chi.append(0.5*np.log(Det_P_chi[i-1][0]))
		#print(len(n_k_phi))
		#print(len(n_k_phi[0]))
	return

# Function to initialize data from a sampling of the fields
# Untested
# Data is organized as fld[sample number][time step]
def data_init_fld_sample(phi, chi, zeta, rt, phi_dot_in=[],chi_dot_in=[],cut=0, sample_max=0):
	phi_temp = []
	chi_temp = []
	zeta_temp = []
	#phi_dot_temp = []
	#chi_dot_temp = []
	if sample_max == 0:
		j_max = db_s_list[rt]
	else:
		j_max = sample_max
	#phi[:] = np.empty(j_max)
	for j in range(1,j_max):
		print(j)
		if rt in range(0,4):		
			#print rt_s_list[rt]
			#print 'rt=',
			#print 'j=',j
			phi_temp = []
			chi_temp = []
			zeta_temp = []
			phi_temp[:] = data_read(lat, rt_s_list[rt], 1, db_size=db_s_list[rt], db_row=j,cut_off=cut)
			chi_temp[:] = data_read(lat, rt_s_list[rt], 2, db_size=db_s_list[rt], db_row=j,cut_off=cut)
			#phi_dot_temp[:] = get_fld_dot(lat, rt, 3, a_scl, db_size=db_s_list[rt], db_row=j,cut_off1=cut)
			#chi_dot_temp[:] = get_fld_dot(lat, rt, 4, a_scl, db_size=db_s_list[rt], db_row=j,cut_off1=cut)
			zeta_temp[:] = data_read(lat, rt_s_list[rt], 5, db_size=db_s_list[rt], db_row=j,cut_off=cut)
			phi.append(phi_temp)
			chi.append(chi_temp)
			#phi_dot_in.append(phi_dot_temp)
			#chi_dot_in.append(chi_dot_temp)
			zeta.append(zeta_temp)
		if rt in range(4,6):
			#zeta_temp1 = []
			phi_temp[:] = data_read(lat, rt_s_list[rt], 1, db_size=db_s_list[rt], db_row=j,cut_off=cut)
			#phi_dot_temp1.append(get_fld_dot(lat, rt, 2, a_scl, db_size=db_s_list[rt], db_row=j,cut_off1=cut))
			zeta_temp[:] = data_read(lat, rt_s_list[rt], 3, db_size=db_s_list[rt], db_row=j,cut_off=cut)
			phi.append(phi_temp)
			#phi_dot_in.append(phi_dot_temp)
			zeta.append(zeta_temp)
	return

def data_init_shorter(a_scl, hub, tau, phi, phi_dot, chi, chi_dot, zeta, zeta_mean, rho_in=[], hub_machine=[],cut=0,index_list=[]):
	
	for i in range(0,len(rt_s_list)):
		a_scl.append(data_read(ener, rt_s_list[i], 1, cut_off=cut))
		#hub.append(get_hub(ener, rt_s_list[i], cut_off1=cut))
		#hub_machine.append(get_hub(ener, rt_s_list[i], cut_off1=cut, units=False))
		#tau.append(data_read(ener, rt_s_list[i], 0, cut_off=cut))
		#zeta_mean.append(data_read(zeta_n, rt_s_list[i], 3, cut_off=cut))
		#rho_in.append(data_read(ener, rt_s_list[i], 2, cut_off=cut))

	for i in range(0,len(rt_s_list)):
		phi_temp2 = []	
		chi_temp2 = []
		#phi_dot_temp2 = []
		#chi_dot_temp2 = []
		#zeta_temp2 = []
		#for j in range(1,n_sample+1):
		for j in range(1,db_s_list[i]):
			print(j)
			if i in range(0,4):
				phi_temp1 = []
				chi_temp1 = []
				#phi_dot_temp1 = []
				#chi_dot_temp1 = []
				#zeta_temp1 = []
				#print rt_s_list[i]
				#print 'i=',i
				#print 'j=',j
				phi_temp1.append(data_read(lat, rt_s_list[i], 1, db_size=db_s_list[i], db_row=j,cut_off=cut))
				chi_temp1.append(data_read(lat, rt_s_list[i], 2, db_size=db_s_list[i], db_row=j,cut_off=cut))
				#phi_dot_temp1.append(get_fld_dot(lat, i, 3, a_scl, db_size=db_s_list[i], db_row=j,cut_off1=cut))
				#chi_dot_temp1.append(get_fld_dot(lat, i , 4, a_scl, db_size=db_s_list[i], db_row=j,cut_off1=cut))
				#zeta_temp1.append(data_read(lat, rt_s_list[i], 5, db_size=db_s_list[i], db_row=j,cut_off=cut))
				phi_temp2.append(phi_temp1)
				chi_temp2.append(chi_temp1)
				#phi_dot_temp2.append(phi_dot_temp1)
				#chi_dot_temp2.append(chi_dot_temp1)
				#zeta_temp2.append(zeta_temp1)
			if i in range(4,6):
				phi_temp1 = []
				#phi_dot_temp1 = []
				#zeta_temp1 = []
				phi_temp1.append(data_read(lat, rt_s_list[i], 1, db_size=db_s_list[i], db_row=j,cut_off=cut))
				#phi_dot_temp1.append(get_fld_dot(lat, i, 2, a_scl, db_size=db_s_list[i], db_row=j,cut_off1=cut))
				#zeta_temp1.append(data_read(lat, rt_s_list[i], 3, db_size=db_s_list[i], db_row=j,cut_off=cut))
				phi_temp2.append(phi_temp1)
				#zeta_temp2.append(zeta_temp1)
		phi.append(phi_temp2)
		chi.append(chi_temp2)
		#phi_dot.append(phi_dot_temp2)
		#chi_dot.append(chi_dot_temp2)
		#zeta.append(zeta_temp2)

	return

# missing variable definitions
def run_append(a_scl, hub, tau, phi, phi_dot, chi, chi_dot, zeta, zeta_mean, rho, hub_mac):
	a_scl_rn.append(a_scl)
	hub_rn.append(hub)
	hub_mac_rn.append(hub_mac)
	tau_rn.append(tau)
	phi_rn.append(phi)
	phi_dot_rn.append(phi_dot)
	chi_rn.append(chi)
	chi_dot_rn.append(chi_dot)
	zeta_rn.append(zeta)
	zeta_mean_rn.append(zeta_mean)
	rho_rn.append(rho)
	return

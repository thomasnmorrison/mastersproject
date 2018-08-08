# dif_plots.py
# Program to create plots from the difference between runs

#### Include packages #####
import matplotlib.pyplot as plt
import math
import mpmath as mp
import numpy as np
import data_in_mod3 as di

##### Data read in #####
# The format is data is kept in a list where each entry is an array storing the data for one run type (eg. with/without dv)
# the ordering can be found in the data_in_mod2 module.
a_scl = []
a_scl_t = np.array([])
tau = []
zeta = []
#for rt in di.rt_list:
#	print rt
#	print len(a_scl)
#	a_scl.append(di.data_read(di.ener, rt, 1))
#	tau.append(di.data_read(di.ener, rt, 0))
#	zeta.append(di.data_read(di.zeta, rt, 3))

phi = []
chi = []
#phi.append(di.data_read(di.spec, di.rt_list[0], 1, db_size=di.db_list[0],scl=di.nl_list[0], root=True))
for i in range(0,len(di.rt_list)):
	a_scl.append(di.data_read(di.ener, di.rt_list[i], 1))
	tau.append(di.data_read(di.ener, di.rt_list[i], 0))
	zeta.append(di.data_read(di.zeta, di.rt_list[i], 3))
	phi.append(di.data_read(di.spec, di.rt_list[i], 1, db_size=di.db_list[i],scl=di.nl_list[i], root=True))
	if i == 0 or i == 3:
		chi.append(np.array([]))
	if i != 0 and i != 3:	# Exception for longitudinal only data
		chi.append(di.data_read(di.spec, di.rt_list[i], 5, db_size=di.db_list[i], db_row=30,scl=di.nl_list[i], root=True))
#print len(tau)
#print len(zeta)

# read in and compute hub

##### Plotting #####
##### Plot Zeta #####
plt.figure(1)
#plt.subplot(2,1,1)
#plt.scatter(tau[0], zeta[1]-zeta[0])
#plt.yscale('log')
#plt.xlim(tau[0][0],tau[0][-1])
#plt.grid(True)
#plt.subplot(2,1,2)
#plt.scatter(tau[2], zeta[3]-zeta[2])
#plt.yscale('log')
#plt.xlim(tau[1][0],tau[1][-1])
#plt.grid(True)
#for i in range(1,len(di.rt_list)/2):
#	plt.subplot(len(di.rt_list)/2,1,i+1)
#	plt.scatter(tau[i], np.absolute(zeta[i+len(di.rt_list)/2]-zeta[i]), s=2)
#	plt.yscale('log')
#	plt.xlim(tau[i][0],tau[i][-1])
#	plt.grid(True)

plt.subplot(1,1,1)
plt.scatter(tau[4], zeta[5]*10**12-zeta[4]*10**12)
#plt.yscale('log')
plt.xlim(tau[4][0],tau[4][-1])
plt.grid(True)
##### Plot Chi #####
plt.figure(2)
plt.subplot(1,1,1)
plt.scatter(tau[4], chi[5]*10**15-chi[4]*10**15)
plt.xlim(tau[1][0],tau[1][-1])
plt.grid(True)
##### Plot Phi #####
#plt.figure(3)
#plt.subplot(2,1,1)
#plt.scatter(tau[0], phi[2]-phi[0])
#plt.xlim(tau[0][0],tau[0][-1])
#plt.grid(True)
#plt.subplot(2,1,2)
#plt.scatter(tau[1], phi[3]-phi[1])
#plt.xlim(tau[1][0],tau[1][-1])
#plt.grid(True)

plt.show()

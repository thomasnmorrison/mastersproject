# sampled_plots.py
# Program to make plots of sampled tragectories

#### Include packages #####
import matplotlib.pyplot as plt
import math
import mpmath as mp
import numpy as np
import data_in_mod3 as di
import plot_def_mod as pd

##### Read-in data #####
a_scl = []
hub = []
tau = []
phi = []
phi_dot = []#does not include mpl factor
chi = []
chi_dot = []#does not include mpl factor
zeta = []
zeta_mean = []

for i in range(0,len(di.rt_s_list)):
	a_scl.append(di.data_read(di.ener, di.rt_s_list[i], 1))
	hub.append(di.get_hub(di.ener, di.rt_s_list[i]))
	tau.append(di.data_read(di.ener, di.rt_s_list[i], 0))
	zeta_mean.append(di.data_read(di.zeta, di.rt_s_list[i], 3))

for i in range(0,len(di.rt_s_list)):
	phi_temp2 = []	
	chi_temp2 = []
	phi_dot_temp2 = []
	chi_dot_temp2 = []
	zeta_temp2 = []
	#for j in range(1,di.n_sample+1):
	for j in range(1,di.db_s_list[i]):
		if i in range(0,4):
			phi_temp1 = []
			chi_temp1 = []
			phi_dot_temp1 = []
			chi_dot_temp1 = []
			zeta_temp1 = []
			#print di.rt_s_list[i]
			#print 'i=',i
			#print 'j=',j
			phi_temp1.append(di.data_read(di.lat, di.rt_s_list[i], 1, db_size=di.db_s_list[i], db_row=j))
			chi_temp1.append(di.data_read(di.lat, di.rt_s_list[i], 2, db_size=di.db_s_list[i], db_row=j))
			phi_dot_temp1.append(di.get_fld_dot(di.lat, i, 3, a_scl, db_size=di.db_s_list[i], db_row=j))
			chi_dot_temp1.append(di.get_fld_dot(di.lat, i , 4, a_scl, db_size=di.db_s_list[i], db_row=j))
			zeta_temp1.append(di.data_read(di.lat, di.rt_s_list[i], 5, db_size=di.db_s_list[i], db_row=j))
			phi_temp2.append(phi_temp1)
			chi_temp2.append(chi_temp1)
			phi_dot_temp2.append(phi_dot_temp1)
			chi_dot_temp2.append(chi_dot_temp1)
			zeta_temp2.append(zeta_temp1)
		if i in range(4,6):
			phi_temp1 = []
			phi_dot_temp1 = []
			zeta_temp1 = []
			phi_temp1.append(di.data_read(di.lat, di.rt_s_list[i], 1, db_size=di.db_s_list[i], db_row=j))
			phi_dot_temp1.append(di.get_fld_dot(di.lat, i, 2, a_scl, db_size=di.db_s_list[i], db_row=j))
			zeta_temp1.append(di.data_read(di.lat, di.rt_s_list[i], 3, db_size=di.db_s_list[i], db_row=j))
			phi_temp2.append(phi_temp1)
			zeta_temp2.append(zeta_temp1)
	phi.append(phi_temp2)
	chi.append(chi_temp2)
	phi_dot.append(phi_dot_temp2)
	chi_dot.append(chi_dot_temp2)
	zeta.append(zeta_temp2)

##### Compute values of interest #####
# Using phi as a time variable it is not necessary to find p like this
end_infl = []
#p = []# formatted as[rt0min,rt0p,rt0max,...]
for i in range(0,len(di.rt_s_list)):
	end_infl.append(di.get_end_infl(a_scl[i],hub[i]))
#	p.append(di.get_p(phi[i],di.phi_p+di.b_1))
#	p.append(di.get_p(phi[i],di.phi_p))
#	p.append(di.get_p(phi[i],di.phi_p-di.b_1))	
	
##### Make plots #####
colours = ['k','b','r','g','c','m','y']
#phi_range = 5
phi_range = np.amin(phi[0][0])

n = 1
### Plot phi ###
pd.plot_phi_hom(n, phi, tau, a_scl, end_infl); n = n+1

### Plot chi ###
pd.plot_chi_hom(n,0, phi, chi, a_scl, end_infl, phi_range); n = n+1

### Plot chi from sampled points ###
pd.plot_chi_sampled(n, 1, phi, chi, phi_range); n = n+1
pd.plot_chi_sampled(n, 3, phi, chi, phi_range); n = n+1

### Plot chi_dot from sampled points ###
pd.plot_chi_dot_sampled(n, 1, phi, chi_dot, phi_range); n = n+1
pd.plot_chi_dot_sampled(n, 3, phi, chi_dot, phi_range); n = n+1

### Plot zeta from sampled points ###
pd.plot_zeta_sampled(n, 1, phi, zeta, phi_range); n = n+1
pd.plot_zeta_sampled(n, 3, phi, zeta, phi_range); n = n+1

### Plot chi difference from sampled points ###
pd.plot_chi_dif_sampled(n, 1, 3, phi, chi, phi_range); n = n+1

### Plot chi_dot difference from sampled points ###
pd.plot_chi_dot_dif_sampled(n, 1, 3, phi, chi_dot, phi_range); n = n+1

### Plot difference in phi from sampled points ###
pd.plot_phi_dif_sampled(n, 1, 3, phi, phi_range); n = n+1

### Plot difference in phi_dot from sampled points ###
pd.plot_phi_dot_dif_sampled(n, 1, 3, phi, phi_dot, phi_range); n = n+1

### Plot zeta difference along sampled trajectories ###
pd.plot_zeta_dif_sampled(n, 1, 3, phi, zeta, phi_range); n = n+1

### Plot mean zeta ###
pd.plot_zeta_mean(n, 1, phi, zeta_mean, phi_range); n = n+1
pd.plot_zeta_mean(n, 3, phi, zeta_mean, phi_range); n = n+1

### Plot horizon ###
pd.plot_horizon(n, 1, phi, a_scl, hub, phi_range); n = n+1





plt.show()



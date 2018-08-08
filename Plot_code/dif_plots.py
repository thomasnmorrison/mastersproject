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
# the ordering can be found in the data_in_mod3 module.
a_scl = []
hub = []
tau = []
zeta = []
phi = []
chi = []

for i in range(0,len(di.rt_list)):
	a_scl.append(di.data_read(di.ener, di.rt_list[i], 1))
	hub.append(di.get_hub(di.ener, di.rt_list[i]))
	tau.append(di.data_read(di.ener, di.rt_list[i], 0))
	zeta.append(di.data_read(di.zeta, di.rt_list[i], 3))
	phi.append(di.data_read(di.spec, di.rt_list[i], 1, db_size=di.db_list[i],scl=di.nl_list[i], root=True))
	if i == 0 or i == 3:
		chi.append(np.array([]))
	if i != 0 and i != 3:	# Exception for longitudinal only data
		chi.append(di.data_read(di.spec, di.rt_list[i], 5, db_size=di.db_list[i], db_row=1,scl=di.nl_list[i], root=True))

##### Read in array for vaules of k #####
k = np.array([])
k = di.data_read(di.spec, di.wdv, 0, cut_off=di.db-1)
# read in data for higher k modes of chi
chi_k = []
for i in range(0,5):
	chi_k.append(di.data_read(di.spec, di.rt_list[4], 5, db_size=di.db_list[4], db_row=i+1,scl=di.nl_list[4], root=True))
	chi_k.append(di.data_read(di.spec, di.rt_list[5], 5, db_size=di.db_list[5], db_row=i+1,scl=di.nl_list[5], root=True))

##### Compute some values of interest #####
end_infl = []
p = []# formatted as[rt0min,rt0p,rt0max,...]
k_p = []
for i in range(0,len(di.rt_list)):
	end_infl.append(di.get_end_infl(a_scl[i],hub[i]))
	p.append(di.get_p(phi[i],di.phi_p+di.b_1))
	p.append(di.get_p(phi[i],di.phi_p))
	p.append(di.get_p(phi[i],di.phi_p-di.b_1))
	# k_p is not properly calculated
	k_p.append(di.get_k_p(a_scl[i],hub[i],end_infl[i],p[3*i+1]))

##### Plotting #####
colors = ['k','b','r','g','c','m','y']
##### Plot Zeta #####
plt.figure(1)
x=[tau[4][0],tau[4][-1]]
y=[-1.2*np.amax(np.absolute(zeta[5]-zeta[4])),0.6*np.amax(np.absolute(zeta[5]-zeta[4]))]
plt.subplot(1,1,1)
zeta_scl = 10**0
plt.scatter(tau[4], (zeta[5]-zeta[4])*zeta_scl, s=2)
#plt.plot([tau[4][p[3*4]],tau[4][p[3*4]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][p[3*4+1]],tau[4][p[3*4+1]]], [y[0],y[1]], linestyle='--', color='b')
#plt.plot([tau[4][p[3*4+2]],tau[4][p[3*4+2]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][end_infl[4]],tau[4][end_infl[4]]], [y[0],y[1]], linestyle='--', color='b')
#plt.yscale('log')
plt.xlim(x[0],x[1])
plt.ylim(y[0],y[1])
plt.grid(True)
plt.title(r'Difference in $\zeta$ with and without $\Delta V$')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\zeta_{\Delta V} - \zeta_{\Delta V=0}$')

##### Plot Chi #####
# to do: plot several modes of chi
#plt.figure(2)
#plt.subplot(1,1,1)
#plt.scatter(tau[4], chi[5]-chi[4], s=2)
#plt.xlim(tau[1][0],tau[1][-1])
#plt.ylim(-1.2*np.amax(np.absolute(chi[5]-chi[4])),1.2*np.amax(np.absolute(chi[5]-chi[4])))
#plt.grid(True)
#plt.title(r'Difference in $\chi$ with and without $\Delta V$')
#plt.xlabel(r'$\tau$')
#plt.ylabel(r'$\chi_{\Delta V} - \chi_{\Delta V=0}$')

plt.figure(3)
x=[tau[4][0],tau[4][-1]]
y=[-1.2*np.amax(np.absolute(chi_k[1]-chi_k[0])),1.2*np.amax(np.absolute(chi_k[1]-chi_k[0]))]
plt.subplot(1,1,1)
for i in range(0,5):
	plt.scatter(tau[4], chi_k[2*i+1]-chi_k[2*i], s=2, color=colors[i%len(colors)], label=r'$k/(a_pH_p) = {0}$'.format(k[i]/(a_scl[4][p[3*4+1]]*hub[4][p[3*4+1]])))
#plt.plot([tau[4][p[3*4]],tau[4][p[3*4]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][p[3*4+1]],tau[4][p[3*4+1]]], [y[0],y[1]], linestyle='--', color='b')
#plt.plot([tau[4][p[3*4+2]],tau[4][p[3*4+2]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][end_infl[4]],tau[4][end_infl[4]]], [y[0],y[1]], linestyle='--', color='b')
plt.xlim(x[0],x[1])
plt.ylim(y[0],y[1])
plt.grid(True)
plt.title(r'Difference in $\chi$ with and without $\Delta V$')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\chi_{\Delta V} - \chi_{\Delta V=0}$')
plt.legend()
##### Plot Phi #####
plt.figure(4)
plt.subplot(2,1,1)
x=[tau[4][0],tau[4][-1]]
y=[-1.2*np.amax(np.absolute(phi[5])),1.2*np.amax(np.absolute(phi[5]))]
plt.scatter(tau[4], phi[5],s=2)
#plt.plot([tau[4][p[3*4]],tau[4][p[3*4]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][p[3*4+1]],tau[4][p[3*4+1]]], [y[0],y[1]], linestyle='--', color='b')
#plt.plot([tau[4][p[3*4+2]],tau[4][p[3*4+2]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][end_infl[4]],tau[4][end_infl[4]]], [y[0],y[1]], linestyle='--', color='b')
plt.xlim(x[0],x[1])
plt.ylim(y[0],y[1])
plt.grid(True)
plt.subplot(2,1,2)
y=[-1.2*np.amax(np.absolute(phi[5]-phi[4])),1.2*np.amax(np.absolute(phi[5]-phi[4]))]
plt.scatter(tau[4], phi[5]-phi[4],s=2)
#plt.plot([tau[4][p[3*4]],tau[4][p[3*4]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][p[3*4+1]],tau[4][p[3*4+1]]], [y[0],y[1]], linestyle='--', color='b')
#plt.plot([tau[4][p[3*4+2]],tau[4][p[3*4+2]]], [y[0],y[1]], linestyle='--', color='b')
plt.plot([tau[4][end_infl[4]],tau[4][end_infl[4]]], [y[0],y[1]], linestyle='--', color='b')
plt.xlim(tau[4][0],tau[4][-1])
plt.ylim(y[0],y[1])
plt.grid(True)

plt.show()

##### trapped_plot4.py #####
##### Script to make plots for trapped/trapped-instability runs focusing on sampled trajectories #####

# to do:
# to do:

#### Include packages ####
import matplotlib.pyplot as plt
import math
import mpmath as mp
import numpy as np
import data_in_mod4 as di

#### Define data lists ####
a_scl = []; hub = []
phi_h = []; chi_h = []
phi = []; chi = []; zeta = []

a_scl_rn = []; hub_rn = [];
phi_h_rn = []; chi_h_rn = []
phi_rn = []; chi_rn = []; zeta_rn = []
ind_p_rn = []


cut_off = 0
rn = [74, 76]
rt = 3
n_sample = 32
cut_off = 4095#0
#n_sample = di.db_s_list[rt]

for i in range(0,len(rn)):
	#### Read in data ####
	a_scl = []; hub = []
	phi_h = []; chi_h = []
	phi = []; chi = []; zeta = []
	di.clear_rt()
	di.init_rt(rn[i])
	di.init_param(rn[i])
	print "run initialized"
	di.data_init_scl_short(a_scl, hub, rt,cut=cut_off)
	print "scl read"
	di.data_init_fld_hom_short(phi_h, chi_h, rt,cut=cut_off)
	di.data_init_fld_sample(phi, chi, zeta, rt, cut=cut_off,sample_max=n_sample+1)
	print "fields read"

	#### Rescale a_scl to phi_p ####
	ind_p = di.get_p(phi_h,di.phi_p[0])
	print "ind_p = ", ind_p
	print "len(a_scl) = ", len(a_scl)
	a_scl_p = a_scl[ind_p]
	print "a_scl_p = ", a_scl_p
	print "phi_h[ind_p] = ", phi_h[ind_p]
	print "hub[ind_p] = ", hub[ind_p]

	#### Multi-run append ####
	if len(rn) != 1:
		a_scl_rn.append(a_scl); hub_rn.append(hub)
		phi_h_rn.append(phi_h); chi_h_rn.append(chi_h)
		phi_rn.append(phi); chi_rn.append(chi); zeta_rn.append(zeta)
		ind_p_rn.append(ind_p)


#di.clear_rt()
#di.init_rt(rn[1])
#di.init_param(rn[1])
#print "run initialized"
#di.data_init_scl_short(a_scl, hub, rt,cut=cut_off)
#print "scl read"
#di.data_init_fld_hom_short(phi_h, chi_h, rt,cut=cut_off)
#di.data_init_fld_sample(phi, chi, zeta, rt, cut=cut_off,sample_max=n_sample+1)
#print "fields read"

#### Rescale a_scl to phi_p ####
#ind_p = di.get_p(phi_h,di.phi_p[0])
#print "ind_p = ", ind_p
#print "len(a_scl) = ", len(a_scl)
#a_scl_p = a_scl[ind_p]
#print "a_scl_p = ", a_scl_p
#print "phi_h[ind_p] = ", phi_h[ind_p]
#print "hub[ind_p] = ", hub[ind_p]

#### Multi-run append ####
#a_scl_rn.append(a_scl); hub_rn.append(hub)
#phi_h_rn.append(phi_h); chi_h_rn.append(chi_h)
#phi_rn.append(phi); chi_rn.append(chi); zeta_rn.append(zeta)
#ind_p_rn.append(ind_p)

#### Make plots ####
colours = ['k','b','r','g','c','m','y']
nfig = 1

#### Difference Plots ####
if len(rn) != 1:
	fig = plt.figure(nfig)
	plt.subplot(3,1,1)
	#x = [np.amin(np.log(a_scl/a_scl_p)), np.amin(np.log(a_scl/a_scl_p))]
	#y = [-1.2*np.amax(np.absolute(np.array(chi))), 1.2*np.amax(np.absolute(np.array(chi)))]
	for j in range(0,n_sample):
		plt.scatter(np.log(a_scl/a_scl_p), np.array(chi_rn[1][j])-np.array(chi_rn[0][j]), s=2, color=colours[j%len(colours)])
	#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.ylabel(r'$\chi_{\Delta V-\Delta V_{trapped}}(M_{Pl})$', fontsize=18)
	#plt.legend()
	plt.grid(True)
	plt.title(r'Sampled Trajectories $(\Delta V-\Delta V_{trapped})$', fontsize=24)

	plt.subplot(3,1,2)
	#y = [-1.2*np.amax(np.absolute(np.array(phi))), 1.2*np.amax(np.absolute(np.array(phi)))]
	for j in range(0,n_sample):
		print j
		plt.scatter(np.log(a_scl/a_scl_p), np.array(phi_rn[1][j])-np.array(phi_rn[0][j]), s=2, color=colours[j%len(colours)])
	#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	#plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
	plt.ylabel(r'$\phi_{\Delta V-\Delta V_{trapped}}(M_{Pl})$', fontsize=18)
	#plt.legend()
	plt.grid(True)

	plt.subplot(3,1,3)
	#y = [-1.2*np.amax(np.absolute(np.array(zeta))), 1.2*np.amax(np.absolute(np.array(zeta)))]
	for j in range(0,n_sample):
		plt.scatter(np.log(a_scl/a_scl_p), np.array(zeta_rn[1][j])-np.array(zeta_rn[0][j]), s=2, color=colours[j%len(colours)])
	#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
	plt.ylabel(r'$\zeta_{\Delta V-\Delta V_{trapped}}$', fontsize=18)
	#plt.legend()
	plt.grid(True)
	nfig = nfig+1

#### Raw Plots ####
fig = plt.figure(nfig)
plt.subplot(3,1,1)
#x = [np.amin(np.log(a_scl/a_scl_p)), np.amin(np.log(a_scl/a_scl_p))]
#y = [-1.2*np.amax(np.absolute(np.array(chi))), 1.2*np.amax(np.absolute(np.array(chi)))]
for j in range(0,n_sample):
	plt.scatter(np.log(a_scl/a_scl_p), chi[j], s=2, color=colours[j%len(colours)])
#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
plt.ylabel(r'$\chi_{\Delta V}(M_{Pl})$', fontsize=18)
#plt.legend()
plt.grid(True)
plt.title(r'Sampled Trajectories $(\Delta V)$', fontsize=24)

plt.subplot(3,1,2)
#y = [-1.2*np.amax(np.absolute(np.array(phi))), 1.2*np.amax(np.absolute(np.array(phi)))]
for j in range(0,n_sample):
	print j
	plt.scatter(np.log(a_scl/a_scl_p), phi[j], s=2, color=colours[j%len(colours)])
#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
#plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
plt.ylabel(r'$\phi_{\Delta V}(M_{Pl})$', fontsize=18)
#plt.legend()
plt.grid(True)

plt.subplot(3,1,3)
#y = [-1.2*np.amax(np.absolute(np.array(zeta))), 1.2*np.amax(np.absolute(np.array(zeta)))]
for j in range(0,n_sample):
	plt.scatter(np.log(a_scl/a_scl_p), zeta[j], s=2, color=colours[j%len(colours)])
#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
plt.ylabel(r'$\zeta_{\Delta V}$', fontsize=18)
#plt.legend()
plt.grid(True)
nfig = nfig+1

#fig = plt.figure(nfig)
#plt.title(r'Sampled Trajectories')
#plt.subplot(3,1,1)
#x = [np.amin(np.log(a_scl/a_scl_p)), np.amin(np.log(a_scl/a_scl_p))]
#y = [-1.2*np.amax(np.absolute(np.array(chi))), 1.2*np.amax(np.absolute(np.array(chi)))]
#for j in range(0,n_sample):
#	plt.scatter(phi_h, chi[j], s=2, color=colours[j%len(colours)])
##plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
#plt.ylabel(r'$\chi_{\Delta V}(M_{Pl})$', fontsize=18)
##plt.legend()
#plt.grid(True)

#plt.subplot(3,1,2)
##y = [-1.2*np.amax(np.absolute(np.array(phi))), 1.2*np.amax(np.absolute(np.array(phi)))]
#for j in range(0,n_sample):
#	print j
#	plt.scatter(phi_h, phi[j], s=2, color=colours[j%len(colours)])
#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
#plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
#plt.ylabel(r'$\phi_{\Delta V}(M_{Pl})$', fontsize=18)
#plt.legend()
#plt.grid(True)

#plt.subplot(3,1,3)
#y = [-1.2*np.amax(np.absolute(np.array(zeta))), 1.2*np.amax(np.absolute(np.array(zeta)))]
#for j in range(0,n_sample):
#	plt.scatter(phi_h, zeta[j], s=2, color=colours[j%len(colours)])
#plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
#plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
#plt.ylabel(r'$\zeta_{\Delta V}$', fontsize=18)
#plt.legend()
#plt.grid(True)

#nfig = nfig+1

plt.show()


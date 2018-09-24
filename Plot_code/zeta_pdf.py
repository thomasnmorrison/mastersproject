##### zeta_pdf.py #####

##### Include packages #####
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import math
import mpmath as mp
import numpy as np
import data_in_mod3 as di
import potential_plots as pp

def get_index(ar_in, search, ascending=True):
	if ascending==True:
		for i in range(0,len(ar_in)):
			if ar_in[i] >= search:
				break
	return i

# Written to only read in one run
cut_off = 1210
cut = cut_off
rn = [42]
rt = 3
#n_sample = 32
n_sample = di.db_s_list[rt]
sample_max=512#32#
n_bin = 64
a_scl = []; hub = []; tau = []; phi = []; phi_dot = []; chi = []; chi_dot = []; zeta = []; zeta_mean = []
##### Long data read in #####
for i in range(0,len(rn)):
	di.clear_rt()
	di.init_rt(rn[i])
	a_scl = []; hub = []; tau = []; phi = []; phi_dot = []; chi = []; chi_dot = []; zeta = []; zeta_mean = []
	di.data_init_zeta_short(a_scl, hub, zeta, zeta_mean, 1, cut=cut_off)	
##### Choose lna snapshots #####
lna = np.arange(0,4.0,0.6)
##### Find indicies for lna snapshots #####
index_slice = []
for i in lna:
	j = get_index(np.log(a_scl[0]), i)
	index_slice.append(j)
print(index_slice)
##### Make histogram of zeta and normalize to PDF #####
zeta_hist = []
print('len(zeta): ',len(zeta))
print(len(zeta[0]))
print(len(zeta[0][0]))
print(len(zeta[0][0][0]))
x = [-np.amax(np.absolute(zeta[0])), np.amax(np.absolute(zeta[0]))]
for i in range(0,len(lna)):
	zeta_h = np.zeros(sample_max)
	for j in range(0,sample_max):
		index_temp = index_slice[i]
		zeta_h[j] = zeta[0][j][0][index_temp]
	H, bin_edge = np.histogram(zeta_h, bins=n_bin, range=(x[0],x[1]))
	print('len(H): ',len(H))
	print(H)
	#H = H/((x[1]-x[0])/n_bin*sample_max)
	H = H*1.0/sample_max/(bin_edge[1]-bin_edge[0])
	zeta_hist.append(H)
##### For each index compute zeta rms, mean, std #####
z_mean = []
z_rms = []
z_std = []
zeta_flat = np.zeros(sample_max)
for i in range(0,cut):
	for j in range(0,sample_max):
		zeta_flat[j] = zeta[0][j][0][i]
	z_mean.append(np.mean(zeta_flat))
	z_std.append(np.std(zeta_flat))
	z_rms.append(math.sqrt(z_mean[i]**2 + z_std[i]**2))
##### Plot zeta rms, mean, std vs lna #####
n=1
fig = plt.figure(n)
plt.subplot(1,1,1)
plt.scatter(np.log(a_scl[0]), z_mean, s=2)
plt.grid(True)
#plt.xlabel(r'$\phi$', fontsize=18)
#plt.ylabel(r'$1/2ln[P_{\phi \phi}P_{\dot{\phi} \dot{\phi}} - P_{\phi \dot{\phi}}^2]$(Machine Units)', fontsize=18)
n=n+1

fig = plt.figure(n)
plt.subplot(1,1,1)
plt.scatter(np.log(a_scl[0]), z_std, s=2)
plt.grid(True)
#plt.xlabel(r'$\phi$', fontsize=18)
#plt.ylabel(r'$1/2ln[P_{\phi \phi}P_{\dot{\phi} \dot{\phi}} - P_{\phi \dot{\phi}}^2]$(Machine Units)', fontsize=18)
n=n+1

fig = plt.figure(n)
plt.suptitle(r'$\zeta_{rms} (\Delta V=0)$', fontsize=20)
plt.subplot(1,1,1)
plt.scatter(np.log(a_scl[0]), z_rms, s=2)
plt.grid(True)
plt.xlim(0,5); plt.ylim(0,0.006)
plt.xlabel(r'$ln(a)$', fontsize=18)
plt.ylabel(r'$\zeta_{rms}$', fontsize=18)
n=n+1

fig = plt.figure(n)
plt.suptitle(r'$\zeta_{rms} (\Delta V=0)$', fontsize=20)
plt.subplot(1,1,1)
plt.scatter(np.log(a_scl[0]), z_rms, s=2)
plt.grid(True)
plt.xlim(0,5); plt.ylim(0,0.00006)
plt.xlabel(r'$ln(a)$', fontsize=18)
plt.ylabel(r'$\zeta_{rms}$', fontsize=18)
n=n+1

##### Plot PDF of zeta at different values of lna #####
fig = plt.figure(n)
plt.suptitle(r'$\zeta$ PDF $(\Delta V=0)$', fontsize=20)
#bins = np.arange(-np.amax(np.absolute(zeta[3])), np.amax(np.absolute(zeta[3])), 2.0*np.amax(np.absolute(zeta[3]))/n_bin)
bins = []
for i in range(0,len(bin_edge)-1):
	bins.append(bin_edge[i])
print(bins)
print(zeta_hist[1])
for i in range(0,6):
	plt.subplot(2,3,i+1)
	plt.bar(bins, zeta_hist[i],width=bin_edge[1]-bin_edge[0])
	#plt.grid(True)
	#plt.xlim(bin_edge[0],bin_edge[-1])
	plt.xlim(-0.015,0.015) 
	plt.ylim(0, 1.2*np.amax(zeta_hist[i]))#plt.ylim(0, 1.2)
	plt.title(r'$ln(a)={0}$'.format(lna[i]), fontsize=12)
	#plt.ylabel(r'$1/2ln[P_{\phi \phi}P_{\dot{\phi} \dot{\phi}} - P_{\phi \dot{\phi}}^2]$(Machine Units)', fontsize=18)
	plt.xticks([-0.01,0.0,0.01])
n=n+1












plt.show()

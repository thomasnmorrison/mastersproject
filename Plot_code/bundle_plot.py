##### bundle_plot.py #####
# to do: make animation in 2d

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

def get_level(lna_in, x, y):
	level = lna_in + 0*x + 0*y
	return level

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

# Written to only read in one run
cut_off = 1210
cut = cut_off
rn = [42]
rt = 3
#n_sample = 32
n_sample = di.db_s_list[rt]
sample_max=256
n_bin = 2**5
index_slice = []
chi_short = []
phi_short = []

##### Testing short data read in #####
di.clear_rt()
di.init_rt(rn[0])
a_scl.append(di.data_read(di.ener, di.rt_s_list[rt], 1, cut_off=cut))
lna = np.arange(0,4.5,0.5)
for i in lna:
	j = get_index(np.log(a_scl[0]), i)
	index_slice.append(j)
#for i in range(1,n_sample):
for i in range(1,sample_max+1):
	print(i)
	phi_temp1 = []
	chi_temp1 = []
	phi_temp1.append(di.data_read(di.lat, di.rt_s_list[rt], 1, db_size=n_sample, db_row=i,cut_off=cut, index_list=index_slice))
	chi_temp1.append(di.data_read(di.lat, di.rt_s_list[rt], 2, db_size=n_sample, db_row=i,cut_off=cut, index_list=index_slice))
	chi_short.append(chi_temp1)
	phi_short.append(phi_temp1)

print len(phi_short)
print len(phi_short[0])
print len(phi_short[0][0])

phi_slice = [] # mean value of phi for a given slice on ln(a)
phi_lim = [] # lower and upper limit on spread of phi at each slice of ln(a)
y = [-1.2*np.amax(np.absolute(chi_short)), 1.2*np.amax(chi_short)]
x = [y[0],y[1]]
bundle_hist = []
for i in range(0,len(lna)):
	#phi_h = phi_short[i]
	#chi_h = chi_short[i]
	#phi_h = np.zeros(n_sample-1)
	#chi_h = np.zeros(n_sample-1)
	phi_h = np.zeros(sample_max)
	chi_h = np.zeros(sample_max)
	#for j in range(0,n_sample-1):
	for j in range(0,sample_max):
		#print(j)
		#print(n_sample)
		phi_h[j] = phi_short[j][0][i]
		chi_h[j] = chi_short[j][0][i]
		#if y[0]>phi_h[j] or y[1]<phi_h[j]:
		#	print('range error: ', phi_h[j])
		#else:
		#	print('range good')
	# Subtract off mean phi
	phi_lim.append([np.amin(phi_h), np.amax(phi_h)])
	phi_mean = np.average(phi_h)
	phi_slice.append(phi_mean)
	#print(np.average(phi_h))
	phi_h = phi_h - phi_mean
	H, xedges, yedges = np.histogram2d(phi_h, chi_h, bins=n_bin, range=[[x[0],x[1]],[y[0],y[1]]]) # May need to do a 1/2 grid offset to centre points in mesh
	bundle_hist.append(H)
x_mesh, y_mesh = np.meshgrid(xedges, yedges)


##### Long data read in #####
#for i in range(0,len(rn)):
#	di.clear_rt()
#	di.init_rt(rn[i])
#	a_scl = []; hub = []; tau = []; phi = []; phi_dot = []; chi = []; chi_dot = []; zeta = []; zeta_mean = []
#	di.data_init_short(a_scl, hub, tau, phi, phi_dot, chi, chi_dot, zeta, zeta_mean, cut=cut_off)

#lna = np.arange(0,4.5,0.5) # Specify values for ln(a)
#phi_slice = [] # mean value of phi for a given slice on ln(a)
#phi_lim = [] # lower and upper limit on spread of phi at each slice of ln(a)
#y = [-np.amax(np.absolute(chi[3])), np.amax(chi[3])]
#x = [y[0],y[1]]
#bundle_hist = []
#for i in lna:
#	j = get_index(np.log(a_scl[3]), i)
#	phi_h = []
#	chi_h = []
#	for k in range(0,n_sample):
#		phi_h.append(phi[3][k][0][j]) # Probably makes more sense to subtract of the mean value of phi
#		chi_h.append(chi[3][k][0][j])
#	# Subtract off mean phi
#	phi_lim.append([np.amin(phi_h), np.amax(phi_h)])
#	phi_mean = np.average(phi_h)
#	phi_slice.append(phi_mean)
#	phi_h = phi_h - phi_mean
#	H, xedges, yedges = np.histogram2d(phi_h, chi_h, bins=n_bin, range=[[x[0],x[1]],[y[0],y[1]]]) # May need to do a 1/2 grid offset to centre points in mesh
#	bundle_hist.append(H)
#x_mesh, y_mesh = np.meshgrid(xedges, yedges)

##### Declare Plot #####
n = 1
fig = plt.figure(n)
ax = fig.add_subplot(2,1,1, projection='3d')
zlna = np.empty(x_mesh.shape)
ones_grid = np.ones(x_mesh.shape)
cscale = 5

##### Loop over lna and plot surfaces #####
for k in range(0,len(lna)):
	Z = bundle_hist[k]
	zlna = ones_grid*lna[k]
	#bundle_hist_plot = ax.plot_surface(x_mesh, y_mesh, zlna, facecolors=cm.binary(cscale*Z.transpose()/(np.amax(Z))), rstride=1, cstride=1, alpha=0.2)
	bundle_hist_plot = ax.plot_surface(zlna, y_mesh, x_mesh, facecolors=cm.OrRd(cscale*Z.transpose()/(np.amax(Z))), rstride=1, cstride=1, alpha=0.75, linewidth=0)
ax.view_init(elev=27, azim=-100)
ax.set_xlabel(r'$ln(a)$')
ax.set_ylabel(r'$\chi (M_{Pl})$')
ax.set_zlabel(r'$\phi (M_{Pl})$')

##### Create transverse potential plot #####
# Plot phi vs V_{,chi chi} with vertical slices at each lna slice
# Create arrays to plot
df = (phi_slice[0]-phi_slice[-1]+0.2)/(2**10)
f1 = np.arange(phi_slice[-1]-0.1, phi_slice[0]+0.1, df)
DVxx = pp.potential(1, 0, di.a2_list[rn[0]-1], di.b_1_list[rn[0]-1], di.phi_p_list[rn[0]-1], di.d_list[rn[0]-1], f1, f1, 0, 2)
ax2 = fig.add_subplot(2,1,2)
x = [phi_slice[-1], phi_slice[0]]
y = [-1.5*(np.amax(np.absolute(DVxx))), 1.5*(np.amax(np.absolute(DVxx)))]
ax2.plot(f1, DVxx, color='k', linewidth=1)
#for i in range(0,len(phi_slice)):
#print(len(phi_slice))
#print(len(lna))
for i in range(0,len(phi_slice)):
	#print i
	l=phi_lim[i][0]
	r=phi_lim[i][1]
	#print(l,r)
	#print(lna[i])
	ax2.fill_between(f1,y[0],y[1], where=(1-((f1>r)+(f1<l))), facecolor='blue', alpha=0.5)
	ax2.plot([l,l],[y[0],y[1]],linestyle='--', color='b')
	ax2.plot([r,r],[y[0],y[1]],linestyle='--', color='b')
plt.xlim(x[1],x[0]), plt.ylim(y[0],y[1])
plt.xlabel(r'$\phi(M_{Pl})$'); plt.ylabel(r'$\Delta V_{,\chi \chi}$')
plt.grid(True)


#n = n+1
#fig = plt.figure(n)
#bundle_hist[4]=bundle_hist[4][:-1,:-1]
#bundle_hist[0].transpose()
#bundle_c = plt.pcolormesh(x_mesh, y_mesh, bundle_hist[3].transpose(),cmap=cm.jet) # Needed to switch order of x and y for consistancy with np.histogram2d
#plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])

plt.show()

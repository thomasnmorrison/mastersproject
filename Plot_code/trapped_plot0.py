##### trapped_plot0.py #####
##### Script to make plots for trapped/trapped-instability runs focusing on a_scl, hub, and the potential #####

# to do: check hub units on plot
# to do: make nyquist frequency plot

#### Include packages ####
import matplotlib.pyplot as plt
import math
import mpmath as mp
import numpy as np
import data_in_mod4 as di
import potential_mod as pot

# function to calculate the power spectrum and write it to P_fld
def get_spectrum(k, P_fldfld, P_fld, L=True):
	# L: varies the definition of power spectrum to include the factor of L**3 or not
	for i in range(0,len(k[0])):
		if (L==False):
			#P_temp = (di.l_lattice[0]/(2.0*math.pi))**3*4.0*math.pi*k[0][i]**3*P_fldfld[i][0] #out of date format
			P_temp = (di.l_lattice[0]/(2.0*math.pi))**3*4.0*math.pi*k[0][i]**3*P_fldfld[i]
		if (L==True):
			#P_temp = k[0][i]**3*P_fldfld[i][0]/(2.0*math.pi**2) #out of date format
			P_temp = k[0][i]**3*P_fldfld[i]/(2.0*math.pi**2)
		P_fld.append(P_temp)
	return

# function to scale power wiht hub
def get_P_scaled(P_scaled, P_fld, hub):
	#P_scaled_temp = np.empty_like(P_fld)
	P_scaled_temp = []
	for i in range(0, len(P_fld)):
		P_scaled_temp.append(4.0*math.pi**2*P_fld[i]/(np.power(hub, 2)))
		#P_scaled_temp[i] = (4.0*math.pi**2*P_fld[i]/(np.power(hub, 2))
	P_scaled[:] = P_scaled_temp
	return

#### Define data lists ####
a_scl = []; hub = []
phi_h = []; chi_h = []

# Written to only read in one run
cut_off = 0#4095#
#cut = cut_off
rn = [76]
rt = 3

#### Read in data ####
di.clear_rt()
di.init_rt(rn[0])
di.init_param(rn[0])
print "run initialized"
di.data_init_scl_short(a_scl, hub, rt,cut=cut_off)
print "scl read"
a_scl=np.array(a_scl); hub=np.array(hub)
di.data_init_fld_hom_short(phi_h, chi_h, rt,cut=cut_off)
print "hom flds read"

#### Rescale a_scl to phi_p ####
ind_p = di.get_p(phi_h,di.phi_p[0])
print "ind_p = ", ind_p
print "len(a_scl) = ", len(a_scl)
a_scl_p = a_scl[ind_p]
print "a_scl_p = ", a_scl_p
print "phi_h[ind_p] = ", phi_h[ind_p]
print "hub[ind_p] = ", hub[ind_p]

#### Make plots ####
end_infl = di.get_end_infl(np.array(a_scl), np.array(hub))
print "end_infl = ", end_infl
print "at the plots"

nfig=1

# Average phi plot
fig = plt.figure(nfig)
plt.subplot(1,1,1)
plt.title(r'$\langle \phi \rangle$ with $\Delta V$', fontsize=22)
x = [np.amin(np.log(a_scl/a_scl_p)), np.amax(np.log(a_scl/a_scl_p))]
y = [0, 1.2*np.amax(np.absolute(phi_h))]
plt.scatter(np.log(a_scl/a_scl_p), phi_h, s=2)
plt.plot([0,0],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
if (end_infl < len(a_scl)) and (end_infl>0):
	plt.plot([math.log(a_scl[end_infl]/a_scl_p),math.log(a_scl[end_infl]/a_scl_p)],[y[0],y[1]],linestyle='--', color='g',label=r'end infl')
plt.xlim(x[0],x[1])
plt.ylim(y[0],y[1])
plt.grid(True)
plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
plt.ylabel(r'$\langle \phi \rangle (M_{Pl})$', fontsize=18)
plt.legend()
nfig=nfig+1

# Hubble plot
fig = plt.figure(nfig)
plt.subplot(1,1,1)
grid_spacing = float(di.l_lattice[0])/((di.n_lattice[0])**(1./3.))*di.mpl[0]
plt.title(r'$1/aH$', fontsize=22)
x = [np.amin(np.log(a_scl/a_scl_p)), np.amax(np.log(a_scl/a_scl_p))]
y = [0, 2.*np.amax(1./(a_scl*hub))]
plt.scatter(np.log(a_scl/a_scl_p), 1./(a_scl*hub), s=2)
plt.plot([0,0],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
plt.plot([x[0],x[1]],[grid_spacing,grid_spacing],linestyle='--',color='r',label=r'lattice spacing')
plt.plot([x[0],x[1]],[grid_spacing*(di.n_lattice[0])**(1./3.),grid_spacing*(di.n_lattice[0])**(1./3.)],linestyle='--',color='k',label=r'lattice size')
if (end_infl < len(a_scl)) and (end_infl>0):
	plt.plot([math.log(a_scl[end_infl]/a_scl_p),math.log(a_scl[end_infl]/a_scl_p)],[y[0],y[1]],linestyle='--', color='g',label=r'end infl')
plt.plot()
plt.xlim(x[0],x[1])
plt.ylim(y[0],y[1])
plt.grid(True)
plt.xlabel(r'$\mathrm{ln}(a/a_p)$', fontsize=18)
plt.ylabel(r'$(aH)^{-1}(M_{Pl}^{-1})$', fontsize=18)
plt.legend()
nfig=nfig+1

# Potential plot
pot.init_mpl(di.mpl)
pot.init_iopo(2, 1)

x = [np.amin(phi_h), np.amax(phi_h)]
y = [-0.01,0.01]
pot.plot_potential(nfig,x,y,di.m2[0], di.lambda_chi[0], di.g2[0], di.beta2[0], di.phi_p[0], 0,0)
nfig=nfig+1

pot.plot_potential_2d(nfig,x,di.m2[0], di.lambda_chi[0], di.g2[0], di.beta2[0], di.phi_p[0], 0,2)
plt.xlim(x[1],x[0])
nfig=nfig+1

plt.show()


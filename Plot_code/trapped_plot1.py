##### trapped_plot1.py #####
##### Script to make plots for trapped/trapped-instability runs focusing on power spectrum #####

# to do: find how spectrum is organized in therms of indicies
# to do: fix get_spectrum to properly get indicies
# to do: scale a to a_p

#### Include packages ####
import matplotlib.pyplot as plt
import math
import mpmath as mp
import numpy as np
import data_in_mod4 as di
#import plot_def_mod as pd

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
k_spec = []
P_phiphi = []; P_dphidphi = []; P_dphiphi_re = []; P_dphiphi_im = []
Det_P = []
n_k_phi = []
P_chichi = []; P_dchidchi = []; P_dchichi_re = []; P_dchichi_im = []
Det_P_chi = []
n_k_chi = []
P_zz = []; P_dzdz = []; P_dzz_re = []; P_dzz_im = []
Det_P_z = []
n_k_z = []

# Written to only read in one run
#cut_off = 1210
#cut = cut_off
rn = [78]
rt = 3

#### Read in data ####
di.clear_rt()
di.init_rt(rn[0])
di.init_param(rn[0])
print "run initialized"
di.data_init_scl_short(a_scl, hub, rt)
print "scl read"
di.data_init_fld_hom_short(phi_h, chi_h, rt)
print "hom flds read"
di.data_init_spec(k_spec, P_phiphi, P_dphidphi, P_dphiphi_re, P_dphiphi_im, P_chichi, P_dchidchi, P_dchichi_re, P_dchichi_im, P_zz, rt, a_scl)
print "spec read"
print "len(P_fldfld) = ", len(P_phiphi)
print "len(P_fldfld[0]) = ", len(P_phiphi[0])
#print "len(P_fldfld[0][0]) = ", len(P_phiphi[0][0])
#### Calculate power #####
P_phi = []; P_chi = []; P_z = []
get_spectrum(k_spec, P_phiphi, P_phi)
get_spectrum(k_spec, P_chichi, P_chi)
get_spectrum(k_spec, P_zz, P_z)
print "len(P_phi) = ", len(P_phi)
print "len(P_phi[0]) = ", len(P_phi[0])

#### Rescale a_scl to phi_p ####
ind_p = di.get_p(phi_h,di.phi_p[0])
print "ind_p = ", ind_p
print "len(a_scl) = ", len(a_scl)
a_scl_p = a_scl[ind_p]
print "a_scl_p = ", a_scl_p
print "phi_h[ind_p] = ", phi_h[ind_p]
print "hub[ind_p] = ", hub[ind_p]
#### Scaled power #####
# choose a specific index at which to calculated the quantity log_10(4*pi*P_phi/(H**2)) array over values of k
P_phi_scaled = []; P_chi_scaled = []; P_z_scaled = []
get_P_scaled(P_phi_scaled, P_phi, hub)
print "len(P_phi_scaled) = ", len(P_phi_scaled)
print "len(P_phi_scaled[0]) = ", len(P_phi_scaled[0])
P_phi_scaled = np.transpose(P_phi_scaled)
print "len(P_phi_scaled) = ", len(P_phi_scaled)
print "len(P_phi_scaled[0]) = ", len(P_phi_scaled[0])

#### Find k_star #####

#### Make plots ####
# list of a_scl values at which to plot (measured by normalizing a=1 at phi_p)
a_list = [2.2]
ind_list = []
# this has a bug since get_p assumes a decreasing function, use index_finder
for a in a_list:
	#ind_list.append(di.get_p(a_scl, a/a_scl_p))
	print (np.argmin(np.absolute(a_scl-a/a_scl_p)))
	ind_list.append(np.argmin(np.absolute(a_scl-a*a_scl_p)))

end_infl = di.get_end_infl(np.array(a_scl), np.array(hub))
print "end_infl = ", end_infl

print "ind_list = ", ind_list
print "a_scl[ind_list[0]] = ", a_scl[ind_list[0]]
print "hub[ind_list[0]] = ", hub[ind_list[0]]
print "at the plots"
nfig=1
#### Plot Power of phi ####
fig = plt.figure(nfig)
plt.subplot(1,1,1)
plt.scatter(np.log(k_spec), np.log10(P_phi_scaled[ind_list[0]]), s=2)
#plt.xlim()
#plt.ylim()
plt.grid(True)
plt.xlabel(r'$ln(k)$', fontsize=18)
plt.ylabel(r'$log_{10}(\frac{4\pi^2P_{\phi}}{H^2})$', fontsize=18)
nfig=nfig+1

plt.show()







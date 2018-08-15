##### plot_def_mod.py #####

##### Include packages #####
import matplotlib.pyplot as plt
import math
import mpmath as mp
import numpy as np
import data_in_mod3 as di

colours = ['k','b','r','g','c','m','y']
##### Define plots #####

### Plot horizon ###
def plot_horizon(fig, rt, phi, a_scl, hub, phi_range):
# to do: add a line for the Nyquist frequency
# to do: subplot with ln(a)
# to do: plot phi_p	
	grid_spacing = float(di.l_lattice)/((di.nl_list[rt])**(1./3.))
	plt.figure(fig)
	plt.subplot(2,1,1)
	plt.title(r'Horizon')
	x = [phi[rt][0][0][0], phi_range]
	y = [0, 1.5/np.amin(a_scl[rt]*hub[rt])]
	plt.scatter(phi[rt][0], 1/(a_scl[rt]*hub[rt]), s=2)
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.plot([x[0],x[1]],[grid_spacing,grid_spacing],color='b',label='lattice spacing')
	plt.plot([x[0],x[1]],[2*grid_spacing,2*grid_spacing],color='g',label=' twice lattice spacing')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$1/aH(M_{Pl}/mpl)$', fontsize=18)
	plt.legend(loc='upper right')
	plt.grid(True)

	plt.subplot(2,1,2)
	x=[math.log(a_scl[0][0]), math.log(a_scl[0][-1])]
	y = [0, 1.5/np.amin(a_scl[rt]*hub[rt])]
	plt.scatter(np.log(a_scl[0]), 1/(a_scl[rt]*hub[rt]), s=2)
	plt.plot([math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p)]), math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p)])],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p+2*di.b_1)]), math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p+2*di.b_1)])],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p-2*di.b_1)]), math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p-2*di.b_1)])],[y[0],y[1]],linestyle='--', color='r')
	plt.plot([x[0],x[1]],[grid_spacing, grid_spacing],color='b',label='lattice spacing')
	plt.plot([x[0],x[1]],[2*grid_spacing,2*grid_spacing],color='g',label=' twice lattice spacing')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'ln$(a)$', fontsize=18); plt.ylabel(r'$1/aH(M_{Pl}/mpl)$', fontsize=18)
	plt.legend(loc='upper right')
	plt.grid(True)

	return

### Plot mean phi ###
def plot_phi_hom(fig, phi, tau, a_scl, end_infl):
	plt.figure(fig)
	plt.subplot(2,1,1) #plot phi vs tau
	plt.title(r'Homogeneous Calculation of $\phi$')
	x = [tau[0][0], tau[0][-1]]
	y = [np.amin(phi[0][0]), 1.1*np.amax(phi[0][0])]
	plt.scatter(tau[0], phi[0][0], s=2)
	#plt.plot([tau[0][end_infl[0]],tau[0][end_infl[0]]], [y[0],y[1]], linestyle='--', color='b',label='End of inflation')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\tau$', fontsize=18); plt.ylabel(r'$\phi (M_{Pl})$', fontsize=18)
	#plt.legend()
	plt.grid(True)

	plt.subplot(2,1,2) #plot phi vs ln(a)
	x=[math.log(a_scl[0][0]), math.log(a_scl[0][-1])]
	plt.scatter(np.log(a_scl[0]), phi[0][0], s=2)
	#plt.plot([math.log(a_scl[0][end_infl[0]]),math.log(a_scl[0][end_infl[0]])], [y[0],y[1]], linestyle='--', color='b',label='End of inflation')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'ln($a$)', fontsize=18); plt.ylabel(r'$\phi (M_{Pl})$', fontsize=18)
	#plt.legend()
	plt.grid(True)

	return

### Plot mean chi ###
def plot_chi_hom(fig, rt, phi, chi, a_scl, end_infl, phi_range):
	plt.figure(fig)
	rt = 0
	plt.subplot(2,1,1) #plot chi vs tau
	plt.title(r'Homogeneous Calculation of $\chi$')
	x=[math.log(a_scl[rt][0]), math.log(a_scl[rt][-1])]
	y = [0.5*np.amin(chi[rt][0]), 1.5*np.amax(chi[rt][0])]
	plt.scatter(np.log(a_scl[rt]), chi[rt][0], s=2)
	#plt.plot([math.log(a_scl[rt][end_infl[rt]]),math.log(a_scl[rt][end_infl[rt]])], [y[0],y[1]], linestyle='--', color='b',label='End of inflation')
	plt.plot([math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p)]), math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p)])],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p+2*di.b_1)]), math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p+2*di.b_1)])],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p-2*di.b_1)]), math.log(a_scl[0][di.index_finder(phi[rt][0][0],di.phi_p-2*di.b_1)])],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'ln($a$)', fontsize=18); plt.ylabel(r'$\chi (M_{Pl})$', fontsize=18)
	plt.legend()
	plt.grid(True)

	plt.subplot(2,1,2) #plot chi vs phi
	x = [phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	plt.scatter(phi[rt][0], chi[rt][0], s=2)
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	#plt.plot([math.log(a_scl[rt][end_infl[rt]]),math.log(a_scl[rt][end_infl[rt]])], [y[0],y[1]], linestyle='--', color='b',label='End of inflation')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\phi (M_{Pl})$', fontsize=18); plt.ylabel(r'$\chi (M_{Pl})$', fontsize=18)
	plt.legend()
	plt.grid(True)

	return

### Plot sampled chi ###
def plot_chi_sampled(fig, rt, phi, chi, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	if rt == 1:
		plt.title(r'Sampled Trajectories $\chi_{\Delta V=0}$')
	if rt == 3:
		plt.title(r'Sampled Trajectories $\chi_{\Delta V}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	y = [1.2*np.amin(chi[rt]), 1.2*np.amax(chi[rt])]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], chi[rt][j], s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	if rt == 1:
		plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\chi_{\Delta V=0}(M_{Pl})$', fontsize=18)
	if rt == 3:
		plt.xlabel(r'$\phi_{\Delta V}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\chi_{\Delta V}(M_{Pl})$', fontsize=18)
	plt.legend()
	plt.grid(True)

	return

### Plot sampled chi_dot ###
def plot_chi_dot_sampled(fig, rt, phi, chi_dot, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	if rt == 1:
		plt.title(r'Sampled Trajectories $\dot{\chi}_{\Delta V=0}$')
	if rt == 3:
		plt.title(r'Sampled Trajectories $\dot{\chi}_{\Delta V}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	print 'len(chi_dot) = ', len(chi_dot)
	print 'len(chi_dot[rt]) = ', len(chi_dot[rt])
	print 'len(chi_dot[rt][0]) = ', len(chi_dot[rt][0])
	print 'len(chi_dot[rt][0][0]) = ', len(chi_dot[rt][0][0])
	y = [-1.2*np.amax(np.absolute(np.array(chi_dot[rt]))), 1.2*np.amax(np.absolute(np.array(chi_dot[rt])))]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], chi_dot[rt][j], s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	if rt == 1:
		plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\dot{\chi}_{\Delta V=0}(M_{Pl}^2/mpl)$', fontsize=18)
	if rt == 3:
		plt.xlabel(r'$\phi_{\Delta V}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\dot{\chi}_{\Delta V}(M_{Pl}^2/mpl)$', fontsize=18)
	plt.legend()
	plt.grid(True)

	return

### Plot sampled difference in chi ###
def plot_chi_dif_sampled(fig, rt, rt2, phi, chi, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	plt.title(r'Difference in Sampled Trajectories $\chi_{\Delta V}-\chi_{\Delta V=0}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	y = [-1.2*np.amax(np.absolute(np.array(chi[rt2])-np.array(chi[rt]))), 1.2*np.amax(np.absolute(np.array(chi[rt2])-np.array(chi[rt])))]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], np.array(chi[rt2][j])-np.array(chi[rt][j]), s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\chi_{\Delta V}-\chi_{\Delta V=0}(M_{Pl})$', fontsize=18)
	plt.legend()
	plt.grid(True)

	return

### Plot sampled difference in phi
def plot_phi_dif_sampled(fig, rt, rt2, phi, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	plt.title(r'Difference in Sampled Trajectories $\phi_{\Delta V}-\phi_{\Delta V=0}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	y = [-1.5*np.amax(np.absolute(np.array(phi[rt2])-np.array(phi[rt]))), 1.5*np.amax(np.absolute(np.array(phi[rt2])-np.array(phi[rt])))]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], np.array(phi[rt2][j])-np.array(phi[rt][j]), s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\phi_{\Delta V}-\phi_{\Delta V=0}(M_{Pl}^2/mpl)$', fontsize=18)
	plt.legend(loc='upper left')
	plt.grid(True)

### Plot sampled difference in chi_dot ###
def plot_chi_dot_dif_sampled(fig, rt, rt2, phi, chi_dot, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	plt.title(r'Difference in Sampled Trajectories $\dot{\chi}_{\Delta V}-\dot{\chi}_{\Delta V=0}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	y = [-1.2*np.amax(np.absolute(np.array(chi_dot[rt2])-np.array(chi_dot[rt]))), 1.2*np.amax(np.absolute(np.array(chi_dot[rt2])-np.array(chi_dot[rt])))]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], np.array(chi_dot[rt2][j])-np.array(chi_dot[rt][j]), s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\dot{\chi}_{\Delta V}-\dot{\chi}_{\Delta V=0}(M_{Pl}^2/mpl)$', fontsize=18)
	plt.legend()
	plt.grid(True)

	return

### Plot sampled difference in phi_dot ###
def plot_phi_dot_dif_sampled(fig, rt, rt2, phi, phi_dot, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	plt.title(r'Difference in Sampled Trajectories $\dot{\phi}_{\Delta V}-\dot{\phi}_{\Delta V=0}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	y = [-1.2*np.amax(np.absolute(np.array(phi_dot[rt2])-np.array(phi_dot[rt]))), 1.2*np.amax(np.absolute(np.array(phi_dot[rt2])-np.array(phi_dot[rt])))]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], np.array(phi_dot[rt2][j])-np.array(phi_dot[rt][j]), s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\dot{\phi}_{\Delta V}-\dot{\phi}_{\Delta V=0}(M_{Pl}^2/mpl)$', fontsize=18)
	plt.legend()
	plt.grid(True)

	return

### Plot sampled difference in zeta ###
def plot_zeta_dif_sampled(fig, rt, rt2, phi, zeta, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	plt.title(r'Difference in Sampled Trajectories $\zeta_{\Delta V}-\zeta_{\Delta V=0}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	#y = [1.2*np.amin(np.array(zeta[rt2])-np.array(zeta[rt])), 1.2*np.amax(np.array(zeta[rt2])-np.array(zeta[rt]))]
	y = [-1.3*np.amax(np.absolute(np.array(zeta[rt2])-np.array(zeta[rt]))), 1.3*np.amax(np.absolute(np.array(zeta[rt2])-np.array(zeta[rt])))]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], np.array(zeta[rt2][j])-np.array(zeta[rt][j]), s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\zeta_{\Delta V}-\zeta_{\Delta V=0}$', fontsize=18)
	plt.legend(loc='lower left')
	plt.grid(True)

	return

### Plot sampled zeta ###
def plot_zeta_sampled(fig, rt, phi, zeta, phi_range):
	
	plt.figure(fig)
	plt.subplot(1,1,1)
	if rt == 1:
		plt.title(r'Sampled Trajectories $\zeta_{\Delta V=0}$')
	if rt == 3:
		plt.title(r'Sampled Trajectories $\zeta_{\Delta V}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	y = [-1.5*np.amax(np.absolute(zeta[rt])), 1.5*np.amax(np.absolute(zeta[rt]))]
	for j in range(0,di.db_s_list[rt]-1):
		plt.scatter(phi[rt][j], zeta[rt][j], s=2, color=colours[j%len(colours)])
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	if rt == 1:
		plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\zeta_{\Delta V=0}$', fontsize=18)
	if rt == 3:
		plt.xlabel(r'$\phi_{\Delta V}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\zeta_{\Delta V}$', fontsize=18)
	plt.legend(loc='upper left')
	plt.grid(True)

	return

### Plot mean zeta ###
def plot_zeta_mean(fig, rt, phi, zeta_mean, phi_range):
	plt.figure(fig)
	plt.subplot(1,1,1)
	if rt == 1:
		plt.title(r'Lattice Averaged $\zeta_{\Delta V=0}$')
	if rt == 3:
		plt.title(r'Lattice Averaged $\zeta_{\Delta V}$')
	x = [phi[rt][0][0][0], phi_range]#[phi[rt][0][0][0], phi[rt][0][0][end_infl[rt]]]
	y = [1.5*np.amin(zeta_mean[rt]), 1.5*np.amax(zeta_mean[rt])]
	plt.scatter(phi[rt][0], zeta_mean[rt], s=2)
	plt.plot([di.phi_p,di.phi_p],[y[0],y[1]],linestyle='--', color='b',label=r'$\phi_p$')
	plt.plot([di.phi_p+2*di.b_1,di.phi_p+2*di.b_1],[y[0],y[1]],linestyle='--', color='r',label=r'$\phi_p \pm 2b$')
	plt.plot([di.phi_p-2*di.b_1,di.phi_p-2*di.b_1],[y[0],y[1]],linestyle='--', color='r')
	plt.xlim(x[0], x[1]); plt.ylim(y[0], y[1])
	if rt == 1:
		plt.xlabel(r'$\phi_{\Delta V=0}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\zeta_{\Delta V=0}$', fontsize=18)
	if rt == 3:
		plt.xlabel(r'$\phi_{\Delta V}(M_{Pl})$', fontsize=18); plt.ylabel(r'$\zeta_{\Delta V}$', fontsize=18)
	plt.legend()
	plt.grid(True)
	return







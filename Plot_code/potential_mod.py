##### potential_mod.py #####

##### Include packages #####
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import math
import mpmath as mp
import numpy as np
import data_in_mod4 as di

# to do:
# to do:

# Note on variable names:
# Throughout f1 is used to refer to the inflaton and f2 is used to refer to the transverse field
# When running with po = 1 or io  = 2 the following variable names are used (compare potential_mod.f90)
# m2 -> lambda1
# g2 -> a2
# beta2 -> b

##### Define parameter definitions #####
mpl = []
io = []
po = []

##### Define initialization functions #####
def init_mpl(mpl_in):
	mpl[:] = mpl_in
	print "mpl = ",mpl
	return

def init_iopo(io_in, po_in):
	io[:] = []
	po[:] = []
	io.append(io_in)
	po.append(po_in)
	print "io = ",io
	print "po = ",po
	return

##### Define potential functions #####
def v(lambda1, lambda2, f1, f2, units=True):
	print "mpl = ",mpl
	print "io = ",io
	print "po = ",po
	if io[0]==1:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)
		v = 0.25*lambda1*f1**4 + 0.25*lambda2*f2**4
	if io[0]==2:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)
		v = 0.5*lambda1*f1**2 + 0.25*lambda2*f2**4
	return v

def v_phi(lambda1, lambda2, f1, f2, units=True):
	if io[0]==1:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)	
		v_p = lambda1*f1**3
	if io[0]==2:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			#lambda2 = lambda2/(mpl[0]**2)	
		v_p = lambda1*f1
	return v_p

def v_chi(lambda1, lambda2, f1, f2, units=True):
	if io[0]==1:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)	
		v_c = lambda2*f2**3
	if io[0]==2:
		if units == True:
#			lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)	
		v_c = lambda2*f2**3
	return v_c

def v_phiphi(lambda1, lambda2, f1, f2, units=True):
	if io[0]==1:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)			
		v_pp = 3.0*lambda1*f1**2
	if io[0]==2:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			#lambda2 = lambda2/(mpl[0]**2)			
		v_pp = lambda1	
	return v_pp

def v_chichi(lambda1, lambda2, f1, f2, units=True):
	if io[0]==1:
		if units == True:
			lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)	
		v_cc = 3.0*lambda2*f2**2
	if io[0]==2:
		if units == True:
			#lambda1 = lambda1/(mpl[0]**2)
			lambda2 = lambda2/(mpl[0]**2)	
		v_cc = 3.0*lambda2*f2**2
	return v_cc

def amp(a2,b, units=True):
	# Setting units=True will put the ampltude in units of M_Pl for fields in M_Pl
	amplitude = a2*math.exp(0.5)/(b)
	if units == True:
		amplitude = amplitude/(mpl[0]**2)
	return amplitude

def delta_v(a2, b, phi_p, f1, f2):
	if po[0]==0:
		dv = 0
	if po[0]==1:
		a2 = a2/(mpl[0]**2); b = b/(mpl[0]**2)	# Units
		dv = (0.5*a2*(f1-phi_p)**2 - 0.5*b)*f2**2
	if po[0]==2:
		dv = -amp(a2,b)*(f1-phi_p)*np.exp(-0.5*(f1-phi_p)**2/(b**2))*f2**2
	return dv

def delta_v_phi(a2, b, phi_p, f1, f2):
	if po[0]==0:
		dv_p = 0
	if po[0]==1:
		a2 = a2/(mpl[0]**2); b = b/(mpl[0]**2)	# Units
		dv_p = a2*(f1-phi_p)*f2**2
	if po[0]==2:
		dv_p = -amp(a2,b)*(1.0-(f1-phi_p)**2/(b**2))*np.exp(-0.5*(f1-phi_p)**2/(b**2))*f2**2
	return dv_p

def delta_v_chi(a2, b, phi_p, f1, f2):
	if po[0]==0:
		dv_c = 0
	if po[0]==1:
		a2 = a2/(mpl[0]**2); b = b/(mpl[0]**2)	# Units
		dv_c = (a2*(f1-phi_p)**2 - b)*f2
	if po[0]==2:
		dv_c = -2.0*amp(a2,b)*(f1-phi_p)*np.exp(-0.5*(f1-phi_p)**2/(b**2))*f2
	return dv_c

def delta_v_phiphi(a2, b, phi_p, f1, f2):
	if po[0]==0:
		dv_pp = 0
	if po[0]==1:
		a2 = a2/(mpl[0]**2); b = b/(mpl[0]**2)	# Units
		dv_pp = a2*f2**2
	if po[0]==2:
		dv_pp = amp(a2,b)/(b**2)*((f1-phi_p)*(2.0-(f1-phi_p)*(f1-phi_p)/(b**2))+(f1-phi_p))*np.exp(-0.5*(f1-phi_p)**2/(b**2))*f2**2
	return dv_pp

def delta_v_phichi(a2, b, phi_p, f1, f2):
	if po[0]==0:
		dv_pc = 0
	if po[0]==1:
		a2 = a2/(mpl[0]**2); b = b/(mpl[0]**2)	# Units
		dv_pc = 2*a2*(f1-phi_p)*f2
	if po[0]==2:
		dv_pc = -2.0*amp(a2,b)*(1.0-(f1-phi_p)*(f1-phi_p)/(b**2))*np.exp(-0.5*(f1-phi_p)**2/(b**2))*f2
	return dv_pc

def delta_v_chichi(a2, b, phi_p, f1, f2):
	if po[0]==0:
		dv_cc = 0
	if po[0]==1:
		a2 = a2/(mpl[0]**2); b = b/(mpl[0]**2)	# Units
		dv_cc = a2*(f1-phi_p)**2 - b
	if po[0]==2:
		dv_cc = -2.0*amp(a2,b)*(f1-phi_p)*np.exp(-0.5*(f1-phi_p)**2/(b**2))
	return dv_cc

def potential(lambda1, lambda2, a2, b, phi_p, f1, f2, ind1, ind2):
	if ind1 == 0 and ind2 == 0:
		pot = v(lambda1, lambda2, f1, f2) + delta_v(a2, b, phi_p, f1, f2)
	elif ind1 == 1 and ind2 == 0:
		pot = v_phi(lambda1, lambda2, f1, f2) + delta_v_phi(a2, b, phi_p, f1, f2)
	elif ind1 == 0 and ind2 == 1:
		pot = v_chi(lambda1, lambda2, f1, f2) + delta_v_chi(a2, b, phi_p, f1, f2)
	elif ind1 == 2 and ind2 == 0:
		pot = v_phiphi(lambda1, lambda2, f1, f2) + delta_v_phiphi(a2, b, phi_p, f1, f2)
	elif ind1 == 1 and ind2 == 1:
		pot = delta_v_phichi(a2, b, phi_p, f1, f2)
	elif ind1 == 0 and ind2 == 2:
		pot = v_chichi(lambda1, lambda2, f1, f2) + delta_v_chichi(a2, b, phi_p, f1, f2)
	else:
		print 'invalid entry of ind in potential'
		return
	return pot

##### Define Potential Plot #####
def plot_potential(fig,x,y,lambda1, lambda2, a2, b, phi_p, ind1,ind2):
	f1 = np.arange(x[0],x[1],(x[1]-x[0])/2.**7)
	f2 = np.arange(y[0],y[1],(y[1]-y[0])/2.**7)
	f1,f2 = np.meshgrid(f1,f2)
	pot = potential(lambda1, lambda2, a2, b, phi_p, f1, f2, ind1, ind2)
	fig1 = plt.figure(fig)
	ax = fig1.add_subplot(1,1,1, projection='3d')
	pot_surface = ax.plot_surface(f1,f2,pot,cmap=cm.coolwarm, antialiased=True)
	return

def plot_potential_2d(fig,x,lambda1, lambda2, a2, b, phi_p, ind1,ind2):
	f1 = np.arange(x[0],x[1],(x[1]-x[0])/2.**8)
	pot = potential(lambda1, lambda2, a2, b, phi_p, f1, 0, ind1, ind2)
	fig1 = plt.figure(fig)
	ax = fig1.add_subplot(1,1,1)
	ax.plot(f1, pot)
	ax.grid(True)
	return


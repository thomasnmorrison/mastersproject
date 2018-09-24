##### zeta_moment.py #####
##### Script to plot the moments of zeta #####

# to do: MAke 4th moment connected plot

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

##### Read in data #####
# Written to only read in one run
cut_off = 1210
cut = cut_off
rn = [43]
rt = 3
#n_sample = 32
n_sample = di.db_s_list[rt]
sample_max=512#32#
pi = math.atan(1.0)
a_scl = []; hub = []; zeta = []; zeta_mean = []; zeta_moment = []

for i in range(0,len(rn)):
	di.clear_rt()
	di.init_rt(rn[i])
	di.data_init_zeta_short(a_scl, hub, zeta, zeta_mean, rt, zeta_moment_in=zeta_moment, cut=cut_off, moment_max=4)	

##### Initialize lna #####
lna = np.log(a_scl[0])
lnH = np.log(hub[0])

##### Take finite difference derivatives #####
# Make use of numpy.gradient
print('len(zeta_moment) = ', len(zeta_moment))
print('len(zeta_moment[0]) = ', len(zeta_moment[0]))
print('len(zeta_moment[0][0]) = ', len(zeta_moment[0][0]))
print('len(zeta_moment[0][0][0]) = ', len(zeta_moment[0][0][0]))
dzmoment = []
dzmoment2 = []
##### NUMPY INCORRECTLY CALCULATES GRADIENTS #####
dzmoment.append(np.gradient(zeta_moment[0][0][0], lna))
dzmoment[0] = np.delete(dzmoment[0],0)
dzmoment.append(np.gradient(zeta_moment[1][0][0], lna))
dzmoment[1] = np.delete(dzmoment[1],0)
dzmoment.append(np.gradient(zeta_moment[2][0][0], lna))
dzmoment[2] = np.delete(dzmoment[2],0)
dzmoment.append(np.gradient(zeta_moment[3][0][0], lna))
dzmoment[3] = np.delete(dzmoment[3],0)
#dzmoment4_con = np.gradient(zeta_moment[3][0][0]-(zeta_moment[1][0][0])**2, lna)
#dzmoment4_con = np.delete(dzmoment4_con,0)
epsilon = -1.0*np.gradient(lnH, lna)

##### Corrected derivative calculation #####
dzmoment2.append(np.diff(zeta_moment[0][0][0])/np.diff(lna))
dzmoment2.append(np.diff(zeta_moment[1][0][0])/np.diff(lna))
dzmoment2.append(np.diff(zeta_moment[2][0][0])/np.diff(lna))
dzmoment2.append(np.diff(zeta_moment[3][0][0])/np.diff(lna))
dzmoment4_con = (np.diff(zeta_moment[3][0][0]-zeta_moment[1][0][0]*zeta_moment[1][0][0]))/np.diff(lna)
epsilon2 = -1.0*np.diff(lnH)/np.diff(lna)

#print(dzmoment[0])
#print(dzmoment[1])
#print(dzmoment[2])
#print(dzmoment[3])
lna = np.delete(lna,0)
hub = np.delete(hub,0)
lnH = np.delete(lnH,0)
epsilon = np.delete(epsilon,0)
print('H: ', hub)
print('lnH: ', lnH)
print('lna: ', lna)
print('epsilon: ', epsilon)

##### Make Plots #####
n=1

#fig = plt.figure(n)
#for i in range(1,4):
#	plt.subplot(3,1,i)
#	x= [0,6]
#	y= [-0.1*np.amax(np.absolute(dzmoment[i])), 1.2*np.amax(np.absolute(dzmoment[i]))]
#	plt.scatter(lna, dzmoment[i],s=2, color='k')
#	plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
#	plt.grid(True)
	#plt.xlabel(r'$\phi$', fontsize=18)
	#plt.ylabel(r'$1/2ln[P_{\phi \phi}P_{\dot{\phi} \dot{\phi}} - P_{\phi \dot{\phi}}^2]$(Machine Units)', fontsize=18)
#n=n+1

fig = plt.figure(n)
plt.suptitle(r'$\frac{\mathrm{d}\langle \Delta \zeta^2 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=24)
plt.subplot(2,1,1)
x= [0,6]
#y= [-0.1*np.amax(np.absolute(dzmoment[1])), 1.2*np.amax(np.absolute(dzmoment[1]))]
y= [-0.1*np.amax(np.absolute(dzmoment2[1])), 1.2*np.amax(np.absolute(dzmoment2[1]))]
#plt.scatter(lna, dzmoment[1],s=2, color='k', label=r'$\frac{\mathrm{d}\langle \Delta \zeta^2 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$')
#plt.scatter(lna, (hub/(2*math.sqrt(2)*pi*np.sqrt(epsilon)))**2 ,s=2, color='r', label=r'$(\frac{H}{2\pi \sqrt{2} \sqrt{\epsilon} M_{Pl}})^2$')
plt.scatter(lna, dzmoment2[1], s=2, color='k', label=r'$\frac{\mathrm{d}\langle \Delta \zeta^2 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$')
plt.scatter(lna, (hub/(2*math.sqrt(2)*pi*np.sqrt(epsilon2)))**2 ,s=2, color='r', label=r'$(\frac{H}{2\pi \sqrt{2} \sqrt{\epsilon} M_{Pl}})^2$')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.legend(loc='upper right')
#plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
#plt.ylabel(r'$\frac{\mathrm{d}\langle \Delta \zeta^2 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=18)
plt.subplot(2,1,2)
y = [-0.1*np.amax(np.absolute((hub/(2*math.sqrt(2)*pi*np.sqrt(epsilon2)))**2)), 1.2*np.amax(np.absolute((hub/(2*math.sqrt(2)*pi*np.sqrt(epsilon2)))**2))]
#plt.scatter(lna, hub/(np.sqrt(epsilon)) ,s=2, color='r')
plt.scatter(lna, (hub/(2*math.sqrt(2)*pi*np.sqrt(epsilon2)))**2 ,s=2, color='r')
#plt.scatter(lna, dzmoment[1],s=2, color='k')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
plt.ylabel(r'$(\frac{H}{2\pi \sqrt{2} \sqrt{\epsilon} M_{Pl}})^2$', fontsize=18)
n=n+1

fig = plt.figure(n)
plt.suptitle(r'Testing', fontsize=24)
plt.subplot(4,1,1)
x= [0,1]
#y= [-1.2*np.amax(np.absolute(lnH)), 1.2*np.amax(np.absolute(lnH))]
y=[-12.5,-12]
plt.scatter(lna, lnH,s=2, color='k')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
#plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
plt.ylabel(r'$ln(H)$', fontsize=20)
plt.subplot(4,1,2)
#y= [-1.2*np.amax(np.absolute(hub)), 1.2*np.amax(np.absolute(hub))]
y=[0.000003, 0.000006]
plt.scatter(lna, hub,s=2, color='k')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.ylabel(r'$H$', fontsize=20)
plt.subplot(4,1,3)
y= [-1.2*np.amax(np.absolute(epsilon)), 1.2*np.amax(np.absolute(epsilon))]
#y=[-0.01,0.15]
plt.scatter(lna, epsilon, s=2, color='k')
plt.scatter(lna, epsilon2, s=2, color='r')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.ylabel(r'$\epsilon (=-\frac{\mathrm{d ln}H}{\mathrm{d ln}a})$', fontsize=20)
#plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
plt.subplot(4,1,4)
#y= [-1.2*np.amax(np.absolute(np.power(epsilon,-1))), 1.2*np.amax(np.absolute(np.power(epsilon,-1)))]
y=[0,5000]
plt.scatter(lna, np.power(epsilon,-1), s=2, color='k')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.ylabel(r'$\epsilon^{-1}$', fontsize=20)
plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
n=n+1

fig = plt.figure(n)
plt.suptitle(r'$\frac{\mathrm{d}\langle \Delta \zeta^3 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=24)
plt.subplot(1,1,1)
x= [0,6]
y= [-0.1*np.amax(np.absolute(dzmoment2[2])), 1.2*np.amax(np.absolute(dzmoment2[2]))]
plt.scatter(lna, dzmoment2[2],s=2, color='k')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
plt.ylabel(r'$\frac{\mathrm{d}\langle \Delta \zeta^3 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=18)
#plt.subplot(2,1,2)
#y = [-0.1*np.amax(np.absolute((hub/(2*math.sqrt(2)*pi))**2)), 1.2*np.amax(np.absolute((hub/(2*math.sqrt(2)*pi))**2))]
#plt.scatter(lna, (hub/(2*math.sqrt(2)*pi))**2 ,s=2)
#plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
#plt.grid(True)
#plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
#plt.ylabel(r'$(\frac{H}{2\pi \sqrt{2} M_{Pl}})^2$', fontsize=18)
n=n+1

fig = plt.figure(n)
plt.suptitle(r'$\frac{\mathrm{d}\langle \Delta \zeta^4 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=24)
plt.subplot(2,1,1)
x= [0,6]
y= [-0.1*np.amax(np.absolute(dzmoment2[3])), 1.2*np.amax(np.absolute(dzmoment2[3]))]
plt.scatter(lna, dzmoment2[3],s=2, color='k')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
#plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
plt.ylabel(r'$\frac{\mathrm{d}\langle \Delta \zeta^4 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=18)
plt.subplot(2,1,2)
y = [-0.1*np.amax(np.absolute((hub/(2*math.sqrt(2)*pi))**2)), 1.2*np.amax(np.absolute((hub/(2*math.sqrt(2)*pi))**2))]
plt.scatter(lna, (hub/(2*math.sqrt(2)*pi))**2 ,s=2)
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
plt.ylabel(r'$(\frac{H}{2\pi \sqrt{2} M_{Pl}})^2$', fontsize=18)
n=n+1

fig = plt.figure(n)
plt.suptitle(r'$\frac{\mathrm{d}\langle \Delta \zeta^4 \rangle- \langle \Delta \zeta^2 \rangle^2}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=24)
plt.subplot(1,1,1)
x= [0,6]
y= [-0.1*np.amax(np.absolute(dzmoment4_con)), 1.2*np.amax(np.absolute(dzmoment4_con))]
plt.scatter(lna, dzmoment4_con,s=2, color='k', label=r'$\frac{\mathrm{d}\langle \Delta \zeta^4 \rangle-\langle \Delta \zeta^2 \rangle^2}{\mathrm{d}\,\mathrm{ln}(a)}$')
plt.scatter(lna, dzmoment2[3],s=2, color='r', label=r'$\frac{\mathrm{d}\langle \Delta \zeta^4 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$')
plt.xlim(x[0],x[1]); plt.ylim(y[0],y[1])
plt.grid(True)
plt.xlabel(r'$\mathrm{ln}(a)$', fontsize=18)
#plt.ylabel(r'$\frac{\mathrm{d}\langle \Delta \zeta^4 \rangle}{\mathrm{d}\,\mathrm{ln}(a)}$', fontsize=18)
plt.legend(loc='upper right')
n=n+1

plt.show()

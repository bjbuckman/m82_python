
def synch_spectrum_z(los_file, suffix, XX):
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	# import matplotlib
	# matplotlib.use('Agg')
	# import matplotlib.pyplot as plt
	# from matplotlib.colors import LogNorm

	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	DD = 3520 #kpc
	x_integrated = 60 #arcsecs
	
	a=py.open(los_file)
	data=a[0].data
	head=a[0].header
	
	x0=head['CRVAL1']
	y0=head['CRVAL2']
	x0=min(x0,y0)
	
	dx=head['CDELT1']
	numx=head['NAXIS1']
	nu0=head['CRVAL3']
	dnu=head['CDELT3']
	numnu=head['NAXIS3']

	xxf=np.linspace(0,numx-1,numx)
	xxf*=dx
	xxf+=x0

	nuf=np.linspace(0,numnu-1,numnu)
	nuf=dnu**nuf
	nuf*=nu0

	data=np.swapaxes(data,2,1)
	data=np.swapaxes(data,0,1)
	
	#Picking out near XX
	FF=data[abs(xxf-XX)<=x_integrated/2]
	FF=np.average(FF,axis=0)
	FF*= x_integrated/radian2arcsec
	# FF=np.sum(FF,axis=0)
	# FF*= dx/radian2arcsec
	
	np.savetxt(los_file+'_z'+suffix+'.dat', xxf)
	np.savetxt(los_file+'_z'+suffix+'_nu.dat', nuf)
	np.savetxt(los_file+'_z'+suffix+'_spec.dat', FF)
	
	return xxf, nuf, FF
	#########################
	## Plotting
	#########################
	# fig = plt.figure()
	
	# f1, ax = plt.subplots(1,2)
	
	# X, Y = np.meshgrid(nuf,Z)
	# V = np.logspace(-30, -10, 20)
	# dat = FF.T
	# CF = ax[0].contourf(X,Y,dat,norm = LogNorm(), levels=V)
	# ax[0].set_xlim([1e6,1e17])
	# ax[0].set_xscale('log')
	# ax[0].set_xlabel(r'$\nu$ (Hz)')
	# ax[0].set_ylim([-1,1])
	# ax[0].set_ylabel(r'$z$ (kpc)')
	# ax[0].set_title(r'$F(E,z)$')
	 
	# plt.colorbar(CF,ax=ax[0],ticks=V,pad=0)
	
	# data = PP.T
	# X1, Y1 = np.meshgrid(nup,Z)
	# V1 = np.linspace(-4, 1, 21)
	# CF1 = ax[1].contourf(X1,Y1,data, levels=V1)
	# ax[1].set_xlim([1e6,1e17])
	# ax[1].set_xscale('log')
	# ax[1].set_xlabel(r'$\nu$ (Hz)')
	# ax[1].set_ylim([-1,1])
	# ax[1].set_title(r'$p(E,z) = 1-2*\alpha$')
	# plt.colorbar(CF1,ticks=V1,ax=ax[1],pad=0)
	
	# plt.suptitle('Synchroton ('+gp_name+')')
	# plt.savefig(gp_name+'_synchrotron_z.png')

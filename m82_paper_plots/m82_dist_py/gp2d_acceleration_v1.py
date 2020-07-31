import sys

def CR_proton_accel_z(emissfilename,RR):
	''' Calculates acceleration of gas from CR momentum loss
		OUTPUT = Acceleration*gas_n_density
		UNITS  = cm/s^2 g/cm^3 = g/cm^2/s^2 '''
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	zmin = head['CRVAL2']
	log_e_min = head['CRVAL3'] #log10(e_min)
	
	xdel = head['CDELT1']
	zdel = head['CDELT2']
	log_efactor = head['CDELT3'] #factor
	
	xnum = head['NAXIS1']
	znum = head['NAXIS2']
	enum = head['NAXIS3']
	nucnum = head['NAXIS4']
	
	Z = np.linspace(0, znum-1, znum)
	Z *= zdel
	Z += zmin
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	luminosity = np.zeros(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
	
	mp = 1.67e-24 #g
	#Picking nuclei
	# data = scidata[0][0]*0.
	for A in range(1,5):
		iA = 0
		for AAA in AA:
			if AAA == A:
				iA+=1
		dat = scidata[abs(AA)==A]
		
		mass = 938.28*A #MeV/c^2
		beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
		
		#Adding all protons together
		# if iA > 1:
		dat = np.sum(dat, axis=0)
		
		#picking out radius
		dat = np.swapaxes(dat, 0, 2)*beta/c*(1.+mass/energy)/3.
		dat = dat[(R-RR) <= 2*xdel]
		dat = np.average(dat, axis=0)
		
		##
		## dat = I_p
		
		dat *= 4.*math.pi
		if A == 1:
			data = dat
		else:
			data+= dat
		
	#integrate over E
	intIp = (data[:,1:]+data[:,:-1])*np.log(energy[1:]/energy[:-1])/2
	intIp = np.sum(intIp, axis=1)*1e6*ev2ergs
	
	#############################
	# gradient of U_cr
	#############################
	Z_grad = (Z[1:]+Z[:-1])/2
	grad_Ip = (intIp[1:]-intIp[:-1])/(zdel*kpc2cm)
	
	accel = -1*grad_Ip/mp
	
	return Z_grad, accel
	
def CR_proton_U_z(emissfilename,RR):
	''' Calculates acceleration of gas from CR momentum loss
		OUTPUT = Acceleration*gas_density
		UNITS  = cm/s^2 g/cm^3 = g/cm^2/s^2 '''
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	zmin = head['CRVAL2']
	log_e_min = head['CRVAL3'] #log10(e_min)
	
	xdel = head['CDELT1']
	zdel = head['CDELT2']
	log_efactor = head['CDELT3'] #factor
	
	xnum = head['NAXIS1']
	znum = head['NAXIS2']
	enum = head['NAXIS3']
	nucnum = head['NAXIS4']
	
	Z = np.linspace(0, znum-1, znum)
	Z *= zdel
	Z += zmin
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	luminosity = np.zeros(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
	
	#Picking electrons
	
	for ii in range(0,nucnum):
		scidata[ii] = scidata[ii]*AA[ii]
	
	scidata = scidata[AA>0]
	mass = mass[AA>0]
	AA = AA[AA>0]
	
	dat = scidata[:,:,:,(R-RR) < 1.5*xdel]
	dat = np.average(dat, axis=3)
	
	beta = np.zeros(dat.shape[0:2])
	EFAC = np.zeros(dat.shape[0:2])
	for ii in range(0,AA.shape[0]):
		beta[ii,:] = np.sqrt(1.-1./(1.+energy/mass[ii])**2.)
		EFAC[ii,:] = (1.+mass[ii]/energy)
	
	for iz in range(0,znum):
		dat[:,:,iz] = dat[:,:,iz]/beta*EFAC
	
	#Adding all protons together
	dat = np.sum(dat, axis=0)
	# dat = np.swapaxes(dat, 0, 2)#*1/3
	
	
	##
	## dat = I_p
	
	dat *= 4.*math.pi/c
	#dat = dat.T
	
	#integrate over E
	intIp = np.sum(dat, axis=0)*math.log(energy[1]/energy[0])*1e6*ev2ergs
	
	return Z, intIp
	
def CR_electron_accel_z(emissfilename,RR):
	''' Calculates acceleration of gas from CR momentum loss
		OUTPUT = Acceleration*gas_density
		UNITS  = cm/s^2 g/cm^3 = g/cm^2/s^2 '''
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	zmin = head['CRVAL2']
	log_e_min = head['CRVAL3'] #log10(e_min)
	
	xdel = head['CDELT1']
	zdel = head['CDELT2']
	log_efactor = head['CDELT3'] #factor
	
	xnum = head['NAXIS1']
	znum = head['NAXIS2']
	enum = head['NAXIS3']
	nucnum = head['NAXIS4']
	
	Z = np.linspace(0, znum-1, znum)
	Z *= zdel
	Z += zmin
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	luminosity = np.zeros(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
	
	#Picking electrons
	dat = scidata[abs(AA)==0]
	
	AA = 1
	mass = 0.5109989461 #MeV/c^2
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	me = 9.11e-28 #g 1.66e-24 #g
	mp = 1.66e-24 #g
	
	#Adding all protons together
	dat = np.sum(dat, axis=0)
	
	#picking out radius
	dat = np.swapaxes(dat, 0, 2)*beta/c*(1+mass/energy)/3.
	dat = dat[(R-RR) <= 2*xdel]
	dat = np.average(dat, axis=0)
	
	##
	## dat = I_p
	
	dat *= 4.*math.pi
	
	#integrate over E
	intIp = (dat[:,1:]+dat[:,:-1])*np.log(energy[1:]/energy[:-1])/2
	intIp = np.sum(intIp, axis=1)*1e6*ev2ergs
	
	#############################
	# gradient of U_cr
	#############################
	Z_grad = (Z[1:]+Z[:-1])/2
	grad_Ip = (intIp[1:]-intIp[:-1])/(zdel*kpc2cm)
	
	accel = -1*grad_Ip/mp
	
	return Z_grad, accel
	
def CR_electron_U_z(emissfilename,RR):
	''' Calculates acceleration of gas from CR momentum loss
		OUTPUT = Acceleration*gas_density
		UNITS  = cm/s^2 g/cm^3 = g/cm^2/s^2 '''
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	zmin = head['CRVAL2']
	log_e_min = head['CRVAL3'] #log10(e_min)
	
	xdel = head['CDELT1']
	zdel = head['CDELT2']
	log_efactor = head['CDELT3'] #factor
	
	xnum = head['NAXIS1']
	znum = head['NAXIS2']
	enum = head['NAXIS3']
	nucnum = head['NAXIS4']
	
	Z = np.linspace(0, znum-1, znum)
	Z *= zdel
	Z += zmin
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	luminosity = np.zeros(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
	
	#Picking electrons
	dat = scidata[abs(AA)==0]
	
	AA = 1
	mass = 0.5109989461 #MeV/c^2
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	me = 9.11e-28 #g 1.66e-24 #g
	mp = 1.66e-24 #g
	
	#Adding all protons together
	dat = np.sum(dat, axis=0)
	
	#picking out radius
	dat = np.swapaxes(dat, 0, 2)/beta/c*(1+mass/energy)#*1/3
	dat = dat[(R-RR) <= 2*xdel]
	dat = np.average(dat, axis=0)
	
	##
	## dat = I_p
	
	dat *= 4.*math.pi
	#dat = dat.T
	
	#integrate over E
	intIp = np.sum(dat, axis=1)*math.log(energy[1]/energy[0])*1e6*ev2ergs
	
	return Z, intIp
	
def m82_accel_r():
	''' Calculates acceleration as a function of r for m82 
		OUTPUT = Acceleration cm/s^2 '''
	import numpy as np
	
	kpc2cm = 3.08568025e21 #cm/kpc
	G = 6.67e-8 #cm^3/g/s^2
	
	Mass_E = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/mass_dyn/mass_dyn_E.txt').T
	Mass_W = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/mass_dyn/mass_dyn_W.txt').T
	
	Mass = (Mass_E + Mass_W)/2

	MASS = Mass[1]*2.0e33 #g
	r_MASS = Mass[0]*kpc2cm #cm
	
	accel = -1*G*MASS/r_MASS**2
	
	R = np.append(-1*np.flipud(Mass[0]), Mass[0])
	ACCEL = np.append(-1*np.flipud(accel), accel)
	
	return R, ACCEL

def m82_potential():
	''' Calculates gravitational potential as a function of r for m82 
		OUTPUT = Potential ergs g^-1 '''
	import numpy as np
	from scipy import interpolate
	
	kpc2cm = 3.08568025e21 #cm/kpc
	G = 6.67e-8 #cm^3/g/s^2
	
	Mass_E = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/mass_dyn/mass_dyn_E.txt').T
	Mass_W = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/mass_dyn/mass_dyn_W.txt').T
	
	Mass = (Mass_E + Mass_W)/2
	
	z = np.linspace(Mass[0].min(),Mass[0].max(), num = 500)
	f_mass = interpolate.interp1d(Mass[0],Mass[1], kind='cubic')#, fill_value=0)
	
	MASS = f_mass(z)
	MASS*= 2.0e33 #g
	r_MASS = z*kpc2cm #cm
	
	U = -1*G*MASS/r_MASS
	
	R = np.append(-1*np.flipud(z), z)
	UU = np.append(-1*np.flipud(U), U)
	
	return R, UU
	
def m82_accel_all(emissfilename,RR):
	import numpy as np
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
	import gas_density
	import magnetic_field
	
	r1, proton = CR_proton_accel_z(emissfilename,RR)
	r2, grav   = m82_accel_r()
	r3, electron = CR_electron_accel_z(emissfilename,RR)
	r4, magnetic = magnetic_field.magnetic_gradient(1)
	
	zmax = 3.0
	znum = 1000
	z = np.linspace(-1*zmax, zmax, znum)
	
	r0w = 0.200 #kpc
	r1w = 1.0
	z0w = 0.050 #kpc
	z1w = 4.0

	rwindw = 0.05

	v0_rw = 1000.
	v0_zw = 0.
	densityw = 1.e-2
	mass_lossw = 5.0
	
	n_density = np.zeros(len(z))
	#print z
	for iz in range(0, len(z)):
		n_density[iz] = gas_density.m82_wind_radial( 0.0, 0.0, z[iz], r0w, r1w, rwindw, z0w, z1w, v0_rw, densityw, mass_lossw)
	
	BB = np.interp(z, r4, magnetic)
	BB/= n_density
	PP = np.interp(z, r1, proton)
	PP/= n_density
	PE = np.interp(z, r3, electron)
	PE/= n_density
	GR = np.interp(z, r2, grav  )
	
	#PLOTS
	fig = plt.figure()
	plt.plot(z ,PP, color='red')
	plt.plot(z ,PE, color='blue')
	plt.plot(z ,GR, color='green')
	plt.plot(z ,BB, color='brown')
	plt.plot(z ,GR+PP, color='black')
	plt.xlim([-1*zmax,zmax])
	plt.xlabel(r'$z$ (kpc)')
	plt.ylabel(r'$a_z$ (cm s$^{-1}$)')
	
	plt.title('m82 Accel_all')
	plt.savefig('m82_accel_all.png')
	

def CR_electron_N(emissfilename):	
	''' Calculates electron number density
		OUTPUT = Acceleration*gas_density
		UNITS  = cm/s^2 g/cm^3 = g/cm^2/s^2 '''
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	zmin = head['CRVAL2']
	log_e_min = head['CRVAL3'] #log10(e_min)
	
	xdel = head['CDELT1']
	zdel = head['CDELT2']
	log_efactor = head['CDELT3'] #factor
	
	xnum = head['NAXIS1']
	znum = head['NAXIS2']
	enum = head['NAXIS3']
	nucnum = head['NAXIS4']
	
	Z = np.linspace(0, znum-1, znum)
	Z *= zdel
	Z += zmin
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	luminosity = np.zeros(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
	
	#Picking electrons
	dat = scidata[abs(AA)==0]
	
	AA = 1
	mass = 0.5109989461 #MeV/c^2
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	me = 9.11e-28 #g 1.66e-24 #g
	mp = 1.66e-24 #g
	
	#Adding all electrons together
	dat = np.sum(dat, axis=0)
	
	#pick location
	dat = dat[:,0,Z==0.]
	
	##
	## dat = I_p
	
	dat *= 4.*math.pi/c/energy
	#dat = dat.T
	
	#integrate over E
	# intIp = np.sum(dat, axis=1)*math.log(energy[1]/energy[0])*1e6*ev2ergs
	
	return energy, dat
	
def synch_power(emissfilename):
	'''
	integrated_flux(emissfile, fileout, objectdist)
	- Finds total Integrated flux of GALPROP emiss file
	
	Units of output file:
	gamma - [MeV], [MeV^2 cm^-2 s^-1 MeV^-1]
	synch - [Hz ], [erg   cm^-2 s^-1 Hz^-1 ]
	
	Units of return value:
	gamma - [erg s^-1]
	synch - [erg s^-1]
	'''
	import os
	import pyfits as py
	import math
	import numpy as np
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	erg_per_MeV = 1.6022e-6
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data[0]
	
	#Distances in kpc
	xmin = head['CRVAL1']
	zmin = head['CRVAL2']
	log_e_min = head['CRVAL3'] #log10(e_min)
	
	xdel = head['CDELT1']
	zdel = head['CDELT2']
	log_efactor = head['CDELT3'] #factor
	
	xnum = head['NAXIS1']
	znum = head['NAXIS2']
	enum = head['NAXIS3']
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin
	
	Z = zmin+np.array(range(0,znum))*zdel
	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	z_0 = np.argwhere(Z==0)
	# print scidata.shape
	scidata = scidata[:,z_0,0]
	
	#integrate over r
	# r_int = math.pi*2.*R*xdel
	# r_int[0] = math.pi*(xdel/2.)**2.
	# scidata = np.sum(scidata*r_int, axis=1)
	
	# nan_array = np.isnan(scidata)
	#scidata = scidata[~nan_array]
	#energy = energy[~nan_array]
	
	#unit conversion and distance projection
	# if 'synchrotron' in emissfilename:
	# total_emission = np.sum(energy[~nan_array]*scidata[~nan_array])*4*math.pi*math.log(efactor)
	# else:
		# total_emission = np.sum(scidata[~nan_array])*4*math.pi*kpc2cm**3*math.log(efactor)*erg_per_MeV
	
	scidata *= 4*math.pi
	
	# np.savetxt(fileout+'.dat', [energy, scidata], '%e')
	
	return energy, scidata[:,0,0]

def freef_absorp(emissfilename):
	'''
	integrated_flux(emissfile, fileout, objectdist)
	- Finds total Integrated flux of GALPROP emiss file
	
	Units of output file:
	gamma - [MeV], [MeV^2 cm^-2 s^-1 MeV^-1]
	synch - [Hz ], [erg   cm^-2 s^-1 Hz^-1 ]
	
	Units of return value:
	gamma - [erg s^-1]
	synch - [erg s^-1]
	'''
	import os
	import pyfits as py
	import math
	import numpy as np
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	erg_per_MeV = 1.6022e-6
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	zmin = head['CRVAL2']
	log_e_min = head['CRVAL3'] #log10(e_min)
	
	xdel = head['CDELT1']
	zdel = head['CDELT2']
	log_efactor = head['CDELT3'] #factor
	
	xnum = head['NAXIS1']
	znum = head['NAXIS2']
	enum = head['NAXIS3']
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin
	
	Z = zmin+np.array(range(0,znum))*zdel
	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	z_0 = np.argwhere(Z==0)
	# print scidata.shape
	scidata = scidata[:,z_0,0]
	
	#integrate over r
	# r_int = math.pi*2.*R*xdel
	# r_int[0] = math.pi*(xdel/2.)**2.
	# scidata = np.sum(scidata*r_int, axis=1)
	
	# nan_array = np.isnan(scidata)
	#scidata = scidata[~nan_array]
	#energy = energy[~nan_array]
	
	#unit conversion and distance projection
	# if 'synchrotron' in emissfilename:
	# total_emission = np.sum(energy[~nan_array]*scidata[~nan_array])*4*math.pi*math.log(efactor)
	# else:
		# total_emission = np.sum(scidata[~nan_array])*4*math.pi*kpc2cm**3*math.log(efactor)*erg_per_MeV
	
	# scidata *= 4*math.pi
	
	# np.savetxt(fileout+'.dat', [energy, scidata], '%e')
	
	return energy, scidata[:,0,0]
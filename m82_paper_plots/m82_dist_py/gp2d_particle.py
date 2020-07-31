import sys

def particle_luminosity(emissfilename, particle):
	##
	## Takes galprop injection file and calculates luminosity
	##
	## Resturns luminosity (erg s^-1)
	##
	
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	from scipy.interpolate import griddata
	
	A = 1
	if particle == 'electron':
		mass = 0.5109989461 #MeV/c^2
	elif particle == 'positron':
		mass = 0.5109989461 #MeV/c^2
	elif particle == 'proton':
		mass = 	938.28
	elif particle == 'Helium' :
		mass = 3723.7#4.*938.28
		A = 4
	else:
		print 'failed'
	
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	
	energy = emin*efactor**np.array(range(0,enum))
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	scidata = np.sum(scidata, axis=1)*zdel #sum over z
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	scidata = np.sum(scidata*r_int, axis=1)
	scidata = np.sum(scidata*energy**2./beta)
	luminosity = scidata*log_efactor*math.log(10.)*4.*math.pi/c*1.e6*ev2ergs*kpc2cm**3.*A
	
	return luminosity
	
def particle_luminosity_z(emissfilename, particle):
	##
	## Takes galprop injection file and calculates luminosity
	##
	## Resturns luminosity (erg s^-1)
	##
	
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	from scipy.interpolate import griddata
	
	A = 1
	if particle == 'electron':
		mass = 0.5109989461 #MeV/c^2
	elif particle == 'positron':
		mass = 0.5109989461 #MeV/c^2
	elif particle == 'proton':
		mass = 	938.28
	elif particle == 'Helium' :
		mass = 3723.7#4.*938.28
		A = 4
	else:
		print 'failed'
	
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	Z = zmin+np.array(range(0,znum))*zdel
	R = xmin+np.array(range(0,xnum))*xdel
	
	energy = emin*efactor**np.array(range(0,enum))
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	scidata = scidata[:,:,0] #sum over z
	scidata = np.sum(scidata.T*energy**2./beta, axis=1)
	luminosity = scidata*log_efactor*math.log(10.)*4.*math.pi/c*1.e6*A
	
	return Z, luminosity

def nuclei_luminosity(emissfilename):
	##
	## Takes galprop full nuclei file and calculates total energy of all species
	##
	## Resturns total energy array (A, Z, erg)
	##
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	from scipy.interpolate import griddata
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open('nuclei_full_54_'+emissfilename+'.gz')
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
			
	for i in range(0,nucnum):
		beta = np.sqrt(1.-1./(1.+energy/mass[i])**2.)
		dat = np.sum(scidata[i]*zdel, axis=1) #sum over z
		dat = np.sum(dat*r_int, axis=1) #sum over r
		dat = np.sum(dat/beta*log_efactor)*math.log(10.) #sum over e
		if AA[i] == 0: 
			A = 1
		else:
			A = AA[i]
		luminosity[i] = dat*4.*math.pi/c*1.e6*ev2ergs*kpc2cm**3.*A
	return AA, ZZ, luminosity

def electron_spectrum(nuclei_file, r, z):
	##
	## Calculates electron spectrum at position r,z
	##
	## Returns Energy, E_kin*U_kin (species,energy)
	##
	
	import pyfits as py
	import numpy as np
	import math
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(nuclei_file)
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	ir = np.argmin(np.abs(R-r))
	iz = np.argmin(np.abs(Z-z))
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	electron_species = 0
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
			electron_species+= 1
		else:
			mass[i] = AA[i]*938.28
	
	#== Picking electrons
	electron_data = scidata[AA==0]
	
	data = electron_data
	#>> data[species][energy][z][r]
	
	AA = 1 #Set to 1 to use formulae
	mass = 0.5109989461 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Picking position
	data = data[:,:,iz,ir]
	#>> data[species][energy]
	
	#factor
	data*= 4.*math.pi/c/beta
	data*= (1+mass/energy)**2
	
	energy+= mass
	
	return energy, data

def electron_spectrum_kinetic(nuclei_file, r, z):
	##
	## Calculates electron spectrum at position r,z
	##
	## Returns Energy, SE_kin*U_kin (species,energy)
	##
	
	import pyfits as py
	import numpy as np
	import math
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(nuclei_file)
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	ir = np.argmin(np.abs(R-r))
	iz = np.argmin(np.abs(Z-z))
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	electron_species = 0
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
			electron_species+= 1
		else:
			mass[i] = AA[i]*938.28
	
	#== Picking electrons
	electron_data = scidata[AA==0]
	
	data = electron_data
	#>> data[species][energy][z][r]
	
	AA = 1 #Set to 1 to use formulae
	mass = 0.5109989461 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Picking position
	data = data[:,:,iz,ir]
	#>> data[species][energy]
	
	#factor
	data*= 4.*math.pi/c/beta
	
	return energy, data
	
def electron_spectrum_total(nuclei_file):
	##
	## Calculates total electron spectrum
	##
	## Returns Energy, E_kin*U_kin (species,energy)
	##
	
	import pyfits as py
	import numpy as np
	import math
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(nuclei_file)
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	electron_species = 0
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
			electron_species+= 1
		else:
			mass[i] = AA[i]*938.28
	
	#== Picking electrons
	electron_data = scidata[AA==0]
	
	AA = 1 #Set to 1 to use formulae
	mass = 0.5109989461 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	data = electron_data
	#>> data[species][energy][z][r]
	
	#integrate in r
	data = np.sum(data*r_int, axis=3)*xdel*kpc2cm
	#>> data[species][energy][z]
	
	#integrate in z
	data = np.sum(data, axis=2)*zdel*kpc2cm
	#>> data[species][energy]
	
	#factor
	data*= 4.*math.pi/c/beta
	
	return energy, data
	
def proton_spectrum(nuclei_file, r, z):
	##
	## Calculates electron spectrum at position r,z
	##
	## Returns Energy, E_kin*U_kin (species,energy)
	##
	
	import pyfits as py
	import numpy as np
	import math
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(nuclei_file)
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	ir = np.argmin(np.abs(R-r))
	iz = np.argmin(np.abs(Z-z))
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	species = 0
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
			species+= 1
	
	#== Picking protons
	data = scidata[AA==1]
	#>> data[species][energy][z][r]
	
	AA = 1 #Set to 1 to use formulae
	mass = 938.28 #0.5109989461 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Picking position
	data = data[:,:,iz,ir]
	#>> data[species][energy]
	
	#factor
	data*= 4.*math.pi/c/beta
	data*= (1+mass/energy)**2
	
	energy+= mass
	
	return energy, data

def proton_spectrum_kinetic(nuclei_file, r, z):
	##
	## Calculates electron spectrum at position r,z
	##
	## Returns Energy, E_kin*U_kin (species,energy)
	##
	
	import pyfits as py
	import numpy as np
	import math
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(nuclei_file)
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	ir = np.argmin(np.abs(R-r))
	iz = np.argmin(np.abs(Z-z))
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	species = 0
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
			species+= 1
	
	#== Picking protons
	data = scidata[AA==1]
	#>> data[species][energy][z][r]
	
	AA = 1 #Set to 1 to use formulae
	mass = 938.28 #0.5109989461 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Picking position
	data = data[:,:,iz,ir]
	#>> data[species][energy]
	
	#factor
	data*= 4.*math.pi/c/beta
	
	return energy, data
	
def proton_spectrum_total(nuclei_file):
	##
	## Calculates total electron spectrum
	##
	## Returns Energy, E_kin*U_kin (species,energy)
	##
	
	import pyfits as py
	import numpy as np
	import math
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	c = 2.99792458e10 #cm/s
	ev2ergs = 1.60218e-12
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(nuclei_file)
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	
	AA = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	species = 0
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
		else:
			mass[i] = AA[i]*938.28
			species+= 1
	
	#== Picking protons
	data = scidata[AA==1]
	#>> data[species][energy][z][r]
	
	AA = 1 #Set to 1 to use formulae
	mass = 938.28 #0.5109989461 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#integrate in r
	data = np.sum(data*r_int, axis=3)*xdel*kpc2cm
	#>> data[species][energy][z]
	
	#integrate in z
	data = np.sum(data, axis=2)*zdel*kpc2cm
	#>> data[species][energy]
	
	#factor
	data*= 4.*math.pi/c/beta
	
	return energy, data
	
def electron_spectrum_z(emissfilename, RR):
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
	hdulist=py.open('nuclei_full_54_'+emissfilename+'.gz')
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
	
	dat = scidata[AA==0]

	#Picking electrons
	AA = 1 #AA[AA==0]
	mass = 0.5109989461 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Adding e+- together
	dat = np.sum(dat, axis=0)
	
	#picking out radius
	dat = np.swapaxes(dat, 0, 2)/beta*AA/energy**2
	dat = dat[(R-RR) <= 2*xdel]
	dat = np.average(dat, axis=0)
	
	dat *= 4.*math.pi/c#*1.e6*ev2ergs
	dat = dat.T
	
	# np.savetxt(emissfilename+'_electron_spectrum.dat', dat)
	# np.savetxt(emissfilename+'_electron_spectrum_z.dat', Z)
	# np.savetxt(emissfilename+'_electron_spectrum_e.dat', energy)
	
	##############################
	#SPECTRAL INDEX
	##############################
	energy_i = np.zeros(len(energy)-1)
	data = np.zeros([len(energy)-1, len(Z)])
	
	for n in range(0,len(energy)-1):
		energy_i[n] = math.sqrt(energy[n]*energy[n+1])
		loge2e1 = math.log(energy[n+1]/energy[n])
		for j in range(0,len(Z)):
			data[n][j]=math.log(dat[n+1][j]/dat[n][j]+1e-40)/loge2e1
	
	# np.savetxt(emissfilename+'_electron_index.dat', data)
	# np.savetxt(emissfilename+'_electron_index_z.dat', Z)
	# np.savetxt(emissfilename+'_electron_index_e.dat', energy_i)
	
	
	#############################
	#PLOTS
	#############################
	fig = plt.figure()
	
	f1, ax = plt.subplots(1,2)
	#ax = fig.add_subplot(1,1,1)# plt.subplot(1,1)
	
	X, Y = np.meshgrid(energy,Z)
	V = np.logspace(-25, -8, 20)
	
	dat = dat.T
	CF = ax[0].contourf(X,Y,dat,norm = LogNorm(), levels=V)
	ax[0].set_xlim([1,1e8])
	ax[0].set_xscale('log')
	ax[0].set_xlabel(r'$E$ (MeV)')
	ax[0].set_ylabel(r'$z$ (kpc)')
	ax[0].set_title(r'$n(E,z)$')
	
	#cbaxes = fig.add_axes([-0.2, 0, 0.5, 0.1]) 
	plt.colorbar(CF,ax=ax[0],ticks=V,pad=0)#,cax = cbaxes)
	
	data = data.T
	X1, Y1 = np.meshgrid(energy_i,Z)
	V1 = np.linspace(-4, 1, 21)
	CF1 = ax[1].contourf(X1,Y1,data, levels=V1)
	ax[1].set_xlim([1,1e8])
	ax[1].set_xscale('log')
	ax[1].set_xlabel(r'$E$ (MeV)')
	ax[1].set_title(r'$p(E,z)$')
	plt.colorbar(CF1,ticks=V1,ax=ax[1],pad=0)
	
	plt.suptitle('Electrons ('+emissfilename+')')
	plt.savefig(emissfilename+'_electron_z.png')

def proton_spectrum_z(emissfilename, RR):
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
	hdulist=py.open('nuclei_full_54_'+emissfilename+'.gz')
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
	
	dat = scidata[abs(AA)==1]
	# print AA
	# print AA==0
	# print scidata.shape
	# print dat.shape
	#Picking electrons
	AA = 1 #AA[AA==0]
	mass = 938.28 #MeV/c^2 mass[AA==0]
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	# print beta.shape
	# print dat.shape
	#Adding e+- together
	dat = np.sum(dat, axis=0)
	# print dat.shape
	#dat = np.swapaxes(dat, 0,2)/beta*AA
	
	#picking out radius
	dat = np.swapaxes(dat, 0, 2)/beta*AA/energy**2
	dat = dat[(R-RR) <= 2*xdel]
	dat = np.average(dat, axis=0)
	
	dat *= 4.*math.pi/c#*1.e6*ev2ergs
	dat = dat.T
	
	# np.savetxt(emissfilename+'_proton_spectrum.dat', dat)
	# np.savetxt(emissfilename+'_proton_spectrum_z.dat', Z)
	# np.savetxt(emissfilename+'_proton_spectrum_e.dat', energy)
	
	##############################
	#SPECTRAL INDEX
	##############################
	energy_i = np.zeros(len(energy)-1)
	data = np.zeros([len(energy)-1, len(Z)])
	
	for n in range(0,len(energy)-1):
		energy_i[n] = math.sqrt(energy[n]*energy[n+1])
		loge2e1 = math.log(energy[n+1]/energy[n])
		for j in range(0,len(Z)):
			data[n][j]=math.log(dat[n+1][j]/dat[n][j]+1e-40)/loge2e1
	
	# np.savetxt(emissfilename+'_proton_index.dat', data)
	# np.savetxt(emissfilename+'_proton_index_z.dat', Z)
	# np.savetxt(emissfilename+'_proton_index_e.dat', energy_i)
	
	
	#############################
	#PLOTS
	#############################
	fig = plt.figure()
	
	f1, ax = plt.subplots(1,2)
	#ax = fig.add_subplot(1,1,1)# plt.subplot(1,1)
	
	X, Y = np.meshgrid(energy,Z)
	V = np.logspace(-25, -8, 20)
	dat = dat.T
	CF = ax[0].contourf(X,Y,dat,norm = LogNorm(), levels=V)
	ax[0].set_xlim([1,1e8])
	ax[0].set_xscale('log')
	ax[0].set_xlabel(r'$E$ (MeV)')
	ax[0].set_ylabel(r'$z$ (kpc)')
	ax[0].set_title(r'$n(E,z)$')
	
	#cbaxes = fig.add_axes([-0.2, 0, 0.5, 0.1]) 
	plt.colorbar(CF,ax=ax[0],ticks=V,pad=0)#,cax = cbaxes)
	
	data = data.T
	X1, Y1 = np.meshgrid(energy_i,Z)
	V1 = np.linspace(-4, 1, 21)
	CF1 = ax[1].contourf(X1,Y1,data, levels=V1)
	ax[1].set_xlim([1,1e8])
	ax[1].set_xscale('log')
	ax[1].set_xlabel(r'$E$ (MeV)')
	ax[1].set_title(r'$p(E,z)$')
	plt.colorbar(CF1,ticks=V1,ax=ax[1],pad=0)
	
	plt.suptitle('Protons ('+emissfilename+')')
	plt.savefig(emissfilename+'_proton_z.png')
	
def radiation_emissivity(emissfilename, r, z):
	##
	## Calculates emissivity from galprop output at position r,z
	##
	## Returns freq/energy, emissivity (energy)
	## Check galprop file for actual units
	##
	
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
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	ir = np.argmin(np.abs(R-r))
	iz = np.argmin(np.abs(Z-z))
	
	emissivity = scidata[0,:,iz,ir]
	return energy, emissivity

def radiation_emissivity_ic(emissfilename, r, z):
	##
	## Calculates emissivity from galprop IC output at position r,z
	##
	## Returns energy, emissivity (energy)
	## Check galprop file for actual units
	##
	
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
	nucnum = head['NAXIS4']
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*(xnum-1)
	zmax = zmin+zdel*(znum-1)
	
	xsize = xmax-xmin
	zsize = zmax-zmin

	R = xmin+np.array(range(0,xnum))*xdel
	Z = zmin+np.array(range(0,znum))*zdel
	energy = emin*efactor**np.array(range(0,enum))
	
	ir = np.argmin(np.abs(R-r))
	iz = np.argmin(np.abs(Z-z))
	
	emissivity = scidata[0,:,iz,ir]
	return energy, emissivity
	
def synch_spectrum_z(gp_name,XX):
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
	
	DD = 3700 #kpc
	
	#GETTING DATA
	a=py.open(gp_name+'synch_spec.fits')
	data=a[0].data
	head=a[0].header
	x0=head['CRVAL1']
	dx=head['CDELT1']
	numx=head['NAXIS1']
	nu0=head['CRVAL3']
	dnu=head['CDELT3']
	numnu=head['NAXIS3']

	xxp=np.linspace(0,numx-1,numx)
	xxp*=dx
	xxp+=x0

	nup=np.linspace(0,numnu-1,numnu)
	nup=dnu**nup
	nup*=nu0

	zmin = head['CRVAL2']
	zdel = head['CDELT2']
	znum = head['NAXIS2']
	Z = np.linspace(0, znum-1, znum)
	Z *= zdel
	Z += zmin
	Z = Z/radian2arcsec*DD
	
	data=np.swapaxes(data,2,1)
	data=np.swapaxes(data,0,1)

	#Picking out near XX
	PP=data[abs(xxp-XX)<3*dx]
	PP=np.average(PP,axis=0)
	PP = -1*(1-2*(PP+1))

	a=py.open(gp_name+'synch.fits')
	data=a[0].data
	head=a[0].header
	x0=head['CRVAL1']
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
	FF=data[abs(xxf-XX)<3*dx]
	FF=np.average(FF,axis=0)

	#########################
	## Plotting
	#########################
	fig = plt.figure()
	
	f1, ax = plt.subplots(1,2)
	
	X, Y = np.meshgrid(nuf,Z)
	V = np.logspace(-30, -10, 20)
	dat = FF.T
	CF = ax[0].contourf(X,Y,dat,norm = LogNorm(), levels=V)
	ax[0].set_xlim([1e6,1e17])
	ax[0].set_xscale('log')
	ax[0].set_xlabel(r'$\nu$ (Hz)')
	ax[0].set_ylim([-1,1])
	ax[0].set_ylabel(r'$z$ (kpc)')
	ax[0].set_title(r'$F(E,z)$')
	 
	plt.colorbar(CF,ax=ax[0],ticks=V,pad=0)
	
	data = PP.T
	X1, Y1 = np.meshgrid(nup,Z)
	V1 = np.linspace(-4, 1, 21)
	CF1 = ax[1].contourf(X1,Y1,data, levels=V1)
	ax[1].set_xlim([1e6,1e17])
	ax[1].set_xscale('log')
	ax[1].set_xlabel(r'$\nu$ (Hz)')
	ax[1].set_ylim([-1,1])
	ax[1].set_title(r'$p(E,z) = 1-2*\alpha$')
	plt.colorbar(CF1,ticks=V1,ax=ax[1],pad=0)
	
	plt.suptitle('Synchroton ('+gp_name+')')
	plt.savefig(gp_name+'_synchrotron_z.png')

def gamma_spectrum_z(gp_name,XX):
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
	
	DD = 3700 #kpc
	
	#GETTING DATA
	a=py.open(gp_name+'gamma_total_spec.fits')
	data=a[0].data
	head=a[0].header
	x0=head['CRVAL1']
	dx=head['CDELT1']
	numx=head['NAXIS1']
	nu0=head['CRVAL3']
	dnu=head['CDELT3']
	numnu=head['NAXIS3']

	xxp=np.linspace(0,numx-1,numx)
	xxp*=dx
	xxp+=x0

	nup=np.linspace(0,numnu-1,numnu)
	nup=dnu**nup
	nup*=nu0

	zmin = head['CRVAL2']
	zdel = head['CDELT2']
	znum = head['NAXIS2']
	Z = np.linspace(0, znum-1, znum)
	Z *= zdel
	Z += zmin
	Z = Z/radian2arcsec*DD
	
	data=np.swapaxes(data,2,1)
	data=np.swapaxes(data,0,1)

	#Picking out near XX
	PP=data[abs(xxp-XX)<3*dx]
	PP=np.average(PP,axis=0)
	PP = -1*(1-2*(PP+1))

	a=py.open(gp_name+'gamma_total.fits')
	data=a[0].data
	head=a[0].header
	x0=head['CRVAL1']
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
	FF=data[abs(xxf-XX)<3*dx]
	FF=np.average(FF,axis=0)

	#########################
	## Plotting
	#########################
	fig = plt.figure()
	
	f1, ax = plt.subplots(1,2)
	
	X, Y = np.meshgrid(nuf,Z)
	V = np.logspace(-6, 4, 20)
	dat = FF.T
	CF = ax[0].contourf(X,Y,dat,norm = LogNorm(), levels=V)
	ax[0].set_xlim([1e1,1e7])
	ax[0].set_xscale('log')
	ax[0].set_xlabel(r'$E$ (MeV)')
	ax[0].set_ylim([-1,1])
	ax[0].set_ylabel(r'$z$ (kpc)')
	ax[0].set_title(r'$F(E,z)$')
	 
	plt.colorbar(CF,ax=ax[0],ticks=V,pad=0)
	
	data = PP.T
	X1, Y1 = np.meshgrid(nup,Z)
	V1 = np.linspace(-4, 1, 21)
	CF1 = ax[1].contourf(X1,Y1,data, levels=V1)
	ax[1].set_xlim([1e1,1e7])
	ax[1].set_xscale('log')
	ax[1].set_xlabel(r'$E$ (MeV)')
	ax[1].set_ylim([-1,1])
	ax[1].set_title(r'$p(E,z) = 1-2\alpha$')
	plt.colorbar(CF1,ticks=V1,ax=ax[1],pad=0)
	
	plt.suptitle('Total Gamma ('+gp_name+')')
	plt.savefig(gp_name+'_gamma_z.png')

def electron_energy_density_z(emissfilename,RR):
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
	hdulist=py.open('nuclei_full_54_'+emissfilename+'.gz')
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
	dat = scidata[AA==0]
	
	AA = 1
	mass = 0.5109989461 #MeV/c^2
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Adding e+- together
	dat = np.sum(dat, axis=0)
	
	#picking out radius
	dat = np.swapaxes(dat, 0, 2)/beta
	dat = dat[(R-RR) <= 2*xdel]
	dat = np.average(dat, axis=0)
	
	dat *= 4.*math.pi/c
	#dat = dat.T
	
	#integrate over E
	U_cr = np.sum(dat, axis=1)*math.log(energy[1]/energy[0])*1e6*ev2ergs
	
	#############################
	# gradient of U_cr
	#############################
	Z_grad = (Z[1:]+Z[:-1])/2
	grad_U = (U_cr[1:]-U_cr[:-1])/(zdel*kpc2cm)
	
	
	#############################
	#PLOTS
	#############################
	fig = plt.figure()
	
	f1, ax = plt.subplots(2,1)
	#ax = fig.add_subplot(1,1,1)# plt.subplot(1,1)
	
	ax[0].plot(Z,U_cr)
	ax[0].set_xlabel(r'$z$ (kpc)')
	ax[0].set_ylabel(r'$U_{cr}$ (erg/cm$^3$)')
	#ax[0].set_title(r'Electrons $U_{cr}(z)$')
	
	ax[1].plot(Z_grad,grad_U)
	ax[1].set_xlabel(r'$z$ (kpc)')
	ax[1].set_ylabel(r'$\frac{d}{dz} U_{cr}(z)$')
	
	plt.suptitle('Electrons ('+emissfilename+')')
	plt.savefig(emissfilename+'_electron_energy_density_z.png')
	
	return Z_grad, grad_U

def proton_energy_density_z(emissfilename,RR):
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
	hdulist=py.open('nuclei_full_54_'+emissfilename+'.gz')
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
	dat = scidata[abs(AA)==1]
	
	AA = 1
	mass = 938.28 #MeV/c^2
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Adding e+- together
	dat = np.sum(dat, axis=0)
	
	#picking out radius
	dat = np.swapaxes(dat, 0, 2)/beta*AA
	dat = dat[(R-RR) <= 2*xdel]
	dat = np.average(dat, axis=0)
	
	dat *= 4.*math.pi/c
	#dat = dat.T
	
	#integrate over E
	U_cr = np.sum(dat, axis=1)*math.log(energy[1]/energy[0])*1e6*ev2ergs
	
	#############################
	# gradient of U_cr
	#############################
	Z_grad = (Z[1:]+Z[:-1])/2
	grad_U = (U_cr[1:]-U_cr[:-1])/(zdel*kpc2cm)
	
	
	#############################
	#PLOTS
	#############################
	fig = plt.figure()
	
	f1, ax = plt.subplots(2,1)
	#ax = fig.add_subplot(1,1,1)# plt.subplot(1,1)
	
	ax[0].plot(Z,U_cr)
	ax[0].set_xlabel(r'$z$ (kpc)')
	ax[0].set_ylabel(r'$U_{cr}$ (erg/cm$^3$)')
	#ax[0].set_title(r'Electrons $U_{cr}(z)$')
	
	ax[1].plot(Z_grad,grad_U)
	ax[1].set_xlabel(r'$z$ (kpc)')
	ax[1].set_ylabel(r'$\frac{d}{dz} U_{cr}(z)$')
	
	plt.suptitle('Protons ('+emissfilename+')')
	plt.savefig(emissfilename+'_proton_energy_density_z.png')
	
	return Z_grad, grad_U

def proton_pressure_z_new(emissfilename,RR):
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
	hdulist=py.open('nuclei_full_54_'+emissfilename+'.gz')
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
	dat = scidata[abs(AA)==1]
	
	AA = 1
	mass = 938.28 #MeV/c^2
	beta = np.sqrt(1.-1./(1.+energy/mass)**2.)
	
	#Adding e+- together
	dat = np.sum(dat, axis=0)
	
	#picking out radius
	dat = np.swapaxes(dat, 0, 2)/beta*AA
	dat = dat[(R-RR) <= 2*xdel]
	dat = np.average(dat, axis=0)
	
	dat *= 4.*math.pi/c
	#dat = dat.T
	
	#integrate over E #NOTE: KINETIC ENERGY DENSITY
	U_cr = np.sum(dat, axis=1)*math.log(energy[1]/energy[0])*1e6*ev2ergs
	
	#############################
	# gradient of U_cr
	#############################
	Z_grad = (Z[1:]+Z[:-1])/2
	grad_U = (U_cr[1:]-U_cr[:-1])/(zdel*kpc2cm)
	
	
	#############################
	#PLOTS
	#############################
	fig = plt.figure()
	
	f1, ax = plt.subplots(2,1)
	#ax = fig.add_subplot(1,1,1)# plt.subplot(1,1)
	
	ax[0].plot(Z,U_cr)
	ax[0].set_xlabel(r'$z$ (kpc)')
	ax[0].set_ylabel(r'$U_{cr}$ (erg/cm$^3$)')
	#ax[0].set_title(r'Electrons $U_{cr}(z)$')
	
	ax[1].plot(Z_grad,grad_U)
	ax[1].set_xlabel(r'$z$ (kpc)')
	ax[1].set_ylabel(r'$\frac{d}{dz} U_{cr}(z)$')
	
	plt.suptitle('Protons ('+emissfilename+')')
	plt.savefig(emissfilename+'_proton_energy_density_z_new.png')
	
	return Z_grad, grad_U

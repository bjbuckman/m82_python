### =======================================================
### Finds integrated flux from any emissfile
###	
def integrated_flux(emissfilename, fileout, objectdist):
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
	
	objectdistCM = objectdist*kpc2cm
	
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	
	if 'IC_emiss' in emissfilename:
		scidata=hdulist[0].data[0][3]
	else:
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
	
	R = xmin+np.array(range(0,xnum))*xdel
	energy = emin*efactor**np.array(range(0,enum))
	
	#integrate over z
	scidata = np.sum(scidata, axis=1)*zdel
	
	#integrate over r
	r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2.
	scidata = np.sum(scidata*r_int, axis=1)
	
	nan_array = np.isnan(scidata)
	#scidata = scidata[~nan_array]
	#energy = energy[~nan_array]
	
	#unit conversion and distance projection
	if 'synchrotron' in emissfilename or 'free_free_emiss' in emissfilename:
		total_emission = np.sum(energy[~nan_array]*scidata[~nan_array])*4*math.pi*kpc2cm**3*math.log(efactor)
	else:
		total_emission = np.sum(scidata[~nan_array])*4*math.pi*kpc2cm**3*math.log(efactor)*erg_per_MeV
	
	scidata *= kpc2cm**3/objectdistCM**2
	
	
	np.savetxt(fileout+'.dat', [energy, scidata], '%e')
	
	return total_emission
	
### =======================================================
### From ntegrated flux files, sum all
###	
def sum_flux3(file1, file2, file3, fileout):
	'''
	Sums 3 files together
	'''
	import numpy as np
	a1 = np.loadtxt(file1)
	a2 = np.loadtxt(file2)
	a3 = np.loadtxt(file3)
	
	energy = a1[0]
	total_flux = a1[1]+a2[1]+a3[1]
	
	np.savetxt(fileout+'.dat', [energy, total_flux], '%e')
	
### =======================================================
### From ntegrated flux files, sum all
###	
def sum_flux2(file1, file2, fileout):
	'''
	Sums 2 files together
	'''
	import numpy as np
	a1 = np.loadtxt(file1)
	a2 = np.loadtxt(file2)
	
	energy = a1[0]
	total_flux = a1[1]+a2[1]
	
	np.savetxt(fileout+'.dat', [energy, total_flux], '%e')	
	
	
	
	
	
	

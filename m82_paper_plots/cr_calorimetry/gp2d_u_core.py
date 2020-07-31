import pyfits as py
import numpy as np
import math
import os
import sys
from sys import stdout

# galdef_name = sys.argv[1]
distance = 3520 #kpc

#constants
pc2cm = 3.08568025e18 #cm/pc
kpc2cm = 3.08568025e21 #cm/kpc
radian2arcsec = 648000./math.pi #acsec/radian
Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
c = 2.99792458e10 #cm/s
ev2ergs = 1.60218e-12

def nuclei_energy_core(emissfilename):
	##
	## Takes galprop full nuclei file and calculates total energy inside core of all species
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
	energy = emin*efactor**np.array(range(0,enum))
	
	r_int = math.pi*((R+xdel/2.)**2-(R-xdel/2.)**2)
	# r_int = math.pi*2.*R*xdel
	r_int[0] = math.pi*(xdel/2.)**2
	
	Z = zmin+np.array(range(0,znum))*zdel
	
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
		dat = scidata[i]
		dat = np.sum(dat[:,:,R<=0.2]*r_int[R<=0.2], axis=2) #sum over r
		dat = np.sum(dat[:,Z<=0.05]*zdel, axis=1) #sum over z
		dat = np.sum(dat/beta*(1.+mass[i]/energy))*math.log(energy[1]/energy[0]) #sum over e
		if AA[i] == 0: 
			A = 1
		else:
			A = AA[i]
		luminosity[i] = dat*4.*math.pi/c*1.e6*A/(math.pi*(0.2+xdel/2.)**2*(0.1+zdel))
	return AA, ZZ, luminosity
	
nuclei_dir = '../fits/raw/'

MODE = 11

if MODE == 0:
	GG = 'm82_05aff1018'
	factor = 2.18287458838
elif MODE == 1:
	GG =  'm82_05b2ff1019'
	factor = 1.55132821598
elif MODE == 2:
	GG =  'm82_06bff1022'
	factor = 1.47745544379
elif MODE == 10:
	GG = 'm82_05awff1019'
	factor = 2.2920183178
	factor = 2.18287458838
elif MODE == 11:
	GG =  'm82_05b2wff1020'
	factor = 1.62889462678
	factor = 1.55132821598
elif MODE == 12:
	GG =  'm82_06bff1022'
	factor = 1.47745544379



nuclei_file = nuclei_dir+GG+'/nuclei_full_54_'+GG+'.gz'

result = nuclei_energy_core(nuclei_file)
print GG
print result[0]
print result[1]
print result[2]*factor
print np.sum(result[2]*factor)

# data = np.loadtxt(galdef_name+'analytics.dat')
# print data

# np.savetxt(galdef_name+'analytics.dat',DATA,'%e')
# np.savetxt(galdef_name+'analytics_nuc.dat',DATA2,'%e')


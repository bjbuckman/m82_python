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

def particle_escape(source_file, edge):
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(source_file)
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
	A  = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
			A[i] = 1.
		else:
			mass[i] = AA[i]*938.28
			A[i] = AA[i]
	
	#get data of edge
	# [x_max, all z]
	scidata_edge1 = np.swapaxes(scidata,0,3)
	scidata_edge1 = scidata_edge1[-(1+edge)]
	#int energy
	scidata_edge1 = np.sum(scidata_edge1, axis=0)*math.log(efactor)
	scidata_edge1 = scidata_edge1[edge:-(1+edge)]*zdel*2.*math.pi*R[-(1+edge)]
	scidata_edge1 = np.sum(scidata_edge1, axis=0)*math.pi*1.e6*ev2ergs*kpc2cm**2.*A
	
	# print scidata_edge1
	
	# [all r, min and max z]
	scidata_edge2 = np.swapaxes(scidata,0,2)
	scidata_edge2 = scidata_edge2[-(1+edge)]+scidata_edge2[edge]
	scidata_edge2 = np.sum(scidata_edge2, axis=0)*math.log(efactor)
	scidata_edge2*= r_int
	scidata_edge2 = np.sum(scidata_edge2[:,:-(1+edge)], axis=1)*A
	scidata_edge2*= math.pi*1.e6*ev2ergs*kpc2cm**2.
	
	# print scidata_edge2
	
	return scidata_edge1 + scidata_edge2

def particle_escape_d(source_file, edge):
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(source_file)
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
	
	Z = zmin+np.array(range(0,znum))*zdel
	
	AA = np.zeros(nucnum)
	A  = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
			A[i] = 1.
		else:
			mass[i] = AA[i]*938.28
			A[i] = AA[i]
	
	#get data of edge
	# [x_max, all z]
	scidata_edge1 = np.swapaxes(scidata,0,3)
	scidata_edge1 = scidata_edge1[int(edge/xdel)]
	#int energy
	scidata_edge1 = np.sum(scidata_edge1, axis=0)*math.log(efactor)
	scidata_edge1 = scidata_edge1[abs(Z)<=edge]*zdel*2.*math.pi*R[int(edge/xdel)]
	scidata_edge1 = np.sum(scidata_edge1, axis=0)*math.pi*1.e6*ev2ergs*kpc2cm**2.*A
	
	# print scidata_edge1
	
	# [all r, min and max z]
	scidata_edge2 = np.swapaxes(scidata,0,2)
	scidata_edge2 = scidata_edge2[int((edge-zmin)/zdel)]+scidata_edge2[int((-edge-zmin)/zdel+1)]
	scidata_edge2 = np.sum(scidata_edge2, axis=0)*math.log(efactor)
	scidata_edge2*= r_int
	scidata_edge2 = np.sum(scidata_edge2[:,R<=edge], axis=1)*A
	scidata_edge2*= math.pi*1.e6*ev2ergs*kpc2cm**2.
	
	# print scidata_edge2
	
	return scidata_edge1 + scidata_edge2

def particle_escape_flux(source_file, edge):
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(source_file)
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
	
	Z = zmin+np.array(range(0,znum))*zdel
	
	AA = np.zeros(nucnum)
	A  = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
			A[i] = 1.
		else:
			mass[i] = AA[i]*938.28
			A[i] = AA[i]
	
	#get data of edge
	# [x_max, all z]
	scidata_edge1p = np.swapaxes(scidata,0,3)
	scidata_edge1 = scidata_edge1p[int(edge/xdel)]
	scidata_edge1-= scidata_edge1p[int(edge/xdel-1)]
	#int energy
	scidata_edge1 = np.sum(scidata_edge1, axis=0)*math.log(efactor)
	scidata_edge1 = scidata_edge1[abs(Z)<=edge]*zdel*2.*math.pi*R[int(edge/xdel)]
	scidata_edge1 = np.sum(scidata_edge1, axis=0)*math.pi*1.e6*ev2ergs*kpc2cm**2.*A
	
	# print scidata_edge1
	
	# [all r, min and max z]
	scidata_edge2p = np.swapaxes(scidata,0,2)
	scidata_edge2 = scidata_edge2p[int((edge-zmin)/zdel)]+scidata_edge2p[int((-edge-zmin)/zdel+1)]
	scidata_edge2-= scidata_edge2p[int((edge-zmin)/zdel-1)]+scidata_edge2p[int((-edge-zmin)/zdel+2)]
	scidata_edge2 = np.sum(scidata_edge2, axis=0)*math.log(efactor)
	scidata_edge2*= r_int
	scidata_edge2 = np.sum(scidata_edge2[:,R<=edge], axis=1)*A
	scidata_edge2*= math.pi*1.e6*ev2ergs*kpc2cm**2.
	
	# print scidata_edge2
	
	return scidata_edge1 + scidata_edge2

def particle_escape_flux_v2(source_file, edge):
	#READING EMISSION FITS-------------------------------------------------
	hdulist=py.open(source_file)
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
	
	Z = zmin+np.array(range(0,znum))*zdel
	
	AA = np.zeros(nucnum)
	A  = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	
	for i in range(0,nucnum):
		AA[i] = head['NUCA'+str(i+1).zfill(3)]
		ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
		if AA[i] == 0:
			mass[i] = 0.5109989461 #MeV/c^2
			A[i] = 1.
		else:
			mass[i] = AA[i]*938.28
			A[i] = AA[i]
	
	#get data of edge
	# [x_max, all z]
	scidata_edge1 = np.sum(scidata, axis=1)*math.log(efactor)
	scidata_edge1p = np.swapaxes(scidata_edge1,0,2)
	scidata_edge1u = scidata_edge1p[int(edge/xdel+1)]
	# scidata_edge1  = scidata_edge1p[int(edge/xdel)]
	scidata_edge1d = scidata_edge1p[int(edge/xdel-1)]
	scidata_edge1 = scidata_edge1u-scidata_edge1d
	
	scidata_edge1 = scidata_edge1[abs(Z)<=edge]*zdel*2.*math.pi*R[int(edge/xdel)]
	scidata_edge1 = np.sum(scidata_edge1, axis=0)*math.pi*1.e6*ev2ergs*kpc2cm**2.*A
	
	# print scidata_edge1
	
	# [all r, min and max z]
	scidata_edge2 = np.sum(scidata, axis=1)*math.log(efactor)
	scidata_edge2p = np.swapaxes(scidata_edge2,0,1)
	scidata_edge2u = scidata_edge2p[int((edge-zmin)/zdel+1.)]+scidata_edge2p[int((-edge-zmin)/zdel+0)]
	# scidata_edge2  = scidata_edge2p[int((edge-zmin)/zdel)]+scidata_edge2p[int((-edge-zmin)/zdel+1)]
	scidata_edge2d = scidata_edge2p[int((edge-zmin)/zdel-1)]+scidata_edge2p[int((-edge-zmin)/zdel+2)]
	scidata_edge2 = scidata_edge2u-scidata_edge2d
	
	scidata_edge2*= r_int
	scidata_edge2 = np.sum(scidata_edge2[:,R<=edge], axis=1)*A
	scidata_edge2*= math.pi*1.e6*ev2ergs*kpc2cm**2.
	
	# print scidata_edge2
	
	return scidata_edge1 + scidata_edge2
	
fits_dir = '../fits/raw/'

MODE = 2

if MODE == 0:
	GG = 'm82_05aff1018'
	factor = 2.18287458838
elif MODE == 1:
	GG =  'm82_05b2ff1019'
	factor = 1.55132821598
elif MODE == 2:
	GG =  'm82_06bff1022'
	factor = 1.47745544379

galdef_name = fits_dir+GG+'/'+GG
nuclei_file = fits_dir+GG+'/nuclei_full_54_'+GG+'.gz'
data = np.loadtxt(galdef_name+'analytics.dat')
print data

for D in abs(np.arange(-3.9,-0.1,0.1)):
	result = particle_escape_flux(nuclei_file,D) #- particle_escape_d(galdef_name, i+1)
	print GG
	print 'd = '+str(D)
	print 'He      '+str(result[9]/data[0])
	print 'prim p+ '+str(result[6]/data[1])
	print 'scdy p+ '+str(result[5]/data[12])
	print 'prim e- '+str(result[3]/data[2])
	print 'scdy e- '+str(result[2]/data[4])
	print 'scdy e+ '+str(result[0]/data[5])
# print result




# np.savetxt(galdef_name+'analytics.dat',DATA,'%e')
# np.savetxt(galdef_name+'analytics_nuc.dat',DATA2,'%e')


#import multiprocessing
import numpy as np
import math
##
#  NO LONGER INTERPOLATING for V3, For interpolation use V2
##

### =======================================================
### Gives coordinates to interpolate x axis
###
def gen_pixel_xindex(xmin, dx, nmin, dn, SYMMETRY):
	"""
	Gets x index for interpolation
	"""
	if SYMMETRY == 1:
		num_coords = int( round(abs(nmin/dn)) +1)
		nmin = xmin
	else:
		num_coords = int( round(abs(2*nmin/dn)) +1)
	coords = np.zeros(num_coords)
	
	for n in range(0, num_coords):
		nn = nmin+n*dn
		coords[n] = (nn-xmin)/dx
	
	return coords
	
	
### =======================================================
### Interpolates x axis
###
def gen_pixel_x(xmin, dx, nmin, dn, data, SYMMETRY):
	"""
	FORM OF DATA = data[x][z][y][energy/freq]
	Interpolates x-axis
	"""
	ix_max = data.shape[0]-1
	coords = gen_pixel_xindex(xmin, dx, nmin, dn, SYMMETRY)
	#print coords
	ic = np.around(coords)
	ic = ic.astype(int)
	#dic = coords-ic
	num_coords = ic.shape[0]
	#print ic
	#print dic
	#ix_max = abs(2*xmin/dx)-1
	
	pixel = np.zeros([num_coords, data.shape[1], data.shape[2], data.shape[3]])
	for i in range(0,num_coords):
		if ic[i]<0 or ic[i]>=ix_max:
			pass
		else:
			pixel[i]+= data[ic[i]] #+ dic[i]*(data[ic[i]+1]-data[ic[i]])
	
	print 'Finished x-interp'
	return pixel	
	
	
### =======================================================
### Returns indices for LOS integration assuming no rotation in x
###
def gen_pixel_yindex(y, nmin, dn, ymin, dy, zmin, dz, theta):
	"""
	Rotation along x-axis
	Returns interpolation indices
	"""
	coordsnum = int( round(abs(2*nmin/dn)) +1)
	#print coordsnum
	coords = np.zeros([coordsnum, 2])
	
	for n in range(0, coordsnum):
		nn = nmin+n*dn
		ydat, zdat = orig_yz(y, nn, theta)
		coords[n] = np.array([(ydat-ymin)/dy, (zdat-zmin)/dz])
		#print coords[n]
	
	return coords	

	
### =======================================================
### Generates pixel for all x_im for specified y_im
###
def gen_pixel_y(y, nmin, dn, ymin, dy, zmin, dz, theta, data):
	"""
	FORM OF DATA = data[y][z][energy/freq][x]
	LOS integral for all x of image at specified y of image
	"""
	coords = gen_pixel_yindex(y, nmin, dn, ymin, dy, zmin, dz, theta)
	#print coords
	ic = np.around(coords)
	ic = ic.astype(int)
	#dic = coords-ic
	num_coords = coords.shape[0]
	
	iy_max = data.shape[0]-1
	iz_max = data.shape[1]-1
	
	x_pixel = np.zeros([data.shape[2], data.shape[3]])
	for i in range(0,num_coords):
		if ic[i][0]<0 or ic[i][1]<0 or ic[i][0]>=iy_max or ic[i][1]>=iz_max:
			pass
		else:
			x_pixel_temp = data[ic[i][0]] #+ dic[i][0]*(data[ic[i][0]+1]-data[ic[i][0]])
			x_pixel+= x_pixel_temp[ic[i][1]] #+ dic[i][1]*(x_pixel_temp[ic[i][1]+1]-x_pixel_temp[ic[i][1]])
	
	print 'Finished y='+str(y)
	return x_pixel
	
	
### =======================================================
### Take Y and Z from image and gives original coordinates for data
###
def orig_yz(yim, zim, theta):
	"""
	Rotates about y-z axis
	"""
	ydat =  math.cos(theta)*yim + math.sin(theta)*zim
	zdat = -math.sin(theta)*yim + math.cos(theta)*zim
	return ydat, zdat

	
### =======================================================
### Converts 2D galprop emiss -> 3D galprop emiss
###
def d2_to_d3(rmin, dr, data):
	"""
	data[energy/freq][z][r] -> data[energy/freq][z][y][x]
	"""
	enum = data.shape[0]
	znum = data.shape[1]
	rnum = data.shape[2]
	data = np.swapaxes(data,0,2)
	#data[r,z,e]
	
	xnum = 2*rnum-1
	ynum = xnum
	
	xmin = -1*(rmin+(rnum-1)*dr)
	ymin = xmin
	
	output = np.zeros([enum, znum, ynum, xnum])
	output = np.swapaxes(output,0,3)
	output = np.swapaxes(output,1,2)
	#output[x,y,z,e]
	
	xx = xmin
	for ix in range(0,xnum):
		yy = ymin
		for iy in range(0,ynum):
			rr  = math.sqrt(xx**2+yy**2)
			irr = (rr-rmin)/dr
			ir  = int(math.floor(irr))
			dir = irr-ir
			#print ir
			if ir < rnum-2:
				output[ix][iy] = data[ir]+dir*(data[ir+1]-data[ir]) 
			yy+= dr
		xx+= dr
		#print xx
	
	output = np.swapaxes(output,0,3)
	output = np.swapaxes(output,1,2)
	#output[e,z,y,x]
	return output
	
	
### =======================================================
### Converts 2D galprop emiss -> 3D galprop emiss
###
def d2_to_d3sym(rmin, dr, data):
	"""
	data[energy/freq][z][r] -> data[energy/freq][z][y][x]
	"""
	enum = data.shape[0]
	znum = data.shape[1]
	rnum = data.shape[2]
	data = np.swapaxes(data,0,2)
	#data[r,z,e]
	
	xnum = rnum
	ynum = 2*rnum-1
	
	xmin = rmin
	ymin = -1*(rmin+(rnum-1)*dr)
	
	output = np.zeros([enum, znum, ynum, xnum])
	output = np.swapaxes(output,0,3)
	output = np.swapaxes(output,1,2)
	#output[x,y,z,e]
	
	xx = xmin
	for ix in range(0,xnum):
		yy = ymin
		for iy in range(0,ynum):
			rr  = math.sqrt(xx**2+yy**2)
			irr = (rr-rmin)/dr
			ir  = int(math.floor(irr))
			dir = irr-ir
			#print ir
			if ir < rnum-2:
				output[ix][iy] = data[ir]+dir*(data[ir+1]-data[ir]) 
			yy+= dr
		xx+= dr
		#print xx
	
	output = np.swapaxes(output,0,3)
	output = np.swapaxes(output,1,2)
	#output[e,z,y,x]
	return output
	
	
### =======================================================
### LOS integrator for all galprop emission files
###	
def LOS_integral(emissfilename, fileout, inclination, objectdist, outRes, int_steps, imageR):
	"""
	Line of sight integrator for all emission files
	inclines by inclination (in degrees)
	projects a distance (in kpc)
	outputs in fits file
	"""
	import pyfits as py
	import os
	import sys
	from sys import stdout
	
	#global scidata
	#global SYMMETRY
	
	SYMMETRY = 1
	
	# outRes = 21 #number of pixels = outRes**2
	# int_steps = 100 #number of integration steps
	map_prop = 1.0 #image proportion of galaxy diameter
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	
	inc_rad = inclination*math.pi/180. #inclination in radians
	
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
	
	print 'Data Read'
	
	# Convert 2D -> 3D
	if SYMMETRY == 1:
		scidata = d2_to_d3sym(xmin, xdel, scidata)
	else:
		scidata = d2_to_d3(xmin, xdel, scidata)
	print 'Finished 2d->3D'
	#print np.amax(scidata)
	
	original_mapsize = 2*xsize#math.sqrt((2.*xsize)**2.+zsize**2.) #kpc
	rotate_mapsize = map_prop*original_mapsize #kpc
	
	#OUTPUT IMAGE------------------------------------------------------	
	kpcPix = 2*imageR/(outRes-1) #kpc/pixel_dim
	
	if SYMMETRY == 1:
		image = np.zeros([enum,outRes,(outRes-1)/2+1])
		xminIm = xmin #-imageR#/2. #kpc
	else:
		image = np.zeros([enum,outRes,outRes]) #output map (pixels)
		xminIm = -imageR #/2. #kpc
	
	yminIm = -imageR#/2. #kpc
	zminIm = -rotate_mapsize/2. #kpc
	
	dxIm = kpcPix #kpc
	dyIm = kpcPix #kpc
	dzIm = rotate_mapsize/(int_steps-1) #kpc
	
	cmxPix = dxIm*kpc2cm #cm/pixel
	cmyPix = dyIm*kpc2cm #cm/pixel
	cmzPix = dzIm*kpc2cm #cm/pixel
	
	cmDist = objectdist*kpc2cm
	volPix = cmxPix*cmyPix*cmzPix #cm^3	
	srPix = dxIm*dyIm/objectdist**2. #steradian/pixel
	Pix2sac = srPix*Sr2sqArcSec #arcsec^2/pixel
	sac2Pix = 1./Pix2sac #pixel/arsec^2
	
	# scidata[energy/freq][z][y][x]
	scidata = np.swapaxes(scidata,0,3)
	if SYMMETRY == 1:
		scidata = gen_pixel_x(xmin, xdel, yminIm, dxIm, scidata, SYMMETRY)
	else:
		scidata = gen_pixel_x(-1*xmax, xdel, xminIm, dxIm, scidata, SYMMETRY)
	scidata = np.swapaxes(scidata,0,2)
	#scidata[y,z,x,e]
	
	image = np.swapaxes(image,0,2)
	image = np.swapaxes(image,0,1)
	# image[y][x][energy/freq]
	
	print 'Starting Integration'
	
	integration_type = 0
	if integration_type == 1:
		# os.system("taskset -p 0xff %d" % os.getpid())
		# num_processors = multiprocessing.cpu_count()
		# print num_processors
		# pool = multiprocessing.Pool()
		#with multiprocessing.Pool(num_processors) as pool:
		pool = multiprocessing.Pool()
		for iy in range(0,outRes):
			#print iy
			#pool.apply_async(xxxx, (iy,))
			result = pool.apply_async(gen_pixel_y, args=(yminIm+iy*dyIm, zminIm, dzIm, -1*xmax, xdel, zmin, zdel, inc_rad, scidata))
			image[iy] = result.get()
			#print image[iy]
		
		pool.close()
		# pool.join()
		# for iy in range(0,outRes):
			# image[iy] = result.get()
		
		#print image.shape
		#print scidata.shape
		
	else:
		for iy in range(0,outRes):
			image[iy] = gen_pixel_y(yminIm+iy*dyIm, zminIm, dzIm, -1*xmax, xdel, zmin, zdel, inc_rad, scidata)
			#print np.amax(image[iy])
	
	image*= volPix/cmxPix/cmyPix
	image = np.swapaxes(image,0,1)
	if SYMMETRY == 1:
		image_temp = np.flipud(image[1:])
		image = np.concatenate((image_temp, image))
	image = np.swapaxes(image,0,2)
	# image[e,y,x]
	
	# Energy File
	energy = np.zeros(enum)
	energy_number = np.zeros(enum)

	if os.access(fileout+'_energy', os.F_OK ):  os.remove(fileout+'_energy')
	energy_file = open(fileout+'_energy','w')
	
	for ip in range(0,enum):  # p
		energy[ip]=emin*efactor**ip
		energy_number[ip] = ip+1
		energy_file.write(str(energy_number[ip])+' '+str(energy[ip])+'\n')
	
	energy_file.close()
	
	#Write array to FITS image--------------------------------------------------  
	hduout = py.PrimaryHDU(image)
	
	if 'synchrotron_emiss' in emissfilename:
		hduout.header['UNIT'  ] = 'erg cm^-2 sr^-1 s^-1 Hz^-1' #hpbw'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		#hduout.header['HPBW'  ] = (hpbw, 'Half power beam width')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Frequency')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'frequency=CRVAL3*CDELT3^index')
	else:
		hduout.header['UNIT'  ] = 'MeV^2 cm^-2 sr^-1 s^-1 MeV^-1'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Energy (MeV)')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'energy=CRVAL3*CDELT3^index')
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 

	
### =======================================================
### Generates pixel for all x_im for specified y_im
###
def gen_pixel_yAB(y, nmin, dn, ymin, dy, zmin, dz, theta, data, dataAB):
	"""
	FORM OF DATA = data[y][z][energy/freq][x]
	LOS integral for all x of image at specified y of image
	"""
	kpc2cm = 3.08568025e21 #cm/kpc
	coords = gen_pixel_yindex(y, nmin, dn, ymin, dy, zmin, dz, theta)
	#print coords
	ic = np.around(coords)
	ic = ic.astype(int)
	#dic = coords-ic
	num_coords = coords.shape[0]
	
	iy_max = data.shape[0]-1
	iz_max = data.shape[1]-1
	
	x_pixel = np.zeros([data.shape[2], data.shape[3]])
	for i in range(0,num_coords):
		if ic[i][0]<0 or ic[i][1]<0 or ic[i][0]>=iy_max or ic[i][1]>=iz_max:
			pass
		else:
			x_pixel_temp = data[ic[i][0]] #+ dic[i][0]*(data[ic[i][0]+1]-data[ic[i][0]])
			x_pixel_temp = x_pixel_temp[ic[i][1]] #+ dic[i][1]*(x_pixel_temp[ic[i][1]+1]-x_pixel_temp[ic[i][1]])
			x_pixel_tempAB = dataAB[ic[i][0]] #+ dic[i][0]*(dataAB[ic[i][0]+1]-dataAB[ic[i][0]])
			x_pixel_tempAB = x_pixel_tempAB[ic[i][1]] #+ dic[i][1]*(x_pixel_tempAB[ic[i][1]+1]-x_pixel_tempAB[ic[i][1]])
			alpha = x_pixel_tempAB + 1e-99
			x_pixel_tempAB = np.exp(-1*alpha*dn*kpc2cm)
			x_pixel = x_pixel_temp*(-1)*np.expm1(-1*alpha*dn*kpc2cm)/alpha + x_pixel*x_pixel_tempAB
	
	print 'Finished y='+str(y)
	return x_pixel
	
	
### =======================================================
### LOS integrator for all galprop emission files
###	
def LOS_integral_radio(emissfilename, fileout, inclination, objectdist, outRes, int_steps, imageR, absorbfilename):
	"""
	Line of sight integrator for all emission files
	inclines by inclination (in degrees)
	projects a distance (in kpc)
	outputs in fits file
	"""
	import pyfits as py
	import os
	import sys
	from sys import stdout
	
	SYMMETRY = 1
	map_prop = 1.0 #image proportion of galaxy diameter
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	
	inc_rad = inclination*math.pi/180. #inclination in radians
	
	#READING EMISSION FITS-------------------------------------------------
	hdulistAB=py.open(absorbfilename)
	headAB=hdulistAB[0].header
	scidataAB=hdulistAB[0].data[0]
	
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
	
	print 'Data read'
	
	# Convert 2D -> 3D
	if SYMMETRY == 1:
		scidata = d2_to_d3sym(xmin, xdel, scidata)
		scidataAB = d2_to_d3sym(xmin, xdel, scidataAB)
	else:
		scidata = d2_to_d3(xmin, xdel, scidata)
		scidataAB = d2_to_d3(xmin, xdel, scidataAB)
	print 'Finished 2d->3D'
	
	original_mapsize = 2*xsize#math.sqrt((2.*xsize)**2.+zsize**2.) #kpc
	rotate_mapsize = map_prop*original_mapsize #kpc
	
	#OUTPUT IMAGE------------------------------------------------------	
	kpcPix = 2*imageR/(outRes-1) #kpc/pixel_dim
	
	if SYMMETRY == 1:
		image = np.zeros([enum,outRes,(outRes-1)/2+1])
		xminIm = xmin #-imageR#/2. #kpc
	else:
		image = np.zeros([enum,outRes,outRes]) #output map (pixels)
		xminIm = -imageR #/2. #kpc
	
	yminIm = -imageR#/2. #kpc
	zminIm = -rotate_mapsize/2. #kpc
	
	dxIm = kpcPix #kpc
	dyIm = kpcPix #kpc
	dzIm = rotate_mapsize/(int_steps-1) #kpc
	
	cmxPix = dxIm*kpc2cm #cm/pixel
	cmyPix = dyIm*kpc2cm #cm/pixel
	cmzPix = dzIm*kpc2cm #cm/pixel
	
	cmDist = objectdist*kpc2cm
	volPix = cmxPix*cmyPix*cmzPix #cm^3	
	srPix = dxIm*dyIm/objectdist**2. #steradian/pixel
	Pix2sac = srPix*Sr2sqArcSec #arcsec^2/pixel
	sac2Pix = 1./Pix2sac #pixel/arsec^2
	
	# scidata[energy/freq][z][y][x]
	scidata = np.swapaxes(scidata,0,3)
	if SYMMETRY == 1:
		scidata = gen_pixel_x(xmin, xdel, yminIm, dxIm, scidata, SYMMETRY)
	else:
		scidata = gen_pixel_x(-1*xmax, xdel, xminIm, dxIm, scidata, SYMMETRY)
	scidata = np.swapaxes(scidata,0,2)
	#scidata[y,z,x,e]
	
	#scidataAB[energy/freq][z][y][x]
	scidataAB = np.swapaxes(scidataAB,0,3)
	if SYMMETRY == 1:
		scidataAB = gen_pixel_x(xmin, xdel, yminIm, dxIm, scidataAB)
	else:
		scidataAB = gen_pixel_x(-1*xmax, xdel, xminIm, dxIm, scidataAB)
	scidataAB = np.swapaxes(scidataAB,0,2)
	#scidata[y,z,x,e]
	
	image = np.swapaxes(image,0,2)
	image = np.swapaxes(image,0,1)
	# image[y][x][energy/freq]
	
	print 'Starting Integration'
	
	integration_type = 0
	if integration_type == 1:
		## ATTEMPTED PARALLELIZATION
		# os.system("taskset -p 0xff %d" % os.getpid())
		# num_processors = multiprocessing.cpu_count()
		# print num_processors
		# pool = multiprocessing.Pool()
		#with multiprocessing.Pool(num_processors) as pool:
		pool = multiprocessing.Pool()
		for iy in range(0,outRes):
			#print iy
			#pool.apply_async(xxxx, (iy,))
			result = pool.apply_async(gen_pixel_yAB, args=(yminIm+iy*dyIm, zminIm, dzIm, -1*xmax, xdel, zmin, zdel, inc_rad, scidata, scidataAB))
			image[iy] = result.get()
			#print image[iy]
		
		pool.close()
		# pool.join()
		# for iy in range(0,outRes):
			# image[iy] = result.get()
		
		#print image.shape
		#print scidata.shape
		##
	else:
		for iy in range(0,outRes):
			image[iy] = gen_pixel_yAB(yminIm+iy*dyIm, zminIm, dzIm, -1*xmax, xdel, zmin, zdel, inc_rad, scidata, scidataAB)
			#print np.amax(image[iy])
	
	#image*= volPix/cmxPix/cmyPix
	image = np.swapaxes(image,0,1)
	if SYMMETRY == 1:
		image_temp = np.flipud(image[1:])
		image = np.concatenate((image_temp, image))
	image = np.swapaxes(image,0,2)
	# image[e,y,x]
	
	# Energy File
	energy = np.zeros(enum)
	energy_number = np.zeros(enum)

	if os.access(fileout+'_freq', os.F_OK ):  os.remove(fileout+'_freq')
	energy_file = open(fileout+'_freq','w')
	
	for ip in range(0,enum):  # p
		energy[ip]=emin*efactor**ip
		energy_number[ip] = ip+1
		energy_file.write(str(energy_number[ip])+' '+str(energy[ip])+'\n')
	
	energy_file.close()
	
	#Write array to FITS image--------------------------------------------------  
	hduout = py.PrimaryHDU(image)
	
	if 'synchrotron_emiss' or 'free_free_emiss' in emissfilename:
		hduout.header['UNIT'  ] = 'erg cm^-2 sr^-1 s^-1 Hz^-1' #hpbw'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		#hduout.header['HPBW'  ] = (hpbw, 'Half power beam width')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Frequency')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'frequency=CRVAL3*CDELT3^index')
	else:
		hduout.header['UNIT'  ] = 'MeV^2 cm^-2 sr^-1 s^-1 MeV^-1'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Energy (MeV)')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'energy=CRVAL3*CDELT3^index')
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 


### =======================================================
### Takes all gamma LOS files and combines them for total gamma LOS
###		
def gamma_total_los(file1, file2, file3, fileout):
	"""
	gamma_total_los(file1, file2, file3)
	Adds 3 gamma-ray los of sight fits files
	"""
	import pyfits as py
	import os
	
	#READING LOS FITS-------------------------------------------------
	hdulist1=py.open(file1)
	hdulist2=py.open(file2)
	hdulist3=py.open(file3)

	hdu1=hdulist1[0].header
	data1=hdulist1[0].data
	data2=hdulist2[0].data
	data3=hdulist3[0].data

	gamma_total=data1+data2+data3
	
	hduout = py.PrimaryHDU(gamma_total, header=hdu1)
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 
	

### =======================================================
### Takes gamma LOS and gives us flux
###	
def gamma_flux(los_file,fileout):
	"""
	gamma_flux(los_file, fileout)
		Calculates flux from gamma los file. And Fermi Integrated Flux.
	"""
	import os
	import pyfits as py
	#import matplotlib.pyplot as plt
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	
	#READING FITS-------------------------------------------------
	hdulist=py.open(los_file)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1'] #arcsec
	ymin = head['CRVAL2'] #arcsec
	emin = head['CRVAL3'] #MeV
	
	xdel = head['CDELT1'] #arcsec
	ydel = head['CDELT2'] #arcsec
	efactor = head['CDELT3'] #energy=emin*efactor^i
	
	xnum = head['NAXIS1']
	ynum = head['NAXIS2']
	enum = head['NAXIS3']
	
	inclination = head['INCLN1']
	objectdist = head['OBJDIS']
	objectdistCM = objectdist*kpc2cm
	pixelsr = xdel*ydel/Sr2sqArcSec
	
	fermi = 0.
	lumin_tot = 0.
	
	energy = np.zeros(enum)
	lumin = np.zeros(enum)
	for i in range(0,enum):
		energy[i] = emin*efactor**i
		lumin[i] = np.sum(scidata[i])*pixelsr
		if not math.isnan(lumin[i]): 
			lumin_tot += lumin[i]
			#if energy[i] >= 100. and energy[i] <= 10.**5. : 
			#	fermi += (efactor-1.)*energy[i]*lumin[i]
	
	lumin_tot *= math.log(efactor)*1.6021e-6*4*math.pi*objectdistCM**2 #to ergs/s

	# #PLOTTING DATA
	# AA=plt.loglog(energy, lumin)
	# plt.title('Gamma Flux')
	# plt.xlabel(r'E (MeV)')
	# plt.ylabel(r'MeV$^2$ s$^{-1}$ cm$^{-2}$ MeV$^{-1}$')
	# plt.savefig(fileout)
	# plt.close()
	
	#PRINTING DATA
	np.savetxt(fileout+'.dat', [energy, lumin], '%e')
	#print('Total Gamma-ray Luminosity = '+str(lumin_tot)+' erg s^-1')
	#print('Fermi Integrated flux (0.1GeV-100GeV) = '+str(fermi)+' MeV^2 s^-1 cm^-2')
	#print('gamma_flux COMPLETED')
	return lumin_tot


### =======================================================
### Takes gamma LOS and gives us spectrum per pixel
###	
def gamma_spectrum(los_file,fileout):
	"""
	gamma_spectrum(los_file, fileout)
		Calculates synchrotron spectrum from input LOS file.
	"""
	import os
	import numpy as np
	import pyfits as py
	import math
	#import matplotlib.pyplot as plt
	
	#READING FITS-------------------------------------------------
	hdulist=py.open(los_file)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1'] #arcsec
	ymin = head['CRVAL2'] #arcsec
	emin = head['CRVAL3'] #Hz
	
	xdel = head['CDELT1'] #arcsec
	ydel = head['CDELT2'] #arcsec
	efactor = head['CDELT3'] #energy=emin*efactor^i
	
	xnum = head['NAXIS1']
	ynum = head['NAXIS2']
	enum = head['NAXIS3']
	
	inclination = head['INCLN1']
	objectdist = head['OBJDIS']
	
	energy = np.zeros(enum)
	for i in range(0,enum):
		energy[i] = emin*efactor**i
	
	energy_s = np.zeros(enum-1)	
	index = np.zeros([enum-1,ynum,xnum])
	for n in range(0,enum-1):
		energy_s[n] = (energy[n]+energy[n+1])/2.
		logf2f1 = math.log(energy[n+1]/energy[n])
		for j in range(0,ynum):
			for i in range(0,xnum): 
				index[n][j][i] = math.log(scidata[n+1][j][i]/scidata[n][j][i]*energy[n]**2./energy[n+1]**2.+1e-40)/logf2f1
		
	#Write array to FITS image--------------------------------------------------  
	hduout = py.PrimaryHDU(index)
	hduout.header['UNIT'  ] = 'spectral index a: flux propto E^(a)'
	hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
	hduout.header['INCLN1'] = (inclination[0], 'Inclination angle 1')
	hduout.header['CRVAL1'] = (xmin, 'arcseconds')
	hduout.header['CRVAL2'] = (ymin, 'arcseconds')
	hduout.header['CRVAL3'] = (energy_s[0], 'Minimum Energy (MeV)')
	hduout.header['CDELT1'] = (xdel, 'arcseconds')
	hduout.header['CDELT2'] = (ydel, 'arcseconds')
	hduout.header['CDELT3'] = (efactor, 'energy=CRVAL3*CDELT3^index')
	hduout.header['OEMIN '] = (emin, 'Original Minimun Energy (MeV)')

	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 
	

### =======================================================
### Takes all radio LOS files and combines them for total radio LOS
###		
def radio_total_los(file1, file2, fileout):
	"""
	radio_total_los(file1, file2, file3)
	Adds 2 radio los of sight fits files
	"""
	import pyfits as py
	import os
	
	#READING LOS FITS-------------------------------------------------
	hdulist1=py.open(file1)
	hdulist2=py.open(file2)

	hdu1=hdulist1[0].header
	data1=hdulist1[0].data
	data2=hdulist2[0].data

	gamma_total=data1+data2
	
	hduout = py.PrimaryHDU(gamma_total, header=hdu1)
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout))
	
	
### =======================================================
### Takes synch LOS and gives us flux
###		
def synch_flux(los_file,fileout):
	"""
	synch_flux(los_file, fileout)
		Calculates synchrotron flux from input LOS file.
	"""
	import os
	import pyfits as py
	#import matplotlib.pyplot as plt
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	
	#READING FITS-------------------------------------------------
	hdulist=py.open(los_file)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1'] #arcsec
	ymin = head['CRVAL2'] #arcsec
	emin = head['CRVAL3'] #Hz
	
	xdel = head['CDELT1'] #arcsec
	ydel = head['CDELT2'] #arcsec
	efactor = head['CDELT3'] #frequency=emin*efactor^i
	
	xnum = head['NAXIS1']
	ynum = head['NAXIS2']
	enum = head['NAXIS3']
	
	inclination = head['INCLN1']
	objectdist = head['OBJDIS']
	objectdistCM = objectdist*kpc2cm
	pixelsr = xdel*ydel/Sr2sqArcSec
	
	energy = np.zeros(enum)
	lumin = np.zeros(enum)
	for i in range(0,enum):
		energy[i] = emin*efactor**i
		lumin[i] = np.sum(scidata[i])*pixelsr
		
	lumin_tot = np.sum(lumin*energy)*math.log(efactor)*4*math.pi*objectdistCM**2

	# #PLOTTING DATA
	# AA=plt.loglog(energy, lumin)
	# plt.title('Synch Flux')
	# plt.xlabel(r'$\nu$ (Hz)')
	# plt.ylabel(r'erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$')
	# plt.savefig(fileout)
	# plt.close()
	
	#PRINTING DATA
	np.savetxt(fileout+'.dat', [energy, lumin], '%e')
	#print('Synch total luminosity ='+str(lumin_tot)+'erg s^-1')
	#print('synch_flux COMPLETED')
	return lumin_tot

	
### =======================================================
### Takes synch LOS and gives us spectrum for each pixel
###		
def synch_spectrum(los_file,fileout):
	"""
	synch_spectrum(los_file, fileout)
		Calculates synchrotron spectrum from input LOS file.
	"""
	import os
	import numpy as np
	import pyfits as py
	import math
	#import matplotlib.pyplot as plt
	
	#READING FITS-------------------------------------------------
	hdulist=py.open(los_file)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1'] #arcsec
	ymin = head['CRVAL2'] #arcsec
	nmin = head['CRVAL3'] #Hz
	
	xdel = head['CDELT1'] #arcsec
	ydel = head['CDELT2'] #arcsec
	nfactor = head['CDELT3'] #freq=nmin*nfactor^i
	
	xnum = head['NAXIS1']
	ynum = head['NAXIS2']
	nnum = head['NAXIS3']
	
	inclination = head['INCLN1']
	#hpbw = head['HPBW']
	objectdist = head['OBJDIS']
	
	freq = np.zeros(nnum)
	for i in range(0,nnum):
		freq[i] = nmin*nfactor**i
	
	freq_s = np.zeros(nnum-1)	
	index = np.zeros([nnum-1,ynum,xnum])
	for n in range(0,nnum-1):
		freq_s[n] = (freq[n]+freq[n+1])/2.
		logf2f1 = math.log(freq[n+1]/freq[n])
		for j in range(0,ynum):
			for i in range(0,xnum): 
				index[n][j][i] = math.log(scidata[n+1][j][i]/scidata[n][j][i]*freq[n]/freq[n+1]+1e-40)/logf2f1
		
	#Write array to FITS image--------------------------------------------------  
	hduout = py.PrimaryHDU(index)
	hduout.header['UNIT'  ] = 'spectral index a: photon num density propto nu^(a)'
	hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
	hduout.header['INCLN1'] = (inclination, 'Inclination angle 1')
	#hduout.header['HPBW'  ] = (hpbw, 'Half power beam width')
	hduout.header['CRVAL1'] = (xmin, 'arcseconds')
	hduout.header['CRVAL2'] = (ymin, 'arcseconds')
	hduout.header['CRVAL3'] = (freq_s[0], 'Minimum Frequency (Hz)')
	hduout.header['CDELT1'] = (xdel, 'arcseconds')
	hduout.header['CDELT2'] = (ydel, 'arcseconds')
	hduout.header['CDELT3'] = (nfactor, 'frequency=CRVAL3*CDELT3^index')
	hduout.header['OFREQ '] = (nmin, 'Original minimum frequency (Hz)')

	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 
	#print('synch_spectrum COMPLETED')
	

### =======================================================
### Generates pixel for all x_im for specified y_im
###
def gen_pixel_yAB2(y, nmin, dn, ymin, dy, zmin, dz, theta, data, data2, dataAB):
	"""
	FORM OF DATA = data[y][z][energy/freq][x]
	LOS integral for all x of image at specified y of image
	"""
	kpc2cm = 3.08568025e21 #cm/kpc
	coords = gen_pixel_yindex(y, nmin, dn, ymin, dy, zmin, dz, theta)
	#print coords
	ic = np.around(coords)
	ic = ic.astype(int)
	#dic = coords-ic
	num_coords = coords.shape[0]
	
	iy_max = data.shape[0]-1
	iz_max = data.shape[1]-1
	
	x_pixel = np.zeros([data.shape[2], data.shape[3]])
	x_pixel2 = np.zeros([data.shape[2], data.shape[3]])
	for i in range(0,num_coords):
		if ic[i][0]<0 or ic[i][1]<0 or ic[i][0]>=iy_max or ic[i][1]>=iz_max:
			pass
		else:
			x_pixel_temp = data[ic[i][0]] #+ dic[i][0]*(data[ic[i][0]+1]-data[ic[i][0]])
			x_pixel_temp2 = data2[ic[i][0]] #+ dic[i][0]*(data2[ic[i][0]+1]-data2[ic[i][0]])
			x_pixel_tempAB = dataAB[ic[i][0]] #+ dic[i][0]*(dataAB[ic[i][0]+1]-dataAB[ic[i][0]])
			
			x_pixel_temp = x_pixel_temp[ic[i][1]] #+ dic[i][1]*(x_pixel_temp[ic[i][1]+1]-x_pixel_temp[ic[i][1]])
			x_pixel_temp2 = x_pixel_temp2[ic[i][1]] #+ dic[i][1]*(x_pixel_temp2[ic[i][1]+1]-x_pixel_temp2[ic[i][1]])
			x_pixel_tempAB = x_pixel_tempAB[ic[i][1]] #+ dic[i][1]*(x_pixel_tempAB[ic[i][1]+1]-x_pixel_tempAB[ic[i][1]])
			
			alpha = x_pixel_tempAB + 1e-99
			x_pixel_tempAB = np.exp(-1*alpha*dn*kpc2cm)
			x_pixel = x_pixel_temp*(-1)*np.expm1(-1*alpha*dn*kpc2cm)/alpha + x_pixel*x_pixel_tempAB
			x_pixel2 = x_pixel_temp2*(-1)*np.expm1(-1*alpha*dn*kpc2cm)/alpha + x_pixel2*x_pixel_tempAB
	
	print 'Finished y='+str(y)
	return x_pixel, x_pixel2

	
### =======================================================
### LOS integrator for all galprop emission files
###	
def LOS_integral_radio2(emissfilename, emissfilename2, fileout, fileout2, inclination, objectdist, outRes, int_steps, imageR, absorbfilename):
	"""
	Line of sight integrator for all emission files
	inclines by inclination (in degrees)
	projects a distance (in kpc)
	outputs in fits file
	"""
	import pyfits as py
	import os
	import sys
	from sys import stdout

	SYMMETRY = 1
	map_prop = 1.0 #image proportion of galaxy diameter
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	
	inc_rad = inclination*math.pi/180. #inclination in radians
	
	#READING EMISSION FITS-------------------------------------------------
	hdulistAB=py.open(absorbfilename)
	headAB=hdulistAB[0].header
	scidataAB=hdulistAB[0].data[0]
	
	hdulist2=py.open(emissfilename2)
	head2=hdulist2[0].header
	scidata2=hdulist2[0].data[0]
	
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
	
	print 'Data read'
	
	# Convert 2D -> 3D
	if SYMMETRY == 1:
		scidata = d2_to_d3sym(xmin, xdel, scidata)
		scidata2 = d2_to_d3sym(xmin, xdel, scidata2)
		scidataAB = d2_to_d3sym(xmin, xdel, scidataAB)
	else:
		scidata = d2_to_d3(xmin, xdel, scidata)
		scidata2 = d2_to_d3(xmin, xdel, scidata2)
		scidataAB = d2_to_d3(xmin, xdel, scidataAB)
	print 'Finished 2d->3D'
	
	original_mapsize = 2*xsize#math.sqrt((2.*xsize)**2.+zsize**2.) #kpc
	rotate_mapsize = map_prop*original_mapsize #kpc
	
	#OUTPUT IMAGE------------------------------------------------------	
	kpcPix = 2*imageR/(outRes-1) #kpc/pixel_dim
	
	if SYMMETRY == 1:
		image = np.zeros([enum,outRes,(outRes-1)/2+1])
		image2 = np.zeros([enum,outRes,(outRes-1)/2+1])
		xminIm = xmin #-imageR#/2. #kpc
	else:
		image = np.zeros([enum,outRes,outRes]) #output map (pixels)
		image2 = np.zeros([enum,outRes,outRes]) #output map (pixels)
		xminIm = -imageR #/2. #kpc
	
	yminIm = -imageR#/2. #kpc
	zminIm = -rotate_mapsize/2. #kpc
	
	dxIm = kpcPix #kpc
	dyIm = kpcPix #kpc
	dzIm = rotate_mapsize/(int_steps-1) #kpc
	
	cmxPix = dxIm*kpc2cm #cm/pixel
	cmyPix = dyIm*kpc2cm #cm/pixel
	cmzPix = dzIm*kpc2cm #cm/pixel
	
	cmDist = objectdist*kpc2cm
	volPix = cmxPix*cmyPix*cmzPix #cm^3	
	srPix = dxIm*dyIm/objectdist**2. #steradian/pixel
	Pix2sac = srPix*Sr2sqArcSec #arcsec^2/pixel
	sac2Pix = 1./Pix2sac #pixel/arsec^2
	
	# scidata[energy/freq][z][y][x]
	scidata = np.swapaxes(scidata,0,3)
	scidata2 = np.swapaxes(scidata2,0,3)
	scidataAB = np.swapaxes(scidataAB,0,3)
	if SYMMETRY == 1:
		scidata = gen_pixel_x(xmin, xdel, yminIm, dxIm, scidata, SYMMETRY)
		scidata2 = gen_pixel_x(xmin, xdel, yminIm, dxIm, scidata2, SYMMETRY)
		scidataAB = gen_pixel_x(xmin, xdel, yminIm, dxIm, scidataAB, SYMMETRY)
	else:
		scidata = gen_pixel_x(-1*xmax, xdel, xminIm, dxIm, scidata, SYMMETRY)
		scidata2 = gen_pixel_x(-1*xmax, xdel, xminIm, dxIm, scidata2, SYMMETRY)
		scidataAB = gen_pixel_x(-1*xmax, xdel, xminIm, dxIm, scidataAB, SYMMETRY)
	scidata = np.swapaxes(scidata,0,2)
	scidata2 = np.swapaxes(scidata2,0,2)
	scidataAB = np.swapaxes(scidataAB,0,2)
	#scidata[y,z,x,e]
	
	image = np.swapaxes(image,0,2)
	image2 = np.swapaxes(image2,0,2)
	
	image = np.swapaxes(image,0,1)
	image2 = np.swapaxes(image2,0,1)
	# image[y][x][energy/freq]
	
	print 'Starting Integration'
	
	integration_type = 0
	if integration_type == 1:
		## ATTEMPTED PARELLELIZATION
		# os.system("taskset -p 0xff %d" % os.getpid())
		# num_processors = multiprocessing.cpu_count()
		# print num_processors
		# pool = multiprocessing.Pool()
		#with multiprocessing.Pool(num_processors) as pool:
		pool = multiprocessing.Pool()
		for iy in range(0,outRes):
			#print iy
			#pool.apply_async(xxxx, (iy,))
			result = pool.apply_async(gen_pixel_yAB2, args=(yminIm+iy*dyIm, zminIm, dzIm, -1*xmax, xdel, zmin, zdel, inc_rad, scidata, scidata2, scidataAB))
			image_temp = result.get()
			image[iy] = image_temp[0]
			image2[iy] = image_temp[1]
			#print image[iy]
		
		pool.close()
		# pool.join()
		# for iy in range(0,outRes):
			# image[iy] = result.get()
		
		#print image.shape
		#print scidata.shape
		##
	else:
		for iy in range(0,outRes):
			image_temp = gen_pixel_yAB2(yminIm+iy*dyIm, zminIm, dzIm, -1*xmax, xdel, zmin, zdel, inc_rad, scidata, scidata2, scidataAB)
			image[iy] = image_temp[0]
			image2[iy] = image_temp[1]
			#print np.amax(image[iy])
	
	#image*= volPix/cmxPix/cmyPix
	image = np.swapaxes(image,0,1)
	image2 = np.swapaxes(image2,0,1)
	if SYMMETRY == 1:
		image_temp = np.flipud(image[1:])
		image_temp2 = np.flipud(image2[1:])
		image = np.concatenate((image_temp, image))
		image2 = np.concatenate((image_temp2, image2))
	image = np.swapaxes(image,0,2)
	image2 = np.swapaxes(image2,0,2)
	# image[e,y,x]
	
	# Energy File
	energy = np.zeros(enum)
	energy_number = np.zeros(enum)

	if os.access(fileout+'_freq', os.F_OK ):  os.remove(fileout+'_freq')
	energy_file = open(fileout+'_freq','w')
	
	for ip in range(0,enum):  # p
		energy[ip]=emin*efactor**ip
		energy_number[ip] = ip+1
		energy_file.write(str(energy_number[ip])+' '+str(energy[ip])+'\n')
	
	energy_file.close()
	
	#Write array to FITS image--------------------------------------------------  
	##IMAGE
	hduout = py.PrimaryHDU(image)
	
	if 'synchrotron_emiss' or 'free_free_emiss' in emissfilename:
		hduout.header['UNIT'  ] = 'erg cm^-2 sr^-1 s^-1 Hz^-1' #hpbw'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		#hduout.header['HPBW'  ] = (hpbw, 'Half power beam width')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Frequency')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'frequency=CRVAL3*CDELT3^index')
	else:
		hduout.header['UNIT'  ] = 'MeV^2 cm^-2 sr^-1 s^-1 MeV^-1'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Energy (MeV)')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'energy=CRVAL3*CDELT3^index')
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 
	
	#Write array to FITS image--------------------------------------------------  
	##IMAGE2
	hduout = py.PrimaryHDU(image2)
	
	if 'synchrotron_emiss' or 'free_free_emiss' in emissfilename2:
		hduout.header['UNIT'  ] = 'erg cm^-2 sr^-1 s^-1 Hz^-1' #hpbw'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		#hduout.header['HPBW'  ] = (hpbw, 'Half power beam width')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Frequency')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'frequency=CRVAL3*CDELT3^index')
	else:
		hduout.header['UNIT'  ] = 'MeV^2 cm^-2 sr^-1 s^-1 MeV^-1'
		hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
		hduout.header['INCLN1'] = (inclination, 'Inclination angle')
		hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CRVAL3'] = (emin, 'Minimum Energy (MeV)')
		hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
		hduout.header['CDELT3'] = (efactor, 'energy=CRVAL3*CDELT3^index')
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout2, os.F_OK ):  os.remove(fileout2)
	hduoutlist.writeto(fileout2)   
	print('FITS image output to '+str(fileout2)) 
	
	
	
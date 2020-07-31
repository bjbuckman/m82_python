##
## Functions for the magnetic field
##
def magnetic_field_2d(name, position, parameters):
	import numpy as np
	import math
	###
	###
	if name == 'm82_01':
		##
		## position = np.meshgrid[r_array, z_array]
		## parameters = [300.e-6, 0.2, 0.05, 0.05, -2., 5.e-6, 0.4, 2.0]
		##
		rr = position[0]
		zz = position[1]
		
		# print rr
		# print zz
		
		B0 =       parameters[0] #300.e-6 #microG
		r_scale =  parameters[1] # 0.2 #kpc
		r_scale1 = parameters[2]
		z_scale =  parameters[3] # #kpc 
		z_scale1 = parameters[4] # #kpc 
		Bz_power = parameters[5] #-2.

		B0_large = parameters[6] #5.0e-6 microG
		r_large =  parameters[7] #0.4
		z_large =  parameters[8] #2.0
		r_sun =	   parameters[9]
		
		B = np.zeros(rr.shape)
		
		#large scale B
		B = B0_large*np.exp(-abs(zz)/z_large)*np.exp(-(rr-r_sun)/r_large)
		
		#Truth arrays
		ARR0 = np.logical_and(rr<=r_scale, abs(zz)<=z_scale)
		ARR1 = np.logical_and(rr<=r_scale, abs(zz)> z_scale)
		ARR2 = np.logical_and(rr> r_scale, abs(zz)<=z_scale)
		ARR3 = np.logical_and(rr> r_scale, abs(zz)> z_scale)
		
		#small scale B
		B += B0*ARR0
		B += B0*((abs(zz)-z_scale+z_scale1)/z_scale1)**Bz_power*ARR1
		B += B0*np.exp(-1*(rr-r_scale)/r_scale1)*ARR2
		B += B0*np.exp(-1*(rr-r_scale)/r_scale1)*((abs(zz)-z_scale+z_scale1)/z_scale1)**Bz_power*ARR3
		
	elif name == 'constant+powerlaw_v2':
		##
		## position = np.meshgrid[r_array, z_array]
		## parameters = [300.e-6, 0.2, 0.05, 0.05, -2., 5.e-6, 0.4, 2.0]
		##
		rr = position[0]
		zz = position[1]
		
		rrrr = np.sqrt(rr**2+zz**2)
		
		B0 =       parameters[0] #300.e-6 #microG
		r_scale =  parameters[1] # 0.2 #kpc
		r_scale1 = parameters[2]
		z_scale =  parameters[3] # #kpc 
		z_scale1 = parameters[4] # #kpc 

		r_pow =    parameters[5]
		r_pow1 =   parameters[6]
		b_pow =    parameters[7]
		
		theta = np.arctan(np.abs(zz)/rr)
		theta_s = np.arctan(z_scale/r_scale)
		
		theta_logic = theta>=theta_s
		theta_logic_n = np.logical_not(theta_logic)
		
		r_newcore = np.zeros(rr.shape)
		r_newscale = np.zeros(rr.shape)
		
		r_newcore+= z_scale/np.sin(theta)*theta_logic
		r_newscale+= z_scale1/np.sin(theta)*theta_logic
		
		r_newcore+= r_scale/np.cos(theta)*theta_logic_n
		r_newscale+= r_scale1/np.cos(theta)*theta_logic_n
		
		pow_logic = rrrr>z_scale
		
		B = np.zeros(rr.shape)
		B_pow = np.zeros(rr.shape)
		
		#Truth arrays
		ARR0 = np.logical_and(rr<=r_scale, abs(zz)<=z_scale)
		ARR1 = np.logical_not(ARR0)
		
		#small scale B
		B += B0*ARR0
		B += B0*ARR1*np.exp(-(rrrr-r_newcore)/r_newscale)
		
		B_pow += pow_logic*((abs(rrrr)-r_pow+r_pow1)/r_pow1)**b_pow
		
		B = np.maximum(B,B_pow)
	
	elif name == 'constant+powerlaw_v3':
		##
		## position = np.meshgrid[r_array, z_array]
		## parameters = [300.e-6, 0.2, 0.05, 0.05, -2., 5.e-6, 0.4, 2.0]
		##
		rr = position[0]
		zz = position[1]
		
		r = np.sqrt(rr**2+zz**2)
		
		B0 =       parameters[0] #300.e-6 #microG
		r_scale =  parameters[1] # 0.2 #kpc
		r_scale1 = parameters[2]
		z_scale =  parameters[3] # #kpc 
		z_scale1 = parameters[4] # #kpc 

		r_pow =    parameters[5]
		r_pow1 =   parameters[6]
		b_pow =    parameters[7]
		
		theta = np.arctan2(np.abs(zz),rr)
		
		r_newcore = np.zeros(rr.shape)
		r_newscale = np.zeros(rr.shape)
		
		r_newcore+= np.sqrt( (z_scale*np.sin(theta))**2 + (r_scale*np.cos(theta))**2 ) 
		r_newscale+= np.sqrt( (z_scale1*np.sin(theta))**2 + (r_scale1*np.cos(theta))**2 ) 
		
		pow_logic = r>r_pow
		
		B = np.zeros(rr.shape)
		B_pow = np.zeros(rr.shape)
		
		#Truth arrays
		ARR0 = r<=r_newcore
		ARR1 = np.logical_not(ARR0)
		
		#small scale B
		B += B0*ARR0
		B += B0*ARR1*np.exp(-(r-r_newcore)/r_newscale)
		
		B_pow += B0*pow_logic*((abs(r)-r_pow+r_pow1)/r_pow1)**b_pow
		
		B = np.maximum(B,B_pow)
		
	elif name == 'constant+powerlaw_v4':
		##
		## position = np.meshgrid[r_array, z_array]
		## parameters = [300.e-6, 0.2, 0.05, 0.05, -2., 5.e-6, 0.4, 2.0]
		##
		rr = position[0]
		zz = position[1]
		
		r = np.sqrt(rr**2+zz**2)
		
		B0 =       parameters[0] #300.e-6 #microG
		r_scale =  parameters[1] # 0.2 #kpc
		r_scale1 = parameters[2]
		z_scale =  parameters[3] # #kpc 
		z_scale1 = parameters[4] # #kpc 

		r_pow =    parameters[5]
		r_pow1 =   parameters[6]
		b_pow =    parameters[7]
		
		theta = np.arctan2(np.abs(zz),rr)
		
		r_newcore = np.zeros(rr.shape)
		r_newscale = np.zeros(rr.shape)
		
		r_newcore+=  ( (np.sin(theta)/z_scale )**2 + (np.cos(theta)/r_scale )**2 )**-0.5 
		r_newscale+= ( (np.sin(theta)/z_scale1)**2 + (np.cos(theta)/r_scale1)**2 )**-0.5 
		
		pow_logic = r>r_pow
		
		B = np.zeros(rr.shape)
		B_pow = np.zeros(rr.shape)
		
		#Truth arrays
		ARR0 = r<=r_newcore
		ARR1 = np.logical_not(ARR0)
		
		#small scale B
		B += B0*ARR0
		B += B0*ARR1*np.exp(-(r-r_newcore)/r_newscale)
		
		B_pow += B0*pow_logic*((abs(r)-r_pow+r_pow1)/r_pow1)**b_pow
		
		B = np.maximum(B,B_pow)
			
	
	return B

def magnetic_gradient(model):
	import numpy as np
	import math
	
	RR = 0
	zmax = 3.0
	znum = 1000
	z = np.linspace(-1*zmax, zmax, znum)
	
	position = np.meshgrid([RR], z)
	
	if model == 1:
		magnetic_field_name = 'm82_01'
		position = np.meshgrid([RR], z)
		parameters = [250.e-6, 0.2, 1.0, 0.05, 4.0, -2.0, 5.e-6, 10.0, 2.0, 8.5]
	else:
		magnetic_field_name = 'm82_01'
		position = np.meshgrid([RR], z)
		parameters = [250.e-6, 0.2, 1.0, 0.05, 4.0, -2.0, 5.e-6, 10.0, 2.0, 8.5]

	B = magnetic_field_2d(magnetic_field_name, position, parameters)
	B = np.average(B, axis=1)

	#(zB, U_B)
	U_B = B**2/(8*math.pi)

	#(zGB, grad_U_B)
	zGB = (z[1:]+z[:-1])/2
	grad_U_B = (U_B[1:]-U_B[:-1])/(z[1:]-z[:-1])
	grad_U_B*= -1
	
	return zGB, grad_U_B
	
	
	
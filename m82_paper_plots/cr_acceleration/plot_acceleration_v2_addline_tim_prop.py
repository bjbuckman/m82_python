import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pyfits as py

import imp
gas_density = imp.load_source("gas_density_v4", "../m82_dist_py/gas_density_v4.py")
magnetic_field = imp.load_source("magnetic_field_v4", "../m82_dist_py/magnetic_field_v4.py")
accel = imp.load_source("gp2d_acceleration_v1", "../m82_dist_py/gp2d_acceleration_v1.py")

c = 3.0e10 ## cm/s
mp = 1.66e-24 ## g
ev2ergs = 1.60218e-12
kpc2cm = 3.08568025e21 ## cm/kpc

## z-array
zmin = 0.0
zmax = 4.0
znum = 401
z = np.linspace(0, zmax, znum)

###
### GAS DENSITY FOR MODEL A (1)
###
nHI_factor = 0.95
nHII_factor = 0.05

density = 150.

r0 = 0.200 ## kpc
r1 = 0.20
z0 = 0.050 ## kpc
z1 = 0.05

v0 = 200.
mass_loss = 5.

nHI = np.zeros(znum)
for iz in range(0,znum):
	nHI[iz] = gas_density.m82_wind_radial_v7( 0, 0, z[iz], density, r0, r1, z0, z1, v0, mass_loss)

nHII_1 = nHII_factor*nHI
nHI_1 = nHI_factor*nHI

###
### [z, nHI] [z, nHII]
#########################

###
### GAS DENSITY FOR MODEL B (2)
###
density = 700.

r0 = 0.200 ## kpc
r1 = 0.20
z0 = 0.050 ## kpc
z1 = 0.05

v0 = 200.
mass_loss = 30.

nHI2 = np.zeros(znum)
for iz in range(0,znum):
	nHI2[iz] = gas_density.m82_wind_radial_v7( 0, 0, z[iz], density, r0, r1, z0, z1, v0, mass_loss)

nHII_2 = nHII_factor*nHI2
nHI_2 = nHI_factor*nHI2

###
### [z, nHI] [z, nHII]
#########################

###
### GAS DENSITY FOR MODEL B' (3)
### 700.,0.2,0.2,0.05,0.05,200.,30.,7,8,9  
density = 900.

r0 = 0.200 ## kpc
r1 = 0.20
z0 = 0.050 ## kpc
z1 = 0.05

v0 = 200.
mass_loss = 30.

nHI2 = np.zeros(znum)
for iz in range(0,znum):
	nHI2[iz] = gas_density.m82_wind_radial_v7( 0, 0, z[iz], density, r0, r1, z0, z1, v0, mass_loss)

nHII_3 = nHII_factor*nHI2
nHI_3 = nHI_factor*nHI2

###
### [z, nHI] [z, nHII]
#########################

####################################################
## ##################################################
## MAGNETIC FIELD 1 #################################
## ##################################################
####################################################
###
### MAGNETIC FIELD 1_1
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000150,0.2,0.2,0.05,0.05,0.05,0.2,-1.4,8,9  

B0 = 150.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.1



B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_1_1 = (abs(grad_magnetic))/(mp*(nHI_1+nHII_1))

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 1_2
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000150,0.2,0.2,0.05,0.05,0.05,0.2,-1.4,8,9  

B0 = 150.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -0.9



B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_1_2 = (abs(grad_magnetic))/(mp*(nHI_1+nHII_1))

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 1_3
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000150,0.2,0.2,0.05,0.05,0.05,0.2,-1.4,8,9  

B0 = 150.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_1_3 = (abs(grad_magnetic))/(mp*(nHI_1+nHII_1))

###
### [z, B] [z, grad_magnetic]
##############################################


####################################################
## ##################################################
## MAGNETIC FIELD 2 #################################
## ##################################################
####################################################
###
### MAGNETIC FIELD 2_1
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.0

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_2_1 = (abs(grad_magnetic))/(mp*(nHI_2+nHII_2))

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 2_2
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.3

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_2_2 = (abs(grad_magnetic))/(mp*(nHI_2+nHII_2))

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 2_3
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.2

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_2_3 = (abs(grad_magnetic))/(mp*(nHI_2+nHII_2))

###
### [z, B] [z, grad_magnetic]
##############################################

####################################################
## ##################################################
## MAGNETIC FIELD 3 #################################
## ##################################################
####################################################
###
### MAGNETIC FIELD 3_1
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -0.1

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_3_1 = (abs(grad_magnetic))/(mp*(nHI_3+nHII_3))

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 3_2
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -0.3

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_3_2 = (abs(grad_magnetic))/(mp*(nHI_3+nHII_3))

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 3_3
###

magnetic_model_name = 'constant+powerlaw_v4'
##B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -0.2

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B = np.zeros(znum)
for iz in range(0,znum):
	B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
	
U_magnetic = B**2./(8.*math.pi)

grad_magnetic_z = z[:-1]+(z[1:]-z[:-1])/2
grad_magnetic = (U_magnetic[1:]-U_magnetic[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_magnetic*= -1

grad_magnetic = np.interp(z, grad_magnetic_z, grad_magnetic)

grad_B = (B[1:]-B[:-1])/(z[1:]-z[:-1])/kpc2cm
grad_B = np.interp(z, grad_magnetic_z, grad_B)
magnetic_tension = B/3.**0.5*grad_B/(4.*math.pi)  

magnetic_accel_3_3 = (abs(grad_magnetic))/(mp*(nHI_3+nHII_3))

###
### [z, B] [z, grad_magnetic]
##############################################

# ############################################
# ############################################
# ############################################

###
### Interstellar Radiation Field
###

isrf_file = '../m82_dist_py/ISRF_m82_v4_flux_z.dat'

data = np.loadtxt(isrf_file)

z_isrf = data[0]
isrf_dat = data[1]#/c

kappa = 1000.
# grad_isrf_z = z_isrf[:-1]+(z_isrf[:-1]-z_isrf[1:])/2.
# grad_isrf = (isrf_dat[1:]-isrf_dat[:-1])/(z_isrf[1:]-z_isrf[:-1])/kpc2cm 


grad_isrf_z = np.interp(z, z_isrf, isrf_dat)
# isrf_accel = grad_isrf_z/(mp*(nHI+nHII))*ev2ergs
# isrf_accel2 = grad_isrf_z/(mp*(nHI2+nHII2))*ev2ergs
isrf_accel = grad_isrf_z*kappa#*ev2ergs
isrf_accel2 = grad_isrf_z*kappa#*ev2ergs
# print isrf_accel
# print isrf_accel2

###
### [z, U_isrf] [grad_isrf_z, grad_isrf]
###########################################

# ############################################
# ############################################
# ############################################

###
### M82
###

m82_z, m82_accel = accel.m82_accel_r()
grav_accel = np.interp(z, m82_z, m82_accel)

###
### [z, m82_accel]

# ############################################
# ############################################
# ############################################

###
### PLOT
###

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 24}

matplotlib.rc('font', **font)

cmap = plt.get_cmap('cubehelix')
fig = plt.figure()
	
f1, ax = plt.subplots(1,1)#, figsize=(18,6))

for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(2)
ax.tick_params(which='major',width=1, length=8, labelsize='small')
ax.tick_params(which='minor',width=0.5, length=5)

MSIZE = 5
LWIDTH = 3
ALPHA = 0.65

# ax.plot(z, abs(grav_accel), c=cmap(0), label='Gravity', linestyle='-', linewidth=LWIDTH, alpha=1.0)
# ax.plot(z, abs(isrf_accel), c=cmap(0.2), label='ISRF', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)

color = cmap(0.5)
ax.fill_between(z, abs(grav_accel), 1.2*abs(grav_accel), facecolor=color, edgecolor='k', alpha=0.5, label='Gravity', zorder=20)

color = cmap(0.75)
ax.fill_between(z, abs(isrf_accel), 1.2*abs(isrf_accel), facecolor=color, edgecolor='k', alpha=0.5, label='ISRF', zorder=19)


###
###
###
### PARTICLES
###
###
###
# MODEL = 1
# MODEL = 4

# if MODEL == 0:
	# galdef = 'm82_21adn1214'
# elif MODEL == 1:
	# galdef = 'm82_21bbv2223'
# elif MODEL == 2:
	# galdef = 'm82_21bdv1411'
# elif MODEL == 3:
	# galdef = 'm82_21azg1212'
# elif MODEL == 4:
	# galdef = 'm82_21azb1014'

dir_outfits = '../fits/raw/'
dir_comp = '../fits/'

galdef_name = []

GG = 'm82_05aff1017'
galdef_name.append(GG)

GG =  'm82_05b2ff1020'
galdef_name.append(GG)

GG =  'm82_06bff1020'
galdef_name.append(GG)

galdef_name_1 = []

GG_1 = 'm82_05aff1019'
galdef_name_1.append(GG_1)

GG_1 = 'm82_05b2ff1022'
galdef_name_1.append(GG_1)

GG_1 = 'm82_06bff1022'
galdef_name_1.append(GG_1)

num_graph = len(galdef_name)

chi = np.zeros(num_graph)
chi_fac = np.ones(num_graph)
for n in range(0,num_graph):
	chi_dat = np.loadtxt(dir_comp+galdef_name[n]+'/'+galdef_name[n]+'chi.dat')
	chi_fac[n] = chi_dat[4]
	GAMMA_CHI = 1
	chi[n] = GAMMA_CHI*(chi_dat[1]+chi_dat[2])+chi_dat[3]

NORM_PLOT = 0
if NORM_PLOT == 1:
	FACTOR = np.ones(num_graph_galdef_tot)
else:
	FACTOR = chi_fac	
	
chi_1 = np.zeros(num_graph)
chi_fac_1 = np.ones(num_graph)
for n in range(0,num_graph):
	chi_dat_1 = np.loadtxt(dir_comp+galdef_name_1[n]+'/'+galdef_name_1[n]+'chi.dat')
	chi_fac_1[n] = chi_dat_1[4]
	GAMMA_CHI = 1
	chi_1[n] = GAMMA_CHI*(chi_dat_1[1]+chi_dat_1[2])+chi_dat_1[3]

if NORM_PLOT == 1:
	FACTOR_1 = np.ones(num_graph_galdef_tot)
else:
	FACTOR_1 = chi_fac_1	

# for n in 1:#range(0,num_graph):
n = 1

nuclei_file = dir_outfits+galdef_name[n]+'/'+'nuclei_full_54_'+galdef_name[n]+'.gz'
nuclei_file_1 = dir_outfits+galdef_name_1[n]+'/'+'nuclei_full_54_'+galdef_name_1[n]+'.gz'

###
### PROTONS
###

protons_z, protons_accel = accel.CR_proton_accel_z(nuclei_file, 0.)
protons_accel = np.interp(z, protons_z, protons_accel)*FACTOR[n]

protons_z_1, protons_accel_1 = accel.CR_proton_accel_z(nuclei_file_1, 0.)
protons_accel_1 = np.interp(z, protons_z_1, protons_accel_1)*FACTOR_1[n]

###
### [z, protons_accel]
#########################

###
### ELECTRONS
###

electrons_z, electrons_accel = accel.CR_electron_accel_z(nuclei_file, 0.)
electrons_accel = np.interp(z, electrons_z, electrons_accel)*FACTOR[n]

electrons_z_1, electrons_accel_1 = accel.CR_electron_accel_z(nuclei_file_1, 0.)
electrons_accel_1 = np.interp(z, electrons_z_1, electrons_accel_1)*FACTOR_1[n]

if n == 0:
	protons_accel/= (nHI_1+nHII_1)
	electrons_accel/= (nHI_1+nHII_1)
	
	protons_accel_1/= (nHI_1+nHII_1)
	electrons_accel_1/= (nHI_1+nHII_1)
elif n == 1:
	protons_accel/= (nHI_2+nHII_2)
	electrons_accel/= (nHI_2+nHII_2)
	
	protons_accel_1/= (nHI_2+nHII_2)
	electrons_accel_1/= (nHI_2+nHII_2)
elif n == 2:
	protons_accel/= (nHI_3+nHII_3)
	electrons_accel/= (nHI_3+nHII_3)
	
	protons_accel_1/= (nHI_3+nHII_3)
	electrons_accel_1/= (nHI_3+nHII_3)
	
###
### [z, electrons_accel]
#########################

color = cmap(0.25) #cmap(1-float(n+1)/(num_graph+1.))

if n == 0:
	ax.fill_between(z, abs(magnetic_accel_1_1), abs(magnetic_accel_1_2), facecolor=color, edgecolor=color, alpha=0.2, label='Magnetic Fields')#, hatch='|')
	# ax.plot(z, isrf_accel, c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, abs(magnetic_accel_1), c=color, label='magnetic', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, abs(isrf_accel), c=color, label='isrf', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	ax.plot(z, abs(magnetic_accel_1_3), c=color, linestyle='-.', linewidth=LWIDTH, alpha=1.0)
elif n == 1:
	ax.fill_between(z, abs(magnetic_accel_2_1), abs(magnetic_accel_2_2), facecolor=color, edgecolor='k', alpha=0.5, label='Magnetic Fields', zorder=18)#, hatch='|')
	# ax.plot(z, isrf_accel2, c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, abs(magnetic_accel_2), c=color, label='magnetic', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, abs(isrf_accel2), c=color, label='isrf', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, abs(magnetic_accel_2_3), c=color, linestyle='-.', linewidth=LWIDTH, alpha=1.0)
# elif n == 2:
	# ax.fill_between(z, abs(magnetic_accel_3_1), abs(magnetic_accel_3_2), facecolor=color, edgecolor=color, alpha=0.2)#, hatch='|')
	# # ax.plot(z, isrf_accel2, c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
	# # ax.plot(z, abs(magnetic_accel_2), c=color, label='magnetic', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	# # ax.plot(z, abs(isrf_accel2), c=color, label='isrf', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, abs(magnetic_accel_3_3), c=color, linestyle='-.', linewidth=LWIDTH, alpha=1.0)

color_p = cmap(0.0)
# color_e = cmap(0.60)
ax.fill_between(z, protons_accel, 1.2*protons_accel_1, facecolor=color_p, edgecolor='k', alpha=0.5, label='Cosmic-Rays', zorder=18)
# ax.fill_between(z, electrons_accel, electrons_accel_1, facecolor=color_e, edgecolor=color_e, alpha=1., linewidth=LWIDTH*2., label='CRe')
# ax.plot(z, protons_accel, c=color, label='CR_proton', linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z, electrons_accel, c=color, label='CR_electron', linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z, electrons_accel_1, c=color, label='CR_electron', linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
# print electrons_accel


# ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.0,3.5])
ax.set_ylim([3.e-11,3.e-6])
ax.set_xlabel(r'$z$ [kpc]')
ax.set_ylabel(r'Acceleration [cm s$^{-2}$]')

# ax.text(0.66, 0.75, r'magnetic field', ha='center', va='center', transform=ax.transAxes)
# ax.text(0.66, 0.55, r'gravity', ha='center', va='center', transform=ax.transAxes)
# ax.text(0.66, 0.35, r'cosmic rays', ha='center', va='center', transform=ax.transAxes)

# ax.set_title(r'$0$ Temperature Acceleration')

ax.legend(loc='upper right', prop={'size':14})

plt.tight_layout()
plot_name = 'plot_accel_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break
plt.close()




















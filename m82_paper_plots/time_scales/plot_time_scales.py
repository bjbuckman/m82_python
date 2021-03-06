import sys
import os
import numpy as np
import pyfits as py
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import imp
gas_density = imp.load_source("gas_density_v4", "../m82_dist_py/gas_density_v4.py")
magnetic_field = imp.load_source("magnetic_field_v4", "../m82_dist_py/magnetic_field_v4.py")



c = 3.0e10
ev2ergs = 1.60218e-12
mp = 1.66e-24 #g
kpc2cm = 3.08568025e21 #cm/kpc

# zmin = 0.0
# zmax = 4.0
# znum = 4001
# z = np.linspace(0, zmax, znum)

z = np.array([0.,0.5,1.0,1.5])
znum = len(z)

LABEL = [r'$z=0.0$ kpc', r'$z=0.5$ kpc', r'$z=1.0$ kpc', r'$z=1.5$ kpc']

###
### GAS DENSITY 1
###
nHI_factor = 0.95
nHII_factor = 0.05

density = 150./100.

r0 = 0.200 #kpc
r1 = 0.20
z0 = 0.050 #kpc
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
### GAS DENSITY 2
###
density = 675./100.

r0 = 0.200 #kpc
r1 = 0.20
z0 = 0.050 #kpc
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
### GAS DENSITY 3
###
density = 1000./100.

r0 = 0.200 #kpc
r1 = 0.20
z0 = 0.050 #kpc
z1 = 0.05

v0 = 200.
mass_loss = 30.

nHI3 = np.zeros(znum)
for iz in range(0,znum):
	nHI3[iz] = gas_density.m82_wind_radial_v7( 0, 0, z[iz], density, r0, r1, z0, z1, v0, mass_loss)

nHII_3 = nHII_factor*nHI2
nHI_3 = nHI_factor*nHI2

###
### [z, nHI] [z, nHII]
#########################

###
### MAGNETIC FIELD 1
###

magnetic_model_name = 'constant+powerlaw_v4'

B0 = 150.e-6/100.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.0

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B_1 = np.zeros(znum)
for iz in range(0,znum):
	B_1[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 2
###

magnetic_model_name = 'constant+powerlaw_v4'

B0 = 325.e-6/100.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.2

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B_2 = np.zeros(znum)
for iz in range(0,znum):
	B_2[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)

###
### [z, B] [z, grad_magnetic]
##############################################

###
### MAGNETIC FIELD 3
###

magnetic_model_name = 'constant+powerlaw_v4'

B0 = 325.e-6/100.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -0.2

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B_3 = np.zeros(znum)
for iz in range(0,znum):
	B_3[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)

###
### [z, B] [z, grad_magnetic]
##############################################

###
### Alfven velocities
###

vA_a = B_1*100.e-6/np.sqrt(4.*np.pi*(nHII_1+nHII_1)*100.*mp)#/1.e5
vA_b = B_2*100.e-6/np.sqrt(4.*np.pi*(nHII_2+nHII_2)*100.*mp)#/1.e5
vA_c = B_3*100.e-6/np.sqrt(4.*np.pi*(nHII_3+nHII_3)*100.*mp)#/1.e5

# print vA_a

###
### Wind velocities
###
def wind_total_data_1(r):
	Z = np.array([0.0000e+00,1.0000e-02,2.0000e-02,3.0000e-02,4.0000e-02,5.0000e-02,6.0000e-02,7.0000e-02,8.0000e-02,9.0000e-02,1.0000e-01,1.1000e-01,1.2000e-01,1.3000e-01,1.4000e-01,1.5000e-01,1.6000e-01,1.7000e-01,1.8000e-01,1.9000e-01,2.0000e-01,2.1000e-01,2.2000e-01,2.3000e-01,2.4000e-01,2.5000e-01,2.6000e-01,2.7000e-01,2.8000e-01,2.9000e-01,3.0000e-01,3.1000e-01,3.2000e-01,3.3000e-01,3.4000e-01,3.5000e-01,3.6000e-01,3.7000e-01,3.8000e-01,3.9000e-01,4.0000e-01,4.1000e-01,4.2000e-01,4.3000e-01,4.4000e-01,4.5000e-01,4.6000e-01,4.7000e-01,4.8000e-01,4.9000e-01,5.0000e-01,5.1000e-01,5.2000e-01,5.3000e-01,5.4000e-01,5.5000e-01,5.6000e-01,5.7000e-01,5.8000e-01,5.9000e-01,6.0000e-01,6.1000e-01,6.2000e-01,6.3000e-01,6.4000e-01,6.5000e-01,6.6000e-01,6.7000e-01,6.8000e-01,6.9000e-01,7.0000e-01,7.1000e-01,7.2000e-01,7.3000e-01,7.4000e-01,7.5000e-01,7.6000e-01,7.7000e-01,7.8000e-01,7.9000e-01,8.0000e-01,8.1000e-01,8.2000e-01,8.3000e-01,8.4000e-01,8.5000e-01,8.6000e-01,8.7000e-01,8.8000e-01,8.9000e-01,9.0000e-01,9.1000e-01,9.2000e-01,9.3000e-01,9.4000e-01,9.5000e-01,9.6000e-01,9.7000e-01,9.8000e-01,9.9000e-01,1.0000e+00,1.0100e+00,1.0200e+00,1.0300e+00,1.0400e+00,1.0500e+00,1.0600e+00,1.0700e+00,1.0800e+00,1.0900e+00,1.1000e+00,1.1100e+00,1.1200e+00,1.1300e+00,1.1400e+00,1.1500e+00,1.1600e+00,1.1700e+00,1.1800e+00,1.1900e+00,1.2000e+00,1.2100e+00,1.2200e+00,1.2300e+00,1.2400e+00,1.2500e+00,1.2600e+00,1.2700e+00,1.2800e+00,1.2900e+00,1.3000e+00,1.3100e+00,1.3200e+00,1.3300e+00,1.3400e+00,1.3500e+00,1.3600e+00,1.3700e+00,1.3800e+00,1.3900e+00,1.4000e+00,1.4100e+00,1.4200e+00,1.4300e+00,1.4400e+00,1.4500e+00,1.4600e+00,1.4700e+00,1.4800e+00,1.4900e+00,1.5000e+00,1.5100e+00,1.5200e+00,1.5300e+00,1.5400e+00,1.5500e+00,1.5600e+00,1.5700e+00,1.5800e+00,1.5900e+00,1.6000e+00,1.6100e+00,1.6200e+00,1.6300e+00,1.6400e+00,1.6500e+00,1.6600e+00,1.6700e+00,1.6800e+00,1.6900e+00,1.7000e+00,1.7100e+00,1.7200e+00,1.7300e+00,1.7400e+00,1.7500e+00,1.7600e+00,1.7700e+00,1.7800e+00,1.7900e+00,1.8000e+00,1.8100e+00,1.8200e+00,1.8300e+00,1.8400e+00,1.8500e+00,1.8600e+00,1.8700e+00,1.8800e+00,1.8900e+00,1.9000e+00,1.9100e+00,1.9200e+00,1.9300e+00,1.9400e+00,1.9500e+00,1.9600e+00,1.9700e+00,1.9800e+00,1.9900e+00,2.0000e+00,2.0100e+00,2.0200e+00,2.0300e+00,2.0400e+00,2.0500e+00,2.0600e+00,2.0700e+00,2.0800e+00,2.0900e+00,2.1000e+00,2.1100e+00,2.1200e+00,2.1300e+00,2.1400e+00,2.1500e+00,2.1600e+00,2.1700e+00,2.1800e+00,2.1900e+00,2.2000e+00,2.2100e+00,2.2200e+00,2.2300e+00,2.2400e+00,2.2500e+00,2.2600e+00,2.2700e+00,2.2800e+00,2.2900e+00,2.3000e+00,2.3100e+00,2.3200e+00,2.3300e+00,2.3400e+00,2.3500e+00,2.3600e+00,2.3700e+00,2.3800e+00,2.3900e+00,2.4000e+00,2.4100e+00,2.4200e+00,2.4300e+00,2.4400e+00,2.4500e+00,2.4600e+00,2.4700e+00,2.4800e+00,2.4900e+00,2.5000e+00,2.5100e+00,2.5200e+00,2.5300e+00,2.5400e+00,2.5500e+00,2.5600e+00,2.5700e+00,2.5800e+00,2.5900e+00,2.6000e+00,2.6100e+00,2.6200e+00,2.6300e+00,2.6400e+00,2.6500e+00,2.6600e+00,2.6700e+00,2.6800e+00,2.6900e+00,2.7000e+00,2.7100e+00,2.7200e+00,2.7300e+00,2.7400e+00,2.7500e+00,2.7600e+00,2.7700e+00,2.7800e+00,2.7900e+00,2.8000e+00,2.8100e+00,2.8200e+00,2.8300e+00,2.8400e+00,2.8500e+00,2.8600e+00,2.8700e+00,2.8800e+00,2.8900e+00,2.9000e+00,2.9100e+00,2.9200e+00,2.9300e+00,2.9400e+00,2.9500e+00,2.9600e+00,2.9700e+00,2.9800e+00,2.9900e+00,3.0000e+00,3.0100e+00,3.0200e+00,3.0300e+00,3.0400e+00,3.0500e+00,3.0600e+00,3.0700e+00,3.0800e+00,3.0900e+00,3.1000e+00,3.1100e+00,3.1200e+00,3.1300e+00,3.1400e+00,3.1500e+00,3.1600e+00,3.1700e+00,3.1800e+00,3.1900e+00,3.2000e+00,3.2100e+00,3.2200e+00,3.2300e+00,3.2400e+00,3.2500e+00,3.2600e+00,3.2700e+00,3.2800e+00,3.2900e+00,3.3000e+00,3.3100e+00,3.3200e+00,3.3300e+00,3.3400e+00,3.3500e+00,3.3600e+00,3.3700e+00,3.3800e+00,3.3900e+00,3.4000e+00,3.4100e+00,3.4200e+00,3.4300e+00,3.4400e+00,3.4500e+00,3.4600e+00,3.4700e+00,3.4800e+00,3.4900e+00,3.5000e+00,3.5100e+00,3.5200e+00,3.5300e+00,3.5400e+00,3.5500e+00,3.5600e+00,3.5700e+00,3.5800e+00,3.5900e+00,3.6000e+00,3.6100e+00,3.6200e+00,3.6300e+00,3.6400e+00,3.6500e+00,3.6600e+00,3.6700e+00,3.6800e+00,3.6900e+00,3.7000e+00,3.7100e+00,3.7200e+00,3.7300e+00,3.7400e+00,3.7500e+00,3.7600e+00,3.7700e+00,3.7800e+00,3.7900e+00,3.8000e+00,3.8100e+00,3.8200e+00,3.8300e+00,3.8400e+00,3.8500e+00,3.8600e+00,3.8700e+00,3.8800e+00,3.8900e+00,3.9000e+00,3.9100e+00,3.9200e+00,3.9300e+00,3.9400e+00,3.9500e+00,3.9600e+00,3.9700e+00,3.9800e+00,3.9900e+00,4.0000e+00,4.0100e+00,4.0200e+00,4.0300e+00,4.0400e+00,4.0500e+00,4.0600e+00,4.0700e+00,4.0800e+00,4.0900e+00,4.1000e+00,4.1100e+00,4.1200e+00,4.1300e+00,4.1400e+00,4.1500e+00,4.1600e+00,4.1700e+00,4.1800e+00,4.1900e+00,4.2000e+00,4.2100e+00,4.2200e+00,4.2300e+00,4.2400e+00,4.2500e+00,4.2600e+00,4.2700e+00,4.2800e+00,4.2900e+00,4.3000e+00,4.3100e+00,4.3200e+00,4.3300e+00,4.3400e+00,4.3500e+00,4.3600e+00,4.3700e+00,4.3800e+00,4.3900e+00,4.4000e+00,4.4100e+00,4.4200e+00,4.4300e+00,4.4400e+00,4.4500e+00,4.4600e+00,4.4700e+00,4.4800e+00,4.4900e+00,4.5000e+00,4.5100e+00,4.5200e+00,4.5300e+00,4.5400e+00,4.5500e+00,4.5600e+00,4.5700e+00,4.5800e+00,4.5900e+00,4.6000e+00,4.6100e+00,4.6200e+00,4.6300e+00,4.6400e+00,4.6500e+00,4.6600e+00,4.6700e+00,4.6800e+00,4.6900e+00,4.7000e+00,4.7100e+00,4.7200e+00,4.7300e+00,4.7400e+00,4.7500e+00,4.7600e+00,4.7700e+00,4.7800e+00,4.7900e+00,4.8000e+00,4.8100e+00,4.8200e+00,4.8300e+00,4.8400e+00,4.8500e+00,4.8600e+00,4.8700e+00,4.8800e+00,4.8900e+00,4.9000e+00,4.9100e+00,4.9200e+00,4.9300e+00,4.9400e+00,4.9500e+00,4.9600e+00,4.9700e+00,4.9800e+00,4.9900e+00,5.0000e+00,5.0100e+00,5.0200e+00,5.0300e+00,5.0400e+00,5.0500e+00,5.0600e+00,5.0700e+00,5.0800e+00,5.0900e+00,5.1000e+00,5.1100e+00,5.1200e+00,5.1300e+00,5.1400e+00,5.1500e+00,5.1600e+00,5.1700e+00,5.1800e+00,5.1900e+00,5.2000e+00,5.2100e+00,5.2200e+00,5.2300e+00,5.2400e+00,5.2500e+00,5.2600e+00,5.2700e+00,5.2800e+00,5.2900e+00,5.3000e+00,5.3100e+00,5.3200e+00,5.3300e+00,5.3400e+00,5.3500e+00,5.3600e+00,5.3700e+00,5.3800e+00,5.3900e+00,5.4000e+00,5.4100e+00,5.4200e+00,5.4300e+00,5.4400e+00,5.4500e+00,5.4600e+00,5.4700e+00,5.4800e+00,5.4900e+00,5.5000e+00,5.5100e+00,5.5200e+00,5.5300e+00,5.5400e+00,5.5500e+00,5.5600e+00,5.5700e+00,5.5800e+00,5.5900e+00,5.6000e+00,5.6100e+00,5.6200e+00,5.6300e+00,5.6400e+00,5.6500e+00,5.6600e+00,5.6700e+00,5.6800e+00,5.6900e+00,5.7000e+00,5.7100e+00,5.7200e+00,5.7300e+00,5.7400e+00,5.7500e+00,5.7600e+00,5.7700e+00,5.7800e+00,5.7900e+00,5.8000e+00,5.8100e+00,5.8200e+00,5.8300e+00,5.8400e+00,5.8500e+00,5.8600e+00,5.8700e+00,5.8800e+00,5.8900e+00,5.9000e+00,5.9100e+00,5.9200e+00,5.9300e+00,5.9400e+00,5.9500e+00,5.9600e+00,5.9700e+00,5.9800e+00,5.9900e+00,6.0000e+00,6.0100e+00,6.0200e+00,6.0300e+00,6.0400e+00,6.0500e+00,6.0600e+00,6.0700e+00,6.0800e+00,6.0900e+00,6.1000e+00,6.1100e+00,6.1200e+00,6.1300e+00,6.1400e+00,6.1500e+00,6.1600e+00,6.1700e+00,6.1800e+00,6.1900e+00,6.2000e+00,6.2100e+00,6.2200e+00,6.2300e+00,6.2400e+00,6.2500e+00,6.2600e+00,6.2700e+00,6.2800e+00,6.2900e+00,6.3000e+00,6.3100e+00,6.3200e+00,6.3300e+00,6.3400e+00,6.3500e+00,6.3600e+00,6.3700e+00,6.3800e+00,6.3900e+00,6.4000e+00,6.4100e+00,6.4200e+00,6.4300e+00,6.4400e+00,6.4500e+00,6.4600e+00,6.4700e+00,6.4800e+00,6.4900e+00,6.5000e+00,6.5100e+00,6.5200e+00,6.5300e+00,6.5400e+00,6.5500e+00,6.5600e+00,6.5700e+00,6.5800e+00,6.5900e+00,6.6000e+00,6.6100e+00,6.6200e+00,6.6300e+00,6.6400e+00,6.6500e+00,6.6600e+00,6.6700e+00,6.6800e+00,6.6900e+00,6.7000e+00,6.7100e+00,6.7200e+00,6.7300e+00,6.7400e+00,6.7500e+00,6.7600e+00,6.7700e+00,6.7800e+00,6.7900e+00,6.8000e+00,6.8100e+00,6.8200e+00,6.8300e+00,6.8400e+00,6.8500e+00,6.8600e+00,6.8700e+00,6.8800e+00,6.8900e+00,6.9000e+00,6.9100e+00,6.9200e+00,6.9300e+00,6.9400e+00,6.9500e+00,6.9600e+00,6.9700e+00,6.9800e+00,6.9900e+00,7.0000e+00,7.0100e+00,7.0200e+00,7.0300e+00,7.0400e+00,7.0500e+00,7.0600e+00,7.0700e+00,7.0800e+00,7.0900e+00,7.1000e+00,7.1100e+00,7.1200e+00,7.1300e+00,7.1400e+00,7.1500e+00,7.1600e+00,7.1700e+00,7.1800e+00,7.1900e+00,7.2000e+00,7.2100e+00,7.2200e+00,7.2300e+00,7.2400e+00,7.2500e+00,7.2600e+00,7.2700e+00,7.2800e+00,7.2900e+00,7.3000e+00,7.3100e+00,7.3200e+00,7.3300e+00,7.3400e+00,7.3500e+00,7.3600e+00,7.3700e+00,7.3800e+00,7.3900e+00,7.4000e+00,7.4100e+00,7.4200e+00,7.4300e+00,7.4400e+00,7.4500e+00,7.4600e+00,7.4700e+00,7.4800e+00,7.4900e+00,7.5000e+00,7.5100e+00,7.5200e+00,7.5300e+00,7.5400e+00,7.5500e+00,7.5600e+00,7.5700e+00,7.5800e+00,7.5900e+00,7.6000e+00,7.6100e+00,7.6200e+00,7.6300e+00,7.6400e+00,7.6500e+00,7.6600e+00,7.6700e+00,7.6800e+00,7.6900e+00,7.7000e+00,7.7100e+00,7.7200e+00,7.7300e+00,7.7400e+00,7.7500e+00,7.7600e+00,7.7700e+00,7.7800e+00,7.7900e+00,7.8000e+00,7.8100e+00,7.8200e+00,7.8300e+00,7.8400e+00,7.8500e+00,7.8600e+00,7.8700e+00,7.8800e+00,7.8900e+00,7.9000e+00,7.9100e+00,7.9200e+00,7.9300e+00,7.9400e+00,7.9500e+00,7.9600e+00,7.9700e+00,7.9800e+00,7.9900e+00,8.0000e+00,8.0100e+00,8.0200e+00,8.0300e+00,8.0400e+00,8.0500e+00,8.0600e+00,8.0700e+00,8.0800e+00,8.0900e+00,8.1000e+00,8.1100e+00,8.1200e+00,8.1300e+00,8.1400e+00,8.1500e+00,8.1600e+00,8.1700e+00,8.1800e+00,8.1900e+00,8.2000e+00,8.2100e+00,8.2200e+00,8.2300e+00,8.2400e+00,8.2500e+00,8.2600e+00,8.2700e+00,8.2800e+00,8.2900e+00,8.3000e+00,8.3100e+00,8.3200e+00,8.3300e+00,8.3400e+00,8.3500e+00,8.3600e+00,8.3700e+00,8.3800e+00,8.3900e+00,8.4000e+00,8.4100e+00,8.4200e+00,8.4300e+00,8.4400e+00,8.4500e+00,8.4600e+00,8.4700e+00,8.4800e+00,8.4900e+00,8.5000e+00,8.5100e+00,8.5200e+00,8.5300e+00,8.5400e+00,8.5500e+00,8.5600e+00,8.5700e+00,8.5800e+00,8.5900e+00,8.6000e+00,8.6100e+00,8.6200e+00,8.6300e+00,8.6400e+00,8.6500e+00,8.6600e+00,8.6700e+00,8.6800e+00,8.6900e+00,8.7000e+00,8.7100e+00,8.7200e+00,8.7300e+00,8.7400e+00,8.7500e+00,8.7600e+00,8.7700e+00,8.7800e+00,8.7900e+00,8.8000e+00,8.8100e+00,8.8200e+00,8.8300e+00,8.8400e+00,8.8500e+00,8.8600e+00,8.8700e+00,8.8800e+00,8.8900e+00,8.9000e+00,8.9100e+00,8.9200e+00,8.9300e+00,8.9400e+00,8.9500e+00,8.9600e+00,8.9700e+00,8.9800e+00,8.9900e+00,9.0000e+00,9.0100e+00,9.0200e+00,9.0300e+00,9.0400e+00,9.0500e+00,9.0600e+00,9.0700e+00,9.0800e+00,9.0900e+00,9.1000e+00,9.1100e+00,9.1200e+00,9.1300e+00,9.1400e+00,9.1500e+00,9.1600e+00,9.1700e+00,9.1800e+00,9.1900e+00,9.2000e+00,9.2100e+00,9.2200e+00,9.2300e+00,9.2400e+00,9.2500e+00,9.2600e+00,9.2700e+00,9.2800e+00,9.2900e+00,9.3000e+00,9.3100e+00,9.3200e+00,9.3300e+00,9.3400e+00,9.3500e+00,9.3600e+00,9.3700e+00,9.3800e+00,9.3900e+00,9.4000e+00,9.4100e+00,9.4200e+00,9.4300e+00,9.4400e+00,9.4500e+00,9.4600e+00,9.4700e+00,9.4800e+00,9.4900e+00,9.5000e+00,9.5100e+00,9.5200e+00,9.5300e+00,9.5400e+00,9.5500e+00,9.5600e+00,9.5700e+00,9.5800e+00,9.5900e+00,9.6000e+00,9.6100e+00,9.6200e+00,9.6300e+00,9.6400e+00,9.6500e+00,9.6600e+00,9.6700e+00,9.6800e+00,9.6900e+00,9.7000e+00,9.7100e+00,9.7200e+00,9.7300e+00,9.7400e+00,9.7500e+00,9.7600e+00,9.7700e+00,9.7800e+00,9.7900e+00,9.8000e+00,9.8100e+00,9.8200e+00,9.8300e+00,9.8400e+00,9.8500e+00,9.8600e+00,9.8700e+00,9.8800e+00,9.8900e+00,9.9000e+00,9.9100e+00,9.9200e+00,9.9300e+00,9.9400e+00,9.9500e+00,9.9600e+00,9.9700e+00,9.9800e+00,9.9900e+00,1.0000e+01])
	
	V = np.array([0.,1.9076e-02,1.9076e-02,2.8732e-02,3.8534e-02,4.8537e-02,5.8802e-02,6.9399e-02,8.0413e-02,9.1941e-02,1.0411e-01,1.1707e-01,1.3103e-01,1.4626e-01,1.6315e-01,1.8230e-01,2.0462e-01,2.3175e-01,2.6706e-01,3.1976e-01,4.9210e-01,6.3111e-01,6.7945e-01,7.1334e-01,7.3961e-01,7.6098e-01,7.7890e-01,7.9425e-01,8.0759e-01,8.1932e-01,8.2973e-01,8.3904e-01,8.4743e-01,8.5503e-01,8.6196e-01,8.6829e-01,8.7409e-01,8.7945e-01,8.8441e-01,8.8900e-01,8.9327e-01,8.9725e-01,9.0097e-01,9.0445e-01,9.0772e-01,9.1079e-01,9.1368e-01,9.1641e-01,9.1899e-01,9.2142e-01,9.2373e-01,9.2592e-01,9.2800e-01,9.2998e-01,9.3187e-01,9.3367e-01,9.3537e-01,9.3701e-01,9.3857e-01,9.4006e-01,9.4150e-01,9.4288e-01,9.4419e-01,9.4546e-01,9.4667e-01,9.4783e-01,9.4896e-01,9.5004e-01,9.5108e-01,9.5208e-01,9.5305e-01,9.5397e-01,9.5488e-01,9.5575e-01,9.5659e-01,9.5740e-01,9.5819e-01,9.5895e-01,9.5969e-01,9.6040e-01,9.6109e-01,9.6177e-01,9.6241e-01,9.6305e-01,9.6366e-01,9.6426e-01,9.6484e-01,9.6539e-01,9.6594e-01,9.6648e-01,9.6699e-01,9.6750e-01,9.6799e-01,9.6846e-01,9.6892e-01,9.6938e-01,9.6982e-01,9.7025e-01,9.7067e-01,9.7108e-01,9.7148e-01,9.7187e-01,9.7225e-01,9.7262e-01,9.7299e-01,9.7333e-01,9.7368e-01,9.7403e-01,9.7435e-01,9.7468e-01,9.7499e-01,9.7530e-01,9.7561e-01,9.7590e-01,9.7619e-01,9.7647e-01,9.7676e-01,9.7703e-01,9.7730e-01,9.7756e-01,9.7781e-01,9.7807e-01,9.7831e-01,9.7855e-01,9.7879e-01,9.7902e-01,9.7925e-01,9.7947e-01,9.7969e-01,9.7991e-01,9.8012e-01,9.8033e-01,9.8053e-01,9.8073e-01,9.8092e-01,9.8112e-01,9.8131e-01,9.8149e-01,9.8168e-01,9.8186e-01,9.8203e-01,9.8222e-01,9.8238e-01,9.8255e-01,9.8272e-01,9.8288e-01,9.8303e-01,9.8320e-01,9.8335e-01,9.8350e-01,9.8365e-01,9.8380e-01,9.8394e-01,9.8409e-01,9.8423e-01,9.8437e-01,9.8450e-01,9.8463e-01,9.8477e-01,9.8490e-01,9.8503e-01,9.8515e-01,9.8528e-01,9.8540e-01,9.8552e-01,9.8564e-01,9.8576e-01,9.8588e-01,9.8599e-01,9.8611e-01,9.8621e-01,9.8633e-01,9.8643e-01,9.8654e-01,9.8664e-01,9.8675e-01,9.8685e-01,9.8695e-01,9.8705e-01,9.8715e-01,9.8725e-01,9.8734e-01,9.8744e-01,9.8753e-01,9.8762e-01,9.8771e-01,9.8781e-01,9.8789e-01,9.8798e-01,9.8806e-01,9.8815e-01,9.8823e-01,9.8832e-01,9.8840e-01,9.8848e-01,9.8856e-01,9.8864e-01,9.8872e-01,9.8880e-01,9.8887e-01,9.8895e-01,9.8902e-01,9.8909e-01,9.8917e-01,9.8924e-01,9.8931e-01,9.8938e-01,9.8945e-01,9.8952e-01,9.8959e-01,9.8966e-01,9.8973e-01,9.8979e-01,9.8985e-01,9.8992e-01,9.8999e-01,9.9004e-01,9.9011e-01,9.9016e-01,9.9023e-01,9.9029e-01,9.9035e-01,9.9041e-01,9.9046e-01,9.9053e-01,9.9058e-01,9.9063e-01,9.9069e-01,9.9075e-01,9.9080e-01,9.9086e-01,9.9091e-01,9.9096e-01,9.9102e-01,9.9107e-01,9.9112e-01,9.9117e-01,9.9122e-01,9.9127e-01,9.9133e-01,9.9137e-01,9.9142e-01,9.9147e-01,9.9152e-01,9.9156e-01,9.9161e-01,9.9166e-01,9.9170e-01,9.9174e-01,9.9179e-01,9.9183e-01,9.9188e-01,9.9192e-01,9.9197e-01,9.9200e-01,9.9205e-01,9.9209e-01,9.9214e-01,9.9217e-01,9.9222e-01,9.9226e-01,9.9229e-01,9.9234e-01,9.9238e-01,9.9241e-01,9.9246e-01,9.9249e-01,9.9253e-01,9.9256e-01,9.9261e-01,9.9264e-01,9.9268e-01,9.9271e-01,9.9275e-01,9.9279e-01,9.9283e-01,9.9286e-01,9.9289e-01,9.9293e-01,9.9296e-01,9.9299e-01,9.9303e-01,9.9306e-01,9.9310e-01,9.9313e-01,9.9317e-01,9.9320e-01,9.9323e-01,9.9326e-01,9.9329e-01,9.9332e-01,9.9335e-01,9.9338e-01,9.9341e-01,9.9344e-01,9.9347e-01,9.9350e-01,9.9354e-01,9.9357e-01,9.9360e-01,9.9362e-01,9.9365e-01,9.9368e-01,9.9371e-01,9.9374e-01,9.9377e-01,9.9379e-01,9.9382e-01,9.9384e-01,9.9387e-01,9.9390e-01,9.9393e-01,9.9395e-01,9.9398e-01,9.9401e-01,9.9403e-01,9.9406e-01,9.9408e-01,9.9411e-01,9.9413e-01,9.9415e-01,9.9418e-01,9.9421e-01,9.9423e-01,9.9425e-01,9.9428e-01,9.9431e-01,9.9432e-01,9.9435e-01,9.9438e-01,9.9439e-01,9.9442e-01,9.9444e-01,9.9446e-01,9.9449e-01,9.9451e-01,9.9453e-01,9.9455e-01,9.9458e-01,9.9460e-01,9.9462e-01,9.9464e-01,9.9466e-01,9.9468e-01,9.9470e-01,9.9472e-01,9.9474e-01,9.9477e-01,9.9478e-01,9.9480e-01,9.9483e-01,9.9485e-01,9.9487e-01,9.9489e-01,9.9491e-01,9.9493e-01,9.9494e-01,9.9496e-01,9.9498e-01,9.9500e-01,9.9502e-01,9.9504e-01,9.9505e-01,9.9508e-01,9.9509e-01,9.9511e-01,9.9513e-01,9.9515e-01,9.9517e-01,9.9518e-01,9.9520e-01,9.9522e-01,9.9524e-01,9.9525e-01,9.9527e-01,9.9529e-01,9.9531e-01,9.9532e-01,9.9534e-01,9.9535e-01,9.9537e-01,9.9539e-01,9.9541e-01,9.9542e-01,9.9544e-01,9.9545e-01,9.9547e-01,9.9548e-01,9.9550e-01,9.9551e-01,9.9553e-01,9.9554e-01,9.9556e-01,9.9558e-01,9.9559e-01,9.9561e-01,9.9562e-01,9.9564e-01,9.9565e-01,9.9567e-01,9.9569e-01,9.9570e-01,9.9572e-01,9.9573e-01,9.9574e-01,9.9576e-01,9.9577e-01,9.9578e-01,9.9580e-01,9.9581e-01,9.9582e-01,9.9584e-01,9.9585e-01,9.9587e-01,9.9588e-01,9.9589e-01,9.9591e-01,9.9592e-01,9.9593e-01,9.9595e-01,9.9596e-01,9.9598e-01,9.9599e-01,9.9600e-01,9.9601e-01,9.9602e-01,9.9604e-01,9.9605e-01,9.9606e-01,9.9607e-01,9.9609e-01,9.9610e-01,9.9611e-01,9.9612e-01,9.9613e-01,9.9615e-01,9.9616e-01,9.9617e-01,9.9618e-01,9.9619e-01,9.9620e-01,9.9622e-01,9.9623e-01,9.9625e-01,9.9625e-01,9.9627e-01,9.9628e-01,9.9629e-01,9.9630e-01,9.9631e-01,9.9632e-01,9.9633e-01,9.9635e-01,9.9635e-01,9.9637e-01,9.9637e-01,9.9639e-01,9.9640e-01,9.9641e-01,9.9642e-01,9.9643e-01,9.9644e-01,9.9645e-01,9.9646e-01,9.9647e-01,9.9648e-01,9.9649e-01,9.9650e-01,9.9652e-01,9.9653e-01,9.9653e-01,9.9655e-01,9.9655e-01,9.9657e-01,9.9657e-01,9.9658e-01,9.9659e-01,9.9660e-01,9.9661e-01,9.9662e-01,9.9663e-01,9.9664e-01,9.9665e-01,9.9666e-01,9.9667e-01,9.9668e-01,9.9669e-01,9.9669e-01,9.9670e-01,9.9671e-01,9.9672e-01,9.9673e-01,9.9674e-01,9.9675e-01,9.9676e-01,9.9677e-01,9.9677e-01,9.9679e-01,9.9680e-01,9.9681e-01,9.9682e-01,9.9682e-01,9.9683e-01,9.9684e-01,9.9685e-01,9.9686e-01,9.9687e-01,9.9687e-01,9.9688e-01,9.9689e-01,9.9690e-01,9.9691e-01,9.9691e-01,9.9692e-01,9.9693e-01,9.9694e-01,9.9695e-01,9.9695e-01,9.9696e-01,9.9697e-01,9.9698e-01,9.9699e-01,9.9699e-01,9.9700e-01,9.9701e-01,9.9702e-01,9.9702e-01,9.9703e-01,9.9704e-01,9.9705e-01,9.9705e-01,9.9707e-01,9.9707e-01,9.9708e-01,9.9709e-01,9.9710e-01,9.9710e-01,9.9711e-01,9.9712e-01,9.9713e-01,9.9713e-01,9.9714e-01,9.9715e-01,9.9715e-01,9.9716e-01,9.9717e-01,9.9717e-01,9.9718e-01,9.9719e-01,9.9720e-01,9.9720e-01,9.9721e-01,9.9721e-01,9.9722e-01,9.9723e-01,9.9724e-01,9.9724e-01,9.9725e-01,9.9725e-01,9.9726e-01,9.9727e-01,9.9727e-01,9.9728e-01,9.9729e-01,9.9729e-01,9.9730e-01,9.9731e-01,9.9731e-01,9.9732e-01,9.9732e-01,9.9733e-01,9.9734e-01,9.9735e-01,9.9736e-01,9.9736e-01,9.9737e-01,9.9738e-01,9.9738e-01,9.9739e-01,9.9739e-01,9.9740e-01,9.9741e-01,9.9741e-01,9.9742e-01,9.9742e-01,9.9743e-01,9.9743e-01,9.9744e-01,9.9745e-01,9.9745e-01,9.9746e-01,9.9746e-01,9.9747e-01,9.9747e-01,9.9748e-01,9.9749e-01,9.9749e-01,9.9750e-01,9.9750e-01,9.9751e-01,9.9751e-01,9.9752e-01,9.9752e-01,9.9753e-01,9.9754e-01,9.9754e-01,9.9755e-01,9.9755e-01,9.9756e-01,9.9756e-01,9.9757e-01,9.9757e-01,9.9758e-01,9.9759e-01,9.9759e-01,9.9760e-01,9.9760e-01,9.9760e-01,9.9761e-01,9.9761e-01,9.9763e-01,9.9763e-01,9.9764e-01,9.9764e-01,9.9765e-01,9.9765e-01,9.9766e-01,9.9766e-01,9.9767e-01,9.9767e-01,9.9768e-01,9.9768e-01,9.9769e-01,9.9769e-01,9.9770e-01,9.9770e-01,9.9771e-01,9.9771e-01,9.9772e-01,9.9772e-01,9.9773e-01,9.9773e-01,9.9773e-01,9.9774e-01,9.9774e-01,9.9775e-01,9.9775e-01,9.9776e-01,9.9776e-01,9.9777e-01,9.9777e-01,9.9778e-01,9.9778e-01,9.9779e-01,9.9779e-01,9.9779e-01,9.9780e-01,9.9780e-01,9.9781e-01,9.9781e-01,9.9782e-01,9.9782e-01,9.9783e-01,9.9783e-01,9.9783e-01,9.9784e-01,9.9784e-01,9.9785e-01,9.9785e-01,9.9786e-01,9.9786e-01,9.9787e-01,9.9787e-01,9.9787e-01,9.9788e-01,9.9788e-01,9.9789e-01,9.9789e-01,9.9789e-01,9.9791e-01,9.9791e-01,9.9791e-01,9.9792e-01,9.9792e-01,9.9793e-01,9.9793e-01,9.9793e-01,9.9794e-01,9.9794e-01,9.9795e-01,9.9795e-01,9.9796e-01,9.9796e-01,9.9796e-01,9.9797e-01,9.9797e-01,9.9797e-01,9.9798e-01,9.9798e-01,9.9799e-01,9.9799e-01,9.9799e-01,9.9800e-01,9.9800e-01,9.9801e-01,9.9801e-01,9.9801e-01,9.9802e-01,9.9802e-01,9.9802e-01,9.9803e-01,9.9803e-01,9.9804e-01,9.9804e-01,9.9804e-01,9.9805e-01,9.9805e-01,9.9805e-01,9.9806e-01,9.9806e-01,9.9806e-01,9.9807e-01,9.9807e-01,9.9807e-01,9.9808e-01,9.9808e-01,9.9809e-01,9.9809e-01,9.9809e-01,9.9810e-01,9.9810e-01,9.9810e-01,9.9811e-01,9.9811e-01,9.9811e-01,9.9812e-01,9.9812e-01,9.9812e-01,9.9813e-01,9.9813e-01,9.9813e-01,9.9814e-01,9.9814e-01,9.9814e-01,9.9815e-01,9.9815e-01,9.9815e-01,9.9816e-01,9.9816e-01,9.9816e-01,9.9817e-01,9.9817e-01,9.9817e-01,9.9818e-01,9.9819e-01,9.9819e-01,9.9819e-01,9.9820e-01,9.9820e-01,9.9820e-01,9.9821e-01,9.9821e-01,9.9821e-01,9.9822e-01,9.9822e-01,9.9822e-01,9.9823e-01,9.9823e-01,9.9823e-01,9.9823e-01,9.9824e-01,9.9824e-01,9.9824e-01,9.9825e-01,9.9825e-01,9.9825e-01,9.9826e-01,9.9826e-01,9.9826e-01,9.9826e-01,9.9827e-01,9.9827e-01,9.9827e-01,9.9828e-01,9.9828e-01,9.9828e-01,9.9829e-01,9.9829e-01,9.9829e-01,9.9830e-01,9.9830e-01,9.9830e-01,9.9830e-01,9.9831e-01,9.9831e-01,9.9831e-01,9.9832e-01,9.9832e-01,9.9832e-01,9.9832e-01,9.9833e-01,9.9833e-01,9.9833e-01,9.9833e-01,9.9834e-01,9.9834e-01,9.9834e-01,9.9835e-01,9.9835e-01,9.9835e-01,9.9835e-01,9.9836e-01,9.9836e-01,9.9836e-01,9.9836e-01,9.9837e-01,9.9837e-01,9.9837e-01,9.9838e-01,9.9838e-01,9.9838e-01,9.9838e-01,9.9839e-01,9.9839e-01,9.9839e-01,9.9839e-01,9.9840e-01,9.9840e-01,9.9840e-01,9.9840e-01,9.9841e-01,9.9841e-01,9.9841e-01,9.9842e-01,9.9842e-01,9.9842e-01,9.9842e-01,9.9842e-01,9.9843e-01,9.9843e-01,9.9843e-01,9.9843e-01,9.9844e-01,9.9844e-01,9.9844e-01,9.9845e-01,9.9845e-01,9.9845e-01,9.9845e-01,9.9845e-01,9.9847e-01,9.9847e-01,9.9847e-01,9.9847e-01,9.9847e-01,9.9848e-01,9.9848e-01,9.9848e-01,9.9849e-01,9.9849e-01,9.9849e-01,9.9849e-01,9.9849e-01,9.9850e-01,9.9850e-01,9.9850e-01,9.9850e-01,9.9851e-01,9.9851e-01,9.9851e-01,9.9851e-01,9.9851e-01,9.9852e-01,9.9852e-01,9.9852e-01,9.9852e-01,9.9853e-01,9.9853e-01,9.9853e-01,9.9853e-01,9.9854e-01,9.9854e-01,9.9854e-01,9.9854e-01,9.9854e-01,9.9855e-01,9.9855e-01,9.9855e-01,9.9855e-01,9.9856e-01,9.9856e-01,9.9856e-01,9.9856e-01,9.9856e-01,9.9857e-01,9.9857e-01,9.9857e-01,9.9857e-01,9.9857e-01,9.9858e-01,9.9858e-01,9.9858e-01,9.9858e-01,9.9859e-01,9.9859e-01,9.9859e-01,9.9859e-01,9.9859e-01,9.9860e-01,9.9860e-01,9.9860e-01,9.9860e-01,9.9860e-01,9.9861e-01,9.9861e-01,9.9861e-01,9.9861e-01,9.9861e-01,9.9862e-01,9.9862e-01,9.9862e-01,9.9862e-01,9.9863e-01,9.9863e-01,9.9863e-01,9.9863e-01,9.9863e-01,9.9863e-01,9.9864e-01,9.9864e-01,9.9864e-01,9.9864e-01,9.9864e-01,9.9865e-01,9.9865e-01,9.9865e-01,9.9865e-01,9.9865e-01,9.9866e-01,9.9866e-01,9.9866e-01,9.9866e-01,9.9866e-01,9.9867e-01,9.9867e-01,9.9867e-01,9.9867e-01,9.9867e-01,9.9867e-01,9.9868e-01,9.9868e-01,9.9868e-01,9.9868e-01,9.9869e-01,9.9869e-01,9.9869e-01,9.9869e-01,9.9869e-01,9.9869e-01,9.9870e-01,9.9870e-01,9.9870e-01,9.9870e-01,9.9870e-01,9.9870e-01,9.9871e-01,9.9871e-01,9.9871e-01,9.9871e-01,9.9871e-01,9.9871e-01,9.9872e-01,9.9872e-01,9.9872e-01,9.9872e-01,9.9873e-01,9.9873e-01,9.9873e-01,9.9873e-01,9.9873e-01,9.9873e-01,9.9874e-01,9.9875e-01,9.9875e-01,9.9875e-01,9.9875e-01,9.9875e-01,9.9875e-01,9.9876e-01,9.9876e-01,9.9876e-01,9.9876e-01,9.9876e-01,9.9876e-01,9.9877e-01,9.9877e-01,9.9877e-01,9.9877e-01,9.9877e-01,9.9877e-01,9.9878e-01])
	
	iZ=0
	for iz in range(0,len(Z)-1):
		if Z[iz+1] >r:
			iZ=iz
			break

	V_interp = V[iZ] + (V[iZ+1]-V[iZ])*(r-Z[iZ])/(Z[iZ+1]-Z[iZ])
	return V_interp

v = np.zeros(znum)
for ii in range(0,znum):
	v[ii] = wind_total_data_1(z[ii])

v_a = v*800.*1.e5
v_b = v*1000.*1.e5
v_c = v*1000.*1.e5

###
### Interstellar Radiation Field
###

# isrf_file = 'ISRF_m82_v1_flux_z.txt'

# data = np.loadtxt(isrf_file)

# z_isrf = data[0]
# isrf_dat = data[1]

# grad_isrf_z = z_isrf[:-1]+(z_isrf[:-1]-z_isrf[1:])/2.
# grad_isrf = (isrf_dat[1:]-isrf_dat[:-1])/(z_isrf[1:]-z_isrf[:-1])/kpc2cm 


# grad_isrf_z = np.interp(z, grad_isrf_z, grad_isrf)
# isrf_accel = grad_isrf_z/(mp*(nHI+nHII))
# isrf_accel2 = grad_isrf_z/(mp*(nHI2+nHII2))

isrf_file = '../m82_dist_py/ISRF_m82_v4.fits.gz'

#READING EMISSION FITS-------------------------------------------------
hdulist=py.open(isrf_file)
head=hdulist[0].header
scidata=(hdulist[0].data)

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

Zi = np.linspace(0, znum-1, znum)
Zi *= zdel
Zi += zmin

emin = 10.**log_e_min  #MeV
efactor = 10.**log_efactor

xmax = xmin+xdel*(xnum-1)
zmax = zmin+zdel*(znum-1)

xsize = xmax-xmin
zsize = zmax-zmin

R = xmin+np.array(range(0,xnum))*xdel
energy = emin*efactor**np.array(range(0,enum))

isrf = np.sum(scidata, axis=0)
isrf = np.sum(isrf, axis=0)*log_efactor
isrf = np.swapaxes(isrf, 0,1)
U_isrf_z = isrf[0]

U_isrf = np.interp(z, Zi, U_isrf_z)/1000.
# print U_isrf
###
### [z, U_isrf] [grad_isrf_z, grad_isrf]
###########################################


##
## Diffusion coefficient
##
D_xx_1 = 4.e27
D_xx_2 = 1.e28

###
### PLOT
###

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

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

MODE = 3

E=1.
for m in range(0,len(z)):
	nu = np.logspace(-3., 4., num = 100)
	m=len(z)-1-m
	
	tau_sy_1 = 7.85e5* nu**-0.5* B_1[m]**-1.5
	tau_IC_1 = 1.90e5* nu**-0.5* B_1[m]**0.5*  U_isrf[m]**-1.0
	tau_br_1 = 3.18e5* nu**0* (nHI_1[m]+nHII_1[m])**-1.0
	tau_i0_1 = 3.67e6* nu**0.5*  B_1[m]**-0.5* (nHI_1[m])**-1.0
	tau_i1_1 = 9.04e5* nu**0.5*  B_1[m]**-0.5* (nHII_1[m])**-1.0
	tau_io_1 = (1./tau_i0_1+1./tau_i1_1)**-1.
	tau_total_1 = (1./tau_sy_1 +1./tau_IC_1 +1./tau_br_1 + 1./tau_io_1)**-1.*3.15e7
	# DD_1 = np.sqrt(D_xx_1*((nu/B_1[m])**0.5/4.)**0.31*tau_total_1) +v_1[m]*tau_total_1
	# DD_1_2 = -np.sqrt(D_xx_1*((nu/B_1[m])**0.5/4.)**0.31*tau_total_1) +v_1[m]*tau_total_1
	
	tau_re_1 = np.ones(len(nu))*9./4.*D_xx_1/vA_a[m]**2
	tau_wi_1 = np.ones(len(nu))*0.5*3.18e21/v_a[m]
	
	tau_pi_1 = 8.2e5*E**-0.28*(E+200.)**0.2*((nHI_1[m]+nHII_1[m]))**-1.0
	
	tau_sy_2 = 7.85e5* nu**-0.5* B_2[m]**-1.5
	tau_IC_2 = 1.90e5* nu**-0.5* B_2[m]**0.5*  U_isrf[m]**-1.0
	tau_br_2 = 3.18e5* nu**0* (nHI_2[m]+nHII_2[m])**-1.0
	tau_i0_2 = 3.67e6* nu**0.5*  B_2[m]**-0.5* (nHI_2[m])**-1.0
	tau_i1_2 = 9.04e5* nu**0.5*  B_2[m]**-0.5* (nHII_2[m])**-1.0
	tau_io_2 = (1./tau_i0_2+1./tau_i1_2)**-1.
	tau_total_2 = (1./tau_sy_2 +1./tau_IC_2 +1./tau_br_2 + 1./tau_io_2)**-1.*3.15e7
	# DD_2 = np.sqrt(D_xx_2*((nu/B_2[m])**0.5/4.)**0.31*tau_total_2) +v_2[m]*tau_total_2
	# DD_2_2 = -np.sqrt(D_xx_2*((nu/B_2[m])**0.5/4.)**0.31*tau_total_2) +v_2[m]*tau_total_2
	
	tau_re_2 = np.ones(len(nu))*9./4.*D_xx_2/vA_b[m]**2
	tau_wi_2 = np.ones(len(nu))*0.5*3.18e21/v_b[m]
	
	tau_pi_2 = 8.2e5*E**-0.28*(E+200.)**0.2*((nHI_2[m]+nHII_2[m]))**-1.0
	
	tau_sy_3 = 7.85e5* nu**-0.5* B_3[m]**-1.5
	tau_IC_3 = 1.90e5* nu**-0.5* B_3[m]**0.5*  U_isrf[m]**-1.0
	tau_br_3 = 3.18e5* nu**0* (nHI_3[m]+nHII_3[m])**-1.0
	tau_i0_3 = 3.67e6* nu**0.5*  B_3[m]**-0.5* (nHI_3[m])**-1.0
	tau_i1_3 = 9.04e5* nu**0.5*  B_3[m]**-0.5* (nHII_3[m])**-1.0
	tau_io_3 = (1./tau_i0_3+1./tau_i1_3)**-1.
	tau_total_3 = (1./tau_sy_3 +1./tau_IC_3 +1./tau_br_3 + 1./tau_io_3)**-1.*3.15e7
	# DD_2 = np.sqrt(D_xx_2*((nu/B_2[m])**0.5/4.)**0.31*tau_total_2) +v_2[m]*tau_total_2
	# DD_2_2 = -np.sqrt(D_xx_2*((nu/B_2[m])**0.5/4.)**0.31*tau_total_2) +v_2[m]*tau_total_2
	
	tau_re_3 = np.ones(len(nu))*9./4.*D_xx_2/vA_c[m]**2
	tau_wi_3 = np.ones(len(nu))*0.5*3.18e21/v_c[m]
	
	tau_pi_3 = 8.2e5*E**-0.28*(E+200.)**0.2*((nHI_3[m]+nHII_3[m]))**-1.0
	
	# print DD_1, DD_1_2
	color = cmap((float(m+1))/(len(z)+2))
	nu*=1.e9
	
	if MODE == 1:
		# ax.fill_between(nu, DD_1/kpc2cm, DD_1_2/kpc2cm, facecolor=color, edgecolor=color, alpha=0.65, label=LABEL[m])
		ax.plot(nu, tau_total_1/3.15e7, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
		ax.plot(nu, tau_re_1/3.15e7, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		ax.plot(nu, tau_wi_1/3.15e7, c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
		ax.scatter(1.e11,tau_pi_1, color=color, marker='*', alpha=ALPHA, s=100)
		# ax.plot(nu, DD_1/kpc2cm, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_IC_1, c=color, linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_br_1, c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_io_1, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_i1_1, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	elif MODE == 2:
		# ax.fill_between(nu, DD_2/kpc2cm, DD_2_2/kpc2cm, facecolor=color, edgecolor=color, alpha=0.65, label=LABEL[m])
		ax.plot(nu, tau_total_2/3.15e7, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label=LABEL[m])
		ax.plot(nu, tau_re_2/3.15e7, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		ax.plot(nu, tau_wi_2/3.15e7, c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
		ax.scatter(1.e11,tau_pi_2, color=color, marker='*', alpha=ALPHA, s=100)
		# ax.plot(nu, DD_2/kpc2cm, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_IC_2, c=color, linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_br_2, c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_io_2, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_i1_2, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	elif MODE == 3:
		# ax.fill_between(nu, DD_2/kpc2cm, DD_2_2/kpc2cm, facecolor=color, edgecolor=color, alpha=0.65, label=LABEL[m])
		ax.plot(nu, tau_total_3/3.15e7, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
		ax.plot(nu, tau_re_3/3.15e7, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		ax.plot(nu, tau_wi_3/3.15e7, c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
		ax.scatter(1.e11,tau_pi_3, color=color, marker='*', alpha=ALPHA, s=100)
		# ax.plot(nu, DD_2/kpc2cm, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_IC_2, c=color, linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_br_2, c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_io_2, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(nu, tau_i1_2, c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	elif MODE == 4:
		ax.plot(nu, tau_total_1/3.15e7, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
		ax.plot(nu, tau_total_2/3.15e7, c=color, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
		# ax.fill_between(nu, DD_1/kpc2cm, DD_1_2/kpc2cm, facecolor=color, edgecolor=color, alpha=0.85, label=LABEL[m])
		# ax.fill_between(nu, DD_2/kpc2cm, DD_2_2/kpc2cm, facecolor=color, edgecolor=color, alpha=0.35, label=None)



ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([8.e+6,1.e+13])
# ax.set_ylim([1.e-2,1.e2])
ax.set_xlabel(r'$\nu$ [Hz]')
ax.set_ylabel(r'$\tau$ [yr]')


if MODE == 1:
	ax.text(0.333,0.333,r'Model A', ha='center', va='center', color='k', transform=ax.transAxes, backgroundcolor='white')
	plot_name = 'plot_timescales_a_'
elif MODE == 2:
	ax.text(0.333,0.333,r'Model B', ha='center', va='center', color='k', transform=ax.transAxes, backgroundcolor='white')
	ax.legend(loc='upper right', prop={'size':14})
	plot_name = 'plot_timescales_b_'
elif MODE == 3:
	ax.text(0.333,0.333,r"Model B'", ha='center', va='center', color='k', transform=ax.transAxes, backgroundcolor='white')
	plot_name = 'plot_timescales_bp_'

# color_a = cmap(1.-1./4.)
# color_b = cmap(1.-2./4)
# color_c = cmap(1.-3./4)

# colorI = cmap(0.1)
# colorII = cmap(0.85)


# ax.plot(Z,v_a, c=color_a, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label='A')
# ax.plot(Z,v_b, c=color_b, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label='B')
# ax.plot(Z,v_c, c=color_c, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label="B'")

# ax.plot(z,vA_a, c=color_a, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,vA_b, c=color_b, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,vA_c, c=color_c, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)

# ax.plot(z,nHI_a5, c=color_a, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI_b5, c=color_b, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)

# ax.plot(z,nHI_ta, c=color_a, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI_tb, c=color_b, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)

# # ax.plot(z,nHI10, c=colorI, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
# # ax.plot(z,nHI11, c=colorI, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI20, c=color_a, linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI21, c=color_b, linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)

# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlim([0.,3.5])
# ax.set_ylim([0.01,2000.])
# ax.set_xlabel(r'$z$ [kpc]')
# ax.set_ylabel(r'$V_0$ [km s$^{-1}$]')

ax.minorticks_on()

# ax.legend(loc='upper right', prop={'size':15})

# ##Colored y-axis label!! >>
# ybox1 = TextArea(r'$U_B$', textprops=dict(color=cmap(0.1),ha='center',va='bottom',rotation=90))
# ybox2 = TextArea(' ,', textprops=dict(color='k',ha='center',va='bottom',rotation=90))
# ybox3 = TextArea(r'$U_{\mathrm{ISRF}}$', textprops=dict(color=cmap(0.4),ha='center',va='bottom',rotation=90))
# ybox4 = TextArea(r'  [eV cm$^{-3}$]', textprops=dict(color='k',ha='center',va='bottom',rotation=90))

# ybox = VPacker(children=[ybox4, ybox3, ybox2, ybox1],align="center", pad=0, sep=-3)

# anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0, frameon=False, bbox_to_anchor=(-0.114, 0.18), bbox_transform=ax.transAxes, borderpad=0.)
# ax.add_artist(anchored_ybox)
# ## <<

# ##Colored y-axis label!! >>
# ybox1 = TextArea(r'$n_{\mathrm{HI}}$', textprops=dict(color=cmap(0.7),ha='center',va='bottom',rotation=270))
# ybox2 = TextArea(r'  [cm$^{-3}$]', textprops=dict(color='k',ha='center',va='bottom',rotation=270))

# ybox = VPacker(children=[ybox1, ybox2],align="center", pad=0, sep=-4)

# anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0, frameon=False, bbox_to_anchor=(1.184, 0.32), bbox_transform=ax.transAxes, borderpad=0.)
# ax.add_artist(anchored_ybox)
# ## <<

# ax.set_title(r'Model Distributions')

plt.tight_layout()
# plot_name = 'plot_timescales_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break
plt.close()




















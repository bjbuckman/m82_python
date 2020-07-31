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
import scipy

import imp
magnetic_field = imp.load_source("magnetic_field_v4", "../m82_dist_py/magnetic_field_v4.py")

#constants
pc2cm = 3.08568025e18 #cm/pc
kpc2cm = 3.08568025e21 #cm/kpc
radian2arcsec = 648000./math.pi #acsec/radian
deg2arcsec = 60*60 #arcsec/deg
Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
c = 3.0e10 #cm/2
h = 6.62e-27 #erg s
q=4.8e-10 #esu
m=9.11e-28 #g

R=np.arange(0,4,0.01)
Z=np.arange(0,4,0.01)

Rnum = len(R)
Znum = len(Z)

RR, ZZ = np.meshgrid(R,Z)

###
### MAGNETIC FIELD 1
###

magnetic_model_name = 'constant+powerlaw_v3'

B0 = 150.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.0

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

B_1 = np.zeros([Rnum,Znum])
for iR in range(0,Rnum):
	for iZ in range(0,Znum):
		B_1[iR][iZ] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([R[iR],Z[iZ]]), B_parameters)
		
B_1 = B_1.T/1.e-6
# print B_1[:10,:10]

###
### [z, B] [z, grad_magnetic]
##############################################


######################################
## PLOT DATA
##

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

cmap = plt.get_cmap('cubehelix')
fig = plt.figure()
	
f1, ax = plt.subplots(1,1)

for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(2)
ax.tick_params(which='major',width=1., length=8, labelsize='small', direction='in', color='white')
ax.tick_params(which='minor',width=0.5, length=5, direction='in', color='white')

ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
ax.yaxis.set_minor_locator(plt.MultipleLocator(10))

MSIZE = 5
LWIDTH = 3
WIDTH = 3
ALPHA = 0.7

## PLOT 2
loglow = -1.
loghigh = 3.
V = np.logspace(loglow, loghigh, 17)

loglow = -1.
loghigh = 3.
V1 = np.logspace(loglow, loghigh, 17)

CF0 = ax.contourf(RR,ZZ,B_1,norm = LogNorm(), levels=V, cmap=cmap)


#colorbar
cf1 = plt.colorbar(CF0, ax=ax, ticks=V1, pad=0.0025, fraction=0.04763)
cf1.ax.tick_params(labelsize='small')
cf1.ax.yaxis.set_tick_params(color='w', width=1, length=5)
cf1.outline.set_linewidth(2) 
cf1.set_label(r'$\mu$G', rotation=270, horizontalalignment='center', verticalalignment='bottom')

#axes stuff
ax.set_axis_bgcolor('black')
ax.set_aspect('equal')
ax.set_xlim([0.01,3.99])
ax.set_ylim([0.01,3.99])
# ax.set_xticks([-100, -50, 0, 50, 100])
# ax.set_yticks([-100, -50, 0, 50, 100])
ax.set_xlabel(r'R [kpc]')
ax.set_ylabel(r'z [kpc]')


plot_name = 'm82_mag_map_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break

plt.close()		



























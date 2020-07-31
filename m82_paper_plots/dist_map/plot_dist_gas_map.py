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
gas_density = imp.load_source("gas_density_v4", "../m82_dist_py/gas_density_v4.py")

def DoRotation(xspan, yspan, RotRad=0):
    """Generate a meshgrid and rotate it by RotRad radians."""

    # Clockwise, 2D rotation matrix
    RotMatrix = np.array([[np.cos(RotRad),  np.sin(RotRad)],
                          [-np.sin(RotRad), np.cos(RotRad)]])

    x, y = np.meshgrid(xspan, yspan)
    return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))

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

DISTANCE = 3570

x=np.arange(-1,1,0.01)
y=np.arange(-1,1,0.01)
z=np.arange(-1,1,0.01)

xnum = len(x)
ynum = len(y)
znum = len(z)

y_diff = y[1]-y[0]

###
### GAS DENSITY 1
###
nHI_factor = 1.00
nHII_factor = 0.


density = 300.

r0 = 0.0001 #kpc
r1 = 0.2
z0 = 0.0001 #kpc
z1 = 0.05

v0 = 200.
mass_loss = 30.0

nHI = np.zeros([znum,ynum,xnum])
for ix in range(0,xnum):
	for iy in range(0,ynum):
		for iz in range(0,znum):
			nHI[iz][iy][ix] = gas_density.m82_wind_radial_v7( x[ix], y[iy], z[iz], density, r0, r1, z0, z1, v0, mass_loss)
	print x[ix]

nHI = np.sum(nHI, axis=1)*y_diff*kpc2cm
# nHI = nHI.T
nHII_1 = nHII_factor*nHI
nHI_1 = nHI_factor*nHI
print np.amax(nHI)


x*= radian2arcsec/DISTANCE
z*= radian2arcsec/DISTANCE

ANGLE = -1*20.
ANGLE *= math.pi/180 
X, Z = DoRotation(x, z, RotRad=ANGLE)

###
### [z, nHI] [z, nHII]
#########################


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
loglow = 21.
loghigh = 24.4
V = np.logspace(loglow, loghigh, 29)

loglow = 21.
loghigh = 24.
V1 = np.logspace(loglow, loghigh, 13)

CF0 = ax.contourf(X,Z,nHI_1,norm = LogNorm(), levels=V, cmap=cmap)

# CF1 = ax.contour(RR,ZZ,nHI_1, levels=[density/math.e**2, density/math.e, density-0.1])


#colorbar
cf1 = plt.colorbar(CF0, ax=ax, ticks=V1, pad=0.0025, fraction=0.04763)
cf1.ax.tick_params(labelsize='small')
cf1.ax.yaxis.set_tick_params(color='w', width=1, length=5)
cf1.outline.set_linewidth(2) 
cf1.set_label(r'$n$ [cm$^{-2}$]', rotation=270, horizontalalignment='center', verticalalignment='bottom')

#axes stuff
ax.set_axis_bgcolor('black')
ax.set_aspect('equal')
ax.set_xlim([32.5,-32.5])
ax.set_ylim([-20,20])
# ax.set_xlim([300,-300])
# ax.set_ylim([-300,300])
# ax.set_xticks([-100, -50, 0, 50, 100])
# ax.set_yticks([-100, -50, 0, 50, 100])
ax.set_xlabel('Relative R.A. [\"]')
ax.set_ylabel('Relative Dec. [\"]')


plot_name = 'm82_gas_image_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break

plt.close()		



























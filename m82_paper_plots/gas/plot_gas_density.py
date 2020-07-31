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


c = 3.0e10
ev2ergs = 1.60218e-12
mp = 1.66e-24 #g
kpc2cm = 3.08568025e21 #cm/kpc

zmin = 0.0
zmax = 4.0
znum = 4001
z = np.linspace(0, zmax, znum)

###############################################
###
### ORIGINAL
### 0.05 scale height + n_wind
###

###
### GAS DENSITY 1
###
nHI_factor = 1.
nHII_factor = 0.0

density = 150.

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
density = 700.

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
### 700.,0.2,0.2,0.05,0.05,200.,30.,7,8,9  
density = 900.

r0 = 0.200 #kpc
r1 = 0.20
z0 = 0.050 #kpc
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

###############################################
###
### ORIGINAL
### 0.05 scale height + leroy wind
###

###
### GAS DENSITY
###
nHI = np.zeros(znum)

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.050 #kpc
z1w = 0.20

power_law = -10.

rwindw = 0.2

v0_rw = 1000.
v0_zw = 0.
densityw = 150.
mass_lossw = 10.0

nHI = np.zeros(znum)
for iz in range(0,znum):
	nHI[iz] = gas_density.m82_wind_radial_v3( 0, 0, z[iz], r0w, r1w, rwindw, z0w, z1w, v0_rw, densityw, mass_lossw, power_law)

nHII_a5 = nHII_factor*nHI
nHI_a5 = nHI_factor*nHI

###
### [z, nHI] [z, nHII]
#########################

###
### GAS DENSITY
###
nHI = np.zeros(znum)

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.050 #kpc
z1w = 0.20

power_law = -10.

rwindw = 0.2

v0_rw = 1000.
v0_zw = 0.
densityw = 600.
mass_lossw = 10.0

nHI = np.zeros(znum)
for iz in range(0,znum):
	nHI[iz] = gas_density.m82_wind_radial_v3( 0, 0, z[iz], r0w, r1w, rwindw, z0w, z1w, v0_rw, densityw, mass_lossw, power_law)

nHII_b5 = nHII_factor*nHI
nHI_b5 = nHI_factor*nHI

###
### [z, nHI] [z, nHII]
#########################

###############################################
###
### 0.05 constant core
### todd wind
###

###
### GAS DENSITY
###
nHI = np.zeros(znum)

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.050 #kpc
z1w = 0.20

power_law = -10.

rwindw = 6

v0_rw = 200.
v0_zw = 0.
densityw = 150.
mass_lossw = 20.0

nHI = np.zeros(znum)
for iz in range(0,znum):
	nHI[iz] = gas_density.m82_constant_powerlaw( 0, 0, z[iz], r0w, r1w, rwindw, z0w, z1w, v0_rw, densityw, mass_lossw, power_law)

nHII_ta = nHII_factor*nHI
nHI_ta = nHI_factor*nHI

###
### [z, nHI] [z, nHII]
#########################

###
### GAS DENSITY
###
nHI = np.zeros(znum)

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.050 #kpc
z1w = 0.20

power_law = -10.

rwindw = 6

v0_rw = 200.
v0_zw = 0.
densityw = 600.
mass_lossw = 20.0

nHI = np.zeros(znum)
for iz in range(0,znum):
	nHI[iz] = gas_density.m82_constant_powerlaw( 0, 0, z[iz], r0w, r1w, rwindw, z0w, z1w, v0_rw, densityw, mass_lossw, power_law)

nHII_tb = nHII_factor*nHI
nHI_tb = nHI_factor*nHI

###
### [z, nHI] [z, nHII]
#########################


###
### GAS DENSITY
###

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.20 #kpc
z1w = 0.20

power_law = -10.

rwindw = 0.002

v0_rw = 1000.
v0_zw = 0.
densityw = 0.
mass_lossw = 10.0

nHI10 = np.zeros(znum)
for iz in range(0,znum):
	nHI10[iz] = gas_density.m82_wind_radial_v2a( 0, 0, z[iz], r0w, r1w, rwindw, z0w, z1w, v0_rw, densityw, mass_lossw, power_law)

nHII10 = nHII_factor*nHI
nHI10*= nHI_factor

###
### [z, nHI] [z, nHII]
#########################

###
### GAS DENSITY
###

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.20 #kpc
z1w = 0.20

power_law = -10.

rwindw = 0.002

v0_rw = 1000.
v0_zw = 0.
densityw = 0.
mass_lossw = 30.0

nHI11 = np.zeros(znum)
for iz in range(0,znum):
	nHI11[iz] = gas_density.m82_wind_radial_v2a( 0, 0, z[iz], r0w, r1w, rwindw, z0w, z1w, v0_rw, densityw, mass_lossw, power_law)

nHII11 = nHII_factor*nHI
nHI11*= nHI_factor

###
### [z, nHI] [z, nHII]
#########################

###
### GAS DENSITY
###

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.20 #kpc
z1w = 0.20

# g=1e-7
# dv=20.e5

# # nHI20 = np.zeros(znum)
# nHI20 = 300.*np.exp(-g*z*kpc2cm/dv**2.)

nHI20 = 150.*np.exp(-z/0.05)

nHII20 = nHII_factor*nHI
nHI20*= nHI_factor

###
### [z, nHI] [z, nHII]
#########################

###
### GAS DENSITY
###

gas_model = 1
nHI_factor = 1.0 #0.95
nHII_factor = 0.#05

r0w = 0.200 #kpc
r1w = 1.000
z0w = 0.20 #kpc
z1w = 0.20

# g=1e-7
# dv=40.e5

# # nHI20 = np.zeros(znum)
# nHI21 = 300.*np.exp(-g*z*kpc2cm/dv**2.)
nHI21 = 600.*np.exp(-z/0.05)

nHII21 = nHII_factor*nHI
nHI21*= nHI_factor

###
### [z, nHI] [z, nHII]
#########################
	


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

color_a = cmap(1.-1./4.)
color_b = cmap(1.-2./4)
color_c = cmap(1.-3./4)

colorI = cmap(0.1)
colorII = cmap(0.85)


ax.plot(z,nHI_1, c=color_a, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label='A')
ax.plot(z,nHI_2, c=color_b, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label='B')
ax.plot(z,nHI_3, c=color_c, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label="B'")

# ax.plot(z,nHI_a5, c=color_a, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI_b5, c=color_b, linestyle='-', linewidth=LWIDTH, alpha=ALPHA)

# ax.plot(z,nHI_ta, c=color_a, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI_tb, c=color_b, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)

# # ax.plot(z,nHI10, c=colorI, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
# # ax.plot(z,nHI11, c=colorI, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI20, c=color_a, linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
# ax.plot(z,nHI21, c=color_b, linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)

# ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.0,3.5])
ax.set_ylim([2.e-3,2.e3])
ax.set_xlabel(r'$z$ [kpc]')
ax.set_ylabel(r'$n$ [cm$^{-3}$]')

ax.legend(loc='upper right', prop={'size':15})

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

# ax.legend(loc='lower right', prop={'size':14})

plt.tight_layout()
plot_name = 'plot_gas_z_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break
plt.close()




















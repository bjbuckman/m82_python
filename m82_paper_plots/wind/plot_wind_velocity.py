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

zmin = 0.0
zmax = 4.0
znum = 4001
z = np.linspace(0, zmax, znum)


###
### GAS DENSITY 1
###
nHI_factor = 1.0
nHII_factor = 1.0

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
density = 675.

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
density = 1000.

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

B0 = 150.e-6
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

B0 = 325.e-6
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

B0 = 325.e-6
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

vA_a = B_1/np.sqrt(4.*np.pi*nHII_1*mp)/1.e5
vA_b = B_2/np.sqrt(4.*np.pi*nHII_2*mp)/1.e5
vA_c = B_3/np.sqrt(4.*np.pi*nHII_3*mp)/1.e5

# print vA_a

###
### Wind velocities
###
data = np.loadtxt('ben_1_cc85_200pc_20Msunyr_a1.0_b1.0.dat', delimiter=',')

Z = data[0]
v0 = data[2]
v0[0]=0.

v_a = v0*800.
v_b = v0*1000.
v_c = v0*1000.

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


ax.plot(Z,v_a, c=color_a, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label='A')
ax.plot(Z,v_b, c=color_b, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label='B')
ax.plot(Z,v_c, c=color_c, linestyle='-', linewidth=LWIDTH, alpha=ALPHA, label="B'")

ax.plot(z,vA_a, c=color_a, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
ax.plot(z,vA_b, c=color_b, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
ax.plot(z,vA_c, c=color_c, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)

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
ax.set_xlim([0.,3.5])
ax.set_ylim([0.01,2000.])
ax.set_xlabel(r'$z$ [kpc]')
ax.set_ylabel(r'$V_0$ [km s$^{-1}$]')

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

# ax.legend(loc='lower right', prop={'size':14})

plt.tight_layout()
plot_name = 'plot_wind_r_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break
plt.close()




















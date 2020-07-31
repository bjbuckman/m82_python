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
magnetic_field = imp.load_source("magnetic_field_v4", "../m82_dist_py/magnetic_field_v4.py")
accel = imp.load_source("gp2d_acceleration_v1", "../m82_dist_py/gp2d_acceleration_v1.py")

c = 3.0e10
mp = 1.66e-24 #g
ev2ergs = 1.60218e-12
kpc2cm = 3.08568025e21 #cm/kpc

zmin = 0.0
zmax = 4.0
znum = 401
z = np.linspace(0, zmax, znum)

####################################################
# ##################################################
# MAGNETIC FIELD 1 #################################
# ##################################################
####################################################
###
### MAGNETIC FIELD 1_3
###

magnetic_model_name = 'constant+powerlaw_v4'

#B_field_parameters   = 0.000150,0.2,0.2,0.05,0.05,0.05,0.2,-1.4,8,9  

B0 = 150.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

# B = np.zeros(znum)
# for iz in range(0,znum):
B_1 = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,0.]), B_parameters)
	

###
### [z, B] [z, grad_magnetic]
##############################################


####################################################
# ##################################################
# MAGNETIC FIELD 3 #################################
# ##################################################
####################################################
###
### MAGNETIC FIELD 3_3
###

magnetic_model_name = 'constant+powerlaw_v4'
#B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -0.2

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

# B = np.zeros(znum)
# for iz in range(0,znum):
B_3 = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,0.]), B_parameters)
	

###
### [z, B] [z, grad_magnetic]
##############################################


####################################################
# ##################################################
# MAGNETIC FIELD 2 #################################
# ##################################################
####################################################
###
### MAGNETIC FIELD 2_3
###

magnetic_model_name = 'constant+powerlaw_v4'
#B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

B0 = 325.e-6
r_scale = 0.2
r_scale1 = 0.2
z_scale = 0.05
z_scale1 = 0.05

r_pow = 0.05
r_pow1 = 0.2
b_pow = -1.2

B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

# B = np.zeros(znum)
# for iz in range(0,znum):
B_2 = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,0.]), B_parameters)
	

###
### [z, B] [z, grad_magnetic]
##############################################

# ############################################
# ############################################
# ############################################



# ############################################
# ############################################
# ############################################

# ###
# ### PLOT
# ###

# font = {'family' : 'serif',
        # 'weight' : 'normal',
        # 'size'   : 18}

# matplotlib.rc('font', **font)

# cmap = plt.get_cmap('cubehelix')
# fig = plt.figure()
	
# f1, ax = plt.subplots(1,1)#, figsize=(18,6))

# for axis in ['top','bottom','left','right']:
	# ax.spines[axis].set_linewidth(2)
# ax.tick_params(which='major',width=1, length=8, labelsize='small')
# ax.tick_params(which='minor',width=0.5, length=5)

# MSIZE = 5
# LWIDTH = 3
# ALPHA = 0.65


###
###
###
### PARTICLES
###
###
###

dir_outfits = '../fits/raw/'
dir_comp = '../fits/'

galdef_name = []

GG = 'm82_05aff1018'
galdef_name.append(GG)

GG =  'm82_05b2ff1021'
galdef_name.append(GG)

GG =  'm82_06bff1021'
galdef_name.append(GG)

# galdef_name_1 = []

# GG_1 = 'm82_05aff1019'
# galdef_name_1.append(GG_1)

# GG_1 = 'm82_05b2ff1022'
# galdef_name_1.append(GG_1)

# GG_1 = 'm82_06bff1022'
# galdef_name_1.append(GG_1)

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
	
# chi_1 = np.zeros(num_graph)
# chi_fac_1 = np.ones(num_graph)
# for n in range(0,num_graph):
	# chi_dat_1 = np.loadtxt(dir_comp+galdef_name_1[n]+'/'+galdef_name_1[n]+'chi.dat')
	# chi_fac_1[n] = chi_dat_1[4]
	# GAMMA_CHI = 1
	# chi_1[n] = GAMMA_CHI*(chi_dat_1[1]+chi_dat_1[2])+chi_dat_1[3]

# if NORM_PLOT == 1:
	# FACTOR_1 = np.ones(num_graph_galdef_tot)
# else:
	# FACTOR_1 = chi_fac_1	

LABEL = ['A','B',"B'"]
for n in range(0,num_graph):
	
	nuclei_file = dir_outfits+galdef_name[n]+'/'+'nuclei_full_54_'+galdef_name[n]+'.gz'
	# nuclei_file_1 = dir_outfits+galdef_name_1[n]+'/'+'nuclei_full_54_'+galdef_name_1[n]+'.gz'

	synch_file = dir_outfits+galdef_name[n]+'/'+'synchrotron_emiss_54_'+galdef_name[n]+'.gz'
	freef_file = dir_outfits+galdef_name[n]+'/'+'free_free_absorb_54_'+galdef_name[n]+'.gz'
	###
	### ELECTRONS
	###

	# en, electrons_n = accel.CR_electron_N(nuclei_file)
	# electrons_n*= FACTOR[n]
	
	# p = np.log(electrons_n[1:]/electrons_n[:-1])/np.log(en[1:]/en[:-1])
	
	nu, synch_p = accel.synch_power(synch_file)
	
	# print synch_p.shape
	
	alpha_1 = np.log(synch_p[1:]/synch_p[:-1])/np.log(nu[1:]/nu[:-1])
	alpha = (alpha_1[1:]+alpha_1[:-1])/2.
	p = 1.-2.*alpha
	
	nu = nu[1:-1]
	synch_p = synch_p[1:-1]
	
	absorb = (p+2.)*c**2/8./np.pi/nu**2*synch_p
	
	nu_f, absorp_f = accel.synch_power(freef_file)
	absorp_f = absorp_f[1:-1]
	
	
	print absorb/absorp_f
	print nu
	
		
	###
	### [z, electrons_accel]
	#########################
	
	# color = cmap(1-float(n+1)/(num_graph+1.))
	
	# if n == 0:
		# ax.fill_between(z, abs(magnetic_accel_1_1), abs(magnetic_accel_1_2), facecolor=color, edgecolor=color, alpha=0.2)#, hatch='|')
		# # ax.plot(z, isrf_accel, c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# # ax.plot(z, abs(magnetic_accel_1), c=color, label='magnetic', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		# # ax.plot(z, abs(isrf_accel), c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(z, abs(magnetic_accel_1_3), c=color, linestyle='-.', linewidth=LWIDTH, alpha=1.0)
	# elif n == 1:
		# ax.fill_between(z, abs(magnetic_accel_2_1), abs(magnetic_accel_2_2), facecolor=color, edgecolor=color, alpha=0.2)#, hatch='|')
		# # ax.plot(z, isrf_accel2, c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# # ax.plot(z, abs(magnetic_accel_2), c=color, label='magnetic', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		# # ax.plot(z, abs(isrf_accel2), c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(z, abs(magnetic_accel_2_3), c=color, linestyle='-.', linewidth=LWIDTH, alpha=1.0)
	# elif n == 2:
		# ax.fill_between(z, abs(magnetic_accel_3_1), abs(magnetic_accel_3_2), facecolor=color, edgecolor=color, alpha=0.2)#, hatch='|')
		# # ax.plot(z, isrf_accel2, c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# # ax.plot(z, abs(magnetic_accel_2), c=color, label='magnetic', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
		# # ax.plot(z, abs(isrf_accel2), c=color, label='isrf', linestyle='-.', linewidth=LWIDTH, alpha=ALPHA)
		# ax.plot(z, abs(magnetic_accel_3_3), c=color, linestyle='-.', linewidth=LWIDTH, alpha=1.0)
	
	# ax.fill_between(z, protons_accel, protons_accel_1, facecolor=color, edgecolor=color, alpha=1.0, linewidth=LWIDTH*2., label=LABEL[n])
	# # ax.fill_between(z, electrons_accel, electrons_accel_1, facecolor=color, edgecolor='black', alpha=ALPHA, hatch='|')
	# # ax.plot(z, protons_accel, c=color, label='CR_proton', linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, electrons_accel, c=color, label=None, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	# ax.plot(z, electrons_accel_1, c=color, label=None, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	# # print electrons_accel

# ax.plot(z,U_isrf, c='black', linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	
# # ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlim([0.0,3.5])
# ax.set_ylim([3.e-2,3.e3])
# ax.set_xlabel(r'$z$ [kpc]')
# ax.set_ylabel(r'Energy Density [eV cm$^{-3}$]')

# # ax.text(0.66, 0.75, r'magnetic field', ha='center', va='center', transform=ax.transAxes)
# # ax.text(0.66, 0.5, r'isrf', ha='center', va='center', transform=ax.transAxes)
# # ax.text(0.66, 0.25, r'cosmic rays', ha='center', va='center', transform=ax.transAxes)

# # ax.set_title(r'$0$ Temperature Acceleration')

# # ax.legend(loc='upper right', prop={'size':14})
# ax.legend(loc='lower right', prop={'size':15})

# plt.tight_layout()
# plot_name = 'plot_dist_'
# for n in range(0,100):
	# if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		# continue
	# else:
		# plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		# break
# plt.close()




















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
# gas_density = imp.load_source("gas_density_v4", "../m82_dist_py/gas_density_v4.py")
# magnetic_field = imp.load_source("magnetic_field_v4", "../m82_dist_py/magnetic_field_v4.py")
# accel = imp.load_source("gp2d_acceleration_v1", "../m82_dist_py/gp2d_acceleration_v1.py")
particle = imp.load_source("gp2d_particle", "../m82_dist_py/gp2d_particle.py")

c = 3.0e10
mp = 1.66e-24 #g
ev2ergs = 1.60218e-12
kpc2cm = 3.08568025e21 #cm/kpc

zmin = 0.0
zmax = 4.0
znum = 401
z = np.linspace(0, zmax, znum)

# ############################################
# ############################################
# ############################################

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

for n in range(0,num_graph):
	LABEL = ['A','B',"B'"]
	
	nuclei_file = dir_outfits+galdef_name[n]+'/'+'nuclei_full_54_'+galdef_name[n]+'.gz'
	nuclei_file_1 = dir_outfits+galdef_name_1[n]+'/'+'nuclei_full_54_'+galdef_name_1[n]+'.gz'
	
	pelec_file = dir_outfits+galdef_name[n]+'/'+'primary_electrons_source_function_54_'+galdef_name[n]+'.gz'
	selec_file = dir_outfits+galdef_name[n]+'/'+'secondary_electrons_source_function_54_'+galdef_name[n]+'.gz'
	sposi_file = dir_outfits+galdef_name[n]+'/'+'secondary_positrons_source_function_54_'+galdef_name[n]+'.gz'
	
	pelec_file_1 = dir_outfits+galdef_name_1[n]+'/'+'primary_electrons_source_function_54_'+galdef_name_1[n]+'.gz'
	selec_file_1 = dir_outfits+galdef_name_1[n]+'/'+'secondary_electrons_source_function_54_'+galdef_name_1[n]+'.gz'
	sposi_file_1 = dir_outfits+galdef_name_1[n]+'/'+'secondary_positrons_source_function_54_'+galdef_name_1[n]+'.gz'
	
	###
	### LEPTONS
	###
	
	pe_z, pe_lum = particle.particle_luminosity_z(pelec_file, 'electron')
	pe_lum*= FACTOR[n]

	pe_z_1, pe_lum_1 = particle.particle_luminosity_z(pelec_file_1, 'electron')
	pe_lum_1*= FACTOR_1[n]

	se_z, se_lum = particle.particle_luminosity_z(selec_file, 'electron')
	se_lum*= FACTOR[n]
	
	se_z_1, se_lum_1 = particle.particle_luminosity_z(selec_file_1, 'electron')
	se_lum_1*= FACTOR_1[n]
	
	sp_z, sp_lum = particle.particle_luminosity_z(sposi_file, 'electron')
	sp_lum*= FACTOR[n]
	
	sp_z_1, sp_lum_1 = particle.particle_luminosity_z(sposi_file_1, 'electron')
	sp_lum_1*= FACTOR_1[n]

		
	###
	### [z, electrons_accel]
	#########################
	
	color = cmap(1-float(n+1)/(num_graph+1.))
	
	# ax.fill_between(pe_z, pe_lum, pe_lum_1, facecolor=color, edgecolor=color, alpha=0.7, linewidth=LWIDTH*2., label=LABEL[n])
	# ax.fill_between(se_z, se_lum, se_lum_1, facecolor=color, edgecolor=color, alpha=0.5, linewidth=LWIDTH*2., label=None)
	# ax.fill_between(sp_z, sp_lum, sp_lum_1, facecolor=color, edgecolor=color, alpha=0.3, linewidth=LWIDTH*2., label=None)
	
	ax.plot(pe_z, pe_lum*pe_z**2, alpha=ALPHA, c=color, linestyle='-', linewidth=LWIDTH, label=LABEL[n])
	ax.plot(se_z, se_lum*se_z**2, alpha=ALPHA, c=color, linestyle=':', linewidth=LWIDTH)
	ax.plot(sp_z, sp_lum*sp_z**2, alpha=ALPHA, c=color, linestyle='-.', linewidth=LWIDTH)

	
# ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.0,3.5])
ax.set_ylim([1.e-18,1.e-13])
ax.set_xlabel(r'$z$ [kpc]')
ax.set_ylabel(r'Injection Luminosity [eV cm$^{-3}$ s$^{-1}$]')

# ax.text(0.66, 0.75, r'magnetic field', ha='center', va='center', transform=ax.transAxes)
# ax.text(0.66, 0.5, r'isrf', ha='center', va='center', transform=ax.transAxes)
# ax.text(0.66, 0.25, r'cosmic rays', ha='center', va='center', transform=ax.transAxes)

# ax.set_title(r'$0$ Temperature Acceleration')

# ax.legend(loc='upper right', prop={'size':14})
ax.legend(loc='upper right', prop={'size':15})

plt.tight_layout()
plot_name = 'plot_injection_e_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break
plt.close()




















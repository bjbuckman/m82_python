import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import os

import imp
particle = imp.load_source("gp2d_particle", "../m82_dist_py/gp2d_particle.py")

SL_SCALE = 0

E2FLAG = 1

DISTANCE = 3520

MODEL = 2

if MODEL == 0:
	galdef = 'm82_05aff1018'
elif MODEL == 1:
	galdef = 'm82_05b2ff1021'
elif MODEL == 2:
	galdef = 'm82_06bff1021'

dir_outfits = '../fits/raw/'
dir_comp = '../fits/'

galdef_name = []

GG = galdef
galdef_name.append(GG)

# GG = 'm82_21bf1022'
# galdef_name.append(GG)

# LABEL_ALL = [r'A',r'B']
na_fac = 1.0

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
	
	
#constants
pc2cm = 3.08568025e18 #cm/pc
kpc2cm = 3.08568025e21 #cm/kpc
radian2arcsec = 648000./math.pi #acsec/radian
Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
c = 3.0e10 #cm/s
h = 6.62e-27 #erg s
q=4.8e-10 #esu
m=9.11e-28 #g

DISTANCE = 3520 #kpc
D_cm = DISTANCE*kpc2cm

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
ax.tick_params(which='major',width=1, length=8, labelsize='small')
ax.tick_params(which='minor',width=0.5, length=5)

MSIZE = 5
LWIDTH = 3

Z = np.array([0.0, 0.2, 0.5, 3.0])#, 1.5])

for n in range(0,num_graph):

	LABEL = [r'$z=0.0\,$kpc',r'$z=0.2\,$kpc',r'$z=0.5\,$kpc',r'$z=3.0\,$kpc']#,r'$z=1.5$ kpc']
	
	for m in range(0,len(Z)):
	
		color = cmap((float(m+1))/(len(Z)+2))
		
		exponent = 0
		###
		###
		### protons
		###
		###
		###
		
		# energy, p_tot = particle.proton_spectrum_total(dir_outfits+galdef_name[n]+'/'+'nuclei_full_54_'+galdef_name[n]+'.gz')
		
		energy, p_cor = particle.proton_spectrum_kinetic(dir_outfits+galdef_name[n]+'/'+'nuclei_full_54_'+galdef_name[n]+'.gz', 0., Z[m])
		
		p_species = p_cor.shape[0]
		
		# p_tot*= energy**exponent
		p_cor*= energy**exponent
		
		p_tot = np.sum(p_cor, axis=0)
		
		LINE = ['-.','--',':',':']
		ax.loglog(energy, np.sum(p_cor, axis=0)*FACTOR[n], label=r''+LABEL[m], c=color, linestyle='-', linewidth=LWIDTH, alpha=0.9)
		for i in range(0,p_species):
			ax.loglog(energy, p_cor[i]*FACTOR[n], label=None, c=color, linestyle=LINE[i], linewidth=LWIDTH, alpha=0.75)
			
		alpha = np.log(p_tot[1:]/p_tot[:-1])/np.log(energy[1:]/energy[:-1])
		alpha = (alpha[:-1]+alpha[1:])/2.
		e_alpha = energy[1:-1]
		
		e_vals = np.array([1.e3,1.e4,1.e5,1.e6])
		alpha_func = np.interp(e_vals,e_alpha,alpha)
		print Z[m], e_vals, alpha_func
		
ax.set_ylabel(r'$E \, \frac{dU}{dE} $ [MeV$^2$ cm$^{-3}$ MeV$^{-1}$]')
ax.set_xlabel(r'$E$ [MeV]')

ax.set_xlim([1.e1,1.0001e8])
ax.set_ylim([2.e-10, 1.e-4])

if MODEL == 0:
	ax.text(0.666,0.333,r'CRp$-$Model A', ha='center', va='center', color='k', transform=ax.transAxes, backgroundcolor='white')
	plot_name = 'plot_particle_p_A'
elif MODEL == 1:
	ax.text(0.666,0.333,r'CRp$-$Model B', ha='center', va='center', color='k', transform=ax.transAxes, backgroundcolor='white')
	plot_name = 'plot_particle_p_B'
elif MODEL == 2:
	ax.text(0.666,0.333,r"CRp$-$Model B'", ha='center', va='center', color='k', transform=ax.transAxes, backgroundcolor='white')
	plot_name = 'plot_particle_p_Bp'

# ax.text(1.e5,1.e-8,r'CRp$-$Model B', backgroundcolor='white')

# ax.legend(loc='upper right', prop={'size':15})

plt.tight_layout()
# plot_name = 'plot_particle_p_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.png', bbox_inches='tight', dpi=400)
		break
plt.close()


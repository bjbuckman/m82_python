import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import os

SL_SCALE = 1

E2FLAG = 1

DISTANCE = 3520

GALDEF = 'm82_12'
CHANGED = 'integrated_gamma'
num_params = 4
num_galdef = 3#np.array([2, 2, 2, 2])
num_galdef_tot = 3 #int(np.prod(num_galdef))

VAR_PLOT = ['$B_{0}$ ', '$\\rho_0=$ ', '$r_{s}=$', '$z_{s}=$']

BEGIN = np.array([10, 10, 10, 10])
num_start = BEGIN

# num_graph_start = np.array([11, 10, 10, 10])
num_graph_galdef = 3#np.array([1, 1, 1, 2])
num_graph_galdef_tot = 3#int(np.prod(num_graph_galdef))

N=np.zeros(num_params).astype(int)
NN=np.zeros([num_graph_galdef_tot, num_params]).astype(int)
galdef_name = []
galdef_name_fill = []
galdef_dir = '/mnt/c/Users/psyko/Physics/proj_gal/1904/m82_paper_plots/fits/'

GG = 'm82_05aff1017'
galdef_name.append(galdef_dir+GG+'/'+GG)

GG = 'm82_05b2ff1019'
galdef_name.append(galdef_dir+GG+'/'+GG)

GG = 'm82_06bff1020'
galdef_name.append(galdef_dir+GG+'/'+GG)

GG2 = 'm82_05aff1019'
galdef_name_fill.append(galdef_dir+GG2+'/'+GG2)

GG2 = 'm82_05b2ff1021'
galdef_name_fill.append(galdef_dir+GG2+'/'+GG2)

GG2 = 'm82_06bff1022'
galdef_name_fill.append(galdef_dir+GG2+'/'+GG2)


LABEL_ALL = [r'A',r'B',r"B'"]
na_fac = 1.0

#THESE SHOULD BE THE SAME
num_graph = len(galdef_name)
num_graph_fill = len(galdef_name_fill)				

chi = np.zeros(num_graph)
chi_fac = np.ones(num_graph)
for n in range(0,num_graph):
	chi_dat = np.loadtxt(galdef_name[n]+'chi.dat')
	chi_fac[n] = chi_dat[4]
	GAMMA_CHI = 1
	chi[n] = GAMMA_CHI*(chi_dat[1]+chi_dat[2])+chi_dat[3]
	
chi_fill = np.zeros(num_graph_fill)
chi_fac_fill = np.ones(num_graph_fill)
for n in range(0,num_graph_fill):
	chi_dat_fill = np.loadtxt(galdef_name_fill[n]+'chi.dat')
	chi_fac_fill[n] = chi_dat_fill[4]
	GAMMA_CHI = 1
	chi_fill[n] = GAMMA_CHI*(chi_dat_fill[1]+chi_dat_fill[2])+chi_dat_fill[3]
				
				
new_val = np.array([[  2.50000000e-04,   0.00000000e+00],
       [  5.00000000e+02,   0.00000000e+00],
       [  1.00000000e+00,   0.00000000e+00],
       [  1.00000000e+00,   2.00000000e+00]])


NORM_PLOT = 0
if NORM_PLOT == 1:
	FACTOR = np.ones(num_graph_galdef_tot)
	FACTOR_fill = np.ones(num_graph_galdef_tot)
else:
	FACTOR = chi_fac
	FACTOR_fill = chi_fac_fill
	
	
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
gs = gridspec.GridSpec(3, 1)
gs.update(wspace=0.0, hspace=0.0)
plt.subplots_adjust(hspace=0.001)

ax1 = plt.subplot(gs[0:2,0])
ax_e = plt.subplot(gs[2,0])#, sharex = ax1)

for axis in ['top','bottom','left','right']:
	ax1.spines[axis].set_linewidth(2)
	ax_e.spines[axis].set_linewidth(2)
ax1.tick_params(which='major',width=1, length=8, labelsize='small')
ax1.tick_params(which='minor',width=0.5, length=5)
ax_e.tick_params(which='major',width=1, length=8, labelsize='small')
ax_e.tick_params(which='minor',width=0.5, length=5)

###
###
### GAMMA-RAY
###
###
###

MSIZE = 10
LWIDTH = 3

#tim
fermi_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_spectrum_t2_larger_roi.txt').T
emin = fermi_data[0]
emax = fermi_data[1]
fermi_energy = np.sqrt(emin*emax)
# ferr_d = (fermi_data[2]-fermi_data[3])*fermi_energy**2/(emax-emin)
# ferr_u = (fermi_data[2]+fermi_data[3])*fermi_energy**2/(emax-emin)
ferr_d = (fermi_data[3])*fermi_energy**2/(emax-emin)
ferr_u = (fermi_data[3])*fermi_energy**2/(emax-emin)
ax1.errorbar(fermi_energy, fermi_energy*fermi_data[2]/np.log(emax/emin), xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label=None, c='black', markersize=MSIZE, alpha=0.5)

# #3fgl x middle energy
# fermi_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_fermi.dat').T
# emin = fermi_data[1]
# emax = fermi_data[2]
# fermi_energy = np.sqrt(emin*emax)
# ferr_d = fermi_data[3]*fermi_energy**2/(emax-emin)
# ferr_u = fermi_data[4]*fermi_energy**2/(emax-emin)
# ax1.errorbar(fermi_energy, fermi_energy*fermi_data[0]/np.log(emax/emin), xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label=None, c='blue', markersize=MSIZE, alpha=0.5)

#3fgl x -2.2 spectrum
# neg_gamma = 2.2
# fermi_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_fermi.dat').T
# emin = fermi_data[1]
# emax = fermi_data[2]
# fermi_energy = np.sqrt(emin*emax)
# fermi_factor = -neg_gamma/(emin**neg_gamma-emax**neg_gamma)*np.sqrt(emin*emax)**(neg_gamma+1.)
# ferr_d = fermi_data[3]*fermi_factor
# ferr_u = fermi_data[4]*fermi_factor
# ax1.errorbar(fermi_energy, fermi_factor*fermi_data[0], xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label=None, c='red', markersize=MSIZE, alpha=0.5)


verit_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_veritas.dat').T
verit_upper = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_veritas_upper.dat')
#verit_data[0] *= 1000
verit_data[1:] *= verit_data[0]
verit_upper[1] *= verit_upper[0]


HATCH = ['-', '|','/']
	
for n in range(0,num_graph):
	
	LABEL = LABEL_ALL[n]
	
	###
	###
	### GAMMA-RAY
	###
	###
	###
	
	#model 1
	gamma_tot = np.loadtxt(galdef_name[n]+'gamma_flux.dat')
	energy = gamma_tot[0]
	
	gamma_tot_mask = np.isfinite(gamma_tot[1])
	gamma_tot = gamma_tot[:,gamma_tot_mask]
	gamma_tot = np.interp(energy, gamma_tot[0], gamma_tot[1])*FACTOR[n]
	
	bremss_tot = np.loadtxt(galdef_name[n]+'bremss_flux.dat')
	bremss_tot_mask = np.isfinite(bremss_tot[1])
	bremss_tot = bremss_tot[:,bremss_tot_mask]
	bremss_tot = np.interp(energy, bremss_tot[0], bremss_tot[1])*FACTOR[n]
	
	pion_tot = np.loadtxt(galdef_name[n]+'pion_flux.dat')
	pion_tot_mask = np.isfinite(pion_tot[1])
	pion_tot = pion_tot[:,pion_tot_mask]
	pion_tot = np.interp(energy, pion_tot[0], pion_tot[1])*FACTOR[n]
	
	ic_tot = np.loadtxt(galdef_name[n]+'IC_flux.dat')
	ic_tot_mask = np.isfinite(ic_tot[1])
	ic_tot = ic_tot[:,ic_tot_mask]
	ic_tot = np.interp(energy, ic_tot[0], ic_tot[1])*FACTOR[n]
	
	fermi_data = np.loadtxt(galdef_name[n]+'fermi.dat').T
	emin = fermi_data[0]
	emax = fermi_data[1]
	fermi_energy = np.sqrt(emin*emax)
	ferr_d = fermi_data[3]
	ferr_u = fermi_data[4]
	
	#model 2
	gamma_tot_fill = np.loadtxt(galdef_name_fill[n]+'gamma_flux.dat')
	gamma_tot_mask_fill = np.isfinite(gamma_tot_fill[1])
	gamma_tot_fill = gamma_tot_fill[:,gamma_tot_mask_fill]
	gamma_tot_fill = np.interp(energy, gamma_tot_fill[0], gamma_tot_fill[1])*FACTOR_fill[n]
	
	bremss_tot_fill = np.loadtxt(galdef_name_fill[n]+'bremss_flux.dat')
	bremss_tot_mask_fill = np.isfinite(bremss_tot_fill[1])
	bremss_tot_fill = bremss_tot_fill[:,bremss_tot_mask_fill]
	bremss_tot_fill = np.interp(energy, bremss_tot_fill[0], bremss_tot_fill[1])*FACTOR_fill[n]
	
	pion_tot_fill = np.loadtxt(galdef_name_fill[n]+'pion_flux.dat')
	pion_tot_mask_fill = np.isfinite(pion_tot_fill[1])
	pion_tot_fill = pion_tot_fill[:,pion_tot_mask_fill]
	pion_tot_fill = np.interp(energy, pion_tot_fill[0], pion_tot_fill[1])*FACTOR_fill[n]
	
	ic_tot_fill = np.loadtxt(galdef_name_fill[n]+'IC_flux.dat')
	ic_tot_mask_fill = np.isfinite(ic_tot_fill[1])
	ic_tot_fill = ic_tot_fill[:,ic_tot_mask_fill]
	ic_tot_fill = np.interp(energy, ic_tot_fill[0], ic_tot_fill[1])*FACTOR_fill[n]
	
	fermi_data_fill = np.loadtxt(galdef_name_fill[n]+'fermi.dat').T
	emin_fill = fermi_data_fill[0]
	emax_fill = fermi_data_fill[1]
	fermi_energy_fill = np.sqrt(emin*emax)
	ferr_d_fill = fermi_data_fill[3]
	ferr_u_fill = fermi_data_fill[4]
	
	#color
	color = cmap(1-float(n+1)/(num_graph+1.))
	ALPHA = 0.6
	
	ax1.loglog(energy,bremss_tot, color=color, alpha=ALPHA, linestyle=':')
	ax1.loglog(energy,bremss_tot_fill, color=color, alpha=ALPHA, linestyle=':')
	
	ax1.loglog(energy,pion_tot, color=color, alpha=ALPHA, linestyle='--')
	ax1.loglog(energy,pion_tot_fill, color=color, alpha=ALPHA, linestyle='--')
	
	ax1.loglog(energy,ic_tot, color=color, alpha=ALPHA, linestyle='-.')
	ax1.loglog(energy,ic_tot_fill, color=color, alpha=ALPHA, linestyle='-.')
	
	ax1.fill_between(energy, gamma_tot, gamma_tot_fill, facecolor=color, edgecolor=color, alpha=0.65, label=LABEL)
	
	ax1.loglog(energy, gamma_tot/1.e6, c=color, linestyle='-', linewidth=LWIDTH, alpha=0.65)
	
	alpha = np.log(gamma_tot[1:]/gamma_tot[:-1])/np.log(energy[1:]/energy[:-1])
	alpha = (alpha[:-1]+alpha[1:])/2.
	e_alpha = energy[1:-1]
	
	e_vals = np.array([1.e3,1.e4,1.e5,1.e6])
	alpha_func = np.interp(e_vals,e_alpha,alpha)
	# print e_alpha, alpha
	print e_vals,alpha_func
	
	# if n==1:
		# ax1.errorbar(fermi_energy, fermi_data[2], xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label=None, c='green', markersize=MSIZE, alpha=0.5)
	
	gamma_ex1 = np.interp(fermi_energy, energy, gamma_tot)
	gamma_er1 = (fermi_data[2] - gamma_ex1)/gamma_ex1
	gamma_error_d1 = (ferr_d + fermi_data[2] - gamma_ex1)/gamma_ex1
	gamma_error_u1 = (ferr_u + fermi_data[2] - gamma_ex1)/gamma_ex1
	
	gamma_ex2 = np.interp(fermi_energy_fill, energy, gamma_tot_fill)
	gamma_er2 = (fermi_data_fill[2] - gamma_ex2)/gamma_ex2
	gamma_error_d2 = (ferr_d_fill + fermi_data_fill[2] - gamma_ex2)/gamma_ex2
	gamma_error_u2 = (ferr_u_fill + fermi_data_fill[2] - gamma_ex2)/gamma_ex2
	
	ax_e.fill_between(fermi_energy, gamma_er1, gamma_er2, facecolor=color, edgecolor=color, alpha=0.65, hatch=HATCH[n])
	
	err1 = np.maximum(gamma_error_u1, gamma_error_u2)
	err2 = np.minimum(gamma_error_d1, gamma_error_d2)
	ax_e.fill_between(fermi_energy, err1, err2, facecolor=color, edgecolor='k', alpha=0.1)
	
	gamma_ex1 = np.interp(verit_data[0], energy, gamma_tot)
	gamma_er1 = (verit_data[1] - gamma_ex1)/gamma_ex1  
	gamma_error_d1 = (verit_data[2]-gamma_ex1)/gamma_ex1
	gamma_error_u1 = (verit_data[3]-gamma_ex1)/gamma_ex1
	
	gamma_ex2 = np.interp(verit_data[0], energy, gamma_tot_fill)
	gamma_er2 = (verit_data[1] - gamma_ex2)/gamma_ex2  
	gamma_error_d2 = (verit_data[2]-gamma_ex2)/gamma_ex2
	gamma_error_u2 = (verit_data[3]-gamma_ex2)/gamma_ex2
	
	ax_e.fill_between(verit_data[0], gamma_er1, gamma_er2, facecolor=color, edgecolor=color, alpha=0.65, hatch=HATCH[n])
	
	err1 = np.maximum(gamma_error_u1, gamma_error_u2)
	err2 = np.minimum(gamma_error_d1, gamma_error_d2)
	ax_e.fill_between(verit_data[0], err1, err2, facecolor=color, edgecolor='k', alpha=0.1)
	
	# ax_e.errorbar(verit_data[0], gamma_er, yerr=[abs(gamma_error_dd), abs(gamma_error_uu)], markersize=MSIZE, c=color, fmt='s', alpha=0.5)


ax_e.plot((1e0,1e15),(0, 0), 'k', linestyle=':')
ax_e.set_ylim([-1., 1])
ax_e.set_xlim([7.e+1,1.e7])
ax_e.set_xscale('log')
ax_e.set_xlabel(r'$E$ [MeV]')
ax_e.set_ylabel(r'$\frac{\Delta E^2 \Phi_E}{E^2 \Phi_E}$')

#PLOTTING DATA
ax1.errorbar(verit_data[0], verit_data[1], yerr=[abs(verit_data[2]-verit_data[1]), abs(verit_data[3]-verit_data[1])], fmt='s', markersize=MSIZE, label=None, c='black', alpha=0.5)
ax1.errorbar(verit_upper[0], verit_upper[1], yerr=verit_upper[1]/1.5, uplims=True, c='black', markersize=MSIZE, alpha=0.5)

ax1.set_ylabel(r'$E^2 \Phi_E$  [MeV cm$^{-2}$ s$^{-1}$]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([7.e+1,1.e7])
ax1.set_ylim([3.e-8, 6.e-6])
# ax1.set_title(r'M82 Integrated $\gamma$-ray Spectrum')	
	
ax1.set_xticklabels([])
ax_e.set_yticklabels(['',-0.5,0.,0.5])

ax1.legend(loc='upper right', prop={'size':15})

plt.tight_layout()
plot_name = 'plot_'+CHANGED+'_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.png', bbox_inches='tight', dpi=400)
		break
#plt.savefig('plot_gamma_'+CHANGED, dpi=400)
plt.close()




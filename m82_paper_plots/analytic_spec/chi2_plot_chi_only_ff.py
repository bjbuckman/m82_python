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
from scipy import interpolate

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

MMODE = 'bp'

## Spectra
USE_DATA = 0
if USE_DATA == 1:
	## Using DATA
	radio_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_radio/klein1988.dat')
	radio_data = radio_data[radio_data[:,0].argsort()]
	radio_data = radio_data.T

	nu = radio_data[0]

	radio_flux = radio_data[1]
	
	alpha = np.zeros(len(nu))
	alpha[1:] += np.log(radio_flux[1:]/radio_flux[:-1])/np.log(nu[1:]/nu[:-1])
	alpha[:-1]+= np.log(radio_flux[1:]/radio_flux[:-1])/np.log(nu[1:]/nu[:-1])
	alpha/= 2.
	# alpha[1:-1] = np.log(radio_flux[2:]/radio_flux[:-2])/np.log(nu[2:]/nu[:-2])
	alpha = 2.*alpha-1.
	
	radio_flux = radio_flux[nu<5.e10]
	alpha = alpha[nu<5.e10]
	nu = nu[nu<5.e10]
	
	radio_flux = radio_flux[nu>1.e8]
	alpha = alpha[nu>1.e8]
	nu = nu[nu>1.e8]
	
	nu/= 1.e9

	#tim
	fermi_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_spectrum_t2_larger_roi.txt').T
	emin = fermi_data[0]
	emax = fermi_data[1]
	fermi_energy = np.sqrt(emin*emax)
	ferr_d = (fermi_data[3])*fermi_energy**2/(emax-emin)
	ferr_u = (fermi_data[3])*fermi_energy**2/(emax-emin)
	fermi_flux = fermi_energy*fermi_data[2]/np.log(emax/emin)
	# ax1.errorbar(fermi_energy, fermi_energy*fermi_data[2]/np.log(emax/emin), xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label=None, c='black', markersize=MSIZE, alpha=0.5)

	energy = fermi_energy
	energy/=1.e3
	
	gamma_tot = fermi_flux
	
	gamma_tot = gamma_tot[energy>1.e3]
	energy = energy[energy>1.e3]
	
	gamma_tot = gamma_tot[energy<1.e7]
	energy = energy[energy<1.e7]
	
	gamma_flux = gamma_tot
	
else:
	## USING GALPROP MODEL
	galdef_dir = '../fits/'
	if MMODE == 'a':
		GG = 'm82_05aff1018'
	elif MMODE == 'b':
		GG = 'm82_05b2ff1020'
	elif MMODE == 'bp':
		GG = 'm82_06bff1021'
	
	galdef_name = galdef_dir+GG+'/'+GG
	
	chi_dat = np.loadtxt(galdef_name+'chi.dat')
	chi_fac = chi_dat[4]
	FACTOR = chi_fac

	#Radio
	synch = np.loadtxt(galdef_name+'synch_flux_na.dat')
	freef = np.loadtxt(galdef_name+'freefree_flux.dat')
	
	flux_synch = synch[1]*FACTOR*1e23
	flux_freef = freef[1]*FACTOR*1e23
	
	radio_flux = flux_synch
	nu = synch[0]
	
	alpha = np.zeros(len(nu))
	alpha[1:] += np.log(radio_flux[1:]/radio_flux[:-1])/np.log(nu[1:]/nu[:-1])
	alpha[:-1]+= np.log(radio_flux[1:]/radio_flux[:-1])/np.log(nu[1:]/nu[:-1])
	alpha/= 2.
	alpha = 2.*alpha-1.
	
	radio_flux = radio_flux[nu<5.e10]
	alpha = alpha[nu<5.e10]
	nu = nu[nu<5.e10]
	
	radio_flux = radio_flux[nu>1.e9]
	alpha = alpha[nu>1.e9]
	nu = nu[nu>1.e9]
	
	nu/= 1.e9
	
	#Gamma-rays
	gamma_tot = np.loadtxt(galdef_name+'gamma_flux.dat')
	energy = gamma_tot[0]
	
	gamma_tot_mask = np.isfinite(gamma_tot[1])
	gamma_tot = gamma_tot[:,gamma_tot_mask]
	gamma_tot = np.interp(energy, gamma_tot[0], gamma_tot[1])*FACTOR
	
	gamma_tot = gamma_tot[energy>1.e3]
	energy = energy[energy>1.e3]
	
	gamma_tot = gamma_tot[energy<1.e7]
	energy = energy[energy<1.e7]
	
	energy/= 1.e3
	gamma_flux = gamma_tot
	
radio_flux/= 1.
gamma_flux/= 1.e-6


numnu = len(nu)
numga = len(energy)

B = np.linspace(10., 1000., num=100)
B/= 100.
numbb = len(B)

def E_nu(nu, B):
	E = (nu/B)**0.5*1.6556
	return E
	
def F_1_calc(E):
	F_0 = -E_max**(alpha_0+1.)
	F = (alpha_0+1.)/(1.+(1./E)**(alpha_0+1.)*F_0)
	return F

def n_B(E, B, a, U, V, H):
	F_1 = F_1_calc(E)
	n = (E**2.*(B**2./t_s+U/t_i)*(1.-(F_1-(a+1.))) -E*V/H/t_w*(F_1-(a+1.)) )*t_o
	n/= 1.+(F_1-(a+1.))*(1.+E*t_o/t_b)
	if n < 0.:
		n = np.nan
	return n

def tau_synch(E, B):
	tau_s = t_s/B**2./E
	return tau_s
	
def tau_elect(n, E, B, U, V, H):
	t_inv = n*(1./E/t_o + 1./t_b) +E*(B**2./t_s+U/t_i)+V/H/t_w
	return 1./t_inv
	
def g_synch(a):
	g = 27.*3.**0.5/(p-a+1.)*2.**(-3.*(5.+a-p)/2.)*np.pi**((2.+a-p)/2.)
	g*= gamma((5.+(p-a))/4.)*gamma((3.*(p-a)-1.)/12.)
	g*= gamma((19.+3.*(p-a))/12.)/gamma((7.+(p-a))/4.)
	return g

######################################
## PLOT DATA
##

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

cmap = plt.get_cmap('cubehelix_r')

fig = plt.figure()
f1, ax = plt.subplots(1,1)

for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(2)
ax.tick_params(which='major',width=1., length=8, labelsize='small', direction='in', color='black')
ax.tick_params(which='minor',width=0.5, length=5, direction='in', color='black')

MSIZE = 5
LWIDTH = 3
WIDTH = 3


###
###
### SPECTRAL JUNK
###
###
print numnu
if MMODE == 'a':
	LABEL = [r'Standard', r'$U_\mathrm{ISRF} / 2$', r'$V_0\times 10$', r'$p_0+0.1$']
	wind = [100., 1000.]
	ax.text(1.e3, 3.e2, r'$V_0 = 100$ km s$^{-1}$', color='k')
	ax.text(1.5e2, 4.e3, 'Model A', color='k')
elif MMODE == 'b':
	LABEL = [r'Standard', r'$U_\mathrm{ISRF} / 2$', r'$V_0 / 10$', r'$p_0+0.1$']
	wind = [1000., 100.]
	ax.text(1.e3, 3.e2, r'$V_0 = 1000$ km s$^{-1}$', color='k')
	ax.text(1.5e2, 4.e3, 'Model B', color='k')
elif MMODE == 'bp':
	LABEL = [r'Standard', r'$U_\mathrm{ISRF} / 2$', r'$V_0 / 10$', r'$p_0+0.1$']
	wind = [1000., 100.]
	ax.text(1.e3, 3.e2, r'$V_0 = 1000$ km s$^{-1}$', color='k')
	ax.text(1.5e2, 4.e3, "Model B'", color='k')

Uisrf = [1000., 500.]
aaaa = [-2.2, -2.1]
color1 = [cmap(0.75), cmap(0.60), cmap(0.45), cmap(0.30)]

for ii in range(0,numnu):
	B = np.linspace(1.0e1,1.e4,num=4000)/100.
	B_1 = np.linspace(1.0e1,1.e4,num=4000)/100.

	U = Uisrf[0]
	U/= 1000.
	
	V = wind[0]
	V/= 100.
	
	H = 0.1
	H/= 0.1

	alpha_0 = aaaa[0]
	AL = alpha[ii]
	
	EN = 1.6556*(nu[ii]/B)**0.5
	EN_1 = 1.6556*(nu[ii]/B_1)**0.5
	E_max = 1.e+6

	t_b = 3.17e5
	t_s = 1.3e6
	t_i = 3.15e5
	t_o = 2.22e6
	t_o_1 = 2.22e6 #7.
	t_o_2 = 5.46e5 #2.8#3.
	t_w = 9.8e5

	F_0 = -E_max**(alpha_0+1.)
	def F_1_calc(E):
		return (alpha_0+1.)/(1.+(1./E)**(alpha_0+1.)*F_0)

	energy_losses = 0
	fac = 1.0
	for energy_losses in (3,4):
		if energy_losses == 0:
			#brems + synch
			F_1 = F_1_calc(B)
			n_0 = ((F_1 - (AL+1))**(-1.) -1)*B**2/t_s*t_b*(nu/B)**(.5)
			n_0*=100*fac
			# beta_0 = 1.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 1:
			#brems + synch + IC
			F_1 = F_1_calc(B)
			n_1 = ((F_1 - (AL+1))**(-1.) -1)*(B**2/t_s + U/t_i)*t_b*(nu/B)**(.5)
			n_1*=100*fac
			# beta_1 = 2./(1.+U/B**2.) - 0.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 2:
			#ioniz + synch
			F_1 = F_1_calc(B)
			n_2 = (2.*(F_1-(AL+1)+1)**(-1.)-1)*B**2.*(nu/B)*t_o/t_s
			n_2*=100*fac
			# beta_2 = 1. + F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1)+1)
		elif energy_losses == 3:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN)
			n_3 = (EN**2.*(B**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN*V/H/t_w*(F_1-(AL+1.)) )*t_o_1/(1.+(F_1-(AL+1.))*(1.+EN*t_o_1/t_b))
			n_3*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.)-1) + (1.+F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.))*(B/nu*t_b/t_o_1 -1.))/(1.-B/nu*t_b/t_o_1*(F_1-(AL+1.)+1.)/(F_1-(AL+1.)))
		elif energy_losses == 4:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN_1)
			n_4 = (EN_1**2.*(B_1**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN_1*V/H/t_w*(F_1-(AL+1.)) )*t_o_2/(1.+(F_1-(AL+1.))*(1.+EN_1*t_o_2/t_b))
			n_4*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1)-1) + (1+F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1))*(B_1/nu*t_b/t_o_2 -1))/(1-B_1/nu*t_b/t_o_2*(F_1-(AL+1)+1)/(F_1-(AL+1)))

	B*= 100.
	B_1*=100.

	# ax.plot(B, n_0, c=cmap(0.2), label='bremss+synch')
	# ax.plot(B, n_1, c=cmap(0.4), label='bremss+synch+IC')
	# ax.plot(B, n_2, c=cmap(0.6), label='ionize+synch')
	ax.plot(B, (n_3), c=color1[0], alpha=(0.9*(ii+2)/(numnu+2)), linewidth=LWIDTH, linestyle='--')
	# ax.plot(B_1, (n_4), c=color1[0], alpha=(0.8*(ii+2)/(numnu+1)), linewidth=LWIDTH, linestyle=':')

	# print B
	# print n_3
	# print n_4
		
for ii in range(0,numnu):
	B = np.linspace(1.0e1,1.e4,num=4000)/100.
	B_1 = np.linspace(1.0e1,1.e4,num=4000)/100.

	U = Uisrf[1]
	U/= 1000.
	
	V = wind[0]
	V/= 100.
	
	H = 0.1
	H/= 0.1

	alpha_0 = aaaa[0]
	AL = alpha[ii]
	
	EN = 1.6556*(nu[ii]/B)**0.5
	EN_1 = 1.6556*(nu[ii]/B_1)**0.5
	E_max = 1.e+6

	t_b = 3.17e5
	t_s = 1.3e6
	t_i = 3.15e5
	t_o = 2.22e6
	t_o_1 = 2.22e6 #7.
	t_o_2 = 5.46e5 #2.8#3.
	t_w = 9.8e5

	F_0 = -E_max**(alpha_0+1.)
	def F_1_calc(E):
		return (alpha_0+1.)/(1.+(1./E)**(alpha_0+1.)*F_0)

	energy_losses = 0
	fac = 1.0
	for energy_losses in (3,4):
		if energy_losses == 0:
			#brems + synch
			F_1 = F_1_calc(B)
			n_0 = ((F_1 - (AL+1))**(-1.) -1)*B**2/t_s*t_b*(nu/B)**(.5)
			n_0*=100*fac
			# beta_0 = 1.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 1:
			#brems + synch + IC
			F_1 = F_1_calc(B)
			n_1 = ((F_1 - (AL+1))**(-1.) -1)*(B**2/t_s + U/t_i)*t_b*(nu/B)**(.5)
			n_1*=100*fac
			# beta_1 = 2./(1.+U/B**2.) - 0.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 2:
			#ioniz + synch
			F_1 = F_1_calc(B)
			n_2 = (2.*(F_1-(AL+1)+1)**(-1.)-1)*B**2.*(nu/B)*t_o/t_s
			n_2*=100*fac
			# beta_2 = 1. + F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1)+1)
		elif energy_losses == 3:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN)
			n_3 = (EN**2.*(B**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN*V/H/t_w*(F_1-(AL+1.)) )*t_o_1/(1.+(F_1-(AL+1.))*(1.+EN*t_o_1/t_b))
			n_3*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.)-1) + (1.+F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.))*(B/nu*t_b/t_o_1 -1.))/(1.-B/nu*t_b/t_o_1*(F_1-(AL+1.)+1.)/(F_1-(AL+1.)))
		elif energy_losses == 4:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN_1)
			n_4 = (EN_1**2.*(B_1**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN_1*V/H/t_w*(F_1-(AL+1.)) )*t_o_2/(1.+(F_1-(AL+1.))*(1.+EN_1*t_o_2/t_b))
			n_4*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1)-1) + (1+F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1))*(B_1/nu*t_b/t_o_2 -1))/(1-B_1/nu*t_b/t_o_2*(F_1-(AL+1)+1)/(F_1-(AL+1)))

	B*= 100.
	B_1*=100.

	# ax.plot(B, n_0, c=cmap(0.2), label='bremss+synch')
	# ax.plot(B, n_1, c=cmap(0.4), label='bremss+synch+IC')
	# ax.plot(B, n_2, c=cmap(0.6), label='ionize+synch')
	ax.plot(B, (n_3), c=color1[1], alpha=(0.8*(ii+2)/(numnu+1)), linewidth=LWIDTH, linestyle='--')
	# ax.plot(B_1, (n_4), c=color1[1], alpha=(0.8*(ii+2)/(numnu+1)), linewidth=LWIDTH, linestyle=':')

	# print B
	# print n_3
	# print n_4

for ii in range(0,numnu):
	B = np.linspace(1.0e1,1.e4,num=4000)/100.
	B_1 = np.linspace(1.0e1,1.e4,num=4000)/100.

	U = Uisrf[0]
	U/= 1000.
	
	V = wind[1]
	V/= 100.
	
	H = 0.1
	H/= 0.1

	alpha_0 = aaaa[0]
	AL = alpha[ii]
	
	EN = 1.6556*(nu[ii]/B)**0.5
	EN_1 = 1.6556*(nu[ii]/B_1)**0.5
	E_max = 1.e+6

	t_b = 3.17e5
	t_s = 1.3e6
	t_i = 3.15e5
	t_o = 2.22e6
	t_o_1 = 2.22e6 #7.
	t_o_2 = 5.46e5 #2.8#3.
	t_w = 9.8e5

	F_0 = -E_max**(alpha_0+1.)
	def F_1_calc(E):
		return (alpha_0+1.)/(1.+(1./E)**(alpha_0+1.)*F_0)

	energy_losses = 0
	fac = 1.0
	for energy_losses in (3,4):
		if energy_losses == 0:
			#brems + synch
			F_1 = F_1_calc(B)
			n_0 = ((F_1 - (AL+1))**(-1.) -1)*B**2/t_s*t_b*(nu/B)**(.5)
			n_0*=100*fac
			# beta_0 = 1.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 1:
			#brems + synch + IC
			F_1 = F_1_calc(B)
			n_1 = ((F_1 - (AL+1))**(-1.) -1)*(B**2/t_s + U/t_i)*t_b*(nu/B)**(.5)
			n_1*=100*fac
			# beta_1 = 2./(1.+U/B**2.) - 0.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 2:
			#ioniz + synch
			F_1 = F_1_calc(B)
			n_2 = (2.*(F_1-(AL+1)+1)**(-1.)-1)*B**2.*(nu/B)*t_o/t_s
			n_2*=100*fac
			# beta_2 = 1. + F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1)+1)
		elif energy_losses == 3:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN)
			n_3 = (EN**2.*(B**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN*V/H/t_w*(F_1-(AL+1.)) )*t_o_1/(1.+(F_1-(AL+1.))*(1.+EN*t_o_1/t_b))
			n_3*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.)-1) + (1.+F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.))*(B/nu*t_b/t_o_1 -1.))/(1.-B/nu*t_b/t_o_1*(F_1-(AL+1.)+1.)/(F_1-(AL+1.)))
		elif energy_losses == 4:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN_1)
			n_4 = (EN_1**2.*(B_1**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN_1*V/H/t_w*(F_1-(AL+1.)) )*t_o_2/(1.+(F_1-(AL+1.))*(1.+EN_1*t_o_2/t_b))
			n_4*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1)-1) + (1+F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1))*(B_1/nu*t_b/t_o_2 -1))/(1-B_1/nu*t_b/t_o_2*(F_1-(AL+1)+1)/(F_1-(AL+1)))

	B*= 100.
	B_1*=100.

	# ax.plot(B, n_0, c=cmap(0.2), label='bremss+synch')
	# ax.plot(B, n_1, c=cmap(0.4), label='bremss+synch+IC')
	# ax.plot(B, n_2, c=cmap(0.6), label='ionize+synch')
	ax.plot(B, (n_3), c=color1[2], alpha=(0.8*(ii+2)/(numnu+1)), linewidth=LWIDTH, linestyle='--')
	# ax.plot(B_1, (n_4), c=color1[2], alpha=(0.8*(ii+2)/(numnu+1)), linewidth=LWIDTH, linestyle=':')

	# print B
	# print n_3
	# print n_4
	
for ii in range(0,numnu):
	B = np.linspace(1.0e1,1.e4,num=4000)/100.
	B_1 = np.linspace(1.0e1,1.e4,num=4000)/100.

	U = Uisrf[0]
	U/= 1000.
	
	V = wind[0]
	V/= 100.
	
	H = 0.1
	H/= 0.1

	alpha_0 = aaaa[1]
	AL = alpha[ii]
	
	EN = 1.6556*(nu[ii]/B)**0.5
	EN_1 = 1.6556*(nu[ii]/B_1)**0.5
	E_max = 1.e+6

	t_b = 3.17e5
	t_s = 1.3e6
	t_i = 3.15e5
	t_o = 2.22e6
	t_o_1 = 2.22e6 #7.
	t_o_2 = 5.46e5 #2.8#3.
	t_w = 9.8e5

	F_0 = -E_max**(alpha_0+1.)
	def F_1_calc(E):
		return (alpha_0+1.)/(1.+(1./E)**(alpha_0+1.)*F_0)

	energy_losses = 0
	fac = 1.0
	for energy_losses in (3,4):
		if energy_losses == 0:
			#brems + synch
			F_1 = F_1_calc(B)
			n_0 = ((F_1 - (AL+1))**(-1.) -1)*B**2/t_s*t_b*(nu/B)**(.5)
			n_0*=100*fac
			# beta_0 = 1.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 1:
			#brems + synch + IC
			F_1 = F_1_calc(B)
			n_1 = ((F_1 - (AL+1))**(-1.) -1)*(B**2/t_s + U/t_i)*t_b*(nu/B)**(.5)
			n_1*=100*fac
			# beta_1 = 2./(1.+U/B**2.) - 0.5 + 0.5*F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1))
		elif energy_losses == 2:
			#ioniz + synch
			F_1 = F_1_calc(B)
			n_2 = (2.*(F_1-(AL+1)+1)**(-1.)-1)*B**2.*(nu/B)*t_o/t_s
			n_2*=100*fac
			# beta_2 = 1. + F_1*(F_1-(alpha_0+1))/(F_1-(AL+1)-1)/(F_1-(AL+1)+1)
		elif energy_losses == 3:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN)
			n_3 = (EN**2.*(B**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN*V/H/t_w*(F_1-(AL+1.)) )*t_o_1/(1.+(F_1-(AL+1.))*(1.+EN*t_o_1/t_b))
			n_3*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.)-1) + (1.+F_1/2.*(F_1-(alpha_0+1.))/(F_1-(AL+1.))*(B/nu*t_b/t_o_1 -1.))/(1.-B/nu*t_b/t_o_1*(F_1-(AL+1.)+1.)/(F_1-(AL+1.)))
		elif energy_losses == 4:
			#ioniz + brems + synch + IC
			F_1 = F_1_calc(EN_1)
			n_4 = (EN_1**2.*(B_1**2./t_s+U/t_i)*(1.-(F_1-(AL+1.))) -EN_1*V/H/t_w*(F_1-(AL+1.)) )*t_o_2/(1.+(F_1-(AL+1.))*(1.+EN_1*t_o_2/t_b))
			n_4*=100*fac
			# beta_3 = -1. + 2./(1.+U/B**2.) + F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1)-1) + (1+F_1/2*(F_1-(alpha_0+1))/(F_1-(AL+1))*(B_1/nu*t_b/t_o_2 -1))/(1-B_1/nu*t_b/t_o_2*(F_1-(AL+1)+1)/(F_1-(AL+1)))

	B*= 100.
	B_1*=100.

	# ax.plot(B, n_0, c=cmap(0.2), label='bremss+synch')
	# ax.plot(B, n_1, c=cmap(0.4), label='bremss+synch+IC')
	# ax.plot(B, n_2, c=cmap(0.6), label='ionize+synch')
	ax.plot(B, (n_3), c=color1[3], alpha=(0.8*(ii+2)/(numnu+1)), linewidth=LWIDTH, linestyle='--')
	# ax.plot(B_1, (n_4), c=color1[3], alpha=(0.8*(ii+2)/(numnu+1)), linewidth=LWIDTH, linestyle=':')

	# print B
	# print n_3
	# print n_4

for ii in range(0,4):
	ax.plot(B, n_3*1.e-10, c=color1[ii], alpha=0.8, linewidth=LWIDTH, linestyle='--', label=LABEL[ii])
ax.legend(loc='lower right', prop={'size':15})
	
#axes stuff
ax.set_axis_bgcolor('white')
ax.set_aspect('equal')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([9.e1,1.e4])
ax.set_ylim([5.e1,1.e4])
ax.set_xlabel(r'$B_0$ [$\mu$G]')
ax.set_ylabel(r'$n_0$ [cm$^{-3}$]')

ax.plot(150,150,'o',color=cmap(0.25), markersize=15, marker='h', linewidth=4, markeredgecolor='red', markeredgewidth=2)#, edgecolor='black')
ax.plot(325,675,'o',color=cmap(0.5), markersize=15, marker='h', linewidth=4, markeredgecolor='red', markeredgewidth=2)#, edgecolor='black') 
ax.plot(325,1000,'o',color=cmap(0.75), markersize=15, marker='h', linewidth=4, markeredgecolor='red', markeredgewidth=2)#, edgecolor='black') 


plot_name = 'm82_chi_analytic_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.png', bbox_inches='tight', dpi=400)
		break

plt.close()	





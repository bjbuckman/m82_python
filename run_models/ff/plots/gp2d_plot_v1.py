import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import os

SL_SCALE = 0

E2FLAG = 1

DISTANCE = 3520

dir_in = '../gp_out'
GALDEF = 'm82_06'
CHANGED = 'dxx'
num_params = 1
num_galdef = 15
num_galdef_tot = int(np.prod(num_galdef))

VAR_PLOT = '$D_{xx}=$ '

BEGIN = 10
num_start = BEGIN

num_graph_start = 1
num_graph_galdef = 1
num_graph_galdef_tot = int(np.prod(num_graph_galdef))

chi = np.zeros(num_graph_galdef_tot)
chi_fac = np.zeros(num_graph_galdef_tot)
chi_ff_fac = np.zeros(num_graph_galdef_tot)
N=np.zeros(num_params).astype(int)
NN=np.zeros([num_graph_galdef_tot, num_params]).astype(int)
galdef_name = []
for n in range(0,num_graph_galdef_tot):
	#Getting suffix for files
	NN[n] = N
	N_index=N+num_graph_start
	N_suffix = ''
	for s in range(0,num_params):
		N_suffix += str(N_index[s])
	
	GNAME = GALDEF+CHANGED+N_suffix
	galdef_name.append(GNAME)
	
	chi_dat = np.loadtxt(dir_in+GNAME+'/'+GNAME+'chi.dat')
	chi_fac[n] = chi_dat[4]
	# chi_ff_fac[n] = chi_dat[13]
	# print chi_dat[13]
	GAMMA_CHI = 20
	chi[n] = GAMMA_CHI*(chi_dat[1]+chi_dat[2])+chi_dat[3]
	
	#Update N
	N[num_params-1] += 1
	for i in (num_params-1 -np.array(range(0,num_params))):
		if N[i] >= num_graph_galdef[i]:
			N[i] = 0
			if i >= 1:
				N[i-1] += 1
				

new_val = 1.5**np.array(range(-11,-11+num_galdef))*8.32e28

#FACTOR=1#.e4
# chi = np.loadtxt(GALDEF+CHANGED+'chi.dat')

NORM_PLOT = 1
if NORM_PLOT == 1:
	FACTOR = np.ones(num_graph_galdef_tot)
else:
	FACTOR = chi_fac

print chi_fac
print chi_ff_fac	
	
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

fig = plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(6, 3, width_ratios=[1, 1, 1]) 
gs.update(wspace=0.0, hspace=0.0)

cmap = plt.get_cmap('jet')

plt.subplots_adjust(hspace=0.001)

ax2 = plt.subplot(gs[0:2,0])
ax2x = plt.subplot(gs[2:3,0], sharex = ax2)
ax1 = plt.subplot(gs[3:5,0])
ax1x = plt.subplot(gs[5:6,0], sharex = ax1)

ax01 = plt.subplot(gs[0:2,1])
ax02 = plt.subplot(gs[0:2,2])
ax11 = plt.subplot(gs[2:4,1])
ax12 = plt.subplot(gs[2:4,2])
ax21 = plt.subplot(gs[4:6,1])
ax22 = plt.subplot(gs[4:6,2])

AX92 = ax21
AX92_o = ax22
AX22 = ax11
AX22_o = ax12
AX6 = ax02
AX3 = ax01

###
###
### GAMMA-RAY
###
###
###

fermi_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_fermi.dat').T
verit_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_veritas.dat').T
verit_upper = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_veritas_upper.dat')
#verit_data[0] *= 1000
verit_data[1:] *= verit_data[0]
verit_upper[1] *= verit_upper[0]

alpha = -2.21
emin = fermi_data[1]
emax = fermi_data[2]
ferr_d = fermi_data[3]*(alpha+1)/(emax**(alpha+1)-emin**(alpha+1))*(emin*emax)**((alpha+2)/2)
ferr_u = fermi_data[4]*(alpha+1)/(emax**(alpha+1)-emin**(alpha+1))*(emin*emax)**((alpha+2)/2)

fermi_energy = np.zeros(fermi_data.shape[1])
fermi_energy = np.sqrt(emin*emax)

e2dnde = fermi_data[0]*(alpha+1)/(emax**(alpha+1)-emin**(alpha+1))*(emin*emax)**((alpha+2)/2)
#print fermi_energy
#print e2dnde

#PLOTTING DATA
#ax1.errorbar(fermi_energy, e2dnde, xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label=r'FERMI $\alpha=-2.2$', c='black')
#plt.errorbar(fermi_energy, fermi_data[0]*fermi_energy, xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(fermi_data[3])*fermi_energy, fermi_data[4]*fermi_energy],fmt='x', c='black')
ax1.errorbar(verit_data[0], verit_data[1], yerr=[abs(verit_data[2]-verit_data[1]), abs(verit_data[3]-verit_data[1])], fmt='s', markersize=3, label='VERITAS', c='black')
ax1.errorbar(verit_upper[0], verit_upper[1], yerr=verit_upper[1]/1.5, uplims=True, c='black', markersize=3)
#ax1.set_title(r'$\gamma$-ray Flux for '+str(CHANGED))
ax1.set_xlabel(r'E (MeV)')
ax1.set_ylabel(r'$E \, F_{\gamma}(E)$  (MeV$^2$ s$^{-1}$ cm$^{-2}$ MeV$^{-1}$)')
ax1.set_xlim([1.e+0,1.e7])
ax1.set_ylim([1.e-10, 1e-4])

ax1x.plot((1e0,1e15),(0, 0), 'k', linestyle=':')
ax1x.set_ylim([-1., 1])

###
###
### RADIO
###
###
###

radio_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_radio/radio_m82.dat')
radio_data = radio_data.T

#PLOTTING DATA
ax2.errorbar(radio_data[0], radio_data[1], yerr=radio_data[2], fmt='o', markersize=3, c='black')
# ax2.set_title('Radio Flux for '+str(CHANGED))
# ax2.set_xlabel(r'$\nu$ (Hz)')
ax2.set_ylabel(r'$F_{radio}(\nu)$ (Jy)')
ax2.set_xlim([1.e+6,1.e+15])
ax2.set_ylim([1.e-3, 1.e2])

c=3e10

r_lambda=92
r_freq=c/r_lambda
ax2.plot((r_freq,r_freq),(1e-30, 1e30), 'k', linestyle=':')

r_lambda=22
r_freq=c/r_lambda
ax2.plot((r_freq,r_freq),(1e-30, 1e30), 'k', linestyle=':')

r_lambda=6
r_freq=c/r_lambda
ax2.plot((r_freq,r_freq),(1e-30, 1e30), 'k', linestyle=':')

r_lambda=3
r_freq=c/r_lambda
ax2.plot((r_freq,r_freq),(1e-30, 1e30), 'k', linestyle=':')

ax2x.plot((1e0,1e15),(0, 0), 'k', linestyle=':')

ax2x.set_xlabel(r'$\nu$ (Hz)')
# ax2x.set_ylabel(r'$F_{data}/F_{radio}$')
ax2x.set_ylim([-1., 1])

###
###
### RADIO FUNCTION OF Z
###
###
###
	
## 92cm
L92cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_N.dat').T
L92cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_S.dat').T


AX92.errorbar(L92cm_N[0],L92cm_N[1],yerr=[abs(L92cm_N[2]-L92cm_N[1]),abs(L92cm_N[3]-L92cm_N[1])],fmt='o',markersize=3,color='black')
AX92.errorbar(L92cm_S[0],L92cm_S[1],yerr=[abs(L92cm_S[2]-L92cm_S[1]),abs(L92cm_S[3]-L92cm_S[1])],fmt='o',markersize=3,color='grey')
AX92.set_yscale('log')
AX92.set_xlim([0,3.5])
AX92.set_ylim([1e-4,1e+2])
# AX22.set_xlabel('kpc')
# AX22.set_ylabel('Jy')
AX92.text(0.75, 0.8, r'$\lambda = 92$ cm', ha='center', va='center', transform=AX92.transAxes)

## 92cm OFFSET
L92cm_NW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_NW.dat').T
L92cm_SW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_SW.dat').T
L92cm_NE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_NE.dat').T
L92cm_SE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_SE.dat').T


AX92_o.errorbar(L92cm_NW[0],L92cm_NW[1],yerr=[abs(L92cm_NW[2]-L92cm_NW[1]),abs(L92cm_NW[3]-L92cm_NW[1])],fmt='o',markersize=3,color='black')
AX92_o.errorbar(L92cm_SW[0],L92cm_SW[1],yerr=[abs(L92cm_SW[2]-L92cm_SW[1]),abs(L92cm_SW[3]-L92cm_SW[1])],fmt='o',markersize=3,color='grey')
AX92_o.errorbar(L92cm_NE[0],L92cm_NE[1],yerr=[abs(L92cm_NE[2]-L92cm_NE[1]),abs(L92cm_NE[3]-L92cm_NE[1])],fmt='o',markersize=3,color='black')
AX92_o.errorbar(L92cm_SE[0],L92cm_SE[1],yerr=[abs(L92cm_SE[2]-L92cm_SE[1]),abs(L92cm_SE[3]-L92cm_SE[1])],fmt='o',markersize=3,color='grey')
AX92_o.set_yscale('log')
AX92_o.set_xlim([0,3.5])
AX92_o.set_ylim([1e-4,1e+2])
# AX22.set_xlabel('kpc')
# AX22.set_ylabel('Jy')
AX92_o.text(0.75, 0.8, r'$\lambda = 92$ cm''\n''$60$" offset', ha='center', va='center', transform=AX92_o.transAxes)

## 22cm
L22cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_N.dat').T
L22cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_S.dat').T


AX22.errorbar(L22cm_N[0],L22cm_N[1],yerr=[abs(L22cm_N[2]-L22cm_N[1]),abs(L22cm_N[3]-L22cm_N[1])],fmt='o',markersize=3,color='black')
AX22.errorbar(L22cm_S[0],L22cm_S[1],yerr=[abs(L22cm_S[2]-L22cm_S[1]),abs(L22cm_S[3]-L22cm_S[1])],fmt='o',markersize=3,color='grey')
AX22.set_yscale('log')
AX22.set_xlim([0,3.5])
AX22.set_ylim([1e-4,1e+2])
# AX22.set_xlabel('kpc')
# AX22.set_ylabel('Jy')
AX22.text(0.75, 0.8, r'$\lambda = 22$ cm', ha='center', va='center', transform=AX22.transAxes)

## 22cm OFFSET
L22cm_NW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_NW.dat').T
L22cm_SW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_SW.dat').T
L22cm_NE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_NE.dat').T
L22cm_SE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_SE.dat').T


AX22_o.errorbar(L22cm_NW[0],L22cm_NW[1],yerr=[abs(L22cm_NW[2]-L22cm_NW[1]),abs(L22cm_NW[3]-L22cm_NW[1])],fmt='o',markersize=3,color='black')
AX22_o.errorbar(L22cm_SW[0],L22cm_SW[1],yerr=[abs(L22cm_SW[2]-L22cm_SW[1]),abs(L22cm_SW[3]-L22cm_SW[1])],fmt='o',markersize=3,color='grey')
AX22_o.errorbar(L22cm_NE[0],L22cm_NE[1],yerr=[abs(L22cm_NE[2]-L22cm_NE[1]),abs(L22cm_NE[3]-L22cm_NE[1])],fmt='o',markersize=3,color='black')
AX22_o.errorbar(L22cm_SE[0],L22cm_SE[1],yerr=[abs(L22cm_SE[2]-L22cm_SE[1]),abs(L22cm_SE[3]-L22cm_SE[1])],fmt='o',markersize=3,color='grey')
AX22_o.set_yscale('log')
AX22_o.set_xlim([0,3.5])
AX22_o.set_ylim([1e-4,1e+2])
# AX22.set_xlabel('kpc')
# AX22.set_ylabel('Jy')
AX22_o.text(0.75, 0.8, r'$\lambda = 22$ cm''\n''$60$" offset', ha='center', va='center', transform=AX22_o.transAxes)


## 6cm
L6cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_N.dat').T
L6cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_S.dat').T

AX6.errorbar(L6cm_N[0],L6cm_N[1],yerr=[abs(L6cm_N[2]-L6cm_N[1]),abs(L6cm_N[3]-L6cm_N[1])],fmt='o',markersize=3,color='black')
AX6.errorbar(L6cm_S[0],L6cm_S[1],yerr=[abs(L6cm_S[2]-L6cm_S[1]),abs(L6cm_S[3]-L6cm_S[1])],fmt='o',markersize=3,color='grey')
AX6.set_yscale('log')
AX6.set_xlim([0,3.5])
AX6.set_ylim([1e-4,1e+2])
# AX6.set_xlabel('kpc')
# AX6.set_ylabel('Jy')
AX6.text(0.75, 0.8, r'$\lambda = 6$ cm', ha='center', va='center', transform=AX6.transAxes)

## 3cm
L3cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_N.dat').T
L3cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_S.dat').T

AX3.errorbar(L3cm_N[0],L3cm_N[1],yerr=[abs(L3cm_N[2]-L3cm_N[1]),abs(L3cm_N[3]-L3cm_N[1])],fmt='o',markersize=3,color='black')
AX3.errorbar(L3cm_S[0],L3cm_S[1],yerr=[abs(L3cm_S[2]-L3cm_S[1]),abs(L3cm_S[3]-L3cm_S[1])],fmt='o',markersize=3,color='grey')
AX3.set_yscale('log')
AX3.set_xlim([0,3.5])
AX3.set_ylim([1e-4,1e+2])
# AX3.set_xlabel('kpc')
# aAX3.set_ylabel('Jy')
AX3.text(0.75, 0.8, r'$\lambda = 3$ cm', ha='center', va='center', transform=AX3.transAxes)


## FREE-FREE
free_free_dir = '/mnt/c/Users/psyko/Physics/proj_gal/1712/bv/plots/ff_fit/m82_16bhf101312/m82_16bhf101312'
free_free_file = '/mnt/c/Users/psyko/Physics/proj_gal/1712/bv/plots/ff_fit/m82_16bhf101312/m82_16bhf101312freefree_flux.dat'
# freef = np.loadtxt(free_free_file)
# freefL = np.loadtxt(free_free_file)
	
for n in range(0,num_graph_galdef_tot):
	
	LABEL = ''
	for m in range(0,num_params):
		LABEL += VAR_PLOT[m]+str('{:.2e}'.format(float(new_val[m][NN[n][m]])))+' '
	LABEL += '$\chi^2=$'+str('{:.1e}'.format(float(chi[n])))
	
	###
	###
	### GAMMA-RAY
	###
	###
	###
	
	gamma_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'gamma_flux.dat')
	gamma_tot_mask = np.isfinite(gamma_tot[1])
	
	if E2FLAG == 0:
		gamma_tot[1] /= gamma_tot[0]
	
	bremss_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'bremss_flux.dat')
	pion_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'pion_flux.dat')
	ic_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'IC_flux.dat')
	color = cmap(float(n)/num_graph_galdef_tot)
	
	# print galdef_name[n]
	# print ic_tot.shape
	
	fermi_data = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'fermi.dat').T
	emin = fermi_data[0]
	emax = fermi_data[1]
	fermi_energy = np.sqrt(emin*emax)
	ferr_d = fermi_data[3]
	ferr_u = fermi_data[4]
	
	#print fermi_data[2]
	
	LABEL = ''
	for m in range(0,num_params):
		LABEL += VAR_PLOT[m]+str('{:.2e}'.format(float(new_val[m][NN[n][m]])))+' '
	LABEL += '$\chi^2=$'+str('{:.1e}'.format(float(chi[n])))
	
	ax1.loglog(gamma_tot[0][gamma_tot_mask], gamma_tot[1][gamma_tot_mask]*FACTOR[n], label=None, c=color, linestyle='-')
	ax1.loglog(bremss_tot[0], bremss_tot[1]*FACTOR[n], label=None, c=color, linestyle=':')
	ax1.loglog(pion_tot[0], pion_tot[1]*FACTOR[n], label=None, c=color, linestyle='--')
	ax1.loglog(ic_tot[0], ic_tot[1]*FACTOR[n], label=None, c=color, linestyle='-.')
	
	ax1.errorbar(fermi_energy, fermi_data[2], xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label='FERMI', c='black', markersize=3)
	
	gamma_ex = np.interp(fermi_energy, gamma_tot[0][gamma_tot_mask], gamma_tot[1][gamma_tot_mask]*FACTOR[n])
	gamma_er = (fermi_data[2] - gamma_ex)/gamma_ex  
	gamma_error_dd = ferr_d/gamma_ex
	gamma_error_uu = ferr_u/gamma_ex
	ax1x.errorbar(fermi_energy, gamma_er, xerr=[abs(emin-fermi_energy),emax-fermi_energy],yerr=[abs(gamma_error_dd), gamma_error_uu], fmt='o', markersize=3, c=color)
	
	gamma_ex = np.interp(verit_data[0], gamma_tot[0][gamma_tot_mask], gamma_tot[1][gamma_tot_mask]*FACTOR[n])
	gamma_er = (verit_data[1] - gamma_ex)/gamma_ex  
	gamma_error_dd = (verit_data[2]-verit_data[1])/gamma_ex
	gamma_error_uu = (verit_data[3]-verit_data[1])/gamma_ex
	ax1x.errorbar(verit_data[0], gamma_er, yerr=[abs(gamma_error_dd), abs(gamma_error_uu)], fmt='o', markersize=3, c=color)
	
	# yerr=[abs(verit_data[2]-verit_data[1]), abs(verit_data[3]-verit_data[1])]
	
	###
	###
	### RADIO
	###
	###
	###
	
	if SL_SCALE == 1:
		synch = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux.dat')
		freef = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux.dat')
		# synchS = np.loadtxt(galdef_name[n]+'synchS_flux.dat')
		synchL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL_flux.dat')
		# freefS = np.loadtxt(galdef_name[n]+'freefreeS_flux.dat')
		freefL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefreeL_flux.dat')
		
	else:
		synch = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux.dat')
		freef = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux.dat')
		# synchS = np.loadtxt(galdef_name[n]+'synchS_flux.dat')
		synchL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux.dat')
		# freefS = np.loadtxt(galdef_name[n]+'freefreeS_flux.dat')
		freefL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux.dat')
		# synch_na = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux_na.dat')
		# freef_na = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux_na.dat')
		
	color = cmap(float(n)/num_graph_galdef_tot)
	# freef[1]*=chi_ff_fac[n]
	# freefL[1]*=chi_ff_fac[n]
	
	LABEL = ''
	for m in range(0,num_params):
		LABEL += VAR_PLOT[m]+str('{:.2e}'.format(float(new_val[m][NN[n][m]])))+' '
	LABEL += '$\chi^2=$'+str('{:.1e}'.format(float(chi[n])))
	
	NFC=1.
	# ax2.loglog(synch_na[0], synch_na[1]*FACTOR[n]*1e23, label=None, c=color, linestyle=':')
	# ax2.loglog(freef_na[0], freef_na[1]*FACTOR[n]*1e23*NFC, label=None, c=color, linestyle=':')
	# ax2.loglog(synchS[0], synchS[1]*FACTOR[n]*1e23, label=None, c=color, linestyle=':')
	ax2.loglog(synchL[0], synchL[1]*FACTOR[n]*1e23, label=None, c=color, linestyle='--')
	# ax2.loglog(freefS[0], freefS[1]*FACTOR[n]*1e23, label=None, c=color, linestyle=':')
	ax2.loglog(freefL[0], freefL[1]*FACTOR[n]*1e23*NFC, label=None, c=color, linestyle='-.')
	# ax2.loglog(synch[0], synch[1]*FACTOR[n]*1e23, label=r''+LABEL, c=color)
	ax2.loglog(synch[0], (synch[1]*FACTOR[n]+freef[1]*FACTOR[n]*NFC)*1e23, label=r''+LABEL, c=color)
	
	
	radio_ex = np.interp(radio_data[0], synch[0], (synch[1]*FACTOR[n]+freef[1]*FACTOR[n]*NFC)*1e23)
	radio_er = (radio_data[1] - radio_ex)/radio_ex  
	radio_error = radio_data[2]/radio_ex
	ax2x.errorbar(radio_data[0], radio_er, yerr=radio_error, fmt='o', markersize=3, c=color)
	
	
	###
	###
	### RADIO FUNCTION OF Z
	###
	###
	###
	
	if SL_SCALE == 1:
		#load data
		zS = np.loadtxt(free_free_dir+'freefreeL.fits_z0.dat')
		nuS = np.loadtxt(free_free_dir+'freefreeL.fits_z0_nu.dat')
		specS = np.loadtxt(free_free_dir+'freefreeL.fits_z0_spec.dat')*chi_ff_fac[n]
		specS*=FACTOR[n]*NFC
		
		zL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL.fits_z0.dat')
		nuL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL.fits_z0_nu.dat')
		specL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL.fits_z0_spec.dat')
		specL*=FACTOR[n]
		
		zS_o = np.loadtxt(free_free_dir+'freefreeL.fits_zoff.dat')
		nuS_o = np.loadtxt(free_free_dir+'freefreeL.fits_zoff_nu.dat')
		specS_o = np.loadtxt(free_free_dir+'freefreeL.fits_zoff_spec.dat')*chi_ff_fac[n]
		specS_o*=FACTOR[n]*NFC
		
		zL_o = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL.fits_zoff.dat')
		nuL_o = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL.fits_zoff_nu.dat')
		specL_o = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL.fits_zoff_spec.dat')
		specL_o*=FACTOR[n]
		
		# zL = np.loadtxt(galdef_name[n]+'synchL_na.fits_z0.dat')
		# nuL = np.loadtxt(galdef_name[n]+'synchL_na.fits_z0_nu.dat')
		# specL = np.loadtxt(galdef_name[n]+'synchL_na.fits_z0_spec.dat')
		# specL*=FACTOR[n]
		
		# zL_o = np.loadtxt(galdef_name[n]+'synchL_na.fits_zoff.dat')
		# nuL_o = np.loadtxt(galdef_name[n]+'synchL_na.fits_zoff_nu.dat')
		# specL_o = np.loadtxt(galdef_name[n]+'synchL_na.fits_zoff_spec.dat')
		# specL_o*=FACTOR[n]
	else:
		#load data
		zS = np.loadtxt(free_free_dir+'freefree.fits_z0.dat')
		nuS = np.loadtxt(free_free_dir+'freefree.fits_z0_nu.dat')
		specS = np.loadtxt(free_free_dir+'freefree.fits_z0_spec.dat')*chi_ff_fac[n]
		specS*=FACTOR[n]*NFC
		
		zL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch.fits_z0.dat')
		nuL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch.fits_z0_nu.dat')
		specL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch.fits_z0_spec.dat')
		specL*=FACTOR[n]
		
		zS_o = np.loadtxt(free_free_dir+'freefree.fits_zoff.dat')
		nuS_o = np.loadtxt(free_free_dir+'freefree.fits_zoff_nu.dat')
		specS_o = np.loadtxt(free_free_dir+'freefree.fits_zoff_spec.dat')*chi_ff_fac[n]
		specS_o*=FACTOR[n]*NFC
		
		zL_o = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch.fits_zoff.dat')
		nuL_o = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch.fits_zoff_nu.dat')
		specL_o = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch.fits_zoff_spec.dat')
		specL_o*=FACTOR[n]
		
	dzS = zS[1]-zS[0]
	dzL = zL[1]-zL[0]
	
	#interpolate data
	data_freq = np.array([c/92, c/22, c/6, c/3])
	y_int = np.array([30, 15, 7, 4.5])
	
	for j in range(0,len(data_freq)):
		NU = data_freq[j] #Hz 
	
		#INTERP GALPROP
		# S
		LOW = specS[ nuS <NU]
		nu_low = nuS[ nuS <NU]
		LOW = LOW[-1]
		nu_low = nu_low[-1]

		HIGH = specS[ nuS >NU]
		nu_high = nuS[ nuS >NU]
		HIGH = HIGH[0]
		nu_high = nu_high[0]

		#specS_temp = LOW + (HIGH-LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low)
		specS_temp = LOW * np.exp(np.log(HIGH/LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low))
		specS_temp*= 1e23
		
		y_integrated = y_int[j] #arcsec
	
		data_z = zS
		data = specS_temp
		YY_z = zS 
		numpoints = data_z.shape[0]
		SPEC = np.zeros([2,numpoints]) 
		for i in range(0,numpoints):
			SPEC[0][i] = YY_z[i]
			if np.any(abs(data_z-YY_z[i])<=y_integrated/2):
				data_temp = data[ abs(data_z-YY_z[i])<=y_integrated/2 ]
				SPEC[1][i] = np.average(data_temp, axis=0)
			else:
				SPEC[1][i] = np.NaN
		SPEC[0]*=DISTANCE/radian2arcsec
		SPEC[1]*= y_integrated/radian2arcsec
		specS_temp = SPEC 
		
		# L
		LOW = specL[ nuL <NU]
		nu_low = nuL[ nuL <NU]
		LOW = LOW[-1]
		nu_low = nu_low[-1]

		HIGH = specL[ nuL >NU]
		nu_high = nuL[ nuL >NU]
		HIGH = HIGH[0]
		nu_high = nu_high[0]

		#specL_temp = LOW + (HIGH-LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low)
		specL_temp = LOW * np.exp(np.log(HIGH/LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low))
		specL_temp*= 1e23
		
		data_z = zL
		data = specL_temp 
		YY_z = zL
		numpoints = data_z.shape[0]
		SPEC = np.zeros([2,numpoints]) 
		for i in range(0,numpoints):
			SPEC[0][i] = YY_z[i]
			if np.any(abs(data_z-YY_z[i])<=y_integrated/2):
				data_temp = data[ abs(data_z-YY_z[i])<=y_integrated/2 ]
				SPEC[1][i] = np.average(data_temp, axis=0)
			else:
				SPEC[1][i] = np.NaN
		SPEC[0]*=DISTANCE/radian2arcsec
		SPEC[1]*= y_integrated/radian2arcsec		
		specL_temp = SPEC 
		
		# S OFFSET
		LOW = specS_o[ nuS <NU]
		nu_low = nuS_o[ nuS <NU]
		LOW = LOW[-1]
		nu_low = nu_low[-1]

		HIGH = specS_o[ nuS >NU]
		nu_high = nuS_o[ nuS >NU]
		HIGH = HIGH[0]
		nu_high = nu_high[0]

		#specS_temp = LOW + (HIGH-LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low)
		specS_o_temp = LOW * np.exp(np.log(HIGH/LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low))
		specS_o_temp*= 1e23
		
		y_integrated = y_int[j] #arcsec
	
		data_z = zS_o
		data = specS_o_temp
		YY_z = zS_o 
		numpoints = data_z.shape[0]
		SPEC = np.zeros([2,numpoints]) 
		for i in range(0,numpoints):
			SPEC[0][i] = YY_z[i]
			if np.any(abs(data_z-YY_z[i])<=y_integrated/2):
				data_temp = data[ abs(data_z-YY_z[i])<=y_integrated/2 ]
				SPEC[1][i] = np.average(data_temp, axis=0)
			else:
				SPEC[1][i] = np.NaN
		SPEC[0]*=DISTANCE/radian2arcsec
		SPEC[1]*= y_integrated/radian2arcsec
		specS_o_temp = SPEC 
		
		# L OFFSET
		LOW = specL_o[ nuL <NU]
		nu_low = nuL_o[ nuL <NU]
		LOW = LOW[-1]
		nu_low = nu_low[-1]

		HIGH = specL_o[ nuL >NU]
		nu_high = nuL_o[ nuL >NU]
		HIGH = HIGH[0]
		nu_high = nu_high[0]

		#specL_o_temp = LOW + (HIGH-LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low)
		specL_o_temp = LOW * np.exp(np.log(HIGH/LOW)*math.log( NU /nu_low)/math.log(nu_high/nu_low))
		specL_o_temp*= 1e23
		
		data_z = zL_o
		data = specL_o_temp 
		YY_z = zL_o
		numpoints = data_z.shape[0]
		SPEC = np.zeros([2,numpoints]) 
		for i in range(0,numpoints):
			SPEC[0][i] = YY_z[i]
			if np.any(abs(data_z-YY_z[i])<=y_integrated/2):
				data_temp = data[ abs(data_z-YY_z[i])<=y_integrated/2 ]
				SPEC[1][i] = np.average(data_temp, axis=0)
			else:
				SPEC[1][i] = np.NaN
		SPEC[0]*=DISTANCE/radian2arcsec
		SPEC[1]*= y_integrated/radian2arcsec		
		specL_o_temp = SPEC 
		
		if j == 0:
			specS_92cm = specS_temp
			specL_92cm = specL_temp
			specS_o_92cm = specS_o_temp
			specL_o_92cm = specL_o_temp
		elif j == 1:
			specS_22cm = specS_temp
			specL_22cm = specL_temp
			specS_o_22cm = specS_o_temp
			specL_o_22cm = specL_o_temp
		elif j == 2:
			specS_6cm = specS_temp
			specL_6cm = specL_temp
		elif j == 3:
			specS_3cm = specS_temp
			specL_3cm = specL_temp
		
		##########
		## Plot
	color = cmap(float(n)/num_graph_galdef_tot)
	
	AX92.plot(specL_92cm[0], specL_92cm[1], c=color, linestyle='--')#, marker="o")
	AX92.plot(specS_92cm[0], specS_92cm[1], c=color, linestyle='-.')#, marker="o")
	AX92.plot(specS_92cm[0], specS_92cm[1]+specL_92cm[1], c=color, linestyle='-')#, marker="o")
	
	AX92_o.plot(specL_o_92cm[0], specL_o_92cm[1], c=color, linestyle='--')#, marker="o")
	AX92_o.plot(specS_o_92cm[0], specS_o_92cm[1], c=color, linestyle='-.')#, marker="o")
	AX92_o.plot(specS_o_92cm[0], specS_o_92cm[1]+specL_o_92cm[1], c=color, linestyle='-')#, marker="o")
	
	AX22.plot(specL_22cm[0], specL_22cm[1], c=color, linestyle='--')#, marker="o")
	AX22.plot(specS_22cm[0], specS_22cm[1], c=color, linestyle='-.')#, marker="o")
	AX22.plot(specS_22cm[0], specS_22cm[1]+specL_22cm[1], c=color, linestyle='-')#, marker="o")
	
	AX22_o.plot(specL_o_22cm[0], specL_o_22cm[1], c=color, linestyle='--')#, marker="o")
	AX22_o.plot(specS_o_22cm[0], specS_o_22cm[1], c=color, linestyle='-.')#, marker="o")
	AX22_o.plot(specS_o_22cm[0], specS_o_22cm[1]+specL_o_22cm[1], c=color, linestyle='-')#, marker="o")
	
	AX6.plot(specL_6cm[0], specL_6cm[1], c=color, linestyle='--')#, marker="o")
	AX6.plot(specS_6cm[0], specS_6cm[1], c=color, linestyle='-.')#, marker="o")
	AX6.plot(specS_6cm[0], specS_6cm[1]+specL_6cm[1], c=color, linestyle='-')#, marker="o")
	
	AX3.plot(specL_3cm[0], specL_3cm[1], c=color, linestyle='--')#, marker="o")
	AX3.plot(specS_3cm[0], specS_3cm[1], c=color, linestyle='-.')#, marker="o")
	AX3.plot(specS_3cm[0], specS_3cm[1]+specL_3cm[1], c=color, linestyle='-')#, marker="o")


# plt.suptitle(r'M82 Synchrotron $F(z)$ (Adebahr et al. 2013)')
# ax01.set_ylabel(r'$F(z)$ (Jy)')
# ax11.set_ylabel(r'$F(z)$ (Jy)')
# ax21.set_ylabel(r'$F(z)$ (Jy)')
# ax02.set_ylabel(r'$F(z)$ (Jy)')
# ax12.set_ylabel(r'$F(z)$ (Jy)')
# ax22.set_ylabel(r'$F(z)$ (Jy)')

# ax01.set_xlabel(r'$z$ (kpc)')
# ax02.set_xlabel(r'$z$ (kpc)')
# ax11.set_xlabel(r'$z$ (kpc)')
# ax12.set_xlabel(r'$z$ (kpc)')
ax21.set_xlabel(r'$z$ (kpc)')
ax22.set_xlabel(r'$z$ (kpc)')

# Turn off tick labels
ax01.set_yticklabels([])
ax11.set_yticklabels([])
ax21.set_yticklabels([])

ax1.text(0.75, 0.8, r'Integrated $\gamma$-ray', ha='center', va='center', transform=ax1.transAxes)
ax2.text(0.75, 0.8, r'Integrated Radio', ha='center', va='center', transform=ax2.transAxes)
ax2.legend(loc='upper left', bbox_to_anchor=(0., -2.05), prop={'size':9})


plot_text = '''BLAH BLAH BLAH '''

fig.text(0.725,0,plot_text,va='top')

plt.tight_layout()
plot_name = 'plot_'+CHANGED+'_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2), bbox_inches='tight', dpi=400)
		break
#plt.savefig('plot_gamma_'+CHANGED, dpi=400)
plt.close()

###
###
### Plot radio and gamma together
###
###

plot_gamma_radio = 1

if plot_gamma_radio == 1:
	h = 6.6260755e-27 #erg s
	ev2ergs = 1.60218e-12
	mev2ev = 1.e6

	fig = plt.figure()
	ax = plt.subplot(111)
	#f, (ax2, ax1) = plt.subplots(1, 2, figsize=(16,9))
	
	###
	###
	### GAMMA-RAY
	###
	###
	###

	# fermi_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_fermi.dat').T
	verit_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_veritas.dat').T
	verit_upper = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_gamma/m82_veritas_upper.dat')
	#verit_data[0] *= 1000
	verit_data[1:] *= verit_data[0]
	verit_upper[1] *= verit_upper[0]

	# alpha = -2.21
	# emin = fermi_data[1]
	# emax = fermi_data[2]
	# ferr_d = fermi_data[3]*(alpha+1)/(emax**(alpha+1)-emin**(alpha+1))*(emin*emax)**((alpha+2)/2)
	# ferr_u = fermi_data[4]*(alpha+1)/(emax**(alpha+1)-emin**(alpha+1))*(emin*emax)**((alpha+2)/2)

	# fermi_energy = np.zeros(fermi_data.shape[1])
	# fermi_energy = np.sqrt(emin*emax)

	# e2dnde = fermi_data[0]*(alpha+1)/(emax**(alpha+1)-emin**(alpha+1))*(emin*emax)**((alpha+2)/2)
	#print fermi_energy
	#print e2dnde

	#PLOTTING DATA
	#ax1.errorbar(fermi_energy, e2dnde, xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label=r'FERMI $\alpha=-2.2$', c='black')
	#plt.errorbar(fermi_energy, fermi_data[0]*fermi_energy, xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(fermi_data[3])*fermi_energy, fermi_data[4]*fermi_energy],fmt='x', c='black')
	ax.errorbar(verit_data[0]*mev2ev, verit_data[1]*mev2ev, yerr=[abs(verit_data[2]-verit_data[1])*mev2ev, abs(verit_data[3]-verit_data[1])*mev2ev], fmt='s', markersize=3, label='VERITAS', c='black')
	ax.errorbar(verit_upper[0]*mev2ev, verit_upper[1]*mev2ev, yerr=verit_upper[1]/1.5*mev2ev, uplims=True, c='black', markersize=3)
	#ax1.set_title(r'$\gamma$-ray Flux for '+str(CHANGED))
	# ax.set_xlabel(r'E (MeV)')
	# ax.set_ylabel(r'$E \, F_{\gamma}(E)$  (MeV$^2$ s$^{-1}$ cm$^{-2}$ MeV$^{-1}$)')
	# ax.set_xlim([1.e+0,1.e7])
	# ax.set_ylim([1.e-10, 1e-4])

	# ax1x.plot((1e0,1e15),(0, 0), 'k', linestyle=':')
	# ax1x.set_ylim([-1., 1])

	###
	###
	### RADIO
	###
	###
	###
	
	radio_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_radio/radio_m82.dat')
	radio_data = radio_data.T
	
	radio_data[1]*= radio_data[0]*1.e-23/ev2ergs
	radio_data[2]*= radio_data[0]*1.e-23/ev2ergs
	radio_data[0]*= h/ev2ergs
	
	
	#PLOTTING DATA
	ax.errorbar(radio_data[0], radio_data[1], yerr=radio_data[2], fmt='o', markersize=3, c='black')
	# ax2.set_title('Radio Flux for '+str(CHANGED))
	# ax2.set_xlabel(r'$\nu$ (Hz)')
	# ax2.set_ylabel(r'$F_{radio}(\nu)$ (Jy)')
	# ax2.set_xlim([1.e+6,1.e+15])
	# ax2.set_ylim([1.e-3, 1.e2])
	
	
	for n in range(0,num_graph_galdef_tot):
	
		LABEL = ''
		for m in range(0,num_params):
			LABEL += VAR_PLOT[m]+str('{:.2e}'.format(float(new_val[m][NN[n][m]])))+' '
		LABEL += '$\chi^2=$'+str('{:.1e}'.format(float(chi[n])))
		
		###
		###
		### GAMMA-RAY
		###
		###
		###
		
		gamma_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'gamma_flux.dat')*mev2ev
		gamma_tot_mask = np.isfinite(gamma_tot[1])
		
		if E2FLAG == 0:
			gamma_tot[1] /= gamma_tot[0]
		
		bremss_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'bremss_flux.dat')*mev2ev
		pion_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'pion_flux.dat')*mev2ev
		ic_tot = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'IC_flux.dat')*mev2ev
		color = cmap(float(n)/num_graph_galdef_tot)
		
		fermi_data = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'fermi.dat').T
		emin = fermi_data[0]*mev2ev
		emax = fermi_data[1]*mev2ev
		fermi_energy = np.sqrt(emin*emax)
		ferr_d = fermi_data[3]*mev2ev
		ferr_u = fermi_data[4]*mev2ev
		
		ax.loglog(gamma_tot[0][gamma_tot_mask], gamma_tot[1][gamma_tot_mask]*FACTOR[n], label=None, c=color, linestyle='-')
		ax.loglog(bremss_tot[0], bremss_tot[1]*FACTOR[n], label=None, c=color, linestyle=':')
		ax.loglog(pion_tot[0], pion_tot[1]*FACTOR[n], label=None, c=color, linestyle='--')
		ax.loglog(ic_tot[0], ic_tot[1]*FACTOR[n], label=None, c=color, linestyle='-.')
		
		ax.errorbar(fermi_energy, fermi_data[2]*mev2ev, xerr=[abs(emin-fermi_energy),emax-fermi_energy], yerr=[abs(ferr_d), ferr_u],fmt='^', label='FERMI', c='black', markersize=3)
		
		# gamma_ex = np.interp(fermi_energy, gamma_tot[0][gamma_tot_mask], gamma_tot[1][gamma_tot_mask]*FACTOR[n])
		# gamma_er = (fermi_data[2] - gamma_ex)/gamma_ex  
		# gamma_error_dd = ferr_d/gamma_ex
		# gamma_error_uu = ferr_u/gamma_ex
		# ax1x.errorbar(fermi_energy, gamma_er, xerr=[abs(emin-fermi_energy),emax-fermi_energy],yerr=[abs(gamma_error_dd), gamma_error_uu], fmt='o', markersize=3, c=color)
		
		# gamma_ex = np.interp(verit_data[0], gamma_tot[0][gamma_tot_mask], gamma_tot[1][gamma_tot_mask]*FACTOR[n])
		# gamma_er = (verit_data[1] - gamma_ex)/gamma_ex  
		# gamma_error_dd = (verit_data[2]-verit_data[1])/gamma_ex
		# gamma_error_uu = (verit_data[3]-verit_data[1])/gamma_ex
		# axx.errorbar(verit_data[0], gamma_er, yerr=[abs(gamma_error_dd), abs(gamma_error_uu)], fmt='o', markersize=3, c=color)
	
		###
		###
		### RADIO
		###
		###
		###
		
		if SL_SCALE == 1:
			synch = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux.dat')
			freef = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux.dat')
			# synchS = np.loadtxt(galdef_name[n]+'synchS_flux.dat')
			synchL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synchL_flux.dat')
			# freefS = np.loadtxt(galdef_name[n]+'freefreeS_flux.dat')
			freefL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefreeL_flux.dat')
			
		else:
			synch = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux.dat')
			freef = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux.dat')
			# synchS = np.loadtxt(galdef_name[n]+'synchS_flux.dat')
			synchL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux.dat')
			# freefS = np.loadtxt(galdef_name[n]+'freefreeS_flux.dat')
			freefL = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux.dat')
			# synch_na = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'synch_flux_na.dat')
			# freef_na = np.loadtxt(dir_in+galdef_name[n]+'/'+galdef_name[n]+'freefree_flux_na.dat')
			
		color = cmap(float(n)/num_graph_galdef_tot)
		
		NFC=1.
		
		synch[1]*= synch[0]/ev2ergs
		synch[0]*= h/ev2ergs
		synch[1]*= FACTOR[n]
		
		freef[1]*= freef[0]/ev2ergs
		freef[0]*= h/ev2ergs
		freef[1]*= FACTOR[n]*chi_ff_fac[n]
		
		synchL[1]*= synchL[0]/ev2ergs
		synchL[0]*= h/ev2ergs
		synchL[1]*= FACTOR[n]
		
		freefL[1]*= freefL[0]/ev2ergs
		freefL[0]*= h/ev2ergs
		freefL[1]*= FACTOR[n]*chi_ff_fac[n]
		
		# ax2.loglog(synch_na[0], synch_na[1]*FACTOR[n]*1e23, label=None, c=color, linestyle=':')
		# ax2.loglog(freef_na[0], freef_na[1]*FACTOR[n]*1e23*NFC, label=None, c=color, linestyle=':')
		# ax2.loglog(synchS[0], synchS[1]*FACTOR[n]*1e23, label=None, c=color, linestyle=':')
		ax.loglog(synchL[0], synchL[1], label=None, c=color, linestyle='--')
		# ax2.loglog(freefS[0], freefS[1]*FACTOR[n]*1e23, label=None, c=color, linestyle=':')
		ax.loglog(freefL[0], freefL[1], label=None, c=color, linestyle='-.')
		# ax2.loglog(synch[0], synch[1]*FACTOR[n]*1e23, label=r''+LABEL, c=color)
		ax.loglog(synch[0], (synch[1]+freef[1]), label=r''+LABEL, c=color)
		
		
		# radio_ex = np.interp(radio_data[0], synch[0], (synch[1]*FACTOR[n]+freef[1]*FACTOR[n]*NFC)*1e23)
		# radio_er = (radio_data[1] - radio_ex)/radio_ex  
		# radio_error = radio_data[2]/radio_ex
		# ax2x.errorbar(radio_data[0], radio_er, yerr=radio_error, fmt='o', markersize=3, c=color)
	
	ax.set_title('GP CR Emission for '+str(CHANGED))
	ax.set_xlabel(r'$E$ (eV)')
	ax.set_ylabel(r'$E^2 \frac{d F(E)}{dE}$ (eV s$^{-1}$ cm$^{-2}$)')
	ax.grid()
	
	#ax.set_xlim([1.e+7,1.e+11])
	ax.set_ylim([1.e-5, 1.e1])
	
#plt.legend(loc=4)
plot_name = 'plot_total_'+CHANGED
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2), dpi=400)
		break
#plt.savefig('plot_gamma_'+CHANGED, dpi=400)
plt.close()


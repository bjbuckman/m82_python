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

GG =  'm82_06bff1020'
galdef_name.append(GG)

galdef_name_1 = []

GG_1 = 'm82_05aff1019'
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

MSIZE = 5
LWIDTH = 3
ALPHA_DAT = 0.5

fig = plt.figure()
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1]) 
gs.update(wspace=0.0, hspace=0.0)

cmap = plt.get_cmap('cubehelix')

plt.subplots_adjust(hspace=0.001)

ax00 = plt.subplot(gs[0,0])
ax01 = plt.subplot(gs[0,1])
ax10 = plt.subplot(gs[1,0])
ax11 = plt.subplot(gs[1,1])

for axis in ['top','bottom','left','right']:
	ax00.spines[axis].set_linewidth(2)
	ax01.spines[axis].set_linewidth(2)
	ax10.spines[axis].set_linewidth(2)
	ax11.spines[axis].set_linewidth(2)

ax00.tick_params(which='major',width=1, length=8, labelsize='small')
ax00.tick_params(which='minor',width=0.5, length=5)
ax01.tick_params(which='major',width=1, length=8, labelsize='small')
ax01.tick_params(which='minor',width=0.5, length=5)
ax10.tick_params(which='major',width=1, length=8, labelsize='small')
ax10.tick_params(which='minor',width=0.5, length=5)
ax11.tick_params(which='major',width=1, length=8, labelsize='small')
ax11.tick_params(which='minor',width=0.5, length=5)

AX92 = ax11
AX22 = ax10
AX6 = ax01
AX3 = ax00



###
###
### RADIO FUNCTION OF Z
###
###
###

XLIM = [-3.5,3.5]

## 92cm
L92cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_N.dat').T
L92cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_S.dat').T
L92_cm_beam = 43.1/radian2arcsec*DISTANCE
AX92.errorbar([2.33], [6.e-1],xerr=L92_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)

AX92.errorbar(L92cm_N[0],L92cm_N[1],yerr=[abs(L92cm_N[2]-L92cm_N[1]),abs(L92cm_N[3]-L92cm_N[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX92.errorbar(-L92cm_S[0],L92cm_S[1],yerr=[abs(L92cm_S[2]-L92cm_S[1]),abs(L92cm_S[3]-L92cm_S[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX92.set_yscale('log')
AX92.set_xlim(XLIM)
AX92.set_ylim([1e-4,2e+1])
AX92.text(0.8, 0.8, r'$\lambda = 92\,\mathrm{cm}$', ha='center', va='center', transform=AX92.transAxes)



# ## 92cm OFFSET
# L92cm_NW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_NW.dat').T
# L92cm_SW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_SW.dat').T
# L92cm_NE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_NE.dat').T
# L92cm_SE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_SE.dat').T


# AX92_o.errorbar(L92cm_NW[0],L92cm_NW[1],yerr=[abs(L92cm_NW[2]-L92cm_NW[1]),abs(L92cm_NW[3]-L92cm_NW[1])],fmt='o',color='black')
# AX92_o.errorbar(L92cm_SW[0],L92cm_SW[1],yerr=[abs(L92cm_SW[2]-L92cm_SW[1]),abs(L92cm_SW[3]-L92cm_SW[1])],fmt='o',color='grey')
# AX92_o.errorbar(L92cm_NE[0],L92cm_NE[1],yerr=[abs(L92cm_NE[2]-L92cm_NE[1]),abs(L92cm_NE[3]-L92cm_NE[1])],fmt='o',color='black')
# AX92_o.errorbar(L92cm_SE[0],L92cm_SE[1],yerr=[abs(L92cm_SE[2]-L92cm_SE[1]),abs(L92cm_SE[3]-L92cm_SE[1])],fmt='o',color='grey')
# AX92_o.set_yscale('log')
# AX92_o.set_xlim([0,3.5])
# AX92_o.set_ylim([1e-4,1e+2])
# # AX22.set_xlabel('kpc')
# # AX22.set_ylabel('Jy')
# AX92_o.text(0.75, 0.8, r'$\lambda = 92$ cm''\n''$60$" offset', ha='center', va='center', transform=AX92_o.transAxes)

## 22cm
L22cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_N.dat').T
L22cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_S.dat').T
L22_cm_beam = 12.7/radian2arcsec*DISTANCE
AX22.errorbar([2.33], [6.e-1],xerr=L22_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)

AX22.errorbar(L22cm_N[0],L22cm_N[1],yerr=[abs(L22cm_N[2]-L22cm_N[1]),abs(L22cm_N[3]-L22cm_N[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX22.errorbar(-L22cm_S[0],L22cm_S[1],yerr=[abs(L22cm_S[2]-L22cm_S[1]),abs(L22cm_S[3]-L22cm_S[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX22.set_yscale('log')
AX22.set_xlim(XLIM)
AX22.set_ylim([1e-4,2e+1])
# AX22.set_xlabel('kpc')
# AX22.set_ylabel('Jy')
AX22.text(0.8, 0.8, r'$\lambda = 22\,\mathrm{cm}$', ha='center', va='center', transform=AX22.transAxes)

# ## 22cm OFFSET
# L22cm_NW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_NW.dat').T
# L22cm_SW = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_SW.dat').T
# L22cm_NE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_NE.dat').T
# L22cm_SE = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_SE.dat').T


# AX22_o.errorbar(L22cm_NW[0],L22cm_NW[1],yerr=[abs(L22cm_NW[2]-L22cm_NW[1]),abs(L22cm_NW[3]-L22cm_NW[1])],fmt='o',color='black')
# AX22_o.errorbar(L22cm_SW[0],L22cm_SW[1],yerr=[abs(L22cm_SW[2]-L22cm_SW[1]),abs(L22cm_SW[3]-L22cm_SW[1])],fmt='o',color='grey')
# AX22_o.errorbar(L22cm_NE[0],L22cm_NE[1],yerr=[abs(L22cm_NE[2]-L22cm_NE[1]),abs(L22cm_NE[3]-L22cm_NE[1])],fmt='o',color='black')
# AX22_o.errorbar(L22cm_SE[0],L22cm_SE[1],yerr=[abs(L22cm_SE[2]-L22cm_SE[1]),abs(L22cm_SE[3]-L22cm_SE[1])],fmt='o',color='grey')
# AX22_o.set_yscale('log')
# AX22_o.set_xlim([0,3.5])
# AX22_o.set_ylim([1e-4,1e+2])
# # AX22.set_xlabel('kpc')
# # AX22.set_ylabel('Jy')
# AX22_o.text(0.75, 0.8, r'$\lambda = 22$ cm''\n''$60$" offset', ha='center', va='center', transform=AX22_o.transAxes)


## 6cm
L6cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_N.dat').T
L6cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_S.dat').T
L6_cm_beam = 12.5/radian2arcsec*DISTANCE
AX6.errorbar([2.33], [6.e-1],xerr=L6_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)

AX6.errorbar(L6cm_N[0],L6cm_N[1],yerr=[abs(L6cm_N[2]-L6cm_N[1]),abs(L6cm_N[3]-L6cm_N[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX6.errorbar(-L6cm_S[0],L6cm_S[1],yerr=[abs(L6cm_S[2]-L6cm_S[1]),abs(L6cm_S[3]-L6cm_S[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX6.set_yscale('log')
AX6.set_xlim(XLIM)
AX6.set_ylim([1e-4,2e+1])
# AX6.set_xlabel('kpc')
# AX6.set_ylabel('Jy')
AX6.text(0.8, 0.8, r'$\lambda = 6\,\mathrm{cm}$', ha='center', va='center', transform=AX6.transAxes)

## 3cm
L3cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_N.dat').T
L3cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_S.dat').T
L3_cm_beam = 7.6/radian2arcsec*DISTANCE
AX3.errorbar([2.33], [6.e-1],xerr=L3_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)

AX3.errorbar(L3cm_N[0],L3cm_N[1],yerr=[abs(L3cm_N[2]-L3cm_N[1]),abs(L3cm_N[3]-L3cm_N[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX3.errorbar(-L3cm_S[0],L3cm_S[1],yerr=[abs(L3cm_S[2]-L3cm_S[1]),abs(L3cm_S[3]-L3cm_S[1])],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX3.set_yscale('log')
AX3.set_xlim(XLIM)
AX3.set_ylim([1e-4,2e+1])
# AX3.set_xlabel('kpc')
# aAX3.set_ylabel('Jy')
AX3.text(0.8, 0.8, r'$\lambda = 3\,\mathrm{cm}$', ha='center', va='center', transform=AX3.transAxes)

def calc_bins(j,y_integrated, zS, specS): 
	
	#INTERP	GALPROP
	# S
	specS_temp = specS
	
	data_z = zS
	data = specS_temp 
	numpoints = data_z.shape[0]
	SPEC = np.zeros([2,numpoints])
	SPEC[0] = data_z
	SPEC[0]*=DISTANCE/radian2arcsec
	SPEC[1] = data[j]
	SPEC[1]*= y_integrated
	specS_temp = SPEC 
	return SPEC


for n in range(0,num_graph):
	LABEL = ['A','B']
	
	gp_file = dir_comp+galdef_name[n]+'/'+galdef_name[n]
	gp_file_1 = dir_comp+galdef_name_1[n]+'/'+galdef_name_1[n]
	
	#load data
	zF = np.loadtxt(gp_file+'_z_freef.dat')
	specF = np.loadtxt(gp_file+'_z0_nu_freef.dat')
	specF*=FACTOR[n]*1.e23
	
	zS = np.loadtxt(gp_file+'_z_synch.dat')
	specS = np.loadtxt(gp_file+'_z0_nu_synch.dat')
	specS*=FACTOR[n]*1.e23
	
	zF_1 = np.loadtxt(gp_file_1+'_z_freef.dat')
	specF_1 = np.loadtxt(gp_file_1+'_z0_nu_freef.dat')
	specF_1*=FACTOR[n]*1.e23
	
	zS_1 = np.loadtxt(gp_file_1+'_z_synch.dat')
	specS_1 = np.loadtxt(gp_file_1+'_z0_nu_synch.dat')
	specS_1*=FACTOR_1[n]*1.e23
	
	#interpolate data
	NU = np.array([c/92, c/22, c/6, c/3])
	y_int = np.array([30, 15, 7, 4.5])
	
	for j in range(0,len(NU)):
		if j == 0:
			specS_92cm = calc_bins(j, y_int[j], zS, specS)
			specF_92cm = calc_bins(j, y_int[j], zF, specF)
			specS_92cm_1 = calc_bins(j, y_int[j], zS_1, specS_1)
			specF_92cm_1 = calc_bins(j, y_int[j], zF_1, specF_1)
		elif j == 1:
			specS_22cm = calc_bins(j, y_int[j], zS, specS)
			specF_22cm = calc_bins(j, y_int[j], zF, specF)
			specS_22cm_1 = calc_bins(j, y_int[j], zS_1, specS_1)
			specF_22cm_1 = calc_bins(j, y_int[j], zF_1, specF_1)
		elif j == 2:
			specS_6cm = calc_bins(j, y_int[j], zS, specS)
			specF_6cm = calc_bins(j, y_int[j], zF, specF)
			specS_6cm_1 = calc_bins(j, y_int[j], zS_1, specS_1)
			specF_6cm_1 = calc_bins(j, y_int[j], zF_1, specF_1)
		elif j == 3:
			specS_3cm = calc_bins(j, y_int[j], zS, specS)
			specF_3cm = calc_bins(j, y_int[j], zF, specF)
			specS_3cm_1 = calc_bins(j, y_int[j], zS_1, specS_1)
			specF_3cm_1 = calc_bins(j, y_int[j], zF_1, specF_1)
		
	##########
	## Plot
	color = cmap(1-float(n+2)/(num_graph+2))
	
	ALPHA = 0.4
	
	AX92.plot(specF_92cm[0], specF_92cm[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX92.plot(specS_92cm[0], specS_92cm[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX92.plot(specF_92cm_1[0], specF_92cm_1[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX92.plot(specS_92cm_1[0], specS_92cm_1[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX92.fill_between(specS_92cm[0], specS_92cm[1]+specF_92cm[1], specS_92cm_1[1]+specF_92cm_1[1], facecolor=color, edgecolor=color, linestyle='-', linewidth=LWIDTH, alpha=0.75, label=LABEL[n])
	
	AX22.plot(specF_22cm[0], specF_22cm[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX22.plot(specS_22cm[0], specS_22cm[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX22.plot(specF_22cm_1[0], specF_22cm_1[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX22.plot(specS_22cm_1[0], specS_22cm_1[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX22.fill_between(specS_22cm[0], specS_22cm[1]+specF_22cm[1], specS_22cm_1[1]+specF_22cm_1[1], facecolor=color, edgecolor=color, linestyle='-', linewidth=LWIDTH, alpha=0.75, label=LABEL[n])
	
	AX6.plot(specF_6cm[0], specF_6cm[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX6.plot(specS_6cm[0], specS_6cm[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX6.plot(specF_6cm_1[0], specF_6cm_1[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX6.plot(specS_6cm_1[0], specS_6cm_1[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX6.fill_between(specS_6cm[0], specS_6cm[1]+specF_6cm[1], specS_6cm_1[1]+specF_6cm_1[1], facecolor=color, edgecolor=color, linestyle='-', linewidth=LWIDTH, alpha=0.75, label=LABEL[n])
	
	AX3.plot(specF_3cm[0], specF_3cm[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX3.plot(specS_3cm[0], specS_3cm[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX3.plot(specF_3cm_1[0], specF_3cm_1[1], c=color, linestyle=':', linewidth=LWIDTH, alpha=ALPHA)
	AX3.plot(specS_3cm_1[0], specS_3cm_1[1], c=color, linestyle='--', linewidth=LWIDTH, alpha=ALPHA)
	AX3.fill_between(specS_3cm[0], specS_3cm[1]+specF_3cm[1], specS_3cm_1[1]+specF_3cm_1[1], facecolor=color, edgecolor=color, linestyle='-', linewidth=LWIDTH, alpha=0.75, label=LABEL[n])


# plt.suptitle(r'M82 Synchrotron $F(z)$ (Adebahr et al. 2013)')
ax00.set_ylabel(r'$F_\lambda (z)$ [Jy]')
ax10.set_ylabel(r'$F_\lambda (z)$ [Jy]')
# ax21.set_ylabel(r'$F(z)$ (Jy)')
# ax02.set_ylabel(r'$F(z)$ (Jy)')
# ax12.set_ylabel(r'$F(z)$ (Jy)')
# ax22.set_ylabel(r'$F(z)$ (Jy)')

# ax01.set_xlabel(r'$z$ (kpc)')
# ax02.set_xlabel(r'$z$ (kpc)')
# ax11.set_xlabel(r'$z$ (kpc)')
# ax12.set_xlabel(r'$z$ (kpc)')
ax10.set_xlabel(r'$z$ [kpc]')
ax11.set_xlabel(r'$z$ [kpc]')

# Turn off tick labels
ax01.set_yticklabels([])
ax11.set_yticklabels([])

ax00.set_yticklabels(['','',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',r'$10^{0}$',r'$10^{1}$'])
ax10.set_yticklabels(['','',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',r'$10^{0}$',r'$10^{1}$'])

ax00.set_xticklabels([])
ax01.set_xticklabels([])

# ax10.set_xticklabels([0,'',1,'',2,'',3])
# ax11.set_xticklabels([0,'',1,'',2,'',3])

# ax1.text(0.75, 0.8, r'Integrated $\gamma$-ray', ha='center', va='center', transform=ax1.transAxes)
# ax2.text(0.75, 0.8, r'Integrated Radio', ha='center', va='center', transform=ax2.transAxes)
AX3.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), prop={'size':15})


# plot_text = '''BLAH BLAH BLAH '''

# fig.text(0.725,0,plot_text,va='top')

plt.tight_layout()
plot_name = 'plot_ext_radio_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break
#plt.savefig('plot_gamma_'+CHANGED, dpi=400)
plt.close()


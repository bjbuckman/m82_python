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

fig = plt.figure(figsize=(12,9))
gs = gridspec.GridSpec(6, 2, width_ratios=[1, 1]) 
gs.update(wspace=0.0, hspace=0.0)
plt.subplots_adjust(hspace=0.001)

cmap = plt.get_cmap('cubehelix')

ax00 = plt.subplot(gs[0:2,0])
ax01 = plt.subplot(gs[0:2,1])
ax10 = plt.subplot(gs[3:5,0])
ax11 = plt.subplot(gs[3:5,1])

ax00e = plt.subplot(gs[2,0])
ax01e = plt.subplot(gs[2,1])
ax10e = plt.subplot(gs[5,0])
ax11e = plt.subplot(gs[5,1])

for axis in ['top','bottom','left','right']:
	ax00.spines[axis].set_linewidth(2)
	ax01.spines[axis].set_linewidth(2)
	ax10.spines[axis].set_linewidth(2)
	ax11.spines[axis].set_linewidth(2)
	
	ax00e.spines[axis].set_linewidth(2)
	ax01e.spines[axis].set_linewidth(2)
	ax10e.spines[axis].set_linewidth(2)
	ax11e.spines[axis].set_linewidth(2)

ax00.tick_params(which='major',width=1, length=8, labelsize='small')
ax00.tick_params(which='minor',width=0.5, length=5)
ax01.tick_params(which='major',width=1, length=8, labelsize='small')
ax01.tick_params(which='minor',width=0.5, length=5)
ax10.tick_params(which='major',width=1, length=8, labelsize='small')
ax10.tick_params(which='minor',width=0.5, length=5)
ax11.tick_params(which='major',width=1, length=8, labelsize='small')
ax11.tick_params(which='minor',width=0.5, length=5)

ax00e.tick_params(which='major',width=1, length=8, labelsize='small')
ax00e.tick_params(which='minor',width=0.5, length=5)
ax01e.tick_params(which='major',width=1, length=8, labelsize='small')
ax01e.tick_params(which='minor',width=0.5, length=5)
ax10e.tick_params(which='major',width=1, length=8, labelsize='small')
ax10e.tick_params(which='minor',width=0.5, length=5)
ax11e.tick_params(which='major',width=1, length=8, labelsize='small')
ax11e.tick_params(which='minor',width=0.5, length=5)

AX92 = ax11
AX22 = ax10
AX6 = ax01
AX3 = ax00

AX92e = ax11e
AX22e = ax10e
AX6e = ax01e
AX3e = ax00e

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

L92cm_S[0]*=-1
L92cm = np.append(L92cm_N[0:2], L92cm_S[0:2], axis=-1)
L92N_error = (abs(L92cm_N[2]-L92cm_N[1]) + abs(L92cm_N[3]-L92cm_N[1]))/2.
L92S_error = (abs(L92cm_S[2]-L92cm_S[1]) + abs(L92cm_S[3]-L92cm_S[1]))/2.
L92_error = np.append(L92N_error, L92S_error)
L92cm = np.append(L92cm, [L92_error], axis=0)
L92cm = L92cm[:,L92cm[0].argsort()] #Single ordered array

L92_cm_beam = 43.1/radian2arcsec*DISTANCE
AX92.errorbar([2.33], [6.e-1],xerr=L92_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1) #plot beam width

AX92.errorbar(L92cm[0],L92cm[1],yerr=L92cm[2],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT) #plot data points
AX92.set_yscale('log')
AX92.set_xlim(XLIM)
AX92.set_ylim([1e-4,2e+1])
AX92.text(0.8, 0.8, r'$\lambda = 92\,\mathrm{cm}$', ha='center', va='center', transform=AX92.transAxes)
AX92e.set_xlim(XLIM)
AX92e.set_ylim([-2,1])
AX92e.plot((-4,4),(0, 0), 'k', linestyle=':')

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

L22cm_S[0]*=-1
L22cm = np.append(L22cm_N[0:2], L22cm_S[0:2], axis=-1)
L22N_error = (abs(L22cm_N[2]-L22cm_N[1]) + abs(L22cm_N[3]-L22cm_N[1]))/2.
L22S_error = (abs(L22cm_S[2]-L22cm_S[1]) + abs(L22cm_S[3]-L22cm_S[1]))/2.
L22_error = np.append(L22N_error, L22S_error)
L22cm = np.append(L22cm, [L22_error], axis=0)
L22cm = L22cm[:,L22cm[0].argsort()] #Single ordered array

L22_cm_beam = 12.7/radian2arcsec*DISTANCE
AX22.errorbar([2.33], [6.e-1],xerr=L22_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)

AX22.errorbar(L22cm[0],L22cm[1],yerr=L22cm[2],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX22.set_yscale('log')
AX22.set_xlim(XLIM)
AX22.set_ylim([1e-4,2e+1])
# AX22.set_xlabel('kpc')
# AX22.set_ylabel('Jy')
AX22.text(0.8, 0.8, r'$\lambda = 22\,\mathrm{cm}$', ha='center', va='center', transform=AX22.transAxes)
AX22e.set_xlim(XLIM)
AX22e.set_ylim([-2,1])
AX22e.plot((-4,4),(0, 0), 'k', linestyle=':')

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

L6cm_S[0]*=-1
L6cm = np.append(L6cm_N[0:2], L6cm_S[0:2], axis=-1)
L6N_error = (abs(L6cm_N[2]-L6cm_N[1]) + abs(L6cm_N[3]-L6cm_N[1]))/2.
L6S_error = (abs(L6cm_S[2]-L6cm_S[1]) + abs(L6cm_S[3]-L6cm_S[1]))/2.
L6_error = np.append(L6N_error, L6S_error)
L6cm = np.append(L6cm, [L6_error], axis=0)
L6cm = L6cm[:,L6cm[0].argsort()] #Single ordered array

L6_cm_beam = 12.5/radian2arcsec*DISTANCE
AX6.errorbar([2.33], [6.e-2],xerr=L6_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)

AX6.errorbar(L6cm[0],L6cm[1],yerr=L6cm[2],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX6.set_yscale('log')
AX6.set_xlim(XLIM)
AX6.set_ylim([1e-5,2e+0])
# AX6.set_xlabel('kpc')
# AX6.set_ylabel('Jy')
AX6.text(0.8, 0.8, r'$\lambda = 6\,\mathrm{cm}$', ha='center', va='center', transform=AX6.transAxes)
AX6e.set_xlim(XLIM)
AX6e.set_ylim([-2,1])
AX6e.plot((-4,4),(0, 0), 'k', linestyle=':')

## 3cm
L3cm_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_N.dat').T
L3cm_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_S.dat').T

L3cm_S[0]*=-1
L3cm = np.append(L3cm_N[0:2], L3cm_S[0:2], axis=-1)
L3N_error = (abs(L3cm_N[2]-L3cm_N[1]) + abs(L3cm_N[3]-L3cm_N[1]))/2.
L3S_error = (abs(L3cm_S[2]-L3cm_S[1]) + abs(L3cm_S[3]-L3cm_S[1]))/2.
L3_error = np.append(L3N_error, L3S_error)
L3cm = np.append(L3cm, [L3_error], axis=0)
L3cm = L3cm[:,L3cm[0].argsort()] #Single ordered array

L3_cm_beam = 7.6/radian2arcsec*DISTANCE
AX3.errorbar([2.33], [6.e-2],xerr=L3_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)

AX3.errorbar(L3cm[0],L3cm[1],yerr=L3cm[2],fmt='o',color='black', markersize=MSIZE, alpha=ALPHA_DAT)
AX3.set_yscale('log')
AX3.set_xlim(XLIM)
AX3.set_ylim([1e-5,2e+0])
# AX3.set_xlabel('kpc')
# aAX3.set_ylabel('Jy')
AX3.text(0.8, 0.8, r'$\lambda = 3\,\mathrm{cm}$', ha='center', va='center', transform=AX3.transAxes)
AX3e.set_xlim(XLIM)
AX3e.set_ylim([-2,1])
AX3e.plot((-4,4),(0, 0), 'k', linestyle=':')

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

HATCH = ['-', '|','/']
for n in range(0,num_graph):
	LABEL = ['A','B',"B'"]
	
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
	color = cmap(1-float(n+1)/(num_graph+1.))
	
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
	
	##########
	## Error calculation & plot
	
	#92cm
	L92_ex1 = np.interp(L92cm[0], specS_92cm[0], specS_92cm[1]+specF_92cm[1])
	L92_er1 = -(L92cm[1]-L92_ex1)/L92_ex1
	L92_error1 = L92cm[2]/L92_ex1
	
	L92_ex2 = np.interp(L92cm[0], specS_92cm_1[0], specS_92cm_1[1]+specF_92cm_1[1])
	L92_er2 = -(L92cm[1]-L92_ex2)/L92_ex2
	L92_error2 = L92cm[2]/L92_ex2
	
	AX92e.fill_between(L92cm[0], L92_er1, L92_er2, facecolor=color, edgecolor=color, alpha=0.65, hatch=HATCH[n])
	
	L92_err1 = np.maximum(L92_er1+L92_error1, L92_er2+L92_error2)
	L92_err2 = np.minimum(L92_er1-L92_error1, L92_er2-L92_error2)
	
	AX92e.fill_between(L92cm[0], L92_err1, L92_err2, facecolor=color, edgecolor='k', alpha=0.1)
	
	#22cm
	L22_ex1 = np.interp(L22cm[0], specS_22cm[0], specS_22cm[1]+specF_22cm[1])
	L22_er1 = -(L22cm[1]-L22_ex1)/L22_ex1
	L22_error1 = L22cm[2]/L22_ex1
	
	L22_ex2 = np.interp(L22cm[0], specS_22cm_1[0], specS_22cm_1[1]+specF_22cm_1[1])
	L22_er2 = -(L22cm[1]-L22_ex2)/L22_ex2
	L22_error2 = L22cm[2]/L22_ex2
	
	AX22e.fill_between(L22cm[0], L22_er1, L22_er2, facecolor=color, edgecolor=color, alpha=0.65, hatch=HATCH[n])
	
	L22_err1 = np.maximum(L22_er1+L22_error1, L22_er2+L22_error2)
	L22_err2 = np.minimum(L22_er1-L22_error1, L22_er2-L22_error2)
	
	AX22e.fill_between(L22cm[0], L22_err1, L22_err2, facecolor=color, edgecolor='k', alpha=0.1)
	
	#6cm
	L6_ex1 = np.interp(L6cm[0], specS_6cm[0], specS_6cm[1]+specF_6cm[1])
	L6_er1 = -(L6cm[1]-L6_ex1)/L6_ex1
	L6_error1 = L6cm[2]/L6_ex1
	
	L6_ex2 = np.interp(L6cm[0], specS_6cm_1[0], specS_6cm_1[1]+specF_6cm_1[1])
	L6_er2 = -(L6cm[1]-L6_ex2)/L6_ex2
	L6_error2 = L6cm[2]/L6_ex2
	
	AX6e.fill_between(L6cm[0], L6_er1, L6_er2, facecolor=color, edgecolor=color, alpha=0.65, hatch=HATCH[n])
	
	L6_err1 = np.maximum(L6_er1+L6_error1, L6_er2+L6_error2)
	L6_err2 = np.minimum(L6_er1-L6_error1, L6_er2-L6_error2)
	
	AX6e.fill_between(L6cm[0], L6_err1, L6_err2, facecolor=color, edgecolor='k', alpha=0.1)
	
	#3cm
	L3_ex1 = np.interp(L3cm[0], specS_3cm[0], specS_3cm[1]+specF_3cm[1])
	L3_er1 = -(L3cm[1]-L3_ex1)/L3_ex1
	L3_error1 = L3cm[2]/L3_ex1
	
	L3_ex2 = np.interp(L3cm[0], specS_3cm_1[0], specS_3cm_1[1]+specF_3cm_1[1])
	L3_er2 = -(L3cm[1]-L3_ex2)/L3_ex2
	L3_error2 = L3cm[2]/L3_ex2
	
	AX3e.fill_between(L3cm[0], L3_er1, L3_er2, facecolor=color, edgecolor=color, alpha=0.65, hatch=HATCH[n])
	
	L3_err1 = np.maximum(L3_er1+L3_error1, L3_er2+L3_error2)
	L3_err2 = np.minimum(L3_er1-L3_error1, L3_er2-L3_error2)
	
	AX3e.fill_between(L3cm[0], L3_err1, L3_err2, facecolor=color, edgecolor='k', alpha=0.1)


# plt.suptitle(r'M82 Synchrotron $F(z)$ (Adebahr et al. 2013)')
ax00.set_ylabel(r'$F_\lambda (z)$ [Jy]')
ax10.set_ylabel(r'$F_\lambda (z)$ [Jy]')

ax00e.set_ylabel(r'$\frac{\Delta F_{\lambda}}{F_{\lambda}}$')
ax10e.set_ylabel(r'$\frac{\Delta F_{\lambda}}{F_{\lambda}}$')

# ax21.set_ylabel(r'$F(z)$ (Jy)')
# ax02.set_ylabel(r'$F(z)$ (Jy)')
# ax12.set_ylabel(r'$F(z)$ (Jy)')
# ax22.set_ylabel(r'$F(z)$ (Jy)')

# ax01.set_xlabel(r'$z$ (kpc)')
# ax02.set_xlabel(r'$z$ (kpc)')
# ax11.set_xlabel(r'$z$ (kpc)')
# ax12.set_xlabel(r'$z$ (kpc)')
ax10e.set_xlabel(r'$z$ [kpc]')
ax11e.set_xlabel(r'$z$ [kpc]')

# Turn off tick labels
ax01.set_yticklabels([])
ax11.set_yticklabels([])

ax01e.set_yticklabels([])
ax11e.set_yticklabels([])

ax00.set_yticklabels(['','',r'10$^{-4}$',r'10$^{-3}$',r'10$^{-2}$',r'10$^{-1}$',r'$10^0$'])
ax10.set_yticklabels(['','',r'10$^{-3}$',r'10$^{-2}$',r'10$^{-1}$',r'10$^{0}$',r'10$^{1}$'])

ax00.set_xticklabels([])
ax01.set_xticklabels([])

ax00e.set_xticklabels([])
ax01e.set_xticklabels([])

ax00e.set_yticks([-1,0,1])
ax10e.set_yticks([-1,0,1])
ax01e.set_yticks([-1,0,1])
ax11e.set_yticks([-1,0,1])

# ax00e.set_yticklabels(['','',r'$-$1','',0,'',1,''])
# ax10e.set_yticklabels(['','',r'$-$1','',0,'',1,''])

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
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.png', bbox_inches='tight', dpi=400)
		break
#plt.savefig('plot_gamma_'+CHANGED, dpi=400)
plt.close()


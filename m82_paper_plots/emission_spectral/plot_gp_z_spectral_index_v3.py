import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pyfits as py
from scipy import interpolate

def DoRotation(xspan, yspan, RotRad=0):
    """Generate a meshgrid and rotate it by RotRad radians."""

    # Clockwise, 2D rotation matrix
    RotMatrix = np.array([[np.cos(RotRad),  np.sin(RotRad)],
                          [-np.sin(RotRad), np.cos(RotRad)]])

    x, y = np.meshgrid(xspan, yspan)
    return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))

SL_SCALE = 0
E2FLAG = 1

DISTANCE = 3520

#Organizing files

dir_outfits = '../fits/raw/'
dir_comp = '../fits/'

dir_in = dir_comp

galdef_name = []
galdef_name_fill = []
galdef_name_mid = []
galdef_name_mid_raw = []
galdef_dir = '/mnt/c/Users/psyko/Physics/proj_gal/1904/m82_paper_plots/fits/'
galdef_raw = '/mnt/c/Users/psyko/Physics/proj_gal/1904/m82_paper_plots/fits/raw/'

GG = 'm82_05aff1017'
galdef_name.append(galdef_dir+GG+'/'+GG)

GG = 'm82_05b2ff1020'
galdef_name.append(galdef_dir+GG+'/'+GG)

GG = 'm82_06bff1020'
galdef_name.append(galdef_dir+GG+'/'+GG)

GG2 = 'm82_05aff1019'
galdef_name_fill.append(galdef_dir+GG2+'/'+GG2)

GG2 = 'm82_05b2ff1022'
galdef_name_fill.append(galdef_dir+GG2+'/'+GG2)

GG2 = 'm82_06bff1022'
galdef_name_fill.append(galdef_dir+GG2+'/'+GG2)

	# MID
GGm = 'm82_05aff1018'
galdef_name_mid.append(galdef_dir+GGm+'/'+GGm)

GGm = 'm82_05b2ff1021'
galdef_name_mid.append(galdef_dir+GGm+'/'+GGm)

GGm = 'm82_06bff1021'
galdef_name_mid.append(galdef_dir+GGm+'/'+GGm)

GGmr = 'm82_05aff1018'
galdef_name_mid_raw.append(galdef_raw+GGmr+'/'+GGmr)

GGmr = 'm82_05b2ff1021'
galdef_name_mid_raw.append(galdef_raw+GGmr+'/'+GGmr)

GGmr = 'm82_06bff1021'
galdef_name_mid_raw.append(galdef_raw+GGmr+'/'+GGmr)

LABEL_ALL = [r'A',r'B',r"B'"]
na_fac = 1.0

#Getting chi factors
num_graph = len(galdef_name)
chi = np.zeros(num_graph)
chi_fac = np.ones(num_graph)
for n in range(0,num_graph):
	chi_dat = np.loadtxt(galdef_name[n]+'chi.dat')
	chi_fac[n] = chi_dat[4]
	GAMMA_CHI = 1
	chi[n] = GAMMA_CHI*(chi_dat[1]+chi_dat[2])+chi_dat[3]
	
num_graph_fill = len(galdef_name_fill)
chi_fill = np.zeros(num_graph_fill)
chi_fac_fill = np.ones(num_graph_fill)
for n in range(0,num_graph_fill):
	chi_dat_fill = np.loadtxt(galdef_name_fill[n]+'chi.dat')
	chi_fac_fill[n] = chi_dat_fill[4]
	GAMMA_CHI = 1
	chi_fill[n] = GAMMA_CHI*(chi_dat_fill[1]+chi_dat_fill[2])+chi_dat_fill[3]
	
num_graph_mid = len(galdef_name_mid)
chi_mid = np.zeros(num_graph_mid)
chi_fac_mid = np.ones(num_graph_mid)
for n in range(0,num_graph_mid):
	chi_dat_mid = np.loadtxt(galdef_name_mid[n]+'chi.dat')
	chi_fac_mid[n] = chi_dat_mid[4]
	GAMMA_CHI = 1
	chi_mid[n] = GAMMA_CHI*(chi_dat_mid[1]+chi_dat_mid[2])+chi_dat_mid[3]
				

NORM_PLOT = 0
if NORM_PLOT == 1:
	FACTOR = np.ones(num_graph)
	FACTOR_fill = np.ones(num_graph_fill)
	FACTOR_mid = np.ones(num_graph_mid)
else:
	FACTOR = chi_fac
	FACTOR_fill = chi_fac_fill
	FACTOR_mid = chi_fac_mid

#constants
pc2cm = 3.08568025e18 #cm/pc
kpc2cm = 3.08568025e21 #cm/kpc
radian2arcsec = 648000./math.pi #acsec/radian
Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
c = 3.0e10 #cm/s
h = 6.62e-27 #erg s
q=4.8e-10 #esu
mp = 1.66e-24 #g
ev2ergs = 1.60218e-12

DISTANCE = 3520 #kpc
D_cm = DISTANCE*kpc2cm

###
### PLOT
###

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)

MSIZE = 5
LWIDTH = 3

fig = plt.figure()
f, ax = plt.subplots(1,1) 

for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(2)
ax.tick_params(which='major',width=1, length=8, labelsize='small')
ax.tick_params(which='minor',width=0.5, length=5)

cmap = plt.get_cmap('cubehelix')

color_dat= [cmap(0.85), cmap(0.6)] #0.65
WHICH_MODES = [2]
HATCHES_DAT = ['','','']
HATCHES = ['|','-','/']

for s in range(0,len(WHICH_MODES)):

	MODE = WHICH_MODES[s]
	color = color_dat[s]
	if MODE == 0:
		# 22cm - 92cm
		low_lambda = 92. #cm
		low_as = 30. #"
		low_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_N.dat').T
		low_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/92cm/92cm_S.dat').T
		high_lambda = 22.
		high_as = 15. #"
		high_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_N.dat').T
		high_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_S.dat').T
	elif MODE == 1:
		# 22cm - 6cm
		low_lambda = 22.
		low_as = 15. #"
		low_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_N.dat').T
		low_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/22cm/22cm_S.dat').T
		high_lambda = 6.
		high_as = 7. #"
		high_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_N_v0.dat').T
		high_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_S_v0.dat').T
	elif MODE == 2:
		# 6cm - 3cm
		low_lambda = 6.
		low_as = 7. #"
		low_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_N_v0.dat').T
		low_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/6cm_S_v0.dat').T
		high_lambda = 3.
		high_as = 4.5 #"
		high_N = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_N.dat').T
		high_S = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/3cm/3cm_S.dat').T

	low_f = c/low_lambda
	high_f = c/high_lambda

	low_S[0] = -low_S[0]
	low_S = low_S[:,::-1]
	low_all = np.append(low_S, low_N, axis=-1)
	low_all[1:] = np.log(low_all[1:]/low_as)
	
	high_S[0] = -high_S[0]
	high_S = high_S[:,::-1]
	high_all = np.append(high_S, high_N, axis=-1)
	high_all[1:] = np.log(high_all[1:]/high_as)
	
	print high_all
	print low_all
	low_func = interpolate.interp1d(low_all[0], low_all, axis=-1, kind='cubic')#, fill_value='extrapolate')
	high_func = interpolate.interp1d(high_all[0], high_all, axis=-1, kind='cubic')#, fill_value='extrapolate')

	z_min = max(min(low_all[0]), min(high_all[0]))
	z_max = min(max(low_all[0]), max(high_all[0]))

	z_num = 1001
	z = np.linspace(z_min, z_max, z_num)

	low = (low_func(z))#+low_func(-z))/2.
	high = (high_func(z))#+high_func(-z))/2.
	
	print low
	print high
	# print (high-low)
	# print np.log(high_f/low_f)

	alpha = (high[1]-low[1])/np.log(high_f/low_f)
	alpha_u = (high[3]-low[2])/np.log(high_f/low_f)
	alpha_l = (high[2]-low[3])/np.log(high_f/low_f)
	alpha_u2 = (high[2]-low[2])/np.log(high_f/low_f)
	alpha_l2 = (high[3]-low[3])/np.log(high_f/low_f)

	print alpha
	print alpha_u
	print alpha_l

	ax.plot(z, alpha, color=color, linestyle='-', alpha=0.9, linewidth=LWIDTH, zorder=0)
	ax.fill_between(z, alpha_u, alpha_l, facecolor=color, edgecolor=color, alpha=0.7, hatch=HATCHES_DAT[s], linewidth=LWIDTH)
	
	# ax.plot(-z, alpha, color=color, linestyle='-', alpha=0.3, linewidth=LWIDTH)
	# ax.fill_between(-z, alpha_u, alpha_l, facecolor=color, edgecolor=color, alpha=0.3, hatch=HATCHES_DAT[s], linewidth=LWIDTH)


	for n in range(0,num_graph):
		for r in range(0,2):
			LABEL = [r'A',r'B',r"B'"]
			
			###
			###
			### RADIO FUNCTION OF Z
			###
			###
			###
			if r == 1:
				#load data
				zS = np.loadtxt(galdef_name_fill[n]+'_z_freef.dat')
				specS = np.loadtxt(galdef_name_fill[n]+'_z0_nu_freef.dat')
				specS*=FACTOR[n]
				
				zL = np.loadtxt(galdef_name_fill[n]+'_z_synch.dat')
				specL = np.loadtxt(galdef_name_fill[n]+'_z0_nu_synch.dat')
				specL*=FACTOR[n]
				
				zS_o = np.loadtxt(galdef_name_fill[n]+'_z_freef.dat')
				specS_o = np.loadtxt(galdef_name_fill[n]+'_zoff_nu_freef.dat')
				specS_o*=FACTOR[n]
				
				zL_o = np.loadtxt(galdef_name_fill[n]+'_z_synch.dat')
				specL_o = np.loadtxt(galdef_name_fill[n]+'_zoff_nu_synch.dat')
				specL_o*=FACTOR[n]
			else:
				#load data
				zS = np.loadtxt(galdef_name[n]+'_z_freef.dat')
				specS = np.loadtxt(galdef_name[n]+'_z0_nu_freef.dat')
				specS*=FACTOR[n]
				
				zL = np.loadtxt(galdef_name[n]+'_z_synch.dat')
				specL = np.loadtxt(galdef_name[n]+'_z0_nu_synch.dat')
				specL*=FACTOR[n]
				
				zS_o = np.loadtxt(galdef_name[n]+'_z_freef.dat')
				specS_o = np.loadtxt(galdef_name[n]+'_zoff_nu_freef.dat')
				specS_o*=FACTOR[n]
				
				zL_o = np.loadtxt(galdef_name[n]+'_z_synch.dat')
				specL_o = np.loadtxt(galdef_name[n]+'_zoff_nu_synch.dat')
				specL_o*=FACTOR[n]
				
			dzS = zS[1]-zS[0]
			dzL = zL[1]-zL[0]
			
			#interpolate data
			data_freq = np.array([c/92, c/22, c/6, c/3])
			y_int = np.array([30, 15, 7, 4.5])
			
			for j in range(0,len(data_freq)):
				NU = data_freq[j] #Hz 
			
				y_integrated = y_int[j]
				#INTERP GALPROP
				# S
				specS_temp = specS*1.e23
				
				data_z = zS
				data = specS_temp 
				numpoints = data_z.shape[0]
				SPEC = np.zeros([2,numpoints])
				SPEC[0] = data_z
				SPEC[0]*=DISTANCE/radian2arcsec
				SPEC[1] = data[j]
				SPEC[1]*= y_integrated
				specS_temp = SPEC 
				
				# L
				specL_temp = specL*1.e23
				
				data_z = zL
				data = specL_temp 
				numpoints = data_z.shape[0]
				SPEC = np.zeros([2,numpoints])
				SPEC[0] = data_z
				SPEC[0]*=DISTANCE/radian2arcsec
				SPEC[1] = data[j]
				SPEC[1]*= y_integrated
				specL_temp = SPEC 
			
				if j == 0 and r == 0:
					specS_92cm = specS_temp
					specL_92cm = specL_temp
					# specS_o_92cm = specS_o_temp
					# specL_o_92cm = specL_o_temp
				elif j == 1 and r == 0:
					specS_22cm = specS_temp
					specL_22cm = specL_temp
					# specS_o_22cm = specS_o_temp
					# specL_o_22cm = specL_o_temp
				elif j == 2 and r == 0:
					specS_6cm = specS_temp
					specL_6cm = specL_temp
				elif j == 3 and r == 0:
					specS_3cm = specS_temp
					specL_3cm = specL_temp
				elif j == 0 and r == 1:
					specS_92cm2 = specS_temp
					specL_92cm2 = specL_temp
					# specS_o_92cm = specS_o_temp
					# specL_o_92cm = specL_o_temp
				elif j == 1 and r == 1:
					specS_22cm2 = specS_temp
					specL_22cm2 = specL_temp
					# specS_o_22cm = specS_o_temp
					# specL_o_22cm = specL_o_temp
				elif j == 2 and r == 1:
					specS_6cm2 = specS_temp
					specL_6cm2 = specL_temp
				elif j == 3 and r == 1:
					specS_3cm2 = specS_temp
					specL_3cm2 = specL_temp
			
		##########
		## Plot
		color = cmap(1-float(n+1)/(num_graph+1.))
		
		ALPHA = 0.60
		
		if MODE == 0:
			spec_22_92 = np.log((specS_22cm[1]+specL_22cm[1])/(specS_92cm[1]+specL_92cm[1])*30./15.)/np.log(92./22.)
			spec_22_92_2 = np.log((specS_22cm2[1]+specL_22cm2[1])/(specS_92cm2[1]+specL_92cm2[1])*30./15.)/np.log(92./22.)
			ax.fill_between(specS_92cm[0], spec_22_92, spec_22_92_2, facecolor=color, edgecolor=color, alpha=ALPHA, hatch=HATCHES[s], linewidth=LWIDTH)
			
			ax.text(0.8, 0.66, r'$92\,\mathrm{cm-}22\,\mathrm{cm}$', ha='center', va='center', transform=ax.transAxes)
			
			L92_cm_beam = 43.1/radian2arcsec*DISTANCE
			ax.errorbar([1.8], [-0.35],xerr=L92_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)
			
			ax.set_xlim([-3,3])
			ax.set_ylim([-2.5,1.0])
			
		elif MODE == 1:
			spec_6_22 = np.log((specS_6cm[1]+specL_6cm[1])/(specS_22cm[1]+specL_22cm[1])*15./7.)/np.log(22./6.)
			spec_6_22_2 = np.log((specS_6cm2[1]+specL_6cm2[1])/(specS_22cm2[1]+specL_22cm2[1])*15./7.)/np.log(22./6.)
			ax.fill_between(specS_92cm[0], spec_6_22, spec_6_22_2, facecolor=color, edgecolor=color, alpha=ALPHA, hatch=HATCHES[s], linewidth=LWIDTH)
			
			ax.text(0.8, 0.66, r'$22\,\mathrm{cm-}6\,\mathrm{cm}$', ha='center', va='center', transform=ax.transAxes)
			
			L22_cm_beam = 12.7/radian2arcsec*DISTANCE
			ax.errorbar([1.8], [-0.95],xerr=L22_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)
			
			ax.set_xlim([-3,3])
			ax.set_ylim([-2.5,0.0])
			
		elif MODE == 2:
			spec_3_6 = np.log((specS_3cm[1]+specL_3cm[1])/(specS_6cm[1]+specL_6cm[1])*7./4.5)/np.log(6./2.)
			spec_3_6_2 = np.log((specS_3cm2[1]+specL_3cm2[1])/(specS_6cm2[1]+specL_6cm2[1])*7./4.5)/np.log(6./2.)
			ax.fill_between(specS_92cm[0], spec_3_6, spec_3_6_2, facecolor=color, edgecolor=color, alpha=ALPHA, hatch=HATCHES[s], linewidth=LWIDTH, label=LABEL[n])
			ax.legend(loc='upper left', prop={'size':15})
			
			ax.text(0.8, 0.66, r'$6\,\mathrm{cm-}3\,\mathrm{cm}$', ha='center', va='center', transform=ax.transAxes)
			
			L6_cm_beam = 12.5/radian2arcsec*DISTANCE
			ax.errorbar([1.8], [-0.625],xerr=L6_cm_beam/2., color='black',linewidth=3, marker='.', markersize=1)
			
			ax.set_xlim([-3,3])
			ax.set_ylim([-2.5,0.5])

			

ax.set_xlabel(r'$z$ [kpc]')
ax.set_ylabel(r'$\alpha$')

plt.tight_layout()
plot_name = 'plot_z_spectral_index_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.png', bbox_inches='tight', dpi=400)
		break
plt.close()




















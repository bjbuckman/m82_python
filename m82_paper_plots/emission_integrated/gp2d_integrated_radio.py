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
CHANGED = 'integrated_radio'
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

ax2 = plt.subplot(gs[0:2,0])
ax_e = plt.subplot(gs[2,0])#, sharex = ax1)

for axis in ['top','bottom','left','right']:
	ax2.spines[axis].set_linewidth(2)
	ax_e.spines[axis].set_linewidth(2)
ax2.tick_params(which='major',width=1, length=8, labelsize='small')
ax2.tick_params(which='minor',width=0.5, length=5)
ax_e.tick_params(which='major',width=1, length=8, labelsize='small')
ax_e.tick_params(which='minor',width=0.5, length=5)


###
###
### RADIO
###
###
###

MSIZE = 5
LWIDTH = 3

radio_data = np.loadtxt('/mnt/c/Users/psyko/Physics/proj_gal/m82_data/integrated_radio/radio_m82.dat')
radio_data = radio_data[radio_data[:,0].argsort()]
radio_data = radio_data.T
# print radio_data[0]
# radio_data = np.sort(radio_data,axis=0)

c=3.e10

HATCH = ['-', '|','/']

for n in range(0,num_graph):
	
	LABEL = LABEL_ALL[n]
	
	###
	###
	### RADIO
	###
	###
	###
	
	synch = np.loadtxt(galdef_name[n]+'synch_flux.dat')
	freef = np.loadtxt(galdef_name[n]+'freefree_flux.dat')
	synchL = np.loadtxt(galdef_name[n]+'synch_flux_na.dat')
	freefL = np.loadtxt(galdef_name[n]+'freefree_flux_na.dat')
	
	fill_synch = np.loadtxt(galdef_name_fill[n]+'synch_flux.dat')
	fill_freef = np.loadtxt(galdef_name_fill[n]+'freefree_flux.dat')
	fill_synchL = np.loadtxt(galdef_name_fill[n]+'synch_flux_na.dat')
	fill_freefL = np.loadtxt(galdef_name_fill[n]+'freefree_flux_na.dat')
	
	color = cmap(1-float(n+1)/(num_graph+1.))
	ALPHA=0.6
	
	#synchrotron unabsorbed
	flux_synch1 = synchL[1]*FACTOR[n]*1e23*na_fac
	flux_synch2 = fill_synchL[1]*FACTOR_fill[n]*1e23
	
	ax2.loglog(synchL[0], flux_synch1, color=color, alpha=ALPHA, linestyle='--')
	ax2.loglog(synchL[0], flux_synch2, color=color, alpha=ALPHA, linestyle='--')
	# ax2.fill_between(synchL[0], flux_synch1, flux_synch2, facecolor=color, edgecolor=color, alpha=0.2, hatch='.')
	# ax2.loglog(synchL[0], synchL[1]*FACTOR[n]*1e23*na_fac, c=color, linestyle='--', label=None, linewidth=LWIDTH, alpha=0.75)
	
	#free-free
	flux_freef1 = freef[1]*FACTOR[n]*1e23
	flux_freef2 = fill_freef[1]*FACTOR_fill[n]*1e23
	
	ax2.loglog(synchL[0], flux_freef1, color=color, alpha=ALPHA, linestyle=':')
	ax2.loglog(synchL[0], flux_freef2, color=color, alpha=ALPHA, linestyle=':')
	# ax2.fill_between(freef[0], flux_freef1, flux_freef2, facecolor=color, edgecolor=color, alpha=0.2, hatch='/')
	# ax2.loglog(freef[0], freef[1]*FACTOR[n]*1e23, c=color, linestyle='-.', label=None, linewidth=LWIDTH, alpha=0.75)
	
	#total
	flux_total1 = (synch[1]*FACTOR[n]+freef[1]*FACTOR[n])*1e23
	flux_total2 = (fill_synch[1]*FACTOR_fill[n]+fill_freef[1]*FACTOR_fill[n])*1e23
	ax2.fill_between(synch[0], flux_total1, flux_total2, facecolor=color, edgecolor=color, alpha=0.65)
	# ax2.loglog(synch[0], flux_total1, label=LABEL, c=color, linewidth=LWIDTH, alpha=0.975)
	# ax2.loglog(synch[0], flux_total2, label=LABEL, c=color, linewidth=LWIDTH, alpha=0.975)

	radio_ex1 = np.interp(radio_data[0], synch[0], flux_total1)
	radio_er1 = -(radio_data[1] - radio_ex1)/radio_ex1
	radio_error1 = radio_data[2]/radio_ex1
	
	radio_ex2 = np.interp(radio_data[0], synch[0], flux_total2)
	radio_er2 = -(radio_data[1] - radio_ex2)/radio_ex2
	radio_error2 = radio_data[2]/radio_ex2
	
	ax_e.fill_between(radio_data[0], radio_er1, radio_er2, facecolor=color, edgecolor=color, alpha=0.65, hatch=HATCH[n])
	
	err1 = np.maximum(radio_er1+radio_error1, radio_er2+radio_error2)
	err2 = np.minimum(radio_er1-radio_error1, radio_er2-radio_error2)
	ax_e.fill_between(radio_data[0], err1, err2, facecolor=color, edgecolor='k', alpha=0.1)
	# ax_e.errorbar(radio_data[0], radio_er1, yerr=radio_error1, fmt='o', markersize=MSIZE, c=color, alpha=0.5)
	
	

ax_e.plot((1e0,1e15),(0, 0), 'k', linestyle=':')
ax_e.set_ylim([-.2, .2])
ax_e.set_xlim([8.e+7,1.e+11])
ax_e.set_xscale('log')
ax_e.set_xlabel(r'$\nu$ [Hz]')
ax_e.set_ylabel(r'$\frac{\Delta F_{\nu}}{F_{\nu}}$')

#PLOTTING DATA
ax2.errorbar(radio_data[0], radio_data[1], yerr=radio_data[2], fmt='o', markersize=MSIZE, c='black', alpha=0.5)
#ax2.set_title('Radio Flux for '+str(CHANGED))
# ax2.set_xlabel(r'$\nu$ [Hz]')
ax2.set_xscale('log')
ax2.set_ylabel(r'$F_{\nu}$ [Jy]')
ax2.set_yscale('log')
ax2.set_xlim([8.e+7,1.e+11])
ax2.set_ylim([3.e-1, 3.e1])
# ax2.set_title(r'M82 Integrated Radio Spectrum')
	
ax2.set_xticklabels([])
ax_e.set_yticklabels(['','',-0.1,'',0.,'',0.1,'',''])

ax2.legend(loc='lower left', prop={'size':15})

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




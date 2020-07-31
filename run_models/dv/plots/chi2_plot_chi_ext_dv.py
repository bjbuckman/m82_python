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

# ffval = 3
gamma_fac = 1.
radio_fac = 1.

chi_files = ['m82_02b1010chi2.fits.gz']
inj_files = ['m82_02b1010inj.fits.gz']

# chi_files = ['m82_21a1010chi2.fits.gz', 'm82_21a2020chi2.fits.gz']
# inj_files = ['m82_21a1010inj.fits.gz', 'm82_21a2020inj.fits.gz']

# chi_files = ['m82_21b1010chi2.fits.gz', 'm82_21b2020chi2.fits.gz']
# inj_files = ['m82_21b1010inj.fits.gz', 'm82_21b2020inj.fits.gz']

#chi_plot can only have 3 values!!
chi_plot = [4,6,8,9]

for x in range(0,len(chi_plot)):
	for n in range(0,len(chi_files)):
		
		#######################################
		## READ DATA
		##
		## 
		#######################################
		fitsfile = chi_files[n]

		hdulist = py.open(fitsfile)
		head=hdulist[0].header
		scidata = hdulist[0].data

		min1 = head['CRVAL1']
		min2 = head['CRVAL2']
		# min3 = head['CRVAL3'] 
		min4 = head['CRVAL4']

		del1 = head['CDELT1']
		del2 = head['CDELT2']
		# del3 = head['CDELT3'] 
		del4 = head['CDELT4']

		num1 = head['NAXIS1']
		num2 = head['NAXIS2']
		# num3 = head['NAXIS3']
		# num4 = head['NAXIS4']
		
		var1 = min1*10.**(np.arange(num1)*del1)
		var2 = min2+np.arange(num2)*del2
		# var1 = min1*10.**(np.arange(num1)*del1)*1.e6
		# var2 = min2*10.**(np.arange(num2)*del2)
		# var3 = min3*10.**(np.arange(num3)*del3)
		# var4 = min4*10.**(np.arange(num4)*del4)

		grid_val_var1, grid_val_var2 = np.meshgrid(var1,var2)
		# scidata = np.swapaxes(scidata, 1, 2)
		# scidata = np.swapaxes(scidata, 2, 3)
		scidata[1]*= radio_fac
		scidata[2]*= gamma_fac

		# print head
		#######################################
		## 
		## Injection
		fitsfile = inj_files[n]

		hdulist = py.open(fitsfile)
		head=hdulist[0].header
		scidata1 = hdulist[0].data

		num4_inj = head['NAXIS3']
		
		scidata1 = np.swapaxes(scidata1, 0, 1)
		scidata1 = np.swapaxes(scidata1, 1, 2)
		# scidata1 = np.swapaxes(scidata1, 0, 3)

		##
		##
		## WHICH CHI^2 TO MINIMIZE
		##
		##
		CHI_MODEL = 2
		chi_mode = chi_plot[x]
		
		grid_dat = np.zeros([num2,num1])
		grid_dat_inj = np.zeros([num2,num1,num4_inj]) 
		grid_dat_ff = np.zeros([num2,num1])
		if CHI_MODEL == 1: 
			for ax2 in range(0,num2):
				for ax1 in range(0,num1):
					grid_dat[ax2][ax1] = scidata[1][ax2][ax1]+scidata[2][ax2][ax1]
					grid_dat_inj[ax2][ax1] = scidata1[ax2][ax1]
		elif CHI_MODEL == 2: 
			for ax2 in range(0,num2):
				for ax1 in range(0,num1):
					grid_dat[ax2][ax1] = scidata[chi_mode][ax2][ax1]#+scidata[2][ax2][ax1]
					grid_dat_inj[ax2][ax1] = scidata1[ax2][ax1]
		elif CHI_MODEL == 100: 
			for ax2 in range(0,num2):
				for ax1 in range(0,num1):
					grid_dat[ax2][ax1] = scidata[chi_mode][ax2][ax1]#+scidata[2][ax2][ax1]
					grid_dat_inj[ax2][ax1] = scidata1[ax2][ax1]
		else:
			print 'You have failed me once more. Check CHI_MODEL'
			
		arr_val_var1_temp = grid_val_var1.flatten(order='F')
		arr_val_var2_temp = grid_val_var2.flatten(order='F')
		arr_dat_temp = grid_dat.flatten(order='F')
		grid_dat_inj = np.swapaxes(grid_dat_inj, 0,1)
		arr_dat_inj_temp = grid_dat_inj.reshape(-1, grid_dat_inj.shape[-1], order='F')
		arr_dat_ff_temp = grid_dat_ff.flatten(order='F')
		
		if n == 0:
			arr_var1 = arr_val_var1_temp
			arr_var2 = arr_val_var2_temp
			arr_dat = arr_dat_temp
			arr_dat_inj = arr_dat_inj_temp
			arr_dat_ff = arr_dat_ff_temp
		else:
			arr_var1 = np.concatenate((arr_var1, arr_val_var1_temp), axis=0)
			arr_var2 = np.concatenate((arr_var2, arr_val_var2_temp), axis=0)
			arr_dat = np.concatenate((arr_dat, arr_dat_temp), axis=0)
			arr_dat_inj = np.concatenate((arr_dat_inj, arr_dat_inj_temp), axis=0)
			arr_dat_ff = np.concatenate((arr_dat_ff, arr_dat_ff_temp), axis=0)


	###
	###
	### INTERPOLATION
	###
	###

	plot_min_var1 = np.amin(arr_var1)
	plot_max_var1 = np.amax(arr_var1)

	plot_min_var2 = np.amin(arr_var2)
	plot_max_var2 = np.amax(arr_var2)

	VAR1 = np.logspace(np.log10(plot_min_var1), np.log10(plot_max_var1), num=10)
	VAR2 = np.linspace(plot_min_var2, plot_max_var2, num=10)

	VAR_ARR1, VAR_ARR2 = np.meshgrid(VAR1,VAR2)

	# interpDAT = interpolate.Rbf(arr_var1, arr_var2, np.log10(arr_dat), smooth=1.e0, epsilon=1.e-1)

	# DAT = interpDAT(VAR_ARR1, VAR_ARR2)
	# DAT = 10.**DAT

	INTERP = 'linear'

	DAT = interpolate.griddata((np.log(arr_var1), arr_var2), arr_dat, (np.log(VAR_ARR1), VAR_ARR2), method=INTERP)
	DAT_INJ = interpolate.griddata((np.log(arr_var1), arr_var2), arr_dat_inj, (np.log(VAR_ARR1), VAR_ARR2), method=INTERP)
	# print DAT

	DAT_INJ = np.swapaxes(DAT_INJ,0,2)
	DAT_INJ = np.swapaxes(DAT_INJ,1,2)
	
	if x == 0:
		DATA0 = DAT
		INJE0 = DAT_INJ
	elif x == 1:
		DATA1 = DAT
		INJE1 = DAT_INJ
	elif x == 2:
		DATA2 = DAT
		INJE2 = DAT_INJ
	elif x == 3:
		DATA3 = DAT
		INJE3 = DAT_INJ
	else:
		print 'x not right size'

MAX_CHI = 1.5e2 #min(min(DATA0.max(), DATA1.max(), DATA2.max()), 1.e1)

print DATA0.max(), DATA1.max(), DATA2.max(), DATA3.max()
print DATA0.min(), DATA1.min(), DATA2.min(), DATA3.min()

MIN_CHI = min(DATA0.min(), DATA1.min(), DATA2.min(), DATA3.min())

# DATA0 = (DATA0.flatten()-MIN_CHI)/MAX_CHI
# DATA1 = (DATA1.flatten()-MIN_CHI)/MAX_CHI
# DATA2 = (DATA2.flatten()-MIN_CHI)/MAX_CHI

# DATA0 = (DATA0.flatten())/MAX_CHI
# DATA1 = (DATA1.flatten())/MAX_CHI
# DATA2 = (DATA2.flatten())/MAX_CHI

DATA0 = (DATA0.flatten())
DATA1 = (DATA1.flatten())
DATA2 = (DATA2.flatten())
DATA3 = (DATA3.flatten())

individual = -1

if individual == 0:
	DATA_tot = DATA0
	print DATA_tot
	DATA_tot/= 50
	
	DATA0 = np.zeros(len(DATA0))
	DATA1 = np.ones(len(DATA1))
	DATA2 = np.ones(len(DATA2))
	DATA3 = np.ones(len(DATA3))

elif individual == 1:
	DATA_tot = DATA1
	print DATA_tot
	DATA_tot/= 50
	
	DATA0 = np.ones(len(DATA0))
	DATA1 = np.zeros(len(DATA1))
	DATA2 = np.ones(len(DATA2))
	DATA3 = np.ones(len(DATA3))
	
elif individual == 2:
	DATA_tot = DATA2
	print DATA_tot
	DATA_tot/= 50
	
	DATA0 = np.ones(len(DATA0))
	DATA1 = np.ones(len(DATA1))
	DATA2 = np.zeros(len(DATA2))
	DATA3 = np.ones(len(DATA3))
	
elif individual == 3:
	DATA_tot = abs(DATA3)
	print DATA_tot
	DATA_tot/= 50
	
	DATA0 = np.zeros(len(DATA0))
	DATA1 = np.zeros(len(DATA1))
	DATA2 = np.zeros(len(DATA2))
	DATA3 = np.zeros(len(DATA3))
	
else:
	DATA_tot = DATA0+DATA1+DATA2+DATA3
	print DATA_tot
	DATA_tot/= 300. #85.#e2

	# DATA_tot-=DATA_tot.min()

	DATA_MIN = np.maximum(DATA0,DATA1)
	DATA_MIN = np.maximum(DATA2,DATA_MIN)
	DATA_MIN = np.maximum(DATA3,DATA_MIN)

	DATA0_arg = np.where(DATA0 == DATA_MIN)
	DATA1_arg = np.where(DATA1 == DATA_MIN)
	DATA2_arg = np.where(DATA2 == DATA_MIN)
	DATA3_arg = np.where(DATA3 == DATA_MIN)

	DATA0 = np.ones(len(DATA0))
	DATA1 = np.ones(len(DATA1))
	DATA2 = np.ones(len(DATA2))
	DATA3 = np.ones(len(DATA3))

	DATA0[DATA0_arg] = 0 
	DATA1[DATA1_arg] = 0 
	DATA2[DATA2_arg] = 0 
	DATA3[DATA3_arg] = 0

	DATA0[DATA3_arg] = 0 
	DATA1[DATA3_arg] = 0 
	DATA2[DATA3_arg] = 0 

	DATA0[DATA0>=1] = 1
	DATA1[DATA1>=1] = 1
	DATA2[DATA2>=1] = 1
	DATA3[DATA3>=1] = 1
	DATA_tot[DATA_tot>=1] = 1

# COLOR = np.array([DATA0, DATA1, DATA2, DATA_tot]).T
COLOR = np.array([DATA0, DATA1, DATA2, DATA_tot]).T
COLOR = 1-COLOR
COLOR = COLOR**0.3

DATA_tot = 1-DATA_tot
# DATA_tot = DATA_tot**0.3

# print DATA0, DATA1, DATA2, DATA3
# COLOR = np.rollaxis(COLOR,0,3)
# print COLOR.shape
# print COLOR
		
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
ax.tick_params(which='major',width=1., length=8, labelsize='small', direction='in', color='white')
ax.tick_params(which='minor',width=0.5, length=5, direction='in', color='white')

MSIZE = 5
LWIDTH = 3
WIDTH = 3


## PLOT 1
ax.set_axis_bgcolor('black')

loglow1 = 1
loghigh1 = 5
V = np.logspace(loglow1, loghigh1, 100)

# CF0 = ax.contourf(VAR1, VAR2, DAT, cmap=cmap_gamma)
# CF0 = ax.contourf(VAR1, VAR2, DATA0, cmap=cmap_gamma, levels=V, norm=LogNorm())
# CF0 = ax[0].contourf(var1, var2, grid_dat, cmap=cmap_gamma, levels=V, norm=LogNorm())

VAR_ARR1 = VAR_ARR1.flatten()
VAR_ARR2 = VAR_ARR2.flatten()
print VAR_ARR1.shape, VAR_ARR2.shape, COLOR.shape
ax.scatter(VAR_ARR1, VAR_ARR2, c=COLOR, s=1700*DATA_tot)#750

# ax.set_aspect('equal')
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim([1.e25,1.e30])
ax.set_ylim([0,2000])
# ax[0].set_xticks(np.arange(-80, 81, 40))
# ax[0].set_yticks(np.arange(-80, 81, 40))

# ax.tick_params(axis='x', colors='black', direction='out', which='major', size=5, width=1)
# ax.tick_params(axis='y', colors='black', direction='out', which='major', size=5, width=1)
# ax.tick_params(axis='x', colors='black', direction='out', which='minor', size=5, width=1)
# ax.tick_params(axis='y', colors='black', direction='out', which='minor', size=5, width=1)

ax.set_xlabel(r'$D_{xx}$ [cm$^2$ s$^{-1}$]')
ax.set_ylabel(r'$V_0$ [km s$^{-1}$]')
# ax.set_title(r'Extended Emission $\chi^2$', verticalalignment='bottom')

# loglow1 = 1
# loghigh1 = 5
# V = np.logspace(loglow1, loghigh1, 22)

# cf0 = plt.colorbar(CF0,ax=ax,pad=0,ticks=V, fraction=0.04763)
# cf0.set_label(r'$\chi_{\tt{weighted}}^2$', rotation=270, horizontalalignment='center', verticalalignment='bottom')

# cf0 = plt.colorbar(CF0,ax=ax,pad=0, fraction=0.04763)
# cf0.set_label(r'$\chi_{\tt{weighted}}^2$', rotation=270, horizontalalignment='center', verticalalignment='bottom')

plot_name = 'm82_chi_ext_plot'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break

plt.close()		


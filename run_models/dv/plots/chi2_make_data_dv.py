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
import scipy

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

#######################################
## READ DATA
##

# Folder
gp_dir = '../gp_out/'
out_file = 'chi2.fits.gz'
out_file1 = 'inj.fits.gz'

#CHANGE GALDEF PARAMETERS
init_galdef = 'm82_02b'
galdefname = init_galdef

suffix = 'dv'

##############################################
##############################################
##############################################


MODEL_NUMBER = 1

if MODEL_NUMBER == 1:
	num_params = 2
	num_galdef = np.zeros(num_params).astype(int)
	num_galdef[0] = 9 #
	num_galdef[1] = 11 #
	
	num_galdef_tot = int(np.prod(num_galdef))
	num_start = np.zeros(num_params).astype(int)
	num_start[0] = 10
	num_start[1] = 10
	
	## arr 1
	arr_1_min = 1.e26
	arr_1_max = 1.e30
	arr_1_diff = np.log10(arr_1_max/arr_1_min)/(num_galdef[0]-1)
	
	arr_1 = []
	for ii in range(0,num_galdef[0]):
		arr_1.append(arr_1_min*10.**(ii*arr_1_diff))
	
	## arr 2
	arr_2_min = 0.
	arr_2_max = 2000.
	arr_2_diff = (arr_2_max-arr_2_min)/(num_galdef[1]-1)
	
	arr_2 = []
	for ii in range(0,num_galdef[1]):
		arr_2.append(arr_2_min+ii*arr_2_diff)
		
else:
	print 'You have failed. Fix the MODEL_NUMBER'
	
# new_val.append(mag_arr) #B_0
# new_val.append(gas_arr) #rho_0
# new_val.append(ff_arr) #rho_0
	
###
### MAKE FITS FILE
###

var1_min = arr_1_min
var2_min = arr_2_min


var1_min_com = 'D_xx [cm^2 s^-1]'
var2_min_com = 'wind [km s^-1]'


var1_dif = arr_1_diff
var2_dif = arr_2_diff


var1_dif_com = 'dxx=min*10^(diff*i)'
var2_dif_com = 'wind=min+diff*i'

#
# Version without chi.dat
#

chi_data = np.zeros(num_galdef)
chi_data_radio = np.zeros(num_galdef)
chi_data_gamma = np.zeros(num_galdef)
chi_data_ext = np.zeros(num_galdef)
chi_data_92 = np.zeros(num_galdef)
chi_data_92_o = np.zeros(num_galdef)
chi_data_22 = np.zeros(num_galdef)
chi_data_22_o = np.zeros(num_galdef)
chi_data_6 = np.zeros(num_galdef)
chi_data_3 = np.zeros(num_galdef)

# rev_num_params = np.fliplr(num_params)
N=np.zeros(num_params).astype(int)
for gt in range(0,num_galdef_tot):
	
	#Getting suffix for files
	N_index=N+num_start
	N_suffix = ''
	for s in range(0,num_params):
		N_suffix += str(N_index[s])
	
	#New galdef file
	new_galdefname = galdefname+suffix+N_suffix
	# data = np.loadtxt(gp_dir+new_galdefname+'/').T
	
	if MODEL_NUMBER < 0:
		gp_file = open(gp_dir+new_galdefname+'/'+new_galdefname+'_analytics.txt', 'r')

		# print new_galdefname
		for line in gp_file:
			if line.startswith('chi^2        ='):
				chi_data[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi^2_radio  ='):
				chi_data_radio[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi^2_gamma  ='):
				chi_data_gamma[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi_ext^2    ='):
				chi_data_ext[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi92^2      ='):
				chi_data_92[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi92_o^2    ='):
				chi_data_92_o[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi22^2      ='):
				chi_data_22[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi22_o^2    ='):
				chi_data_22_o[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi6^2       ='):
				chi_data_6[N[0]][N[1]][N[2]] = float(line[15:-2])
			if line.startswith('chi3^2      ='):
				chi_data_3[N[0]][N[1]][N[2]] = float(line[15:-2])
		print new_galdefname+' done'
		gp_file.close()
	else:
		# gp_chi_data = np.loadtxt(gp_dir+new_galdefname+'/'+new_galdefname+'chi.dat')
		# chi_data[N[0]][N[1]] = gp_chi_data[0]
		# chi_data_radio[N[0]][N[1]] = gp_chi_data[1]
		# chi_data_gamma[N[0]][N[1]] = gp_chi_data[2]
		# chi_data_ext[N[0]][N[1]] = gp_chi_data[3]
		# chi_data_92[N[0]][N[1]] = gp_chi_data[4]
		# chi_data_92_o[N[0]][N[1]] = gp_chi_data[5]
		# chi_data_22[N[0]][N[1]] = gp_chi_data[6]
		# chi_data_22_o[N[0]][N[1]] = gp_chi_data[7]
		# chi_data_6[N[0]][N[1]] = gp_chi_data[8]
		# chi_data_3[N[0]][N[1]] = gp_chi_data[9]
		gp_file = open(gp_dir+new_galdefname+'/'+new_galdefname+'_analytics.txt', 'r')

		# print new_galdefname
		for line in gp_file:
			if line.startswith('chi^2        ='):
				chi_data[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi^2_radio  ='):
				chi_data_radio[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi^2_gamma  ='):
				chi_data_gamma[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi_ext^2    ='):
				chi_data_ext[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi92^2      ='):
				chi_data_92[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi92_o^2    ='):
				chi_data_92_o[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi22^2      ='):
				chi_data_22[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi22_o^2    ='):
				chi_data_22_o[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi6^2       ='):
				chi_data_6[N[0]][N[1]] = float(line[15:-2])
			if line.startswith('chi3^2       ='):
				chi_data_3[N[0]][N[1]] = float(line[15:-2])
		print new_galdefname+' done'
		gp_file.close()
		
	#Update N
	N[num_params-1] += 1
	for i in (num_params-1 -np.array(range(0,num_params))):
		if N[i] >= num_galdef[i]:
			N[i] = 0
			if i >= 1:
				N[i-1] += 1
	
chi_data = np.swapaxes(chi_data, 0, 1)
chi_data_radio = np.swapaxes(chi_data_radio, 0, 1)
chi_data_gamma = np.swapaxes(chi_data_gamma, 0, 1)
chi_data_ext = np.swapaxes(chi_data_ext, 0, 1)
chi_data_92 = np.swapaxes(chi_data_92, 0, 1)
chi_data_92_o = np.swapaxes(chi_data_92_o, 0, 1)
chi_data_22 = np.swapaxes(chi_data_22, 0, 1)
chi_data_22_o = np.swapaxes(chi_data_22_o, 0, 1)
chi_data_6 = np.swapaxes(chi_data_6, 0, 1)
chi_data_3 = np.swapaxes(chi_data_3, 0, 1)

image = np.array([chi_data, chi_data_radio, chi_data_gamma, chi_data_ext, chi_data_92, chi_data_92_o, chi_data_22, chi_data_22_o, chi_data_6, chi_data_3])

N_index=num_start
N_suffix = ''
for s in range(0,num_params):
	N_suffix += str(N_index[s])	
fileout = './'+galdefname+N_suffix+out_file
#Write array to FITS image--------------------------------------------------  
hduout = py.PrimaryHDU(image)
hduout.header['UNIT'  ] = 'chi^2'
hduout.header['CRVAL1'] = (var1_min, var1_min_com)
hduout.header['CRVAL2'] = (var2_min, var2_min_com)
# hduout.header['CRVAL3'] = (var3_min, var3_min_com)
hduout.header['CRVAL4'] = (1, '0=tot, 1=rad, 2=gam')
hduout.header['CDELT1'] = (var1_dif, var1_dif_com)
hduout.header['CDELT2'] = (var2_dif, var2_dif_com)
# hduout.header['CDELT3'] = (var3_dif, var3_dif_com)
hduout.header['CDELT4'] = (1, 'none')

hduoutlist = py.HDUList([hduout])
if os.access(fileout, os.F_OK ):  os.remove(fileout)
hduoutlist.writeto(fileout)   
print('FITS image output to '+str(fileout))

#
# INJECTION FITS FILE
#

inj_data = np.zeros(num_galdef)
inj_proton = np.zeros(num_galdef)
inj_electron = np.zeros(num_galdef)

N=np.zeros(num_params).astype(int)
for gt in range(0,num_galdef_tot):
	
	#Getting suffix for files
	N_index=N+num_start
	N_suffix = ''
	for s in range(0,num_params):
		N_suffix += str(N_index[s])
	
	#New galdef file
	new_galdefname = galdefname+suffix+N_suffix
	# data = np.loadtxt(gp_dir+new_galdefname+'/').T
	
	gp_file = open(gp_dir+new_galdefname+'/'+new_galdefname+'_analytics.txt', 'r')

	# print new_galdefname
	for line in gp_file:
		if line.startswith('-Total Injected     ='):
			inj_data[N[0]][N[1]] = float(line[21:-11])
		if line.startswith('Pion+Secondary / P+He           ='):
			inj_proton[N[0]][N[1]] = float(line[33:-2])
		if line.startswith('e+- Emission / All e+-          ='):
			inj_electron[N[0]][N[1]] = float(line[33:-2])
	print new_galdefname+' done'
	gp_file.close()
	
	#Update N
	N[num_params-1] += 1
	for i in (num_params-1 -np.array(range(0,num_params))):
		if N[i] >= num_galdef[i]:
			N[i] = 0
			if i >= 1:
				N[i-1] += 1
	
inj_data = np.swapaxes(inj_data, 0, 1)
inj_proton = np.swapaxes(inj_proton, 0, 1)
inj_electron = np.swapaxes(inj_electron, 0, 1)

SN_eff = 0.1
SN_nrg = 1.e51
inj_data*= 3.15e7/(SN_nrg*SN_eff)

inj_proton*=2.

image = np.array([inj_data, inj_proton, inj_electron])
# print image.shape
N_index=num_start
N_suffix = ''
for s in range(0,num_params):
	N_suffix += str(N_index[s])	
fileout = './'+galdefname+N_suffix+out_file1
#Write array to FITS image--------------------------------------------------  
hduout = py.PrimaryHDU(image)
hduout.header['UNIT'  ] = 'injection SN yr^-1'
hduout.header['effic' ] = (SN_eff, 'efficiency')
hduout.header['SNENER'] = (SN_nrg, 'energy per SN')
hduout.header['CRVAL1'] = (var1_min, var1_min_com)
hduout.header['CRVAL2'] = (var2_min, var2_min_com)
# hduout.header['CRVAL3'] = (var3_min, var3_min_com)
hduout.header['CRVAL4'] = (1, '0=nrg, 1=prot, 2=elec')
hduout.header['CDELT1'] = (var1_dif, var1_dif_com)
hduout.header['CDELT2'] = (var2_dif, var2_dif_com)
# hduout.header['CDELT3'] = (var3_dif, var3_dif_com)
hduout.header['CDELT4'] = (1, 'none')

hduoutlist = py.HDUList([hduout])
if os.access(fileout, os.F_OK ):  os.remove(fileout)
hduoutlist.writeto(fileout)   
print('FITS image output to '+str(fileout))



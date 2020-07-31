
import numpy as np
import math 
import os
import sys


def log_interp(e_arr,f_arr,e_val):
	low = f_arr[e_arr < e_val]
	e_low = e_arr[e_arr < e_val]
	low = low[-1]
	e_low = e_low[-1]
	
	high = f_arr[e_arr > e_val]
	e_high = f_arr[e_arr > e_val]
	high = high[0]
	e_high = e_high[0]
	
	# print(low, high, e_low, e_high, e_val)
	
	return low + (high-low)*math.log(e_val/e_low)/math.log(e_high/e_low)
	

h = 6.6260755e-27 #erg s
ev2ergs = 1.60218e-12
mev2ev = 1.e6

dir_outfits = './fits/raw/'
dir_comp = './fits/'

Emin = 1.e2 #MeV
Emax = 1.e3 #MeV


MODE_arr = np.array([0,1,2])
for ii in range(0,len(MODE_arr)):

	MODE = MODE_arr[ii]
	if MODE == 0:
		galdef_name = 'm82_05aff1018'
		MODEL = 'A'
	elif MODE == 1:
		galdef_name = 'm82_05b2ff1021'
		MODEL = 'B'
	elif MODE == 2:
		galdef_name = 'm82_06bff1021'
		MODEL = 'Bp'

	gamma_file = dir_outfits+galdef_name+'/'+galdef_name+'gamma_flux.dat'
	
	gamma_tot = np.loadtxt(gamma_file)*mev2ev 
	gamma_tot_mask = np.isfinite(gamma_tot[1])
	
	gamma_E = gamma_tot[0,gamma_tot_mask]
	gamma_flux = gamma_tot[1,gamma_tot_mask]*ev2ergs
	
	gflux_up = log_interp(gamma_E, gamma_flux, Emax)
	gflux_dn = log_interp(gamma_E, gamma_flux, Emin)
	
	gflux_tot = (gflux_up + gflux_dn)/2.*math.log(Emax/Emin)
	
	print(MODEL+' Total gamma-ray flux: {:1.2e}'.format(gflux_tot))
	
	
	radio_file = dir_outfits+galdef_name+'/'+galdef_name+'radio_flux.dat'
	
	radio_tot = np.loadtxt(radio_file)
	radio_tot_mask = np.isfinite(radio_tot[1])
	
	radio_freq = radio_tot[0,radio_tot_mask]
	radio_flux = radio_tot[1,radio_tot_mask]*radio_freq
	
	radio_E = radio_freq*h/ev2ergs 
	
	rflux_up = log_interp(radio_E, radio_flux, Emax)
	rflux_dn = log_interp(radio_E, radio_flux, Emin)
	
	rflux_tot = (rflux_up + rflux_dn)/2.*math.log(Emax/Emin)
	
	print(MODEL+' Total radio flux: {:1.2e}'.format(rflux_tot))
	
	print(MODEL+' Total flux: {:1.2e}'.format(rflux_tot+gflux_tot))
	
	
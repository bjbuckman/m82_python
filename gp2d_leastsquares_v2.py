import numpy as np
import math
import sys
import os
from scipy import interpolate

galdef_name = sys.argv[1]

extended_emiss = 1
SL_SCALE = 0

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
#CONVERT GP DATA TO OBSERVED DATA FORM

data_directory = '/users/PCON0003/cond0064/project_m82/m82_data/'
# data_directory = '/mnt/c/Users/psyko/Physics/proj_gal/m82_data/'

#FERMI DATA
# DATA = ph cm^-2 s^-1 
# Energy = MeV
fermi_data = np.loadtxt(data_directory+'integrated_gamma/m82_fermi.dat').T
num_fermi = fermi_data[0].size
emin = fermi_data[1]
emax = fermi_data[2]

#VERITAS DATA
# DATA = cm^-2 s^-1 [E dN/dE]
# Energy = MeV
verit_data = np.loadtxt(data_directory+'integrated_gamma/m82_veritas.dat').T
verit_upper = np.loadtxt(data_directory+'integrated_gamma/m82_veritas_upper.dat')
num_verit = verit_data[0].size
verit_e = verit_data[0]

#RADIO DATA
# DATA = Jy
# Frequency = Hz
radio_data = np.loadtxt(data_directory+'integrated_radio/radio_m82.dat').T
num_radio = radio_data[0].size
radio_f = radio_data[0]


chi_tot = 1.e25
chi_veritas = 1.e25
chi_radio = 1.e25
chi_fermi = 1.e25
small_factor = 1.

MAXPOW = 1000
FACTOR = 1.01**(np.array(range(-MAXPOW,MAXPOW)))


for FAC in FACTOR:
	chi = 0.
	#
	#GAMMA RAY STUFFS!!!
	#
	gamma_tot = np.loadtxt(galdef_name+'gamma_flux.dat')
	gamma_tot_mask = np.isfinite(gamma_tot[1])
	gamma_tot = gamma_tot.T[gamma_tot_mask].T
	gamma_tot[1] *= FAC
	# MeV ergs/cm^2/s
	gp_energy = gamma_tot[0]
	
	#FERMI JUNK!!!
	gp_fermi = gamma_tot[1]/gp_energy
	chi_f = 0.
	for i in range(0,num_fermi):
		gp_data = 0.
		emin_i = np.searchsorted(gp_energy-emin[i],0)
		emax_i = np.searchsorted(gp_energy-emax[i],0)
		#interpolate beginning and end
		e_temp = 10.**(math.log10(gp_fermi[emin_i])+math.log10(gp_fermi[emin_i+1]/gp_fermi[emin_i])*math.log10(emin[i]/gp_energy[emin_i])/math.log10(gp_energy[emin_i+1]/gp_energy[emin_i]))
		gp_data += math.sqrt(e_temp*gp_fermi[emin_i+1])*math.log10(gp_energy[emin_i+1]/emin[i])
		e_temp = 10.**(math.log10(gp_fermi[emax_i])+math.log10(gp_fermi[emax_i+1]/gp_fermi[emax_i])*math.log10(emax[i]/gp_energy[emax_i])/math.log10(gp_energy[emax_i+1]/gp_energy[emax_i]))
		gp_data += math.sqrt(e_temp*gp_fermi[emax_i])*math.log10(emax[i]/gp_energy[emax_i])
		#integrate
		for k in range(emin_i+1, emax_i):
			gp_data += math.sqrt(gp_fermi[k]*gp_fermi[k+1])*math.log10(gp_energy[k+1]/gp_energy[k])
		#print gp_data
		if fermi_data[0][i] > gp_data:
			sigma = abs(fermi_data[3][i])
			chi_f += (gp_data-fermi_data[0][i])**2./sigma**2.
		if fermi_data[0][i] < gp_data:
			sigma = abs(fermi_data[4][i])
			chi_f += (gp_data-fermi_data[0][i])**2./sigma**2.
	chi += chi_f

	#VERITAS JUNK!!!
	gp_verit = gamma_tot[1]/gp_energy	
	chi_v = 0.
	for i in range(0,num_verit):
		gp_data = 0.
		gp_i = np.searchsorted(gp_energy-verit_e[i],0)
		e_temp = 10.**(math.log10(gp_verit[gp_i])+math.log10(gp_verit[gp_i+1]/gp_verit[gp_i])*math.log10(verit_e[i]/gp_energy[gp_i])/math.log10(gp_energy[gp_i+1]/gp_energy[gp_i]))
		gp_data = e_temp #math.sqrt(e_temp*gp_verit[gp_i])*math.log10(verit_e[i]/gp_energy[gp_i])
		if verit_data[1][i] > gp_data:
			sigma = verit_data[2][i]-verit_data[1][i]
			chi_v += (gp_data-verit_data[1][i])**2./sigma**2.
		if verit_data[1][i] < gp_data:
			sigma = verit_data[3][i]-verit_data[1][i]
			chi_v += (gp_data-verit_data[1][i])**2./sigma**2.
	#VERITAS UPPER BOUND
	gp_i = np.searchsorted(gp_energy-verit_e[i],0)
	e_temp = 10.**(math.log10(gp_verit[gp_i])+math.log10(gp_verit[gp_i+1]/gp_verit[gp_i])*math.log10(verit_e[i]/gp_energy[gp_i])/math.log10(gp_energy[gp_i+1]/gp_energy[gp_i]))
	gp_data = e_temp
	if verit_upper[1] < gp_data:
		sigma = verit_upper[1]
		chi_v += (gp_data-verit_upper[1])**2./sigma**2.
	
	chi += chi_v
	
	#
	#RADIO STUFFS
	#
	synch = np.loadtxt(galdef_name+'synch_flux.dat')
	freef = np.loadtxt(galdef_name+'freefree_flux.dat')
	# Hz erg->Jy
	synch[1] *= 1.e23*FAC
	freef[1] *= 1.e23*FAC#*10.
	freq = synch[0]
	gp_synch = synch[1]+freef[1]
	chi_s = 0.
	for i in range(0,num_radio):
		gp_data = 0.
		if radio_f[i] > 1e9 and radio_f[i] < 300e9:
			gp_i = np.searchsorted(freq-radio_f[i],0) 
			if gp_i < freq.size-1:
				f_temp = 10.**(math.log10(gp_synch[gp_i])+math.log10(gp_synch[gp_i+1]/gp_synch[gp_i])*math.log10(radio_f[i]/freq[gp_i])/math.log10(freq[gp_i+1]/freq[gp_i]))
				gp_data = f_temp #math.sqrt(e_temp*gp_verit[gp_i])*math.log10(verit_e[i]/gp_energy[gp_i])
				if radio_data[1][i] > gp_data:
					sigma = radio_data[2][i]
					chi_s += (gp_data-radio_data[1][i])**2./sigma**2.
				if radio_data[1][i] < gp_data:
					sigma = radio_data[2][i]
					chi_s += (gp_data-radio_data[1][i])**2./sigma**2.
	chi += chi_s

	if chi < chi_tot:
		chi_tot = chi
		chi_veritas = chi_v
		chi_radio = chi_s
		chi_fermi = chi_f
		small_factor = FAC

ANALYTICS = np.loadtxt(galdef_name+'analytics.dat')
ANALYTICS *= small_factor

NUCLEI_DAT = np.loadtxt(galdef_name+'analytics_nuc.dat')
A = NUCLEI_DAT[0]
Z = NUCLEI_DAT[1]
NUC = NUCLEI_DAT[2]
NUC *= small_factor

if os.access(galdef_name+"_analytics.txt", os.F_OK ):  os.remove(galdef_name+"_analytics.txt")
outfile = open(galdef_name+"_analytics.txt", "w")
outfile.write('GALDEF = '+galdef_name+'\n')
outfile.write('\n')
outfile.write('INJECTION LUMINOSITIES:\n')
outfile.write('-Primary protons    = %E ergs s^-1\n' % ANALYTICS[1])
outfile.write('-Primary He         = %E ergs s^-1\n' % ANALYTICS[0])
outfile.write('-Primary electrons  = %E ergs s^-1\n' % ANALYTICS[2])
outfile.write('-Total Injected     = %E ergs s^-1\n' % ANALYTICS[3])
outfile.write('\n')
outfile.write('-Secondary protons    = %E ergs s^-1\n' % ANALYTICS[12])
outfile.write('-Secondary electrons  = %E ergs s^-1\n' % ANALYTICS[4])
outfile.write('-Secondary positrons  = %E ergs s^-1\n' % ANALYTICS[5])
outfile.write('\n')
outfile.write('RADIATION LUMINOSITIES:\n')
outfile.write('-Synchrotron  = %E ergs s^-1\n' % ANALYTICS[6])
outfile.write('-free-free    = %E ergs s^-1\n' % ANALYTICS[13])
outfile.write('\n')
outfile.write('-Total Gamma     = %E ergs s^-1\n' % ANALYTICS[7])
outfile.write('-Bremsstrahlung  = %E ergs s^-1\n' % ANALYTICS[8])
outfile.write('-Pion Decay      = %E ergs s^-1\n' % ANALYTICS[9])
outfile.write('-Inverse Compton = %E ergs s^-1\n' % ANALYTICS[10])
outfile.write('\n')
outfile.write('Emmission / Injection           = '+str((ANALYTICS[6]+ANALYTICS[7])/ANALYTICS[3])+'\n')
outfile.write('Emmission+Secondary / Injection = '+str((ANALYTICS[6]+ANALYTICS[7]+ANALYTICS[4]+ANALYTICS[5]+ANALYTICS[12])/ANALYTICS[3])+'\n')
outfile.write('Pion+Secondary / P+He           = '+str((ANALYTICS[12]+ANALYTICS[4]+ANALYTICS[5]+ANALYTICS[9])/(ANALYTICS[1]+ANALYTICS[0]))+'\n')
outfile.write('e+- Emission / All e+-          = '+str((ANALYTICS[6]+ANALYTICS[8]+ANALYTICS[10])/(ANALYTICS[2]+ANALYTICS[4]+ANALYTICS[5]))+'\n')
outfile.write('\n')
outfile.write('chi^2        = '+str(chi_tot)+'\n')
outfile.write('chi^2_radio  = '+str(chi_radio)+'\n')
outfile.write('chi^2_gamma  = '+str(chi_veritas+chi_fermi)+'\n')
outfile.write('chi^2 factor = '+str(small_factor)+'\n')
outfile.write('\n')
outfile.write('ENERGY IN CRS\n')
for j in range(0,len(NUC)):
	outfile.write('A = '+str(A[j]).zfill(2)+'  Z = '+str(Z[j]).zfill(2)+'  E = '+str(NUC[j])+'\n')
#outfile.write('Total energy = '+str(np.sum(nuc))+'\n')
outfile.close()

#NEW FERMI POINTS!!!
gamma_tot = np.loadtxt(galdef_name+'gamma_flux.dat')
gamma_tot_mask = np.isfinite(gamma_tot[1])
gamma_tot = gamma_tot.T[gamma_tot_mask].T
gamma_tot[1] *= small_factor
# MeV ergs/cm^2/s
gp_energy = gamma_tot[0]

#FERMI JUNK!!!

CC = np.zeros([num_fermi,5+1])
for i in range(0,num_fermi):

	gp_fermi = gamma_tot[1]/gp_energy
	gp_data = 0.
	emin_i = np.searchsorted(gp_energy-emin[i],0)
	emax_i = np.searchsorted(gp_energy-emax[i],0)
	#interpolate beginning and end
	e_temp = 10.**(math.log10(gp_fermi[emin_i])+math.log10(gp_fermi[emin_i+1]/gp_fermi[emin_i])*math.log10(emin[i]/gp_energy[emin_i])/math.log10(gp_energy[emin_i+1]/gp_energy[emin_i]))
	gp_data += math.sqrt(e_temp*gp_fermi[emin_i+1])*math.log10(gp_energy[emin_i+1]/emin[i])
	e_temp = 10.**(math.log10(gp_fermi[emax_i])+math.log10(gp_fermi[emax_i+1]/gp_fermi[emax_i])*math.log10(emax[i]/gp_energy[emax_i])/math.log10(gp_energy[emax_i+1]/gp_energy[emax_i]))
	gp_data += math.sqrt(e_temp*gp_fermi[emax_i])*math.log10(emax[i]/gp_energy[emax_i])
	#integrate
	for k in range(emin_i+1, emax_i):
		gp_data += math.sqrt(gp_fermi[k]*gp_fermi[k+1])*math.log10(gp_energy[k+1]/gp_energy[k])
	EdNdE = gp_data
	
	gp_fermi = gamma_tot[1]
	gp_data = 0.
	#interpolate beginning and end
	e_temp = 10.**(math.log10(gp_fermi[emin_i])+math.log10(gp_fermi[emin_i+1]/gp_fermi[emin_i])*math.log10(emin[i]/gp_energy[emin_i])/math.log10(gp_energy[emin_i+1]/gp_energy[emin_i]))
	gp_data += math.sqrt(e_temp*gp_fermi[emin_i+1])*math.log10(gp_energy[emin_i+1]/emin[i])
	e_temp = 10.**(math.log10(gp_fermi[emax_i])+math.log10(gp_fermi[emax_i+1]/gp_fermi[emax_i])*math.log10(emax[i]/gp_energy[emax_i])/math.log10(gp_energy[emax_i+1]/gp_energy[emax_i]))
	gp_data += math.sqrt(e_temp*gp_fermi[emax_i])*math.log10(emax[i]/gp_energy[emax_i])
	#integrate
	for k in range(emin_i+1, emax_i):
		gp_data += math.sqrt(gp_fermi[k]*gp_fermi[k+1])*math.log10(gp_energy[k+1]/gp_energy[k])
	E2dNdE = gp_data
	
	CC[i][0] = emin[i]
	CC[i][1] = emax[i]
	CC[i][2] = fermi_data[0][i]/EdNdE*E2dNdE
	CC[i][3] = fermi_data[3][i]/EdNdE*E2dNdE
	CC[i][4] = fermi_data[4][i]/EdNdE*E2dNdE
	CC[i][5] = 1/EdNdE*E2dNdE
	
	np.savetxt(galdef_name+'fermi.dat', CC, '%e')


if os.access(galdef_name+"_analytics.txt", os.F_OK ):  os.remove(galdef_name+"_analytics.txt")
outfile = open(galdef_name+"_analytics.txt", "w")
outfile.write('GALDEF = '+galdef_name+'\n')
outfile.write('\n')
outfile.write('INJECTION LUMINOSITIES:\n')
outfile.write('-Primary protons    = %E ergs s^-1\n' % ANALYTICS[1])
outfile.write('-Primary He         = %E ergs s^-1\n' % ANALYTICS[0])
outfile.write('-Primary electrons  = %E ergs s^-1\n' % ANALYTICS[2])
outfile.write('-Total Injected     = %E ergs s^-1\n' % ANALYTICS[3])
outfile.write('\n')
outfile.write('-Secondary electrons  = %E ergs s^-1\n' % ANALYTICS[4])
outfile.write('-Secondary positrons  = %E ergs s^-1\n' % ANALYTICS[5])
outfile.write('\n')
outfile.write('RADIATION LUMINOSITIES:\n')
outfile.write('-Synchrotron  = %E ergs s^-1\n' % ANALYTICS[6])
outfile.write('-Total Gamma  = %E ergs s^-1\n' % ANALYTICS[7])
outfile.write('\n')
outfile.write('-Bremsstrahlung  = %E ergs s^-1\n' % ANALYTICS[8])
outfile.write('-Pion Decay      = %E ergs s^-1\n' % ANALYTICS[9])
outfile.write('-Inverse Compton = %E ergs s^-1\n' % ANALYTICS[10])
outfile.write('\n')
outfile.write('Emmission / Injection           = '+str((ANALYTICS[6]+ANALYTICS[7])/ANALYTICS[3])+'\n')
outfile.write('Emmission+Secondary / Injection = '+str((ANALYTICS[6]+ANALYTICS[7]+ANALYTICS[4]+ANALYTICS[5])/ANALYTICS[3])+'\n')
outfile.write('Pion+Secondary / P+He           = '+str((ANALYTICS[4]+ANALYTICS[5]+ANALYTICS[9])/(ANALYTICS[1]+ANALYTICS[0]))+'\n')
outfile.write('e+- Emission / All e+-          = '+str((ANALYTICS[6]+ANALYTICS[8]+ANALYTICS[10])/(ANALYTICS[2]+ANALYTICS[4]+ANALYTICS[5]))+'\n')
outfile.write('\n')
outfile.write('chi^2        = '+str(chi_tot)+'\n')
outfile.write('chi^2_radio  = '+str(chi_radio)+'\n')
outfile.write('chi^2_gamma  = '+str(chi_veritas+chi_fermi)+'\n')
outfile.write('chi^2 factor = '+str(small_factor)+'\n')

	
###
###
### chi^2 for extended emission junk
###
###
if extended_emiss == 1:
	###
	###
	### RADIO FUNCTION OF Z
	###
	###
	###
	
	if SL_SCALE == 1:
		#load data
		zS = np.loadtxt(galdef_name+'freefreeL.fits_z0.dat')
		nuS = np.loadtxt(galdef_name+'freefreeL.fits_z0_nu.dat')
		specS = np.loadtxt(galdef_name+'freefreeL.fits_z0_spec.dat')
		specS*=small_factor#*NFC
		
		zL = np.loadtxt(galdef_name+'synchL.fits_z0.dat')
		nuL = np.loadtxt(galdef_name+'synchL.fits_z0_nu.dat')
		specL = np.loadtxt(galdef_name+'synchL.fits_z0_spec.dat')
		specL*=small_factor
		
		zS_o = np.loadtxt(galdef_name+'freefreeL.fits_zoff.dat')
		nuS_o = np.loadtxt(galdef_name+'freefreeL.fits_zoff_nu.dat')
		specS_o = np.loadtxt(galdef_name+'freefreeL.fits_zoff_spec.dat')
		specS_o*=small_factor#*NFC
		
		zL_o = np.loadtxt(galdef_name+'synchL.fits_zoff.dat')
		nuL_o = np.loadtxt(galdef_name+'synchL.fits_zoff_nu.dat')
		specL_o = np.loadtxt(galdef_name+'synchL.fits_zoff_spec.dat')
		specL_o*=small_factor
		
		# zL = np.loadtxt(galdef_name+'synchL_na.fits_z0.dat')
		# nuL = np.loadtxt(galdef_name+'synchL_na.fits_z0_nu.dat')
		# specL = np.loadtxt(galdef_name+'synchL_na.fits_z0_spec.dat')
		# specL*=small_factor
		
		# zL_o = np.loadtxt(galdef_name+'synchL_na.fits_zoff.dat')
		# nuL_o = np.loadtxt(galdef_name+'synchL_na.fits_zoff_nu.dat')
		# specL_o = np.loadtxt(galdef_name+'synchL_na.fits_zoff_spec.dat')
		# specL_o*=small_factor
	else:
		#load data
		zS = np.loadtxt(galdef_name+'freefree.fits_z0.dat')
		nuS = np.loadtxt(galdef_name+'freefree.fits_z0_nu.dat')
		specS = np.loadtxt(galdef_name+'freefree.fits_z0_spec.dat')
		specS*=small_factor#*NFC
		
		zL = np.loadtxt(galdef_name+'synch.fits_z0.dat')
		nuL = np.loadtxt(galdef_name+'synch.fits_z0_nu.dat')
		specL = np.loadtxt(galdef_name+'synch.fits_z0_spec.dat')
		specL*=small_factor
		
		zS_o = np.loadtxt(galdef_name+'freefree.fits_zoff.dat')
		nuS_o = np.loadtxt(galdef_name+'freefree.fits_zoff_nu.dat')
		specS_o = np.loadtxt(galdef_name+'freefree.fits_zoff_spec.dat')
		specS_o*=small_factor#*NFC
		
		zL_o = np.loadtxt(galdef_name+'synch.fits_zoff.dat')
		nuL_o = np.loadtxt(galdef_name+'synch.fits_zoff_nu.dat')
		specL_o = np.loadtxt(galdef_name+'synch.fits_zoff_spec.dat')
		specL_o*=small_factor
		
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
			specS_92cm = interpolate.interp1d(specS_temp[0],specS_temp[1])
			specL_92cm = interpolate.interp1d(specL_temp[0],specL_temp[1])
			specS_o_92cm = interpolate.interp1d(specS_o_temp[0],specS_o_temp[1])
			specL_o_92cm = interpolate.interp1d(specL_o_temp[0],specL_o_temp[1])#specL_o_temp
		elif j == 1:
			specS_22cm = interpolate.interp1d(specS_temp[0],specS_temp[1])#specS_temp
			specL_22cm = interpolate.interp1d(specL_temp[0],specL_temp[1])
			specS_o_22cm = interpolate.interp1d(specS_o_temp[0],specS_o_temp[1])
			specL_o_22cm = interpolate.interp1d(specL_o_temp[0],specL_o_temp[1])
		elif j == 2:
			specS_6cm = interpolate.interp1d(specS_temp[0],specS_temp[1])
			specL_6cm = interpolate.interp1d(specL_temp[0],specL_temp[1])
		elif j == 3:
			specS_3cm = interpolate.interp1d(specS_temp[0],specS_temp[1])
			specL_3cm = interpolate.interp1d(specL_temp[0],specL_temp[1])
		
	# chi^2 calculation
	# READ IN DATA
	## 92cm
	L92cm_N = np.loadtxt(data_directory+'92cm/92cm_N.dat').T
	L92cm_S = np.loadtxt(data_directory+'92cm/92cm_S.dat').T
	
	L92 = np.append(L92cm_N,L92cm_S,axis=1)
	
	
	L92sd = (abs(L92[2]-L92[1]) + abs(L92[3]-L92[1]))/2. 
	
	gp_L92 = specS_92cm(L92[0]) + specL_92cm(L92[0])
	# print gp_L92.shape, L92.shape, L92sd.shape
	chi92 = np.sum((gp_L92-L92[1])**2/L92sd**2)
	
	
	## 92cm OFFSET
	L92cm_NW = np.loadtxt(data_directory+'92cm/92cm_NW.dat').T
	L92cm_SW = np.loadtxt(data_directory+'92cm/92cm_SW.dat').T
	L92cm_NE = np.loadtxt(data_directory+'92cm/92cm_NE.dat').T
	L92cm_SE = np.loadtxt(data_directory+'92cm/92cm_SE.dat').T
	
	L92_o = np.append(L92cm_NW,L92cm_SW,axis=1)
	L92_o = np.append(L92_o,L92cm_NE,axis=1)
	L92_o = np.append(L92_o,L92cm_SE,axis=1)
	
	L92sd_o = (abs(L92_o[2]-L92_o[1]) + abs(L92_o[3]-L92_o[1]))/2.
	
	gp_L92_o = specS_o_92cm(L92_o[0]) + specL_o_92cm(L92_o[0])
	
	chi92_o = np.sum((gp_L92_o-L92_o[1])**2/L92sd_o**2)
	
	## 22cm
	L22cm_N = np.loadtxt(data_directory+'22cm/22cm_N.dat').T
	L22cm_S = np.loadtxt(data_directory+'22cm/22cm_S.dat').T
	
	L22 = np.append(L22cm_N,L22cm_S,axis=1)
	
	L22sd = (abs(L22[2]-L22[1]) + abs(L22[3]-L22[1]))/2. 
	
	gp_L22 = specS_22cm(L22[0]) + specL_22cm(L22[0])
	
	chi22 = np.sum((gp_L22-L22[1])**2/L22sd**2.)
	
	## 22cm OFFSET
	L22cm_NW = np.loadtxt(data_directory+'22cm/22cm_NW.dat').T
	L22cm_SW = np.loadtxt(data_directory+'22cm/22cm_SW.dat').T
	L22cm_NE = np.loadtxt(data_directory+'22cm/22cm_NE.dat').T
	L22cm_SE = np.loadtxt(data_directory+'22cm/22cm_SE.dat').T
	
	L22_o = np.append(L22cm_NW,L22cm_SW,axis=1)
	L22_o = np.append(L22_o,L22cm_NE,axis=1)
	L22_o = np.append(L22_o,L22cm_SE,axis=1)
	
	L22sd_o = (abs(L22_o[2]-L22_o[1]) + abs(L22_o[3]-L22_o[1]))/2.
	
	gp_L22_o = specS_o_22cm(L22_o[0]) + specL_o_22cm(L22_o[0])
	
	chi22_o = np.sum((gp_L22_o-L22_o[1])**2/L22sd_o**2.)
	
	## 6cm
	L6cm_N = np.loadtxt(data_directory+'6cm/6cm_N.dat').T
	L6cm_S = np.loadtxt(data_directory+'6cm/6cm_S.dat').T
	
	L6 = np.append(L6cm_N,L6cm_S,axis=1)
	
	L6sd = (abs(L6[2]-L6[1]) + abs(L6[3]-L6[1]))/2. 
	
	gp_L6 = specS_6cm(L6[0]) + specL_6cm(L6[0])
	
	chi6 = np.sum((gp_L6-L6[1])**2/L6sd**2.)
	
	## 3cm
	L3cm_N = np.loadtxt(data_directory+'3cm/3cm_N.dat').T
	L3cm_S = np.loadtxt(data_directory+'3cm/3cm_S.dat').T
	
	L3 = np.append(L3cm_N,L3cm_S,axis=1)
	
	L3sd = (abs(L3[2]-L3[1]) + abs(L3[3]-L3[1]))/2. 
	
	gp_L3 = specS_3cm(L3[0]) + specL_3cm(L3[0])
	
	chi3 = np.sum((gp_L3-L3[1])**2/L3sd**2.)
	
	chi_ext = chi92 + chi92_o + chi22 + chi22_o + chi6 + chi3
	chi_all = chi_ext + chi_tot
	save_array = [chi_tot, chi_fermi, chi_veritas, chi_radio, small_factor, chi92, chi92_o, chi22, chi22_o, chi6, chi3, chi_ext, chi_all]
	np.savetxt(galdef_name+'chi.dat', save_array, '%e')
	
	outfile.write('chi92^2      = '+str(chi92)+'\n')
	outfile.write('chi92_o^2    = '+str(chi92_o)+'\n')
	outfile.write('chi22^2      = '+str(chi22)+'\n')
	outfile.write('chi22_o^2    = '+str(chi22_o)+'\n')
	outfile.write('chi6^2       = '+str(chi6)+'\n')
	outfile.write('chi3^2       = '+str(chi3)+'\n')
	outfile.write('chi_ext^2    = '+str(chi_ext)+'\n')
	outfile.write('chi_all^2    = '+str(chi_all)+'\n')
else:
	np.savetxt(galdef_name+'chi.dat', [chi_tot, chi_fermi, chi_veritas, chi_radio, small_factor], '%e')


outfile.write('\n')
outfile.write('ENERGY IN CRS\n')
for j in range(0,len(NUC)):
	outfile.write('A = '+str(A[j]).zfill(2)+'  Z = '+str(Z[j]).zfill(2)+'  E = '+str(NUC[j])+'\n')
#outfile.write('Total energy = '+str(np.sum(nuc))+'\n')
outfile.close()



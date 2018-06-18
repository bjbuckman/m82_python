import gp2d_LOSimage_v3 as gplos
import gp2d_emission_v2 as gpemi
import numpy as np
#import matplotlib.pyplot as plt
import sys

#GALDEF NAME
galdef_name = str(sys.argv[1])

# PARAMETERS
inclination = 80
distance = 3520  	# kpc

##SINGLE SCALE
outRes = 201 		# ^2 pixels
int_steps = 1001 	# LOS integration steps
imageR = 1.0 		# kpc Window of image

##FOR MULTISCALE
SL_SCALE = 0 #SET TO 1 FOR TWO SCALE

outResL = 101 		# ^2 pixels
int_stepsL = 701 	# LOS integration steps
imageRL = 4.0 		# kpc Window of image

outResS = 101 		# ^2 pixels
int_stepsS = 501 	# LOS integration steps
imageRS = 0.80 		# kpc Window of image

##PRINT PARAMETERS
PRINTGAMMA = 0
NO_ABSORPTION = 0
ABSORPTION = 1

#CALCULATING LINE OF SIGHT
if PRINTGAMMA == 1:
	if SL_SCALE ==1:
		pass
		# gplos.LOS_integral('bremss_emiss_54_'+galdef_name+'.gz',galdef_name+'bremss.fits',inclination, distance, outRes, int_steps, imageR)
		# gplos.LOS_integral('pion_decay_emiss_54_'+galdef_name+'.gz',galdef_name+'pion.fits',inclination, distance, outRes, int_steps, imageR)
		# gplos.LOS_integral('IC_emiss_54_'+galdef_name+'.gz',galdef_name+'IC.fits',inclination, distance, outRes, int_steps, imageR)
		# gplos.gamma_total_los(galdef_name+'bremss.fits',galdef_name+'pion.fits',galdef_name+'IC.fits',galdef_name+'gamma_total.fits')
		
		# gplos.gamma_spectrum(galdef_name+'gamma_total.fits',galdef_name+'gamma_total_spec.fits')
		# gplos.gamma_spectrum(galdef_name+'bremss.fits',galdef_name+'bremss_spec.fits')
		# gplos.gamma_spectrum(galdef_name+'pion.fits',galdef_name+'pion_spec.fits')
		# gplos.gamma_spectrum(galdef_name+'IC.fits', galdef_name+'IC_spec.fits')
	else:
		gplos.LOS_integral('bremss_emiss_54_'+galdef_name+'.gz',galdef_name+'bremss.fits',inclination, distance, outRes, int_steps, imageR)
		gplos.LOS_integral('pion_decay_emiss_54_'+galdef_name+'.gz',galdef_name+'pion.fits',inclination, distance, outRes, int_steps, imageR)
		gplos.LOS_integral('IC_emiss_54_'+galdef_name+'.gz',galdef_name+'IC.fits',inclination, distance, outRes, int_steps, imageR)
		gplos.gamma_total_los(galdef_name+'bremss.fits',galdef_name+'pion.fits',galdef_name+'IC.fits',galdef_name+'gamma_total.fits')
		
		gplos.gamma_spectrum(galdef_name+'gamma_total.fits',galdef_name+'gamma_total_spec.fits')
		gplos.gamma_spectrum(galdef_name+'bremss.fits',galdef_name+'bremss_spec.fits')
		gplos.gamma_spectrum(galdef_name+'pion.fits',galdef_name+'pion_spec.fits')
		gplos.gamma_spectrum(galdef_name+'IC.fits', galdef_name+'IC_spec.fits')
	
	
if NO_ABSORPTION == 1:
	if SL_SCALE == 1:
		gplos.LOS_integral('synchrotron_emiss_54_'+galdef_name+'.gz', galdef_name+'synchS_na.fits', inclination, distance, outResS, int_stepsS, imageRS)
		gplos.LOS_integral('free_free_emiss_54_'+galdef_name+'.gz', galdef_name+'freefreeS_na.fits', inclination, distance, outResS, int_stepsS, imageRS)
		gplos.radio_total_los(galdef_name+'synchS_na.fits',galdef_name+'freefreeS_na.fits',galdef_name+'radioS_na.fits')

		gplos.LOS_integral('synchrotron_emiss_54_'+galdef_name+'.gz', galdef_name+'synchL_na.fits', inclination, distance, outResL, int_stepsL, imageRL)
		gplos.LOS_integral('free_free_emiss_54_'+galdef_name+'.gz', galdef_name+'freefreeL_na.fits', inclination, distance, outResL, int_stepsL, imageRL)
		gplos.radio_total_los(galdef_name+'synchL_na.fits',galdef_name+'freefreeL_na.fits',galdef_name+'radioL_na.fits')

		#CALCULATING SPECTRUM IN EACH PIXEL
		gplos.synch_spectrum(galdef_name+'synchS_na.fits',galdef_name+'synchS_na_spec.fits')
		gplos.synch_spectrum(galdef_name+'freefreeS_na.fits',galdef_name+'freefreeS_na_spec.fits')
		gplos.synch_spectrum(galdef_name+'radioS_na.fits',galdef_name+'radioS_na_spec.fits')

		gplos.synch_spectrum(galdef_name+'synchL_na.fits',galdef_name+'synchL_na_spec.fits')
		gplos.synch_spectrum(galdef_name+'freefreeL_na.fits',galdef_name+'freefreeL_na_spec.fits')
		gplos.synch_spectrum(galdef_name+'radioL_na.fits',galdef_name+'radioL_na_spec.fits')
		
	else:
		gplos.LOS_integral('synchrotron_emiss_54_'+galdef_name+'.gz', galdef_name+'synch.fits', inclination, distance, outRes, int_steps, imageR)
		gplos.LOS_integral('free_free_emiss_54_'+galdef_name+'.gz', galdef_name+'freefree.fits', inclination, distance, outRes, int_steps, imageR)
		gplos.radio_total_los(galdef_name+'synch.fits',galdef_name+'freefree.fits',galdef_name+'radio.fits')

		#CALCULATING SPECTRUM IN EACH PIXEL
		gplos.synch_spectrum(galdef_name+'synch.fits',galdef_name+'synch_spec.fits')
		gplos.synch_spectrum(galdef_name+'freefree.fits',galdef_name+'freefree_spec.fits')
		gplos.synch_spectrum(galdef_name+'radio.fits',galdef_name+'radio_spec.fits')
		
	
	
if ABSORPTION == 1:	
	if SL_SCALE == 1:
		gplos.LOS_integral_radio2('synchrotron_emiss_54_'+galdef_name+'.gz', 'free_free_emiss_54_'+galdef_name+'.gz', galdef_name+'synchS.fits', galdef_name+'freefreeS.fits', inclination, distance, outResS, int_stepsS, imageRS, 'free_free_absorb_54_'+galdef_name+'.gz')
		gplos.LOS_integral_radio2('synchrotron_emiss_54_'+galdef_name+'.gz', 'free_free_emiss_54_'+galdef_name+'.gz', galdef_name+'synchL.fits', galdef_name+'freefreeL.fits', inclination, distance, outResL, int_stepsL, imageRL, 'free_free_absorb_54_'+galdef_name+'.gz')

		gplos.radio_total_los(galdef_name+'synchS.fits',galdef_name+'freefreeS.fits',galdef_name+'radioS.fits')
		gplos.radio_total_los(galdef_name+'synchL.fits',galdef_name+'freefreeL.fits',galdef_name+'radioL.fits')

		#CALCULATING SPECTRUM IN EACH PIXEL
		gplos.synch_spectrum(galdef_name+'synchS.fits',galdef_name+'synchS_spec.fits')
		gplos.synch_spectrum(galdef_name+'freefreeS.fits',galdef_name+'freefreeS_spec.fits')
		gplos.synch_spectrum(galdef_name+'radioS.fits',galdef_name+'radioS_spec.fits')

		gplos.synch_spectrum(galdef_name+'synchL.fits',galdef_name+'synchL_spec.fits')
		gplos.synch_spectrum(galdef_name+'freefreeL.fits',galdef_name+'freefreeL_spec.fits')
		gplos.synch_spectrum(galdef_name+'radioL.fits',galdef_name+'radioL_spec.fits')
		
	else:
		gplos.LOS_integral_radio2('synchrotron_emiss_54_'+galdef_name+'.gz', 'free_free_emiss_54_'+galdef_name+'.gz', galdef_name+'synch.fits', galdef_name+'freefree.fits', inclination, distance, outRes, int_steps, imageR, 'free_free_absorb_54_'+galdef_name+'.gz')
		gplos.radio_total_los(galdef_name+'synch.fits',galdef_name+'freefree.fits',galdef_name+'radio.fits')

		#CALCULATING SPECTRUM IN EACH PIXEL
		gplos.synch_spectrum(galdef_name+'synch.fits',galdef_name+'synch_spec.fits')
		gplos.synch_spectrum(galdef_name+'freefree.fits',galdef_name+'freefree_spec.fits')
		gplos.synch_spectrum(galdef_name+'radio.fits',galdef_name+'radio_spec.fits')


#CALCULATING TOTAL FLUX AT EARTH
#gplos.synch_flux(galdef_name+'synch.fits',galdef_name+'synch_flux')
#gplos.gamma_flux(galdef_name+'gamma_total.fits',galdef_name+'gamma_total_flux')
#gplos.gamma_flux(galdef_name+'bremss.fits',galdef_name+'bremss_flux')
#gplos.gamma_flux(galdef_name+'pion.fits',galdef_name+'pion_flux')
#gplos.gamma_flux(galdef_name+'IC.fits', galdef_name+'IC_flux')

##PLOT ALL GAMMA TOGETHER
#gamma_tot = np.loadtxt(galdef_name+'gamma_total_flux.dat')
#bremss_tot = np.loadtxt(galdef_name+'bremss_flux.dat')
#pion_tot = np.loadtxt(galdef_name+'pion_flux.dat')
#ic_tot = np.loadtxt(galdef_name+'IC_flux.dat')

##PLOTTING DATA
#plt.loglog(gamma_tot[0], gamma_tot[1], label='gamma_tot')
#plt.loglog(bremss_tot[0], bremss_tot[1], label='bremss')
#plt.loglog(pion_tot[0], pion_tot[1], label='pion')
#plt.loglog(ic_tot[0], ic_tot[1], label='IC')
#plt.title('Gamma Flux')
#plt.xlabel(r'E (MeV)')
#plt.ylabel(r'MeV$^2$ s$^{-1}$ cm$^{-2}$ MeV$^{-1}$')
#plt.legend(loc=4)
#plt.savefig(galdef_name+'all_gamma_spec')
#plt.close()

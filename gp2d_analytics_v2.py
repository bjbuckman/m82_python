import sys
import gp2d_particle_luminosity_v1 as particle
import gp2d_emission_v2 as gpem
import gp2d_LOSimage_v3 as gplos
import os
import numpy as np

galdef_name = sys.argv[1]
distance = 3520 #kpc

SL_SCALE = 0

##PRINT PARAMETERS
PRINTGAMMA = 0
NO_ABSORPTION = 1 #DON'T NEED FLUX FITS! ALWAYS 1
ABSORPTION = 1    #NEEDS CLUX FITS

if os.path.isfile('primary_Helium_source_function_54_'+galdef_name+'.gz'):
	pHe = particle.particle_luminosity('primary_Helium_source_function_54_'+galdef_name+'.gz','Helium')
else:
	pHe = 0.
pH = particle.particle_luminosity('primary_protons_source_function_54_'+galdef_name+'.gz','proton')
if os.path.isfile('secondary_protons_source_function_54_'+galdef_name+'.gz'):
	psH = particle.particle_luminosity('secondary_protons_source_function_54_'+galdef_name+'.gz','proton')
else:
	psH = 0.
pem = particle.particle_luminosity('primary_electrons_source_function_54_'+galdef_name+'.gz','electron')
sem = particle.particle_luminosity('secondary_electrons_source_function_54_'+galdef_name+'.gz','electron')
sep = particle.particle_luminosity('secondary_positrons_source_function_54_'+galdef_name+'.gz','positron')
total_injection = pHe + pH + pem

##Setting all variables
synch = 0
synchS = 0
synchL = 0
synch_na = 0

freef = 0
freefS = 0
freefL = 0
freef_na = 0

radio = 0
radioS = 0
radioL = 0
radio_na = 0

brems = 0
pions = 0
icomp = 0
gamma = 0
ratio = 0

if SL_SCALE == 1:
	if ABSORPTION == 1:
		synchS = gplos.synch_flux(galdef_name+'synchS.fits',galdef_name+'synchS_flux')
		freefS = gplos.synch_flux(galdef_name+'freefreeS.fits',galdef_name+'freefreeS_flux')
		radioS = gplos.synch_flux(galdef_name+'radioS.fits',galdef_name+'radioS_flux')
		
		synchL = gplos.synch_flux(galdef_name+'synchL.fits',galdef_name+'synchL_flux')
		freefL = gplos.synch_flux(galdef_name+'freefreeL.fits',galdef_name+'freefreeL_flux')
		radioL = gplos.synch_flux(galdef_name+'radioL.fits',galdef_name+'radioL_flux')
else:
	if ABSORPTION == 1:
		synch = gplos.synch_flux(galdef_name+'synch.fits',galdef_name+'synch_flux')
		freef = gplos.synch_flux(galdef_name+'freefree.fits',galdef_name+'freefree_flux')
		radio = gplos.synch_flux(galdef_name+'radio.fits',galdef_name+'radio_flux')

if NO_ABSORPTION == 1:
	synch_na = gpem.integrated_flux('synchrotron_emiss_54_'+galdef_name+'.gz',galdef_name+'synch_flux_na',distance)
	freef_na = gpem.integrated_flux('free_free_emiss_54_'+galdef_name+'.gz',galdef_name+'freefree_flux_na',distance)
	gpem.sum_flux2(galdef_name+'synch_flux.dat',galdef_name+'freefree_flux.dat',galdef_name+'radio_flux_na')
	radio_na = synch_na + freef_na		
		
brems = gpem.integrated_flux('bremss_emiss_54_'+galdef_name+'.gz',galdef_name+'bremss_flux',distance)
pions = gpem.integrated_flux('pion_decay_emiss_54_'+galdef_name+'.gz',galdef_name+'pion_flux',distance)
icomp = gpem.integrated_flux('IC_emiss_54_'+galdef_name+'.gz',galdef_name+'IC_flux',distance)

gpem.sum_flux3(galdef_name+'bremss_flux.dat',galdef_name+'pion_flux.dat',galdef_name+'IC_flux.dat',galdef_name+'gamma_flux')

gamma = brems+pions+icomp

ratio = (gamma+radio_na)/total_injection

AA, ZZ, nuclei = particle.nuclei_luminosity(galdef_name)

print pHe, pH, pem, total_injection, sem, sep, synch_na, gamma, brems, pions, icomp, ratio, psH, freef_na
#if os.access(galdef_name+"_analytics.dat", os.F_OK ):  os.remove(galdef_name+"_analytics.dat")
DATA = np.array([pHe, pH, pem, total_injection, sem, sep, synch_na, gamma, brems, pions, icomp, ratio, psH, freef_na])
DATA2 = np.array([AA,ZZ,nuclei])
np.savetxt(galdef_name+'analytics.dat',DATA,'%e')
np.savetxt(galdef_name+'analytics_nuc.dat',DATA2,'%e')


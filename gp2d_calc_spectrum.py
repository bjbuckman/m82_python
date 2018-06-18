import numpy as np
import os
import sys
import gp2d_LOSanalyze_v1 as los

galdef_name = sys.argv[1]

SL_SCALE = 0

if SL_SCALE == 1:
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'radioS.fits', '0', 0)
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'radioS.fits', 'off',60)
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'radioL_na.fits', '0', 0)
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'radioL_na.fits', 'off', 60)
	
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'synchS.fits', '0', 0)
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'synchS.fits', 'off',60)
	z, nu, data = los.synch_spectrum_z(galdef_name+'synchL.fits', '0', 0)
	z, nu, data = los.synch_spectrum_z(galdef_name+'synchL.fits', 'off', 60)
	
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'freefS.fits', '0', 0)
	# z, nu, data = los.synch_spectrum_z('./'+new_galdefname+'/'+new_galdefname+'freefS.fits', 'off',60)
	z, nu, data = los.synch_spectrum_z(galdef_name+'freefreeL.fits', '0', 0)
	z, nu, data = los.synch_spectrum_z(galdef_name+'freefreeL.fits', 'off', 60)
else:
	z, nu, data = los.synch_spectrum_z(galdef_name+'synch.fits', '0', 0)
	z, nu, data = los.synch_spectrum_z(galdef_name+'synch.fits', 'off', 60)

	z, nu, data = los.synch_spectrum_z(galdef_name+'freefree.fits', '0', 0)
	z, nu, data = los.synch_spectrum_z(galdef_name+'freefree.fits', 'off', 60)
	
import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import imp
magnetic_field = imp.load_source("magnetic_field_v4", "../m82_dist_py/magnetic_field_v4.py")
import pyfits as py

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker


c = 3.0e10
mp = 1.66e-24 #g
ev2ergs = 1.60218e-12

MODE = 2

if MODE == 0:
	###
	### MAGNETIC FIELD 1_3
	###

	magnetic_model_name = 'constant+powerlaw_v4'

	#B_field_parameters   = 0.000150,0.2,0.2,0.05,0.05,0.05,0.2,-1.4,8,9  

	B0 = 150.e-6
	r_scale = 0.2
	r_scale1 = 0.2
	z_scale = 0.05
	z_scale1 = 0.05

	r_pow = 0.05
	r_pow1 = 0.2
	b_pow = -1.

	B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

	# B = np.zeros(znum)
	# for iz in range(0,znum):
		# B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
		
	# U_magnetic = B**2./(8.*math.pi)


	# magnetic_accel_1_3 = abs(U_magnetic)/ev2ergs

	###
	### [z, B] [z, grad_magnetic]
	##############################################
	
elif MODE == 1:
	###
	### MAGNETIC FIELD 2_3
	###

	magnetic_model_name = 'constant+powerlaw_v4'
	#B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

	B0 = 325.e-6
	r_scale = 0.2
	r_scale1 = 0.2
	z_scale = 0.05
	z_scale1 = 0.05

	r_pow = 0.05
	r_pow1 = 0.2
	b_pow = -1.2

	B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

	# B = np.zeros(znum)
	# for iz in range(0,znum):
		# B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
		
	# U_magnetic = B**2./(8.*math.pi)

	# magnetic_accel_2_3 = (abs(U_magnetic))/ev2ergs#/(mp*(nHI2+nHII2))

	###
	### [z, B] [z, grad_magnetic]
	##############################################
	
elif MODE == 2:
	###
	### MAGNETIC FIELD 3_3
	###

	magnetic_model_name = 'constant+powerlaw_v4'
	#B_field_parameters   = 0.000325,0.2,0.2,0.05,0.05,0.05,0.2,-0.4,8,9 

	B0 = 325.e-6
	r_scale = 0.2
	r_scale1 = 0.2
	z_scale = 0.05
	z_scale1 = 0.05

	r_pow = 0.05
	r_pow1 = 0.2
	b_pow = -0.2

	B_parameters = [B0, r_scale, r_scale1, z_scale, z_scale1, r_pow, r_pow1, b_pow]

	# B = np.zeros(znum)
	# for iz in range(0,znum):
		# B[iz] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([0.,z[iz]]), B_parameters)
		
	# U_magnetic = B**2./(8.*math.pi)

	# magnetic_accel_3_3 = (abs(U_magnetic))/ev2ergs#/(mp*(nHI2+nHII2))

	###
	### [z, B] [z, grad_magnetic]
	##############################################
	

###
### U_cr
###
ZZ = 0.

if MODE == 0:
	GG = 'm82_05aff1018'
	factor = 2.18287458838
elif MODE == 1:
	GG =  'm82_05b2ff1019'
	factor = 1.55132821598
elif MODE == 2:
	GG =  'm82_06bff1022'
	factor = 1.47745544379
	

dir = '../fits/'
galdef = GG
nuclei_file = dir+'raw/'+galdef+'/nuclei_full_54_'+galdef+'.gz'

chi_file = dir+galdef+'/'+galdef+'chi.dat'
chi_dat = np.loadtxt(chi_file)
chi_fac = chi_dat[4]

#constants
pc2cm = 3.08568025e18 #cm/pc
kpc2cm = 3.08568025e21 #cm/kpc
radian2arcsec = 648000./math.pi #acsec/radian
Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
c = 2.99792458e10 #cm/s
ev2ergs = 1.60218e-12

#READING EMISSION FITS-------------------------------------------------
hdulist=py.open(nuclei_file)
head=hdulist[0].header
scidata=hdulist[0].data

#Distances in kpc
xmin = head['CRVAL1']
zmin = head['CRVAL2']
log_e_min = head['CRVAL3'] #log10(e_min)

xdel = head['CDELT1']
zdel = head['CDELT2']
log_efactor = head['CDELT3'] #factor

xnum = head['NAXIS1']
znum = head['NAXIS2']
enum = head['NAXIS3']
nucnum = head['NAXIS4']

emin = 10.**log_e_min  #MeV
efactor = 10.**log_efactor

xmax = xmin+xdel*(xnum-1)
zmax = zmin+zdel*(znum-1)

xsize = xmax-xmin
zsize = zmax-zmin

R = xmin+xdel*np.array(range(0,xnum))
Z = zmin+zdel*np.array(range(0,znum))
energy = emin*efactor**np.array(range(0,enum))

r_int = math.pi*2.*R*xdel
# r_int[0] = math.pi*(xdel/2.)**2.

### Magnetic Field
B = np.zeros([znum,xnum])
for iz in range(0,znum):
	for ir in range(0,xnum):
		B[iz][ir] = magnetic_field.magnetic_field_2d(magnetic_model_name, np.array([R[ir],Z[iz]]), B_parameters)

U_magnetic = B**2/(8*math.pi)/ev2ergs

print U_magnetic.shape
print scidata.shape

scidata_n = np.zeros([nucnum,xnum,znum])
for i in range(0,nucnum):
	
	AA = np.zeros(nucnum)
	AAt = np.zeros(nucnum)
	ZZ = np.zeros(nucnum)
	mass = np.ones(nucnum)
	
	AA[i] = head['NUCA'+str(i+1).zfill(3)]
	ZZ[i] = head['NUCZ'+str(i+1).zfill(3)]
	
	AAt[i] = max(AA[i],1) 
	if AA[i] == 0:
		mass[i] = 0.5109989461 #MeV/c^2
	else:
		mass[i] = AA[i]*938.28

	dat = scidata[i]*U_magnetic

	beta = np.sqrt(1.-1./(1.+energy/mass[i])**2.)

	#picking out radius
	dat = np.swapaxes(dat*r_int, 0, 2)*4.*math.pi/beta/c*AAt[i]
	dat = np.sum(dat, axis=2)*math.log(energy[1]/energy[0])*1e6
	scidata_n[i] = dat

print scidata_n.shape

def CR_U_z(emissfilename,ZZZZ, r_core, z_core):
	''' Calculates stuff '''
	U_cr = np.zeros(nucnum)
	for i in range(0,nucnum):
		dat = scidata_n[i]
		dat = dat[R <= r_core]
		dat = np.sum(dat, axis=0)
		dat = dat[abs(Z-ZZZZ) <= z_core]
		U_cr[i] = np.sum(dat, axis=0)*zdel
		# dat/= math.pi*r_core**2.*2.*z_core
	
	return ZZ, AA, U_cr




rr = np.linspace(0.,3.5,num=50)
U_cr = np.zeros([len(rr),10])

for ir in range(0,len(rr)):
	ZZZ, A, U_cr[ir] = CR_U_z(nuclei_file, 0., rr[ir], rr[ir])

# Z, A, U_cr = CR_U_z(nuclei_file, ZZ)
# U_cr*= chi_fac #/ev2ergs


###
### PLOT
###

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

cmap = plt.get_cmap('cubehelix')
fig = plt.figure()
	
f1, ax = plt.subplots(1,1)#, figsize=(18,6))

for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(2)
ax.tick_params(which='major',width=1, length=8, labelsize='small')
ax.tick_params(which='minor',width=0.5, length=5)

MSIZE = 5
LWIDTH = 3
ALPHA = 0.65

U_cr = U_cr.T
# U_cr = abs
# print U_cr

norm = np.amax([U_cr[0], U_cr[2], U_cr[3]])
# print norm
 
ax.plot(rr, U_cr[3]/norm, c=cmap(0.25), label='pri_e-', linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
ax.plot(rr, U_cr[2]/norm, c=cmap(0.50), label='sec_e-', linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
ax.plot(rr, U_cr[0]/norm, c=cmap(0.75), label='sec_e+', linestyle='-', linewidth=LWIDTH, alpha=ALPHA)
ax.plot(rr, (U_cr[0]+U_cr[2])/norm, c='black', label='sec_e+-', linestyle='-', linewidth=LWIDTH, alpha=ALPHA)

blahg = (U_cr[0]+U_cr[2])/norm
print blahg[-1] 

# ax.set_yscale('log')
ax.set_xlim([0.0,3.5])
ax.set_ylim([0.,1.5])
ax.set_xlabel(r'$z$ [kpc]')
ax.set_ylabel(r'$U_B$-weighted $\frac{U_{cr}(<=z)}{\mathrm{max_z}U_{cr}(z)}$ []')

ax.grid(b=True, which='major', color='black', linestyle='-')

ax.legend(loc='upper right', prop={'size':14})

plt.tight_layout()
plot_name = 'plot_avg_synch_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break
plt.close()

# print (U_cr[2]+U_cr[0])/U_cr[3]

# file_name = 'U_CR_'+galdef+'_Bcor_5.dat'
# save_arr = np.array([Z,A,U_cr])
# np.savetxt(file_name, save_arr.T)


















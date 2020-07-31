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

def DoRotation(xspan, yspan, RotRad=0):
    """Generate a meshgrid and rotate it by RotRad radians."""

    # Clockwise, 2D rotation matrix
    RotMatrix = np.array([[np.cos(RotRad),  np.sin(RotRad)],
                          [-np.sin(RotRad), np.cos(RotRad)]])

    x, y = np.meshgrid(xspan, yspan)
    return np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([x, y]))

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

# #22cm Data
# R22_01 = np.loadtxt('./22cm/22cm_01.txt').T
# R22_02 = np.loadtxt('./22cm/22cm_02.txt').T
# R22_03 = np.loadtxt('./22cm/22cm_03.txt').T
# R22_04 = np.loadtxt('./22cm/22cm_04.txt').T
# R22_05 = np.loadtxt('./22cm/22cm_05.txt').T
# R22_06 = np.loadtxt('./22cm/22cm_06.txt').T
# R22_07 = np.loadtxt('./22cm/22cm_07.txt').T
# R22_08 = np.loadtxt('./22cm/22cm_08.txt').T
# R22_09 = np.loadtxt('./22cm/22cm_09.txt').T
# R22_10 = np.loadtxt('./22cm/22cm_10.txt').T
# R22_11 = np.loadtxt('./22cm/22cm_11.txt').T
# R22_12 = np.loadtxt('./22cm/22cm_12.txt').T
# R22_13 = np.loadtxt('./22cm/22cm_13.txt').T
# R22_14 = np.loadtxt('./22cm/22cm_14.txt').T
# R22_15 = np.loadtxt('./22cm/22cm_15.txt').T

# R22_flux = 2**np.linspace(0,14,15)*90.e-6
# R22_flux = R22_flux[::-1]
# print R22_flux

# #Angle translation
# dec_trans = 47 #"
# R22_01[1]-= dec_trans
# R22_02[1]-= dec_trans
# R22_03[1]-= dec_trans
# R22_04[1]-= dec_trans
# R22_05[1]-= dec_trans
# R22_06[1]-= dec_trans
# R22_07[1]-= dec_trans
# R22_08[1]-= dec_trans
# R22_09[1]-= dec_trans
# R22_10[1]-= dec_trans
# R22_11[1]-= dec_trans
# R22_12[1]-= dec_trans
# R22_13[1]-= dec_trans
# R22_14[1]-= dec_trans
# R22_15[1]-= dec_trans

# ra_trans = -8.2
# dec = 69.68*math.pi/180
# ra_factor = 360/24*math.cos(dec)#(1-math.cos(dec))/(math.cos(dec))**(-1)
# R22_01[0] = ra_factor*(R22_01[0]-ra_trans)
# R22_02[0] = ra_factor*(R22_02[0]-ra_trans)
# R22_03[0] = ra_factor*(R22_03[0]-ra_trans)
# R22_04[0] = ra_factor*(R22_04[0]-ra_trans)
# R22_05[0] = ra_factor*(R22_05[0]-ra_trans)
# R22_06[0] = ra_factor*(R22_06[0]-ra_trans)
# R22_07[0] = ra_factor*(R22_07[0]-ra_trans)
# R22_08[0] = ra_factor*(R22_08[0]-ra_trans)
# R22_09[0] = ra_factor*(R22_09[0]-ra_trans)
# R22_10[0] = ra_factor*(R22_10[0]-ra_trans)
# R22_11[0] = ra_factor*(R22_11[0]-ra_trans)
# R22_12[0] = ra_factor*(R22_12[0]-ra_trans)
# R22_13[0] = ra_factor*(R22_13[0]-ra_trans)
# R22_14[0] = ra_factor*(R22_14[0]-ra_trans)
# R22_15[0] = ra_factor*(R22_15[0]-ra_trans)


###
###
### GALPROP MAP
###
###

dir_outfits = '../fits/raw/'
dir_comp = '../fits/'

MODE = 1
if MODE == 0:
	galdef_name = 'm82_05aff1018'
elif MODE == 1:
	galdef_name = 'm82_05b2ff1021'
elif MODE == 2:
	galdef_name = 'm82_06bff1021'
	
MODEE = 0
if MODEE == 0:
	MAP_FREQ = 1.e3 #MeV
elif MODEE == 1:
	MAP_FREQ = 1.e6 #MeV

chi_dat = np.loadtxt(dir_comp+galdef_name+'/'+galdef_name+'chi.dat')
chi_fac = chi_dat[4]
FACTOR = chi_fac	

emissfilename = dir_outfits+galdef_name+'/'+galdef_name+'gamma_total.fits'

##################
### READ DATA
###
hdulist = py.open(emissfilename)
head = hdulist[0].header
scidata = hdulist[0].data

xmin = head['CRVAL1'] #arcsec
ymin = head['CRVAL2'] #arcsec
nmin = head['CRVAL3'] #Hz
xmin = min(xmin,ymin)

xdel = head['CDELT1'] #arcsec
ydel = head['CDELT2'] #arcsec
nfactor = head['CDELT3'] #freq=nmin*nfactor^i

xnum = head['NAXIS1']
ynum = head['NAXIS2']
nnum = head['NAXIS3']

inclination = np.zeros(3)
inclination[0] = head['INCLN1']
objectdist = head['OBJDIS']

freq = nmin*nfactor**np.linspace(0, nnum-1, nnum)
x = np.linspace(0, xnum-1, xnum)*xdel+xmin
y = np.linspace(0, ynum-1, ynum)*ydel+ymin

###
### INTERPOLATING
###
# MAP_FREQ = 1.e6 #MeV

LOW = scidata[freq<MAP_FREQ]
nu_low = freq[freq<MAP_FREQ]
LOW = LOW[-1]
nu_low = nu_low[-1]

HIGH = scidata[freq>MAP_FREQ]
nu_high = freq[freq>MAP_FREQ]
HIGH = HIGH[0]
nu_high = nu_high[0]

gp1 = LOW + (HIGH-LOW)*math.log(MAP_FREQ/nu_low)/math.log(nu_high/nu_low)

### UNITS
gp1*= FACTOR

##
## ROTATION
##

ANGLE = -1*20.
ANGLE *= math.pi/180 
X, Y = DoRotation(x, y, RotRad=ANGLE)

# x_flat = X.flatten(order='F')
# y_flat = Y.flatten(order='F')
# gp1_flat = gp1.flatten(order='F')

# ### BEAM SIZE
# NEW = np.zeros(gp1.shape)

# beamx = 7.6
# beamy = 7.3

# for i in range(0,xnum):
	# for j in range(0,ynum):
		# d_pix2 = ((x_flat-X[i][j])**2/beamx**2+(y_flat-Y[i][j])**2/beamy**2)
		# d_pix2 = np.exp(-d_pix2*4*math.log(2))
		# NEW[i][j] = np.sum(gp1_flat*d_pix2)/np.sum(d_pix2)

# gp1 = NEW*(xdel*ydel/Sr2sqArcSec)*math.pi

######################################
## PLOT DATA
##

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

cmap = plt.get_cmap('cubehelix')
fig = plt.figure()
	
f1, ax = plt.subplots(1,1)

for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(2)
ax.tick_params(which='major',width=1, length=8, labelsize='small', direction='in', color='white')
ax.tick_params(which='minor',width=0.5, length=5, direction='in', color='white')

ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
ax.yaxis.set_minor_locator(plt.MultipleLocator(10))

MSIZE = 5
WIDTH = 3

## PLOT

loglow = -3.
loghigh = 3.
V = np.logspace(loglow, loghigh, 60)
# NN = (np.log10(R22_flux)-loglow)/(loghigh-loglow)

# ax.plot(R22_01[0], R22_01[1], color=cmap(NN[0]), linewidth=WIDTH)
# ax.plot(R22_02[0], R22_02[1], color=cmap(NN[1]), linewidth=WIDTH)
# ax.plot(R22_03[0], R22_03[1], color=cmap(NN[2]), linewidth=WIDTH)
# ax.plot(R22_04[0], R22_04[1], color=cmap(NN[3]), linewidth=WIDTH)
# ax.plot(R22_05[0], R22_05[1], color=cmap(NN[4]), linewidth=WIDTH)
# ax.plot(R22_06[0], R22_06[1], color=cmap(NN[5]), linewidth=WIDTH)
# ax.plot(R22_07[0], R22_07[1], color=cmap(NN[6]), linewidth=WIDTH)
# ax.plot(R22_08[0], R22_08[1], color=cmap(NN[7]), linewidth=WIDTH)
# ax.plot(R22_09[0], R22_09[1], color=cmap(NN[8]), linewidth=WIDTH)
# ax.plot(R22_10[0], R22_10[1], color=cmap(NN[9]), linewidth=WIDTH)
# ax.plot(R22_11[0], R22_11[1], color=cmap(NN[10]), linewidth=WIDTH)
# ax.plot(R22_12[0], R22_12[1], color=cmap(NN[11]), linewidth=WIDTH)
# ax.plot(R22_13[0], R22_13[1], color=cmap(NN[12]), linewidth=WIDTH)
# ax.plot(R22_14[0], R22_14[1], color=cmap(NN[13]), linewidth=WIDTH)
# ax.plot(R22_15[0], R22_15[1], color=cmap(NN[14]), linewidth=WIDTH)

#GALPROP
#TICKS
loglow = -3.
loghigh = 3.
V1 = np.logspace(loglow, loghigh, 7)

#plot
CF1 = ax.contourf(X,Y,gp1,norm = LogNorm(), levels=V, cmap=cmap)

#colorbar
cf1 = plt.colorbar(CF1, ax=ax, ticks=V1, pad=0.0025, fraction=0.04763)
cf1.ax.tick_params(labelsize='small')
cf1.ax.yaxis.set_tick_params(color='w', width=1, length=5)
cf1.outline.set_linewidth(2) 
cf1.set_label(r'MeV$^2$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$', rotation=270, horizontalalignment='center', verticalalignment='bottom')

#axes stuff
ax.set_axis_bgcolor('black')
ax.set_aspect('equal')
ax.set_xlim([100,-100])
ax.set_ylim([-100,100])
ax.set_xticks([-80, -40, 0, 40, 80])
ax.set_yticks([-80, -40, 0, 40, 80])
ax.set_xlabel('Relative R.A. [\"]')
ax.set_ylabel('Relative Dec. [\"]')

if MODEE == 0:
	ax.text(0.95,0.1,r'1 GeV', ha='right', va='top', color='w', transform=ax.transAxes)
	EN = '1Gev'
elif MODEE == 1:
	ax.text(0.95,0.1,r'1 TeV', ha='right', va='top', color='w', transform=ax.transAxes)
	EN = '1TeV'

if MODE == 0:
	ax.text(0.05,0.9,'Model A', ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_gamma'+EN+'_map_A'
elif MODE == 1:
	ax.text(0.05,0.9,'Model B', ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_gamma'+EN+'_map_B'
elif MODE == 2:
	ax.text(0.05,0.9,"Model B'", ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_gamma'+EN+'_map_Bp'

# #beam
# from matplotlib.patches import Ellipse
# beam = Ellipse(xy=[80,-80] , width=7.6, height=7.3)
# ax.add_artist(beam)
# beam.set_facecolor('white')
# beam.set_alpha(0.6)

# plot_name = 'm82_gamma_ray_map_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.pdf'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.pdf', bbox_inches='tight', dpi=400)
		break

plt.close()		



























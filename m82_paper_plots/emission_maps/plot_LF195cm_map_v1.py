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

dir_dat = './m82_LOFAR_154MHz/'
file_begin = 'm82_LOFAR_154MHz_'
file = dir_dat+file_begin

#195cm Data
LF1_01_0 = np.loadtxt(file+'01_0.txt').T
LF1_01_1 = np.loadtxt(file+'01_1.txt').T
LF1_01_2 = np.loadtxt(file+'01_2.txt').T
LF1_02_0 = np.loadtxt(file+'02_0.txt').T
LF1_02_1 = np.loadtxt(file+'02_1.txt').T
LF1_03_0 = np.loadtxt(file+'03_0.txt').T
LF1_04_0 = np.loadtxt(file+'04_0.txt').T
LF1_05_0 = np.loadtxt(file+'05_0.txt').T
LF1_06_0 = np.loadtxt(file+'06_0.txt').T
LF1_06_1 = np.loadtxt(file+'06_1.txt').T
LF1_07_0 = np.loadtxt(file+'07_0.txt').T
LF1_07_1 = np.loadtxt(file+'07_1.txt').T
LF1_07_2 = np.loadtxt(file+'07_2.txt').T
LF1_07_3 = np.loadtxt(file+'07_3.txt').T
LF1_flux_t = np.loadtxt(file+'flux.txt').T
LF1_flux_t = LF1_flux_t[0]*LF1_flux_t[1]*1.e-3
print LF1_flux_t

LF1_flux = np.zeros(14)
LF1_flux[0:3] = LF1_flux_t[0]
LF1_flux[3:5] = LF1_flux_t[1]
LF1_flux[5] = LF1_flux_t[2]
LF1_flux[6] = LF1_flux_t[3]
LF1_flux[7] = LF1_flux_t[4]
LF1_flux[8:10] = LF1_flux_t[5]
LF1_flux[10:14] = LF1_flux_t[6]

print LF1_flux

#Angle translation
# dec_trans = 47 #"
# LF1_01_0[1]-= dec_trans
# LF1_01_1[1]-= dec_trans
# LF1_01_2[1]-= dec_trans
# LF1_02_0[1]-= dec_trans
# LF1_02_1[1]-= dec_trans
# LF1_03_0[1]-= dec_trans
# LF1_04_0[1]-= dec_trans
# LF1_05_0[1]-= dec_trans
# LF1_06_0[1]-= dec_trans
# LF1_06_1[1]-= dec_trans
# LF1_07_0[1]-= dec_trans
# LF1_07_1[1]-= dec_trans
# LF1_07_2[1]-= dec_trans
# LF1_07_3[1]-= dec_trans

# ra_trans = -8.2
# dec = 69.68*math.pi/180
# ra_factor = 360/24*math.cos(dec)#(1-math.cos(dec))/(math.cos(dec))**(-1)
# LF1_01_0[0] = ra_factor*(LF1_01_0[0]-ra_trans)
# LF1_01_1[0] = ra_factor*(LF1_01_1[0]-ra_trans)
# LF1_01_2[0] = ra_factor*(LF1_01_2[0]-ra_trans)
# LF1_02_0[0] = ra_factor*(LF1_02_0[0]-ra_trans)
# LF1_02_1[0] = ra_factor*(LF1_02_1[0]-ra_trans)
# LF1_03_0[0] = ra_factor*(LF1_03_0[0]-ra_trans)
# LF1_04_0[0] = ra_factor*(LF1_04_0[0]-ra_trans)
# LF1_05_0[0] = ra_factor*(LF1_05_0[0]-ra_trans)
# LF1_06_0[0] = ra_factor*(LF1_06_0[0]-ra_trans)
# LF1_06_1[0] = ra_factor*(LF1_06_1[0]-ra_trans)
# LF1_07_0[0] = ra_factor*(LF1_07_0[0]-ra_trans)
# LF1_07_1[0] = ra_factor*(LF1_07_1[0]-ra_trans)
# LF1_07_2[0] = ra_factor*(LF1_07_2[0]-ra_trans)
# LF1_07_3[0] = ra_factor*(LF1_07_3[0]-ra_trans)

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

chi_dat = np.loadtxt(dir_comp+galdef_name+'/'+galdef_name+'chi.dat')
chi_fac = chi_dat[4]
FACTOR = chi_fac	

emissfilename = dir_outfits+galdef_name+'/'+galdef_name+'radio.fits'

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
MAP_FREQ = 0.154e9 #GHz

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
gp1*= 1e23*FACTOR

##
## ROTATION
##

ANGLE = -1*20.
ANGLE *= math.pi/180 
X, Y = DoRotation(x, y, RotRad=ANGLE)

x_flat = X.flatten(order='F')
y_flat = Y.flatten(order='F')
gp1_flat = gp1.flatten(order='F')

### BEAM SIZE
NEW = np.zeros(gp1.shape)

beamx = 4.66/2.
beamy = 3.56/2.

for i in range(0,xnum):
	for j in range(0,ynum):
		d_pix2 = ((x_flat-X[i][j])**2/beamx**2+(y_flat-Y[i][j])**2/beamy**2)
		d_pix2 = np.exp(-d_pix2/2.)
		NEW[i][j] = np.sum(gp1_flat*d_pix2)#/np.sum(d_pix2)

gp1 = NEW*(xdel*ydel/Sr2sqArcSec)*2./np.pi

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
WIDTH = 1

## PLOT 1
# loglow1 = -2
# loghigh1 = 3
# V = np.logspace(loglow1, loghigh1, 81)

# loglow1 = -9
# loghigh1 = 1 
# V1 = np.logspace(loglow1, loghigh1, 21)

# WIDTH = 3 #2.5

#LOFAR 154MHz
# NN = (np.log10(R22_flux)-loglow1)/(loghigh1-loglow1)
# ax[0].plot(LF1_01_0[0], LF1_01_0[1], color=cmap(NN[0]), linewidth=WIDTH)
# ax[0].plot(LF1_01_1[0], LF1_01_1[1], color=cmap(NN[1]), linewidth=WIDTH)
# ax[0].plot(LF1_01_2[0], LF1_01_2[1], color=cmap(NN[2]), linewidth=WIDTH)
# ax[0].plot(LF1_02_0[0], LF1_02_0[1], color=cmap(NN[3]), linewidth=WIDTH)
# ax[0].plot(LF1_02_1[0], LF1_02_1[1], color=cmap(NN[4]), linewidth=WIDTH)
# ax[0].plot(LF1_03_0[0], LF1_03_0[1], color=cmap(NN[5]), linewidth=WIDTH)
# ax[0].plot(LF1_04_0[0], LF1_04_0[1], color=cmap(NN[6]), linewidth=WIDTH)
# ax[0].plot(LF1_05_0[0], LF1_05_0[1], color=cmap(NN[7]), linewidth=WIDTH)
# ax[0].plot(LF1_06_0[0], LF1_06_0[1], color=cmap(NN[8]), linewidth=WIDTH)
# ax[0].plot(LF1_06_1[0], LF1_06_1[1], color=cmap(NN[9]), linewidth=WIDTH)
# ax[0].plot(LF1_07_0[0], LF1_07_0[1], color=cmap(NN[10]), linewidth=WIDTH)
# ax[0].plot(LF1_07_1[0], LF1_07_1[1], color=cmap(NN[11]), linewidth=WIDTH)
# ax[0].plot(LF1_07_2[0], LF1_07_2[1], color=cmap(NN[12]), linewidth=WIDTH)
# ax[0].plot(LF1_07_3[0], LF1_07_3[1], color=cmap(NN[13]), linewidth=WIDTH)
# ax[0].plot(R22_15[0], R22_15[1], color=cmap(NN[14]), linewidth=WIDTH)

# ax.set_axis_bgcolor('black')
# ax.set_aspect('equal')
# ax.set_xlim([-80,80])
# ax.set_ylim([-80,80])
# ax.set_xticks(np.arange(-80, 81, 40))
# ax.set_yticks(np.arange(-80, 81, 40))

# ax.set_xlabel('Relative R.A. [\"]')
# ax.set_ylabel('Relative Dec. [\"]')
# ax.set_title(r'$E^2 \Phi_{\gamma} (E)$, $E=$'+str(gp2_energy)+' GeV', verticalalignment='bottom')

#GALPROP
# CF0 = ax.contourf(X2,Y2,gp2,norm = LogNorm(), levels=V, cmap=cmap_gamma)
# cf0 = plt.colorbar(CF0,ax=ax,pad=0,ticks=V, fraction=0.04763)
# cf0.set_label(r'MeV$^2$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$', rotation=270, horizontalalignment='center', verticalalignment='bottom')

## PLOT 2

loglow = -4.
loghigh = 1.
V = np.logspace(loglow, loghigh, 41)
NN = (np.log10(LF1_flux)-loglow)/(loghigh-loglow)

ax.plot(LF1_01_0[0], LF1_01_0[1], color=cmap(NN[0]), linewidth=WIDTH)
ax.plot(LF1_01_1[0], LF1_01_1[1], color=cmap(NN[1]), linewidth=WIDTH)
ax.plot(LF1_01_2[0], LF1_01_2[1], color=cmap(NN[2]), linewidth=WIDTH)
ax.plot(LF1_02_0[0], LF1_02_0[1], color=cmap(NN[3]), linewidth=WIDTH)
ax.plot(LF1_02_1[0], LF1_02_1[1], color=cmap(NN[4]), linewidth=WIDTH)
ax.plot(LF1_03_0[0], LF1_03_0[1], color=cmap(NN[5]), linewidth=WIDTH)
ax.plot(LF1_04_0[0], LF1_04_0[1], color=cmap(NN[6]), linewidth=WIDTH)
ax.plot(LF1_05_0[0], LF1_05_0[1], color=cmap(NN[7]), linewidth=WIDTH)
ax.plot(LF1_06_0[0], LF1_06_0[1], color=cmap(NN[8]), linewidth=WIDTH)
ax.plot(LF1_06_1[0], LF1_06_1[1], color=cmap(NN[9]), linewidth=WIDTH)
ax.plot(LF1_07_0[0], LF1_07_0[1], color=cmap(NN[10]), linewidth=WIDTH)
ax.plot(LF1_07_1[0], LF1_07_1[1], color=cmap(NN[11]), linewidth=WIDTH)
ax.plot(LF1_07_2[0], LF1_07_2[1], color=cmap(NN[12]), linewidth=WIDTH)
ax.plot(LF1_07_3[0], LF1_07_3[1], color=cmap(NN[13]), linewidth=WIDTH)

#GALPROP
#TICKS
loglow = -4.
loghigh = 1.
V1 = np.logspace(loglow, loghigh, 6)

CF1 = ax.contourf(X,Y,gp1,norm = LogNorm(), levels=V, cmap=cmap)
#colorbar
cf1 = plt.colorbar(CF1, ax=ax, ticks=V1, pad=0.0025, fraction=0.04763)
cf1.ax.tick_params(labelsize='small')
cf1.ax.yaxis.set_tick_params(color='w', width=1, length=5)
cf1.outline.set_linewidth(2) 
cf1.set_label('Jy/beam', rotation=270, horizontalalignment='center', verticalalignment='bottom')

#axes stuff
ax.set_axis_bgcolor('black')
ax.set_aspect('equal')
ax.set_xlim([100,-100])
ax.set_ylim([-100,100])
ax.set_xticks([-80, -40, 0, 40, 80])
ax.set_yticks([-80, -40, 0, 40, 80])
ax.set_xlabel('Relative R.A. [\"]')
ax.set_ylabel('Relative Dec. [\"]')

ax.text(0.95,0.1,r'195 cm', ha='right', va='top', color='w', transform=ax.transAxes)
if MODE == 0:
	ax.text(0.05,0.9,'Model A', ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_154MHz_map_A'
elif MODE == 1:
	ax.text(0.05,0.9,'Model B', ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_154MHz_map_B'
elif MODE == 2:
	ax.text(0.05,0.9,"Model B'", ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_154MHz_map_Bp'

#beam
from matplotlib.patches import Ellipse
beam = Ellipse(xy=[80,-80] , width=4.66, height=3.57)
ax.add_artist(beam)
beam.set_facecolor('white')
beam.set_alpha(0.6)

# plot_name = 'm82_154MHz_map_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.png', bbox_inches='tight', dpi=400)
		break

plt.close()	



























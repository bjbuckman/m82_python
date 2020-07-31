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

data_dir = '/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm-22cm_spectral_index/6-22cm_' 

#data
R01 = np.loadtxt(data_dir+'01.txt').T
R02 = np.loadtxt(data_dir+'02.txt').T
R03 = np.loadtxt(data_dir+'03.txt').T
R04 = np.loadtxt(data_dir+'04.txt').T
R05 = np.loadtxt(data_dir+'05.txt').T

R_flux = -0.3-np.linspace(0,4,5)*0.3
# R_flux = R_flux[::-1]
print R_flux

#Angle translation
dec_trans = 47 #"
R01[1]-= dec_trans
R02[1]-= dec_trans
R03[1]-= dec_trans
R04[1]-= dec_trans
R05[1]-= dec_trans

ra_trans = -8.2
dec = 69.68*math.pi/180
ra_factor = 360/24*math.cos(dec)#(1-math.cos(dec))/(math.cos(dec))**(-1)
R01[0] = ra_factor*(R01[0]-ra_trans)
R02[0] = ra_factor*(R02[0]-ra_trans)
R03[0] = ra_factor*(R03[0]-ra_trans)
R04[0] = ra_factor*(R04[0]-ra_trans)
R05[0] = ra_factor*(R05[0]-ra_trans)

###
###
### GALPROP MAP
###
###

dir_outfits = '../fits/raw/'
dir_comp = '../fits/'

galdef_name = 'm82_03bff1013'

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
MAP_FREQ0 = 5.000e9 #GHz
MAP_FREQ1 = 1.364e9 #GHz

### First Frequency
LOW = scidata[freq<MAP_FREQ0]
nu_low = freq[freq<MAP_FREQ0]
LOW = LOW[-1]
nu_low = nu_low[-1]

HIGH = scidata[freq>MAP_FREQ0]
nu_high = freq[freq>MAP_FREQ0]
HIGH = HIGH[0]
nu_high = nu_high[0]

gp_0 = LOW + (HIGH-LOW)*math.log(MAP_FREQ0/nu_low)/math.log(nu_high/nu_low)

### Second Frequency
LOW = scidata[freq<MAP_FREQ1]
nu_low = freq[freq<MAP_FREQ1]
LOW = LOW[-1]
nu_low = nu_low[-1]

HIGH = scidata[freq>MAP_FREQ1]
nu_high = freq[freq>MAP_FREQ1]
HIGH = HIGH[0]
nu_high = nu_high[0]

gp_1 = LOW + (HIGH-LOW)*math.log(MAP_FREQ1/nu_low)/math.log(nu_high/nu_low)

### UNITS
gp_0*= 10e23*FACTOR
gp_1*= 10e23*FACTOR

##
## ROTATION
##

ANGLE = -1*20.
ANGLE *= math.pi/180 
X, Y = DoRotation(x, y, RotRad=ANGLE)

x_flat = X.flatten(order='F')
y_flat = Y.flatten(order='F')
gp0_flat = gp_0.flatten(order='F')
gp1_flat = gp_1.flatten(order='F')

# ### BEAM SIZE 1

# beamx = 12.5
# beamy = 11.7

# NEW0 = np.zeros(gp_0.shape)
# for i in range(0,xnum):
	# for j in range(0,ynum):
		# d_pix2 = ((x_flat-X[i][j])**2/beamx**2+(y_flat-Y[i][j])**2/beamy**2)
		# d_pix2 = np.exp(-d_pix2*4*math.log(2))
		# NEW0[i][j] = np.sum(gp0_flat*d_pix2)/np.sum(d_pix2)
		
# gp_0 = NEW0
# gp0_flat = gp_0.flatten(order='F')
		
### BEAM SIZE 2

beamx = 12.7
beamy = 11.8

NEW0 = np.zeros(gp_0.shape)
NEW1 = np.zeros(gp_1.shape)
for i in range(0,xnum):
	for j in range(0,ynum):
		d_pix2 = ((x_flat-X[i][j])**2/beamx**2+(y_flat-Y[i][j])**2/beamy**2)
		d_pix2 = np.exp(-d_pix2*4*math.log(2))
		NEW0[i][j] = np.sum(gp0_flat*d_pix2)/np.sum(d_pix2)
		NEW1[i][j] = np.sum(gp1_flat*d_pix2)/np.sum(d_pix2)	

gp_0 = NEW0*(xdel*ydel/Sr2sqArcSec)
gp_1 = NEW1*(xdel*ydel/Sr2sqArcSec)

## Calculate spectral index
alpha = np.log(gp_1/gp_0)/np.log(MAP_FREQ1/MAP_FREQ0)

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
ALPHA = 0.7

## PLOT 2

low_a = -2.4
high_a = -0.
V = np.linspace(low_a, high_a, num=25)
NN = ((R_flux)-low_a)/(high_a-low_a)

ax.plot(R01[0], R01[1], color=cmap(NN[0]), linewidth=WIDTH, alpha=ALPHA)
ax.plot(R02[0], R02[1], color=cmap(NN[1]), linewidth=WIDTH, alpha=ALPHA)
ax.plot(R03[0], R03[1], color=cmap(NN[2]), linewidth=WIDTH, alpha=ALPHA)
ax.plot(R04[0], R04[1], color=cmap(NN[3]), linewidth=WIDTH, alpha=ALPHA)
ax.plot(R05[0], R05[1], color=cmap(NN[4]), linewidth=WIDTH, alpha=ALPHA)

#GALPROP
#TICKS
# loglow = -4.
# loghigh = 1.
# V1 = np.logspace(loglow, loghigh, 21)

low_a = -2.4
high_a = -0.
V = np.linspace(low_a, high_a, num=25)

low_a = -2.4
high_a = -0.
V1 = np.linspace(low_a, high_a, num=7)

CF1 = ax.contourf(X,Y,alpha, cmap=cmap, levels=V)
#colorbar
cf1 = plt.colorbar(CF1, ax=ax, ticks=V1, pad=0.0025, fraction=0.04763)
cf1.ax.tick_params(labelsize='small')
cf1.ax.yaxis.set_tick_params(color='w', width=1, length=5)
cf1.outline.set_linewidth(2) 
cf1.set_label(r'$\alpha$', rotation=270, horizontalalignment='center', verticalalignment='bottom')

#axes stuff
ax.set_axis_bgcolor('white')
ax.set_aspect('equal')
ax.set_xlim([100,-100])
ax.set_ylim([-100,100])
ax.set_xticks([-80, -40, 0, 40, 80])
ax.set_yticks([-80, -40, 0, 40, 80])
ax.set_xlabel('Relative R.A. [\"]')
ax.set_ylabel('Relative Dec. [\"]')

#beam
from matplotlib.patches import Ellipse
beam = Ellipse(xy=[80,-80] , width=12.7, height=11.8)
ax.add_artist(beam)
beam.set_facecolor('white')
beam.set_alpha(0.6)

plot_name = 'm82_spectral_index_map_6-22_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2), bbox_inches='tight', dpi=400)
		break

plt.close()		



























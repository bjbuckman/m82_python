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

# #######################################
# ## READ DATA
# ##
# data_dir = '/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6cm/contour/6cm_' 

# #data
# R01 = np.loadtxt(data_dir+'01.txt').T
# R02 = np.loadtxt(data_dir+'02.txt').T
# R03 = np.loadtxt(data_dir+'03.txt').T
# R04 = np.loadtxt(data_dir+'04.txt').T
# R05 = np.loadtxt(data_dir+'05.txt').T
# R06 = np.loadtxt(data_dir+'06.txt').T
# R07 = np.loadtxt(data_dir+'07.txt').T
# R08 = np.loadtxt(data_dir+'08.txt').T
# R09 = np.loadtxt(data_dir+'09.txt').T
# R10 = np.loadtxt(data_dir+'10.txt').T
# R11 = np.loadtxt(data_dir+'11.txt').T
# R12 = np.loadtxt(data_dir+'12.txt').T
# R13 = np.loadtxt(data_dir+'13.txt').T
# # R14 = np.loadtxt(data_dir+'14.txt').T
# # R15 = np.loadtxt(data_dir+'15.txt').T

# R_flux = 2**np.linspace(0,12,13)*90.e-6
# R_flux = R_flux[::-1]
# print R_flux

# #Angle translation
# dec_trans = 47 #"
# R01[1]-= dec_trans
# R02[1]-= dec_trans
# R03[1]-= dec_trans
# R04[1]-= dec_trans
# R05[1]-= dec_trans
# R06[1]-= dec_trans
# R07[1]-= dec_trans
# R08[1]-= dec_trans
# R09[1]-= dec_trans
# R10[1]-= dec_trans
# R11[1]-= dec_trans
# R12[1]-= dec_trans
# R13[1]-= dec_trans
# # R14[1]-= dec_trans
# # R15[1]-= dec_trans

# ra_trans = -8.2
# dec = 69.68*math.pi/180
# ra_factor = 360/24*math.cos(dec)#(1-math.cos(dec))/(math.cos(dec))**(-1)
# R01[0] = ra_factor*(R01[0]-ra_trans)
# R02[0] = ra_factor*(R02[0]-ra_trans)
# R03[0] = ra_factor*(R03[0]-ra_trans)
# R04[0] = ra_factor*(R04[0]-ra_trans)
# R05[0] = ra_factor*(R05[0]-ra_trans)
# R06[0] = ra_factor*(R06[0]-ra_trans)
# R07[0] = ra_factor*(R07[0]-ra_trans)
# R08[0] = ra_factor*(R08[0]-ra_trans)
# R09[0] = ra_factor*(R09[0]-ra_trans)
# R10[0] = ra_factor*(R10[0]-ra_trans)
# R11[0] = ra_factor*(R11[0]-ra_trans)
# R12[0] = ra_factor*(R12[0]-ra_trans)
# R13[0] = ra_factor*(R13[0]-ra_trans)
# # R14[0] = ra_factor*(R14[0]-ra_trans)
# # R15[0] = ra_factor*(R15[0]-ra_trans)


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
MAP_FREQ = 5.000e9 #GHz

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

beamx = 12.5/2.
beamy = 11.7/2.

for i in range(0,xnum):
	for j in range(0,ynum):
		d_pix2 = ((x_flat-X[i][j])**2/beamx**2+(y_flat-Y[i][j])**2/beamy**2)
		d_pix2 = np.exp(-d_pix2/2.)
		NEW[i][j] = np.sum(gp1_flat*d_pix2)#/np.sum(d_pix2)

gp1 = NEW*(xdel*ydel/Sr2sqArcSec)*2./np.pi

flux_outside = 0
flux_total = 0
for ii in range(0,xnum):
	for ij in range(0,ynum):
		flux_total+= gp1[ii][ij]
		if X[ii][ij]**2+Y[ii][ij]**2>=30.**2:
			flux_outside+=gp1[ii][ij]
print 'flux_outside/flux_total'			
print flux_outside/flux_total

######################################
## LOAD DATA
##

data_file = '/mnt/c/Users/psyko/Physics/proj_gal/m82_data/6.fits'
hdulist = py.open(data_file)
head = hdulist[0].header
scidata = hdulist[0].data

RAmin = head['CRVAL1']
RAdel = head['CDELT1']
RAnum = scidata.shape[0]

DEmin = head['CRVAL2']
DEdel = head['CDELT2']
DEnum = scidata.shape[1]

# M82_RA = 148.969687
# M82_DE = 69.679383

RA = RAdel*np.arange(RAnum)
DE = DEdel*np.arange(DEnum)

RA-= RA.min()/2.
DE-= DE.max()/2.

# RA = RAmin+RAdel*np.arange(-RAnum/2,RAnum/2)
# DE = DEmin+DEdel*np.arange(-DEnum/2,DEnum/2)

# RA-= M82_RA
# DE-= M82_DE

RA*=3600.
DE*=3600.

# print RA,DE

RAC, DEC = np.meshgrid(RA,DE)



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

## PLOT

loglow = -4.
loghigh = 1.
V = np.logspace(loglow, loghigh, 21)
# NN = (np.log10(R_flux)-loglow)/(loghigh-loglow)

CF0 = ax.contour(RAC,DEC,scidata,norm = LogNorm(), levels=V, cmap=cmap)

# ax.plot(R01[0], R01[1], color=cmap(NN[0]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R02[0], R02[1], color=cmap(NN[1]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R03[0], R03[1], color=cmap(NN[2]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R04[0], R04[1], color=cmap(NN[3]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R05[0], R05[1], color=cmap(NN[4]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R06[0], R06[1], color=cmap(NN[5]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R07[0], R07[1], color=cmap(NN[6]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R08[0], R08[1], color=cmap(NN[7]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R09[0], R09[1], color=cmap(NN[8]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R10[0], R10[1], color=cmap(NN[9]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R11[0], R11[1], color=cmap(NN[10]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R12[0], R12[1], color=cmap(NN[11]), linewidth=WIDTH, alpha=ALPHA)
# ax.plot(R13[0], R13[1], color=cmap(NN[12]), linewidth=WIDTH, alpha=ALPHA)
# # ax.plot(R14[0], R14[1], color=cmap(NN[13]), linewidth=WIDTH, alpha=ALPHA)
# # ax.plot(R15[0], R15[1], color=cmap(NN[14]), linewidth=WIDTH, alpha=ALPHA)

#GALPROP
#TICKS
loglow = -4.
loghigh = 1.
V1 = np.logspace(loglow, loghigh, 6)

loglow = -4.
loghigh = 1.
V = np.logspace(loglow, loghigh, 41)

#plot
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

ax.text(0.95,0.1,r'6 cm', ha='right', va='top', color='w', transform=ax.transAxes)
if MODE == 0:
	ax.text(0.05,0.9,'Model A', ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_06cm_map_A'
elif MODE == 1:
	ax.text(0.05,0.9,'Model B', ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_06cm_map_B'
elif MODE == 2:
	ax.text(0.05,0.9,"Model B'", ha='left', va='bottom', color='w', transform=ax.transAxes)
	plot_name = 'm82_06cm_map_Bp'

#beam
from matplotlib.patches import Ellipse
beam = Ellipse(xy=[80,-80] , width=12.5, height=11.7)
ax.add_artist(beam)
beam.set_facecolor('white')
beam.set_alpha(0.6)
beam.set_zorder(90)


# plot_name = 'm82_06cm_map_'
for n in range(0,100):
	if os.path.isfile(plot_name+str(n).zfill(2)+'.png'):
		continue
	else:
		plt.savefig(plot_name+str(n).zfill(2)+'.png', bbox_inches='tight', dpi=400)
		break

plt.close()		



























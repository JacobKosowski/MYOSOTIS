#Input parameters for synthetic observation

# Files
project_name='test_snap_III'
filterfile='Filters/Generic/Bessell/Generic-Bessell.V.dat'
# filestar='Examples/Teststar2.txt'
filestar='Snaps/snap_data.txt'
filecloud= 'Examples/NoCloud'

#Output type
filetype='fits' #Enter 'fits' or 'hdf' for respective image outputs. Enter 'both' to recieve both formats. Note that spectroscopy data can only be written as 'fits'
getnoise = True
stretchted = False # Stretching only applies for png images
stretchfactor = 0.5


# LTRprovided = True


Columndensities='sph'#'sph'
OBtreatment='yes'
EXTmodel='Dmodel'

Rv=3.1
metallicityZ=1.0 #should be 1.0 (for solar) or 0.5 (for LMC). this will affect choosing the evolutionary and atmosphere models

alphai=0.0
bettai=0.0
gammai=0.0

res=0.05#0.05
fovx=50
fovy=50 #arcsec #default 20
distance=200000#200000 #pc default 50000, 1AU=4.8481e-6pc

xpix=round(fovx/res) # Size of x-axis
ypix=round(fovy/res) # Size of y-axis

fwhm=0.11
SNR=0.0
noise2add=6.15883e-12   #noise in the unit of flux/pix2=erg/cm2/s/A/pix2

PSFtype='gaussian'
Adaptiveoptics='no'
seeing=0.8
SR=0.8

# spectroscopy = 'no'     # Spectroscopy output, choose 'yes' or 'no'
# lminspec     = 5000.    # Minimum wavelength [A] should be set within your filter transparency
# lmaxspec     = 6000.    # Maximum wavelength [A] 
# Rspec        = 100 #700      # Spectral resolution (please check your SED library, should not be larger than the resolution of the SEDs)
# velocitydis  = 'no' 		# Applies Doppler shit to spectra based on stellar velocity


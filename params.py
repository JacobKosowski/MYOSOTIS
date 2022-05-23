#Input parameters for synthetic observation

project_name='test_BessV'
filterfile='Filters/Generic/Bessell/Generic-Bessell.V.dat'
filestar='Examples/Teststar2.txt'
# filestar='Snaps/snap_data.txt'
Columndensities='sph'
filecloud= 'Examples/NoCloud'
OBtreatment='yes'
EXTmodel='Dmodel'
Rv=3.1
metallicityZ=1.0 #should be 1.0 (for solar) or 0.5 (for LMC). this will affect choosing the evolutionary and atmosphere models
alphai=0.0
bettai=0.0
gammai=0.0
res=0.05 #0.05
fovx=20. #20
fovy=20. #arcsec
distance=50000 #pc default 50000, 1AU=4.8481e-6pc
fwhm=0.11
SNR=0.0
noise2add=6.15883e-12   #noise in the unit of flux/pix2=erg/cm2/s/A/pix2
spectroscopy = 'no'     # Spectroscopy output, choose 'yes' or 'no'
lminspec     = 5000.    # Minimum wavelength [A] should be set within your filter transparency
lmaxspec     = 6000.    # Maximum wavelength [A] 
Rspec        = 100 #700      # Spectral resolution (please check your SED library, should not be larger than the resolution of the SEDs)
velocitydis  = 'no'
PSFtype='gaussian'
Adaptiveoptics='no'
seeing=0.8
SR=0.8

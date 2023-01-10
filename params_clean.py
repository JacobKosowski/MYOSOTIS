#Input parameters for synthetic observation

#Resources
Ncpus = 16 # Number of CPUS used in Multiprocessing (SED Reading)
Numbacpus = 16 # Number of CPUS used by Numba
numba_parallel = True

# Files
#filterfile='Filters/hst/acs_wfc/HST-ACS_WFC.F814W_81.dat'
#filterfile='Filters/hst/acs_wfc/HST-ACS_WFC.F606W_81.dat'
#filterfile='Filters/hst/acs_wfc/HST-ACS_WFC.F435W_81.dat'

project_name='mcl_IMF75_100k_Rh6_W3'
filterfile='Filters/generic/bessell/Generic-Bessell.V.dat'
filestar='Snaps/IMF/mcl_IMF75_100k_Rh6_W3.txt'
filecloud= 'Cloud/NoCloud'

#Output type
filetype='fits' #Enter 'fits' or 'hdf' for respective image outputs. Enter 'both' to recieve both formats. Note that spectroscopy data can only be written as 'fits'
getnoise = False
stretchted = False # Stretching only applies for png images
stretchfactor = 0.5


LTRprovided = True


Columndensities='sph'#'sph'
OBtreatment='yes'
EXTmodel='Dmodel'

Rv=3.1
metallicityZ=1.0 #should be 1.0 (for solar) or 0.5 (for LMC). this will affect choosing the evolutionary and atmosphere models

alphai=0.0 #Rotation alpha
bettai=0.0 #Rotation beta
gammai=0.0 #Rotation gamma

res=0.8 #arcsec
fovx=4000
fovy=4000 #arcsec #default 20
distance=5500#200000 #pc default 50000, 1AU=4.8481e-6pc

xpix=round(fovx/res) # Size of x-axis
ypix=round(fovy/res) # Size of y-axis


SNR=0.0
noise2add=6.15883e-12   #noise in the unit of flux/pix2=erg/cm2/s/A/pix2

PSFtype='moffat' # 'gaussian' or 'moffat'
fwhm=1 #Full-width half max. Only applies to gaussian psf
alph = 1.0  #Alpha parameter. Only applies to moffat psf
bet = 1.5 #Beta parameter. Only applies to moffat psf. MUST BE > 1

Adaptiveoptics='no'
seeing=0.8
SR=0.8

spectroscopy = 'no'     # Spectroscopy output, choose 'yes' or 'no'
lminspec     = 4500.    # Minimum wavelength [A] should be set within your filter transparency
lmaxspec     = 7500.    # Maximum wavelength [A] 
Rspec        = 3000 #700      # Spectral resolution (please check your SED library, should not be larger than the resolution of the SEDs)
velocitydis  = 'yes' 		# Applies Doppler shit to spectra based on stellar velocity

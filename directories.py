import params_clean as params

# List of directories for SEDs, etc.

#SEDs
foldersed='SEDs/'
sed_vega='SEDs/vegaf'
sed_z1=foldersed+'Z1/'
sed_z0p5=foldersed+'Z0p5/'

#Isochrones
iso_z0p015='Evolutionary/Z0p015.dat'
iso_z0p008='Evolutionary/Z0p008.dat'

if (params.metallicityZ == 1.0):
    isochrones=iso_z0p015
    foldersed=sed_z1
elif (params.metallicityZ == 0.5):
    isochrones=iso_z0p008
    foldersed=sed_z0p5
else: print('!!!metallicityZ should be 1.0 (for solar metallicity) or 0.5 (for LMC)')

#Interstellar Dust Models
draine3p1 = 'Dust/Draine3p1.txt'
draine4 = 'Dust/Draine4.txt'
draine5p5 = 'Dust/Draine5p5.txt'


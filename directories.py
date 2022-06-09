import params_clean as params
import os, shutil
import functions

#######################################################################################################################
# Structures output files

if not os.path.exists(params.project_name):
    os.makedirs(params.project_name)
else:
    shutil.rmtree(params.project_name)
    os.makedirs(params.project_name)


#######################################################################################################################
# List of directories and key files for SEDs, dust models, etc.

outputim     = params.project_name+'/'+params.project_name+'_image' 
outputimnoise=params.project_name+'/'+params.project_name+'_imageNoise'
outputspecFL = params.project_name+'/'+params.project_name+'_cube_spectra.fits'
outputspecL  = params.project_name+'/'+params.project_name+'_Lambda.txt'
outputstarinfo  = params.project_name+'/'+params.project_name+'_star_info.txt'

#SEDs
foldersed='SEDs/'
sed_vega='SEDs/vegaf'
sed_z1=foldersed+'Z1/'
sed_z0p5=foldersed+'Z0p5/'
mergedseds='merged.hdf'

#Isochrones
iso_z0p015='Evolutionary/Z0p015.dat'
iso_z0p008='Evolutionary/Z0p008.dat'

if (params.metallicityZ == 1.0):
    isochrones=iso_z0p015
    # foldersed=sed_z1+'FEA/'
    foldersed=sed_z1+'HDF/'
    # foldersed=sed_z1+'CSV/'
elif (params.metallicityZ == 0.5):
    isochrones=iso_z0p008
    foldersed=sed_z0p5
else: print('!!!metallicityZ should be 1.0 (for solar metallicity) or 0.5 (for LMC)')



#Interstellar Dust Models
draine3p1 = 'Dust/Draine3p1.txt'
draine4 = 'Dust/Draine4.txt'
draine5p5 = 'Dust/Draine5p5.txt'

if (params.Rv == 3.1): 
    Drainemodel=draine3p1
elif (params.Rv == 4.0): 
    Drainemodel=draine4
elif (params.Rv == 5.5): 
    Drainemodel=draine5p5
else: print('For Dmodel, R_V should be 3.1 or 4.0 or 5.5. If you need other Rv values please choose Fmodel')

def display():

    functions.myso_logo('outim')
    print(outputim)
    print('   ')

    if params.spectroscopy=='yes':
        functions.myso_logo('outspec')
        print(outputspecFL)
        print(outputspecL)
        print('   ')
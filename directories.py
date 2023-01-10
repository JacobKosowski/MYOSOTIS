import params_clean as params
import os, shutil
import functions

#######################################################################################################################
# Structures output files

outdirs = "Runs/"+params.project_name

def outputdirs():
    if not os.path.exists(outdirs):
        os.makedirs(outdirs)
    else:
        shutil.rmtree(outdirs)
        os.makedirs(outdirs)

#######################################################################################################################
# List of directories and key files for SEDs, dust models, etc.

outputim        = outdirs+'/'+params.project_name+'_image' 
outputimnoise   = outdirs+'/'+params.project_name+'_imageNoise'
outputspecFL    = outdirs+'/'+params.project_name+'_cube_spectra.fits'
outputspecL     = outdirs+'/'+params.project_name+'_Lambda.txt'
outputstarinfo  = outdirs+'/'+params.project_name+'_star_info.txt'
outputspecinfo  = outdirs+'/'+params.project_name+'_spec_info.txt'

#SEDs
foldersed='SEDs/'
sed_vega='SEDs/vegaf'
sed_z1=foldersed+'Z1/'
sed_z0p5=foldersed+'Z0p5/'
sedlists='Lists/'
mergedseds='merged.hdf'

#Isochrones
iso_z0p015='Evolutionary/Z0p015.dat'
iso_z0p008='Evolutionary/Z0p008.dat'

if (params.metallicityZ == 1.0):
    isochrones=iso_z0p015
    sedlists=sed_z1+sedlists
    foldersed=sed_z1+'FEA/'

    # foldersed=sed_z1+'HDF/'
    # foldersed=sed_z1+'CSV/'
elif (params.metallicityZ == 0.5):
    isochrones=iso_z0p008
    foldersed=sed_z0p5+'FEA/'
    sedlists=sed_z0p5+sedlists
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


def display_output():

    functions.myso_logo('outim')
    print(outputim)
    print('   ')

    if params.spectroscopy=='yes':
        functions.myso_logo('outspec')
        print(outputspecFL)
        print(outputspecL)
        print('   ')

import numpy as np
import timeit
import params_clean as params
import functions
import directories
import constants
import parallel_functions
import profiling
from datetime import datetime
import pandas as pd

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# This version of MYOSOTIS is built for parallelization.
# It assumes that the user provides Lum, Teff, and Log(g) of the input stars.
# It does not support spectroscopy.
# This was made from a copy of myosotis_clean.py
###################################################################################################################################################################################################
# This particular version is designed to support multiprocessing.
# This was made from a copy of myosotis_parallel.py

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Start of runtime code

def main():

    start = timeit.default_timer()
    s_init = timeit.default_timer()
    functions.myso_logo('logo')

    directories.outputdirs()

    #######################################################################################################################
    # Interstellar dust models
    DrainearrLam, DrainearrK = functions.draine_dust_model()

    #######################################################################################################################
    # 
    print('FoV: ',params.fovx,'" x ',params.fovy,'" =',int(params.xpix),'[pix] x',int(params.ypix),'[pix]')
    sceneim=np.full((int(params.ypix),int(params.xpix)),0.0)


    lambdaF,weight,par2vega = functions.filters()


    #######################################################################################################################
    # 
    massstar,logagestar,kzstar,rhostar,loglstar,Teffstar,loggstar=np.loadtxt(params.filestar,usecols=(6,7,8,9,10,11,12),unpack=True)


    if (params.Columndensities == 'sph'):
        masspar=np.loadtxt(params.filecloud,usecols=(6),unpack=True)


    nstar,newx_,newy_,newz_,vxstar,vystar,vzstar,distancestar,newxcloud,newycloud,newzcloud,newhcloud = functions.rotation()
    newx,newy,newz = newx_,newy_,newz_

    profiling.prediction(nstar)

    pc2pixstar=constants.rad2arcsec/distancestar/params.res 


    #######################################################################################################################
    # Reading list of SEDs

    teffsed,loggsed,metallicity,lh,vtur,sedname=np.loadtxt(directories.sedlists+'kseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4', 'col5', 'col6'), 'formats':(float,float,float,float,float,'|S60')},usecols=(0,1,2,3,4,5),unpack=True)
    nseds=len(teffsed)
    sedname=sedname.astype('U64')

    if (params.OBtreatment == 'yes'):
        teffsedOB,loggsedOB,metallicityOB,sednameOB=np.loadtxt(directories.sedlists+'Tseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4'), 'formats':(float,float,float,'|S60')},usecols=(0,1,2,3),unpack=True)
        nsedsOB=len(teffsedOB)
        functions.myso_logo('tlusty')
        sednameOB=sednameOB.astype('U64')


    nfovstars=0
    jj = 0
    for ii in range(nstar):

        if ((abs(newx_[ii]) < (params.xpix/2)-1) and (abs(newy_[ii]) < (params.ypix/2)-1)): # Check if star is within chosen FOV
            nfovstars += 1
            jj+=1
        else:
            jj=jj
            logagestar = np.delete(logagestar,[jj])
            massstar = np.delete(massstar,[jj])
            kzstar = np.delete(kzstar,[jj])
            newx = np.delete(newx,[jj])
            newy = np.delete(newy,[jj])
            newz = np.delete(newz,[jj])
            pc2pixstar = np.delete(pc2pixstar,[jj])

            Teffstar = np.delete(Teffstar,[jj])
            loggstar = np.delete(loggstar,[jj])
            loglstar = np.delete(loglstar,[jj])

    e_init = timeit.default_timer()

    profiling.init_time = e_init-s_init

    if nfovstars==0:
        raise Exception("No stars in FOV")
    else:
        print('Number of valid stars in the FoV: ',nfovstars)
        print('Image size: '+str(params.xpix)+'x'+str(params.ypix),"\n")


    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # Parallelization here

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    s = timeit.default_timer()

    AVstar, readsed = parallel_functions.clouds_and_SEDs(nfovstars, sedname,nseds,teffsed,loggsed,sednameOB,nsedsOB,teffsedOB,loggsedOB, newx,newy,newz, newxcloud,newycloud,newzcloud,masspar,newhcloud,rhostar, pc2pixstar,Teffstar,loggstar)
    
    e = timeit.default_timer()
    profiling.loop1_time = e-s

    print('Loop 1:',"{:.5}".format(profiling.loop1_time/60),'[min]')
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    s = timeit.default_timer()

    lenWF,indxs = parallel_functions.length_wave_flux(nfovstars,readsed)

    wavelengths = np.empty(lenWF)
    fluxes = np.empty(lenWF)

    wavelengths,fluxes = parallel_functions.sed_load(readsed,lenWF,indxs,wavelengths,fluxes,nfovstars)

    e = timeit.default_timer()

    profiling.sedload_time = e-s
    print('Wavelenght/Flux Setup:',"{:.5}".format(profiling.sedload_time/60),'[min]')
    # #-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    s = timeit.default_timer()

    sceneim,fluxstar,mag = parallel_functions.BC_and_PSF(sceneim,nfovstars,newx,newy,wavelengths,fluxes,indxs,lambdaF,weight,AVstar,Teffstar,loglstar,distancestar,par2vega,DrainearrLam,DrainearrK)
    
    e = timeit.default_timer()
    profiling.loop2_time = e-s

    print('Loop 2:',"{:.5}".format(profiling.loop2_time/60),'[min]')
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    
        # "{:.3}".format(T0/60),'[m]',"{:.3}".format(100*T1/T0),"{:.3}".format(100*T2/T0),"{:.3}".format(100*T3/T0),"{:.3}".format(100*T4/T0),"{:.3}".format(100*T5/T0),"{:.3}".format(100*T6/T0),'%')

    s = timeit.default_timer()

    faintestflux, noise2addim, noise = functions.noise_for_image(fluxstar)

    #######################################################################################################################
    # Create output logs+images and print output filenames and execution time
    sceneimFL = []
    functions.create_image(sceneim,noise2addim,sceneimFL)

    functions.log_output(nfovstars,noise,faintestflux,massstar,logagestar,kzstar,Teffstar,loggstar,loglstar,AVstar,mag,newx+params.xpix/2.0,newy+params.ypix/2.0,readsed,directories.outputstarinfo)

    e = timeit.default_timer()
    profiling.data_save_time = e-s


    directories.display_output()


    stop = timeit.default_timer()
    profiling.total_time=stop-start


    profiling.basic_profile()

    profiling.profiling()

    # profiling.profiling_output_vs_nstars("nstar_output_parallel_feather.txt",nfovstars)
    # profiling.profiling_output_vs_imsize("imsize_output_parallel_feather.txt",)

if __name__ == '__main__':
    main()
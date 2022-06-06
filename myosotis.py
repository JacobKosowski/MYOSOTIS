#!/usr/bin/env python
import numpy as np
import timeit
import params_clean as params
import functions
import directories
import constants
import parallel_functions
import numba.typed as nt
from datetime import datetime
import pandas as pd

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# This version of MYOSOTIS is built for parallelization. It assumes that the user provides Lum, Teff, and Log(g) of the input stars. It does not support spectroscopy.
# This was made from a copy of myosotis_clean.py

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Start of runtime code

start = timeit.default_timer()
functions.myso_logo('logo')

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

pc2pixstar=constants.rad2arcsec/distancestar/params.res 



#######################################################################################################################
# Reading list of SEDs

teffsed,loggsed,metallicity,lh,vtur,sedname=np.loadtxt(directories.foldersed+'kseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4', 'col5', 'col6'), 'formats':(float,float,float,float,float,'|S60')},usecols=(0,1,2,3,4,5),unpack=True)
nseds=len(teffsed)
sedname=sedname.astype('U64')

if (params.OBtreatment == 'yes'):
    teffsedOB,loggsedOB,metallicityOB,sednameOB=np.loadtxt(directories.foldersed+'Tseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4'), 'formats':(float,float,float,'|S60')},usecols=(0,1,2,3),unpack=True)
    nsedsOB=len(teffsedOB)
    functions.myso_logo('tlusty')
    sednameOB=sednameOB.astype('U64')


#wavelength,flux=np.loadtxt(foldersed+sedname[10],comments=['fn:', '#'],unpack=True)
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


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Parallelization here

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s = timeit.default_timer()
AVstar, readsed = parallel_functions.clouds_and_SEDs(nfovstars, sedname,nseds,teffsed,loggsed,sednameOB,nsedsOB,teffsedOB,loggsedOB, newx,newy,newz, newxcloud,newycloud,newzcloud,masspar,newhcloud,rhostar, pc2pixstar,Teffstar,loggstar)
e = timeit.default_timer()
print('Loop 1',(e-s)/60,'[min]')
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
s = timeit.default_timer()


indxs = np.array([0])

nWF = 0
for ii in range(nfovstars):
    if 'NextGen' in readsed[ii]:
        nWF+=21311
    elif 'fnew' in readsed[ii]:
        nWF+=1221
    else:
        nWF+=19997
    indxs = np.append(indxs,nWF)

wavelengths = np.empty(nWF)
flux = np.empty(nWF)

T1 = 0
T2 = 0
for ii in range(nfovstars):
    s_ = timeit.default_timer()
    # w,f=np.loadtxt(directories.foldersed+readsed[ii],comments=['fn:', '#'],unpack=True)
    if 'NextGen' in readsed[ii]:
        w,f = np.array(pd.read_csv(directories.foldersed+readsed[ii],comment='#',sep=r'\s+')).T
    else:
        w,f = np.array(pd.read_csv(directories.foldersed+readsed[ii],comment='#',sep='\t').dropna(axis=1)).T

    e_ = timeit.default_timer()
    T1 +=e_-s_

    s_ = timeit.default_timer()
    wavelengths[indxs[ii]:indxs[ii+1]] = w
    flux[indxs[ii]:indxs[ii+1]] = f
    e_ = timeit.default_timer()
    T2 +=e_-s_

print(T1,T2)

e = timeit.default_timer()
print('Wavelenght/Flux Setup',(e-s)/60,'[min]')
# #-----------------------------------------------------------------------------------------------------------------------------------------------------------------
s = timeit.default_timer()
sceneim,fluxstar,mag = parallel_functions.BC_and_PSF(sceneim,nfovstars,newx,newy,wavelengths,flux,indxs,lambdaF,weight,AVstar,Teffstar,loglstar,distancestar,par2vega,DrainearrLam,DrainearrK)
e = timeit.default_timer()
print('Loop 2',(e-s)/60,'[min]')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


print('Number of valid stars in the FoV: ',nfovstars,"\n")
    # "{:.3}".format(T0/60),'[m]',"{:.3}".format(100*T1/T0),"{:.3}".format(100*T2/T0),"{:.3}".format(100*T3/T0),"{:.3}".format(100*T4/T0),"{:.3}".format(100*T5/T0),"{:.3}".format(100*T6/T0),'%')


faintestflux, noise2addim, noise = functions.noise_for_image(fluxstar)

#######################################################################################################################
# Create output logs+images and print output filenames and execution time
sceneimFL = []
functions.create_image(sceneim,noise2addim,sceneimFL)

functions.log_output(nfovstars,noise,faintestflux,massstar,logagestar,kzstar,Teffstar,loggstar,loglstar,AVstar,mag,newx+params.xpix/2.0,newy+params.ypix/2.0,readsed,directories.outputstarinfo)

functions.myso_logo('outim')
print(directories.outputim)
print('   ')
functions.myso_logo('outspec')
print(directories.outputspecFL)
print(directories.outputspecL)
print('   ')
stop = timeit.default_timer()

print('Simulation time using N core(s): ', (stop - start)/60. ,'[min]')
print("Date and Time =", datetime.now().strftime("%d/%m/%Y %H:%M:%S")) 

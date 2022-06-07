#!/usr/bin/env python
import numpy as np
import timeit
import params_clean as params
import functions
import directories
import constants
import parallel_functions
import profile
import numba.typed as nt
from datetime import datetime
import pandas as pd
import multiprocessing
from itertools import repeat

init_time = 0
loop1_time = 0
loop2_time = 0
sedload_time = 0
data_save_time = 0
total_time = 0


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# This version of MYOSOTIS is built for parallelization. It assumes that the user provides Lum, Teff, and Log(g) of the input stars. It does not support spectroscopy.
# This was made from a copy of myosotis_clean.py

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Start of runtime code

start = timeit.default_timer()
s_init = timeit.default_timer()
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

profile.prediction(nstar)

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

e_init = timeit.default_timer()

init_time = e_init-s_init

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Parallelization here

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
s = timeit.default_timer()
AVstar, readsed = parallel_functions.clouds_and_SEDs(nfovstars, sedname,nseds,teffsed,loggsed,sednameOB,nsedsOB,teffsedOB,loggsedOB, newx,newy,newz, newxcloud,newycloud,newzcloud,masspar,newhcloud,rhostar, pc2pixstar,Teffstar,loggstar)
e = timeit.default_timer()
loop1_time = e-s
print('Loop 1:',(e-s)/60,'[min]')
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
s = timeit.default_timer()


indxs = np.empty(nfovstars+1,dtype=int)
indxs[0] = 0

nWF = 0
for ii in range(nfovstars):
    if 'NextGen' in readsed[ii]:
        nWF+=21312
    elif 'fnew' in readsed[ii]:
        nWF+=1221
    else:
        nWF+=19998
    indxs[ii+1] = nWF

wavelengths = np.empty(nWF)
flux = np.empty(nWF)

# def loaddata():
#     Ncpus = 1
#     pool = multiprocessing.Pool(Ncpus)

#     chunksize = int(len(readsed)/Ncpus)
#     file = directories.foldersed+'merged.hdf'


#     DF = pool.starmap(pd.read_hdf, zip(repeat(file),readsed))
#     pool.close()
#     pool.join()
#     return DF
    

# if __name__ == '__main__': 
#     DF = loaddata()

T1 = 0
T2 = 0
for ii in range(nfovstars):
    s_ = timeit.default_timer()
    

    w,f = np.array(pd.read_hdf(directories.foldersed+'merged.hdf',readsed[ii])).T

    # w,f=np.loadtxt(directories.foldersed+readsed[ii],comments=['fn:', '#'],unpack=True)

    # if 'NextGen' in readsed[ii]:
    #     w,f = np.array(pd.read_csv(directories.foldersed+readsed[ii],comment='#',sep=r'\s+')).T
    # else:
    #     w,f = np.array(pd.read_csv(directories.foldersed+readsed[ii],comment='#',sep='\t').dropna(axis=1)).T

    e_ = timeit.default_timer()
    T1 +=e_-s_

    s_ = timeit.default_timer()
    wavelengths[indxs[ii]:indxs[ii+1]] = w
    flux[indxs[ii]:indxs[ii+1]] = f
    e_ = timeit.default_timer()
    T2 +=e_-s_



e = timeit.default_timer()
# print(T1/(e-s),T2/(e-s))
sedload_time = e-s
print('Wavelenght/Flux Setup:',(e-s)/60,'[min]')
# #-----------------------------------------------------------------------------------------------------------------------------------------------------------------
s = timeit.default_timer()
sceneim,fluxstar,mag = parallel_functions.BC_and_PSF(sceneim,nfovstars,newx,newy,wavelengths,flux,indxs,lambdaF,weight,AVstar,Teffstar,loglstar,distancestar,par2vega,DrainearrLam,DrainearrK)
e = timeit.default_timer()
print('Loop 2:',(e-s)/60,'[min]')
loop2_time = e-s
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


print('Number of valid stars in the FoV: ',nfovstars,"\n")
    # "{:.3}".format(T0/60),'[m]',"{:.3}".format(100*T1/T0),"{:.3}".format(100*T2/T0),"{:.3}".format(100*T3/T0),"{:.3}".format(100*T4/T0),"{:.3}".format(100*T5/T0),"{:.3}".format(100*T6/T0),'%')

s = timeit.default_timer()

faintestflux, noise2addim, noise = functions.noise_for_image(fluxstar)

#######################################################################################################################
# Create output logs+images and print output filenames and execution time
sceneimFL = []
functions.create_image(sceneim,noise2addim,sceneimFL)

functions.log_output(nfovstars,noise,faintestflux,massstar,logagestar,kzstar,Teffstar,loggstar,loglstar,AVstar,mag,newx+params.xpix/2.0,newy+params.ypix/2.0,readsed,directories.outputstarinfo)

e = timeit.default_timer()
data_save_time = e-s

functions.myso_logo('outim')
print(directories.outputim)
print('   ')
functions.myso_logo('outspec')
print(directories.outputspecFL)
print(directories.outputspecL)
print('   ')

stop = timeit.default_timer()
total_time=stop-start

print('Simulation time using N core(s): ', (stop - start)/60. ,'[min]')
print("Date and Time =", datetime.now().strftime("%d/%m/%Y %H:%M:%S")) 
print("")
print("init_time:", init_time)
print("loop1_time:",loop1_time)
print("loop2_time:",loop2_time)
print("sedload_time:",sedload_time)
print("data_save_time:",data_save_time)
print("total_time:",total_time)

# run = [nfovstars,total_time,init_time,loop1_time,sedload_time,loop2_time,data_save_time]

# def file_write(filename,data):
#     data = np.array(data)
#     f = open(filename, "a")

#     for r in data:
#         f.write(str(r))
#         f.write("\t")
#     f.write("\n")

# file_write("output.txt",run)

#!/usr/bin/env python
from math import *
import random
import numpy as np
import timeit
import params
import functions
import directories
import constants

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

sceneimFL = []
if (params.spectroscopy == 'yes'):
    sceneimFL,Lspecarr,nspec = functions.spec_wave_out()

#######################################################################################################################
# 

if params.LTRprovided:
    massstar,logagestar,kzstar,rhostar,loglstar,Teffstar,loggstar=np.loadtxt(params.filestar,usecols=(6,7,8,9,10,11,12),unpack=True)
else: 
    massstar,logagestar,kzstar,rhostar=np.loadtxt(params.filestar,usecols=(6,7,8,9),unpack=True)


if (params.Columndensities == 'sph'):
    masspar=np.loadtxt(params.filecloud,usecols=(6),unpack=True)


nstar,newx,newy,newz,vxstar,vystar,vzstar,distancestar,newxcloud,newycloud,newzcloud,newhcloud = functions.rotation()

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
for ii in range(nstar):
    if ((abs(newx[ii]) < (params.xpix/2)-1) and (abs(newy[ii]) < (params.ypix/2)-1)): # Check if star is within chosen FOV
        nfovstars += 1
    else:
        np.delete(logagestar,[ii])
        np.delete(massstar,[ii])
        np.delete(kzstar,[ii])
        np.delete(newx,[ii])
        np.delete(newy,[ii])
        np.delete(newz,[ii])
        np.delete(pc2pixstar,[ii])

        if params.LTRprovided:
            np.delete(Teffstar,[ii])
            np.delete(loggstar,[ii])
            np.delete(loglstar,[ii])

if not params.LTRprovided:
    Teffstar=np.zeros(nfovstars)
    loggstar=np.zeros(nfovstars)
    loglstar=np.zeros(nfovstars)
sedstar=np.zeros(nfovstars,dtype=np.uint64)
fluxstar=np.zeros(nfovstars)
columncloud=np.zeros(nfovstars) #column density of the cloud in front of each star
AVstar=np.zeros(nfovstars)
mag=np.zeros(nfovstars)
readsed=np.full(nfovstars,None)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Add parallelization here
# Each star is fully independent so this loop can be split as many times as needed

#Varibles for loop
#Inputs
#   xpix, ypix,newx, newy, newz, newxcloud, newycloud, newzcloud, masspar, newhcloud, pc2pixstar, rhostar, logagestar, massstar, kzstar, ziso, logageiso, miniiso, mactiso, logliso, logteff, loggiso,
#   Additional inputs********
T0,T1,T2,T3,T4,T5,T6 = 0,0,0,0,0,0,0

s0 = timeit.default_timer()
iso_load_flag = True
nstars=0
for ii in range(nfovstars):
    if True :#(fovstars_check[ii]): # Check if star is within chosen FOV
        nstars += 1
        # break
    s = timeit.default_timer()

    # Density of gas clouds
    if (params.Columndensities == 'sph'):
        columncloud[ii]=functions.COLUMN_DENSITY(newx[ii],newy[ii],newz[ii],newxcloud,newycloud,newzcloud,masspar,newhcloud)
        cloudden=columncloud[ii]*((1.989/1.673534)/(((3.0857/pc2pixstar[ii])**2.))) #*1.0e21 ;convert column density unit [Msun/pix^2] --> [10^21 * Mhydrogen/cm^2]
    else: cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857**2.)))#*1.0e21 ;convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]
    AVstar[ii]=cloudden/2.21 #e21 ;Guver&Ozel2009: The relation between Optical Extinction and Hydrogen column density
    e = timeit.default_timer()
    T1+=(e-s)

    s = timeit.default_timer()

    # If lum, teff, and radius are not given, use isochrones to calculate appropriate values
    if not params.LTRprovided:
        if iso_load_flag:
            functions.myso_logo('iso')
            ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso=np.loadtxt(directories.isochrones,usecols=(0,1,2,3,4,5,6),unpack=True)
            iso_load_flag = False

        Teffstar[ii], loggstar[ii], loglstar[ii], flag = functions.get_iso_params(logagestar[ii],massstar[ii],kzstar[ii],ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso)
        if flag: nstars-=1; continue

    e = timeit.default_timer()
    T2+=(e-s)

    #############################################################
    # Matching and finding seds
    s = timeit.default_timer()
    if ((params.OBtreatment == 'no') or (Teffstar[ii] < 15000.)):

        deltaT=abs(Teffstar[ii]-teffsed)
        deltagarr=np.full(nseds,99.)

        for jj in range(nseds): 
            if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsed[jj])
    #        sedstar[ii]=where(deltagarr eq min(deltagarr))
        sedstar[ii]=functions.FINDCLOSE(min(deltagarr),deltagarr)

        readsed[ii]=sedname[sedstar[ii]]
        wavelength,flux=np.loadtxt(directories.foldersed+sedname[sedstar[ii]],comments=['fn:', '#'],unpack=True)

    elif ((params.OBtreatment == 'yes') and (Teffstar[ii] >= 15000.)):

        deltaT=abs(Teffstar[ii]-teffsedOB)
        deltagarr=np.full(nsedsOB,99.)

        for jj in range(nsedsOB):  
            if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsedOB[jj])
#            sedstar[ii]=where(deltagarr eq min(deltagarr))
        sedstar[ii]=functions.FINDCLOSE(min(deltagarr),deltagarr)

        readsed[ii]=sednameOB[sedstar[ii]]
        wavelength,flux=np.loadtxt(directories.foldersed+sednameOB[sedstar[ii]],comments=['fn:', '#'],unpack=True)

    e = timeit.default_timer()
    T3+=(e-s)

    s = timeit.default_timer()
    bc1=functions.BCcal(wavelength,flux,lambdaF,weight,AVstar[ii],params.Rv,Teffstar[ii],par2vega,params.EXTmodel,DrainearrLam,DrainearrK)
    mag[ii]=constants.Mbolsun-2.5*loglstar[ii]-(bc1)+5.0*np.log10(distancestar[ii]/10.0)
    fluxstar[ii]=10.**(mag[ii]/(-2.5))
    e = timeit.default_timer()
    T4+=(e-s)


    s = timeit.default_timer()

    if (params.PSFtype == 'gaussian'):
        airy1=functions.makeGaussian((params.xpix,params.ypix),params.fwhm/params.res,center=(newx[ii]+params.xpix/2.,newy[ii]+params.ypix/2.))

    airy1=airy1/np.sum(airy1) #to be sure about normaized total flux across FOV 

    if (params.Adaptiveoptics == 'yes'):
        halo=functions.makeGaussian((params.xpix,params.ypix),params.seeing/params.res,center=(newx[ii]+params.xpix/2.,newy[ii]+params.ypix/2.))
        halo=halo/np.sum(halo)
        sceneim += fluxstar[ii]*(params.SR*airy1+(1.0-params.SR)*halo)
    if (params.Adaptiveoptics == 'no'): sceneim += fluxstar[ii]*airy1

    e = timeit.default_timer()
    T5+=(e-s)

    ###############################################################################################
    # Spectroscopy
    # Can be removed

    s = timeit.default_timer()  
    if (params.spectroscopy == 'yes'):
        sceneimFL = functions.spectroscopy(sceneimFL,nspec,wavelength,lambdaF,weight,flux,Lspecarr,AVstar[ii],Teffstar[ii],DrainearrLam,DrainearrK,loglstar[ii],distancestar[ii],newx[ii],newy[ii],vzstar[ii])
    
    e = timeit.default_timer()
    T6+=(e-s)
    # break
#Outputs
#   nstars, sceneim, sceneimFL,
e0 = timeit.default_timer()
T0=(e0-s0)

print('Number of valid stars in the FoV: ',nstars,"\n"
    "{:.3}".format(T0/60),'[m]',"{:.3}".format(100*T1/T0),"{:.3}".format(100*T2/T0),"{:.3}".format(100*T3/T0),"{:.3}".format(100*T4/T0),"{:.3}".format(100*T5/T0),"{:.3}".format(100*T6/T0),'%')


faintestflux, noise2addim, noise = functions.noise_for_image(flux)

#######################################################################################################################
# Create output logs+images and print output filenames and execution time

functions.create_image(sceneim,noise2addim,sceneimFL)

functions.log_output(noise,faintestflux,massstar,logagestar,kzstar,Teffstar,loggstar,loglstar,AVstar,mag,newx+params.xpix/2.0,newy+params.ypix/2.0,readsed,directories.outputstarinfo)

functions.myso_logo('outim')
print(directories.outputim)
print('   ')
functions.myso_logo('outspec')
print(directories.outputspecFL)
print(directories.outputspecL)
print('   ')
stop = timeit.default_timer()

print('Simulation time using 1 core: ', (stop - start)/60. ,'[min]')

#!/usr/bin/env python
from math import *
import random
import numpy as np
import timeit
import params_clean as params
import functions
import directories
import constants

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Start of runtime code

start = timeit.default_timer()
functions.myso_logo('logo')

# Reads in parameters
Columndensities=params.Columndensities
OBtreatment=params.OBtreatment
res=params.res
fovx=params.fovx
fovy=params.fovy
fwhm=params.fwhm
SNR=params.SNR
noise2add=params.noise2add
Rv=params.Rv
distance=params.distance
metallicityZ=params.metallicityZ

spectroscopy = params.spectroscopy     # Spectroscopy output, choose 'yes' or 'no'
lminspec     = params.lminspec    # Minimum wavelength [A] should be set within your filter transparency
lmaxspec     = params.lmaxspec   # Maximum wavelength [A] 
Rspec        = params.Rspec #700      # Spectral resolution (please check your SED library, should not be larger than the resolution of the SEDs)

velocitydis  = params.velocitydis

PSFtype=params.PSFtype
Adaptiveoptics=params.Adaptiveoptics
seeing=params.seeing
SR=params.SR

#######################################################################################################################
# Interstellar dust models

DrainearrLam, DrainearrK=[0.0],[0.0]
if (params.EXTmodel == 'Dmodel'):
    DrainearrLamu,drainealbedo,drainecos,draineC,DrainearrKu,drainecos2=np.loadtxt(directories.Drainemodel,usecols=(0,1,2,3,4,5),unpack=True)
    DrainearrLamu = DrainearrLamu*1.0E+4
    DrainearrKu = DrainearrKu/constants.DraineKappaV
    DrainearrLam, DrainearrK = zip(*sorted(zip(DrainearrLamu,DrainearrKu)))

#######################################################################################################################
# 

lambdaF,weight=np.loadtxt(params.filterfile,unpack=True)

xpix=round(fovx/res)
ypix=round(fovy/res)

print('FoV: ',fovx,'" x ',fovy,'" =',int(xpix),'[pix] x',int(ypix),'[pix]')

sceneim=np.full((int(ypix),int(xpix)),0.0)

wavelengthvega,fluxvega=np.loadtxt(directories.sed_vega,unpack=True)
wavelength_weightvega=np.interp(wavelengthvega,lambdaF,weight)
par2vega=0.0
for i in range(1, (len(wavelengthvega))):
  par2vega += fluxvega[i]*wavelengthvega[i]*wavelength_weightvega[i]*(wavelengthvega[i]-wavelengthvega[i-1])

#######################################################################################################################
# 

Lspecarr=[]
sceneimFL = []
if (spectroscopy == 'yes'):
    functions.myso_logo('spec')
    Lspecarr.append(lminspec)
    while (Lspecarr[-1] <= lmaxspec): Lspecarr.append(Lspecarr[-1]+Lspecarr[-1]/Rspec)
    nspec=len(Lspecarr)
    sceneL=np.zeros(nspec)
    sceneimFL=np.full((nspec,int(ypix),int(xpix)),0.0)
    lun=open(directories.outputspecL,"w")
    lun.write('# Lambda[A] \n')
    for ll in range(nspec): lun.write("%f\n" %(Lspecarr[ll]))  
    lun.close()

#######################################################################################################################
# 

massstar,logagestar,kzstar,rhostar=np.loadtxt(params.filestar,usecols=(6,7,8,9),unpack=True)
if (params.Columndensities == 'sph'): masspar=np.loadtxt(params.filecloud,usecols=(6),unpack=True)

nstar,newx,newy,newz,vxstar,vystar,vzstar,distancestar,newxcloud,newycloud,newzcloud,newhcloud = functions.rotation()

pc2pixstar=206264.806247/distancestar/params.res



#######################################################################################################################
# Reading list of SEDs

teffsed,loggsed,metallicity,lh,vtur,sedname=np.loadtxt(directories.foldersed+'kseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4', 'col5', 'col6'), 'formats':(float,float,float,float,float,'|S60')},usecols=(0,1,2,3,4,5),unpack=True)
nseds=len(teffsed)

if (OBtreatment == 'yes'):
    teffsedOB,loggsedOB,metallicityOB,sednameOB=np.loadtxt(directories.foldersed+'Tseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4'), 'formats':(float,float,float,'|S60')},usecols=(0,1,2,3),unpack=True)
    nsedsOB=len(teffsedOB)
    functions.myso_logo('tlusty')


#wavelength,flux=np.loadtxt(foldersed+sedname[10],comments=['fn:', '#'],unpack=True)

nfovstars=0
for ii in range(nstar):
    if ((abs(newx[ii]) < (xpix/2)-1) and (abs(newy[ii]) < (ypix/2)-1)): # Check if star is within chosen FOV
        nfovstars += 1
    else:
        np.delete(logagestar,[ii])
        np.delete(massstar,[ii])
        np.delete(kzstar,[ii])
        np.delete(newx,[ii])
        np.delete(newy,[ii])
        np.delete(newz,[ii])
        np.delete(pc2pixstar,[ii])

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

iso_load_flag = True
nstars=0
for ii in range(nfovstars):
    if True :#(fovstars_check[ii]): # Check if star is within chosen FOV
        nstars += 1

    # Density of gas clouds
    if (params.Columndensities == 'sph'):
        columncloud[ii]=functions.COLUMN_DENSITY(newx[ii],newy[ii],newz[ii],newxcloud,newycloud,newzcloud,masspar,newhcloud)
        cloudden=columncloud[ii]*((1.989/1.673534)/(((3.0857/pc2pixstar[ii])**2.))) #*1.0e21 ;convert column density unit [Msun/pix^2] --> [10^21 * Mhydrogen/cm^2]
    else: cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857**2.)))#*1.0e21 ;convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]

    #        cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857**2))) #*1.0e21 #convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]
    AVstar[ii]=cloudden/2.21 #e21 ;Guver&Ozel2009: The relation between Optical Extinction and Hydrogen column density
    
    # If lum, teff, and radius are not given, use isochrones to calculate appropriate values
    if not params.LTRprovided:
        if iso_load_flag:
            functions.myso_logo('iso')
            ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso=np.loadtxt(directories.isochrones,usecols=(0,1,2,3,4,5,6),unpack=True)
            iso_load_flag = False

        Teffstar[ii], loggstar[ii], loglstar[ii], flag = functions.get_iso_params(logagestar[ii],massstar[ii],kzstar[ii],ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso)
        if flag: nstars-=1; continue


    deltaT=abs(Teffstar[ii]-teffsed)

    deltagarr=np.full(nseds,99.)


    for jj in range(nseds): 
        if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsed[jj])
#        sedstar[ii]=where(deltagarr eq min(deltagarr))
    sedstar[ii]=functions.FINDCLOSE(min(deltagarr),deltagarr)

    readsed[ii]=sedname[sedstar[ii]].decode('UTF-8')

#        readcol,foldersed+sedname[sedstar[ii]],wavelength,flux,/silent
    wavelength,flux=np.loadtxt(directories.foldersed+sedname[sedstar[ii]].decode('UTF-8'),comments=['fn:', '#'],unpack=True)

    if ((OBtreatment == 'yes') and (Teffstar[ii] >= 15000.)):
        deltaT=abs(Teffstar[ii]-teffsedOB)
        deltagarr=np.full(nsedsOB,99.)
        for jj in range(nsedsOB):  
            if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsedOB[jj])
#            sedstar[ii]=where(deltagarr eq min(deltagarr))
        sedstar[ii]=functions.FINDCLOSE(min(deltagarr),deltagarr)

        readsed[ii]=sednameOB[sedstar[ii]].decode('UTF-8')
#            readcol,foldersed+sednameOB[sedstar[ii]],wavelength,flux,/silent
        wavelength,flux=np.loadtxt(directories.foldersed+sednameOB[sedstar[ii]].decode('UTF-8'),comments=['fn:', '#'],unpack=True)


    bc1=functions.BCcal(wavelength,flux,lambdaF,weight,AVstar[ii],Rv,Teffstar[ii],par2vega,params.EXTmodel,DrainearrLam,DrainearrK)
    mag[ii]=constants.Mbolsun-2.5*loglstar[ii]-(bc1)+5.0*log10(distancestar[ii]/10.0)
    fluxstar[ii]=10.**(mag[ii]/(-2.5))


    if (PSFtype == 'gaussian'):
        airy1=functions.makeGaussian((xpix,ypix),fwhm/res,center=(newx[ii]+xpix/2.,newy[ii]+ypix/2.))

        airy1=airy1/np.sum(airy1) #to be sure about normaized total flux across FOV 

    if (Adaptiveoptics == 'yes'):
        halo=functions.makeGaussian((xpix,ypix),seeing/res,center=(newx[ii]+xpix/2.,newy[ii]+ypix/2.))
        halo=halo/np.sum(halo)
        sceneim += fluxstar[ii]*(SR*airy1+(1.0-SR)*halo)
    if (Adaptiveoptics == 'no'): sceneim += fluxstar[ii]*airy1


    if (spectroscopy == 'yes'):
        for ll in range(2,nspec-3): 
#                linterp,lambda,weight,wavelength,wavelength_weight
            wavelength_weight=np.interp(wavelength,lambdaF,weight)
            bc3=functions.BCcals(wavelength,flux*wavelength_weight,[Lspecarr[ll-1],Lspecarr[ll],Lspecarr[ll+1],Lspecarr[ll+2]],[0.0,1.0,1.0,0.0],AVstar[ii],Rv,Teffstar[ii],params.EXTmodel,DrainearrLam,DrainearrK)
            mag3=constants.Mbolsun-2.5*loglstar[ii]-bc3+5.0*log10(distancestar[ii]/10.0)
            fluxstar3=float(10.**(mag3/(-2.5)))
            if (PSFtype == 'gaussian'):
                airy3=functions.makeGaussian((xpix,ypix),fwhm/res,center=(newx[ii]+xpix/2.,newy[ii]+ypix/2.))

            airy3=airy3/(np.sum(airy3))

            if (velocitydis == 'yes'):
                shiftv=vzstar[ii]*Lspecarr[ll]/3.0E+5 #shift is in A
                lambdashift=Lspecarr[ll]-shiftv        #if vz>0 ==> source comes toward observer ==> lambdashift<lamda0 (blue-shift)
                llchannel=functions.FINDCLOSE(lambdashift,Lspecarr)
                if (Adaptiveoptics == 'yes'): sceneimFL[llchannel,:,:] += fluxstar3*(SR*airy3+ (1.0-SR)*halo)
                if (Adaptiveoptics == 'no' ): sceneimFL[llchannel,:,:] += fluxstar3*airy3

            elif (velocitydis == 'no'):
                if (Adaptiveoptics == 'yes'): sceneimFL[ll,:,:] += fluxstar3*(SR*airy3+ (1.0-SR)*halo)
                if (Adaptiveoptics == 'no' ): sceneimFL[ll,:,:] += fluxstar3*airy3
#Outputs
#   nstars, sceneim, sceneimFL,


print('Number of valid stars in the FoV: ',nstars)


faintestflux, noise2addim, noise = functions.noise_for_image(flux,xpix,ypix)

#######################################################################################################################
# Create output logs+images and print output filenames and execution time

functions.create_image(sceneim,noise2addim,sceneimFL)

functions.log_output(noise,faintestflux,massstar,logagestar,kzstar,Teffstar,loggstar,loglstar,AVstar,mag,newx+xpix/2.0,newy+ypix/2.0,readsed,directories.outputstarinfo)

functions.myso_logo('outim')
print(directories.outputim)
print('   ')
functions.myso_logo('outspec')
print(directories.outputspecFL)
print(directories.outputspecL)
print('   ')
stop = timeit.default_timer()

print('Simulation time using 1 core: ', (stop - start)/60. ,'[min]')

import numpy as np
from math import sqrt, log10, log
import params_clean as params
import directories
import constants
from numba import jit, prange
import scipy

pi = np.pi
parallel = True
nopython = True
#######################################################################################################################
# Functions to be parallelized

# Misc.
@jit(nopython=nopython)
def FINDCLOSE(par1,pararr):
    nclose=np.argsort(np.abs(np.add(-par1,pararr)))
    return nclose[0]

@jit(nopython=nopython,parallel=parallel)
def distNxM(x,y,x0,y0):
    M = np.empty((len(x),len(y)))
    for i in prange(len(x)):
        for j in prange(len(y)):
            M[i][j] = (x[i]-x0)**2+(y[j]-y0)**2
    return M

@jit(nopython=nopython,parallel=parallel)
def div_mat_num(mat,num):
    temp = np.empty((len(mat),len(mat[0])))
    for i in prange(len(mat)):
        for j in prange(len(mat[0])):
            temp[i][j] = mat[i][j]/num
    return temp



# Dust and Gas
@jit(nopython=nopython)
def Kernel_fn(dr, h):
    if (dr/h < 1.0):
        w = (1/(pi*h**3.)) * (1-1.5*(dr/h)**2. + 0.75*(dr/h)**3.)
    elif (dr/h < 2.0):
        w = (1/(pi*h**3.)) * (0.25 * (2.0-(dr/h))**3.)
    else: w = 0.0

    return w
  
@jit(nopython=nopython,parallel=parallel)
def COLUMN_DENSITY(xstar,ystar,zstar,xarr,yarr,zarr,marr,harr):
    nc=len(marr)
    den2dtot=0.0
    for ii in prange(nc):
        if ((sqrt((xstar-xarr[ii])**2+(ystar-yarr[ii])**2) < 2.*harr[ii]) and (zstar < zarr[ii])):
            dr=sqrt((xstar-xarr[ii])**2+(ystar-yarr[ii])**2)
            lc=2*sqrt(4.*harr[ii]*harr[ii]-dr*dr)
            myres=20
            deltalc=lc/myres
            denarr=np.zeros(myres)
            lcarr=np.zeros(myres)
            for jj in range(myres): denarr[jj]= marr[ii]*Kernel_fn(sqrt(dr*dr+(jj*lc/myres-lc/2.)**2.),harr[ii])*deltalc
            den2dtot = den2dtot+ np.sum(denarr)
    return den2dtot

@jit(nopython=nopython)
def extinctions(lam,Av,Rv):
    Alam=0.0  
    xarr=1./lam
    if ((xarr >= 0.3) and (xarr < 1.1)):
      ax1=0.574*(xarr**1.61)
      bx1=-0.527*(xarr**1.61)
      Alam=Av*(ax1+bx1/Rv)
    elif ((xarr >= 1.1) and (xarr < 3.3)):
      y = xarr-1.82
      ax2 = 1.0+0.17699*y-0.50447*(y**2.0)-0.02427*(y**3.0)+0.72085*(y**4.0) +0.01979*(y**5.0)-0.77530*(y**6.0)+0.32999*(y**7.0)
      bx2 = 1.41338*y+2.28305*(y**2.0)+1.07233*(y**3.0)-5.38434*(y**4.0)-0.62251*(y**5.0)+5.30260*(y**6.0) -2.09002*(y**7.0)
      Alam=Av*(ax2+bx2/Rv)
    elif ((xarr >= 3.3) and (xarr < 5.9)):
      ax3 = 1.752 -0.316*xarr-0.104/((xarr-4.67)**2.0 + 0.341)
      bx3 = (-3.090)+1.825*xarr+1.206/((xarr-4.62)**2.0+0.263)
      Alam=Av*(ax3+bx3/Rv)
    elif ((xarr >= 5.9) and (xarr < 8.0)):
      fax4 = -0.04473*(xarr - 5.9)**2.0 - 0.009779*(xarr -5.9)**3.0
      fbx4 = 0.2130*(xarr-5.9)**2.0 + 0.1207*(xarr-5.9)**3.0
      ax4 = 1.752 - 0.316*xarr - (0.104/((xarr - 4.67)**2.0 + 0.341)) +fax4
      bx4 = -3.090 + 1.825*xarr + (1.206/((xarr - 4.62)**2.0 + 0.263)) +fbx4
      Alam=Av*(ax4+bx4/Rv)
    else:  print('Lambda,',lam,', is out of range for Fitzpatrick (1999) model!!!')
    return Alam


# Bolometric Correction
@jit(nopython=nopython,parallel=parallel)
def BCcal(wavelength,flux,lambdaF,weight,AVstar,Rv,Teff,par2vega,EXTmodel,DrainearrLam,DrainearrK): # sed, filter
    #  BC= Mbol_sun -2.5 log10 [4!PI 10pc^2 Sigma Teff^4 / Lsun] + 2.5 .... Girardi + 2002
 
    n_wave=len(wavelength)
    # I put the zero weight at the edges just to solve the problem of edges in linear interpolation
    weight[0]=0.0
    weight[len(weight)-1]=0.0
    wavelength_weight=np.interp(wavelength,lambdaF,weight)
    par1=0.0
    alamarr=np.zeros(len(wavelength)) #fltarr((size(wavelength))[-1])

    if (EXTmodel == 'Dmodel'):
     kappased=np.interp(wavelength,DrainearrLam,DrainearrK)
     alamarr=AVstar*kappased

  
    for i in prange(1, len(wavelength)):
      ltemp=wavelength[i]*1.0e-4 
      if (EXTmodel == 'Fmodel' and wavelength[i] >= 1250. and wavelength[i] <= 33333.):
          alamarr[i]=extinctions(ltemp,AVstar,Rv)
      par1 += wavelength[i]*flux[i]*(10.0**(-0.4*alamarr[i]))*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])
    BCfilter = constants.Mbolsun + constants.bolconstant - 10.*log10(Teff)+2.5*log10(par1/par2vega)
    return BCfilter


# PSF
@jit(nopython=nopython,parallel=parallel)
def gauss(x,y,x0,y0,fwhm):
    gauss = np.empty((len(x),len(y)),dtype=np.float64)

    D = distNxM(x,y,x0,y0)

    for i in prange(len(x)):
        for j in prange(len(y)):

            gauss[i][j] = np.exp(-4*log(2) * (D[i][j] / fwhm**2))
    return gauss

@jit(nopython=nopython)
def makeGaussian(size, fwhm, center):
    """ Make a square gaussian kernel.
    size is the length-array of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size[0], 1)
    y = np.arange(0, size[1], 1)
    if center is None:
        x0 = np.floor_divide(size[0],2)
        y0 = np.floor_divide(size[1],2)
    else:
        x0 = center[0]
        y0 = center[1]

    
    g = gauss(x,y,x0,y0,fwhm)
    return g


#SED
@jit(nopython=nopython,parallel=parallel)
def sed_ind_search(nseds,Teffstar,teffsed,loggstar,loggsed):
    deltaT=np.abs(Teffstar-teffsed)
    deltagarr=np.full(nseds,99.)
    for jj in prange(nseds): 
        if (deltaT[jj] == min(deltaT)): 
            deltagarr[jj] = np.abs(loggstar-loggsed[jj])
    return deltagarr

# Main Loop
@jit(nopython=nopython,parallel=parallel)
def clouds_and_SEDs(nfovstars,sedname,nseds,teffsed,loggsed,sednameOB,nsedsOB,teffsedOB,loggsedOB,newx,newy,newz,newxcloud,newycloud,newzcloud,masspar,newhcloud,rhostar,pc2pixstar,Teffstar,loggstar):

    columncloud=np.empty(nfovstars,dtype=np.float64) #column density of the cloud in front of each star
    AVstar=np.empty(nfovstars,dtype=np.float64)
    readsed=np.empty(nfovstars,dtype='U64')
    sedstar=np.empty(nfovstars,dtype=np.uint64)

    for ii in prange(nfovstars):

        #############################################################
        # Density of gas clouds
        if (params.Columndensities == 'sph'):
            columncloud[ii]=COLUMN_DENSITY(newx[ii],newy[ii],newz[ii],newxcloud,newycloud,newzcloud,masspar,newhcloud)
            cloudden=columncloud[ii]*((1.989/1.673534)/(((3.0857/pc2pixstar[ii])**2.))) #*1.0e21 ;convert column density unit [Msun/pix^2] --> [10^21 * Mhydrogen/cm^2]
        else: cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857**2.)))#*1.0e21 ;convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]
        AVstar[ii]=cloudden/2.21 #e21 ;Guver&Ozel2009: The relation between Optical Extinction and Hydrogen column density

        #############################################################
        # Matching and finding seds

        if ((params.OBtreatment == 'no') or (Teffstar[ii] < 15000.)):

            deltagarr = sed_ind_search(nseds,Teffstar[ii],teffsed,loggstar[ii],loggsed)

            sedstar[ii]=FINDCLOSE(min(deltagarr),deltagarr)
            readsed[ii]=sedname[sedstar[ii]]

        elif ((params.OBtreatment == 'yes') and (Teffstar[ii] >= 15000.)):

            deltagarr = sed_ind_search(nsedsOB,Teffstar[ii],teffsedOB,loggstar[ii],loggsedOB)

            sedstar[ii]=FINDCLOSE(min(deltagarr),deltagarr)

            readsed[ii]=sednameOB[sedstar[ii]]

    return AVstar, readsed

@jit(nopython=nopython,parallel=parallel)
def BC_and_PSF(sceneim,nfovstars,newx,newy,wavelength,flux,indxs,lambdaF,weight,AVstar,Teffstar,loglstar,distancestar,par2vega,DrainearrLam,DrainearrK):

    fluxstar=np.empty(nfovstars,dtype=np.float64)
    mag=np.empty(nfovstars,dtype=np.float64)

    for ii in prange(nfovstars):
        WL = wavelength[indxs[ii]:indxs[ii+1]]
        FL = flux[indxs[ii]:indxs[ii+1]]

        bc1=BCcal(WL,FL,lambdaF,weight,AVstar[ii],params.Rv,Teffstar[ii],par2vega,params.EXTmodel,DrainearrLam,DrainearrK)
        mag[ii]=constants.Mbolsun-2.5*loglstar[ii]-(bc1)+5.0*log10(distancestar[ii]/10.0)
        fluxstar[ii]=10.**(mag[ii]/(-2.5))

        if (params.PSFtype == 'gaussian'):
            airy1 = makeGaussian(np.array([params.xpix,params.ypix]),params.fwhm/params.res,center=[newx[ii]+params.xpix/2.,newy[ii]+params.ypix/2.])
        airy1 = div_mat_num(airy1,np.sum(airy1)) #to be sure about normaized total flux across FOV 

        if (params.Adaptiveoptics == 'yes'):
            halo = makeGaussian(np.array([params.xpix,params.ypix]),params.seeing/params.res,center=[newx[ii]+params.xpix/2.,newy[ii]+params.ypix/2.])
            halo = div_mat_num(halo,np.sum(halo))
            sceneim += fluxstar[ii]*(params.SR*airy1+(1.0-params.SR)*halo)

        if (params.Adaptiveoptics == 'no'):
            sceneim += airy1*fluxstar[ii]

    return sceneim,fluxstar,mag
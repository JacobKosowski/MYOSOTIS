import numpy as np
from math import sqrt, log, log10
from scipy.linalg import expm
from PIL import Image
from astropy.io import fits
import pandas as pd
import constants
import params_clean as params
import directories

pi = np.pi

#######################################################################################################################
# Functions

# Misc.
def myso_logo(wh):
    if (wh == 'logo'):        
        print(' ================================================================================= ')
        print('|   ===============   Make   Your   Own Synthetic   ObservaTIonS  =============   |')
        print('| =============================================================================== |')
        print('||                                                                               ||')
        print('|| |%|       |%| |%|   |%| |%||%||%| |%||%||%| |%||%||%| |%||%||%| |%| |%||%||%| ||')
        print('|| |%||%| |%||%| |%|   |%| |%|   |%| |%|       |%|   |%|    |%|    |%| |%|       ||')
        print('|| |%|  |%|  |%|    |%|    |%|   |%| |%||%||%| |%|   |%|    |%|    |%| |%||%||%| ||')
        print('|| |%|  |%|  |%|    |%|    |%|   |%|       |%| |%|   |%|    |%|    |%|       |%| ||')
        print('|| |%|  |%|  |%|    |%|    |%||%||%| |%||%||%| |%||%||%|    |%|    |%| |%||%||%  ||')
        print('||                                                                               ||')
        print('| =============================================================================== |')
        print('||        MMMMMMMMMMMMMMMMMWN00KX0OOO0KNMMMMMMMMMMMMMMMMMMMM      ||')
        print('||        MMMMMMMMMMMMMMMMW0dloddolcccld0NMMMMMMMMMMMMMMMMMM      ||')
        print('||        MMMMMMMMMMMMMMMNklcloooooollccckWMMMMMMMMMMMMMMMMM      ||')
        print('||        MMMMMMMMMMMMMMMKo::clooooollccclOWMMMMMMMMMMMMMMMM      ||')
        print('||        MMMMMMMMMMMMMMWOl:ccclllllllccccoKMMMMMMMMMMMMMMMM      ||')
        print('||        MMMMMMMMMMMMMMNxccccccclllllcccclkNMMMMMWWNXNNWMMM      ||')
        print('||        MMMMMMMWWWWWWWNxcccc::::ccccc::ccxNMWNKOkxdoodxONM      ||')
        print('||        MMMWX0OkkkxxxkOxcc::::;;::cc:::clONX0kxdxxxdoolcxN      ||')
        print('||        MMXkllloooooooollloc::;;;:::::ldkOOkxkkkkkxxxddodK      ||')
        print('||        MNd::cclloodddoolcoxd:,;:c:;:dOxddddddxdddddddoooO      ||')
        print('||        Xxlccclllloooooll::dOkodkkxoxko:clooolccccccccc:cO      ||')
        print('||        kllllccccc:cccclllldkkxdllooxxollcc::::cccccccclkN      ||')
        print('||        klllcccc:::::::::cokxdc,...llool;,;;::::::::cldKWM      ||')
        print('||        Xklcccc:::;;;;;;;;ldool;...llolddol:,,,;;::lxKWMMM      ||')
        print('||        MWX0xlc::;;;,,,;cdkxdxdolcloxkxlllc;,,:oxk0XWMMMMM      ||')
        print('||        MMMMWNK0Okxol;;::;;:lxolxkOxdkkoc;,,,,c0WMMMMMMMMM      ||')
        print('||        MMMMMMMMMMWNxc:::;:cllccx00klldddolc:::l0WMMMMMMMM      ||')
        print('||        MMMMMMMMMMWOlccccccllcccccllodddddolcccco0WMMMMMMM      ||')
        print('||        MMMMMMMMMM0occ::clccclcccclodddddddoollccoKMMMMMMM      ||')
        print('||        MMMMMMMMMWkccccccccccclolcdkxdxxxdddooolllxNMMMMMM      ||')
        print('||        MMMMMMMMMNxcc::ccccclloxO0XWX0Oxxddddoooolo0WMMMMM      ||')
        print('||        MMMMMMMMMNx::::ccclllokXWMMMMMWXKOxddooolllkWMMMMM      ||')
        print('||        MMMMMMMMMMXxl::cllclxKWMMMMMMMMMMWNK0OxdodONMMMMMM      ||')
        print('||        MMMMMMMMMMMWNOlccokKWMMMMMMMMMMMMMMMMMWNNWMMMMMMMM      ||')
        print('||        MMMMMMMMMMMMMW0xkKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM      ||')
        print(' =======================     FORGET ME NOT     ===================  ')

    if (wh == 'star'): 
        print(' =============================================================== ')
        print('|               READING YOUR STELLAR SOURCES FILE               |')
        print(' =============================================================== ')
        
    if (wh == 'cloud'):
        print(' =============================================================== ')
        print('|                   READING YOUR CLOUD  FILE                    |')
        print(' =============================================================== ')

    if (wh  == 'filter'):
        print(' =============================================================== ')
        print('|                         READING FILTER                        |')
        print(' =============================================================== ')
 
    if (wh  == 'iso'):
        print(' =============================================================== ')
        print('|                  READING EVOLUTIONARY MODELS                  |')
        print(' =============================================================== ')

    if (wh  == 'spec'):
        print(' =============================================================== ')
        print('|                  SPECTRA CHANNEL is CREATED                   |')
        print(' =============================================================== ')
 
    if (wh  == 'tlusty'):
        print(' =============================================================== ')
        print('|                USING TLUSTY for MASSIVE STARs                 |')
        print(' =============================================================== ')

    if (wh  == 'ao'):
        print(' =============================================================== ')
        print('|              ADOPTIVE OPTICS is being APPLIED                 |')
        print(' =============================================================== ')
 
    if (wh  == 'outim'):
        print(' =============================================================== ')
        print('|            OUTPUT SYNTHETIC IMAGE is WRITTEN in :             |')
        print(' =============================================================== ')
 
    if (wh  == 'outspec'):
        print(' =============================================================== ')
        print('|          OUTPUT SYNTHETIC SPECTRA CUBE is WRITTEN in :        |')
        print(' =============================================================== ')

def FINDCLOSE(par1,pararr):
    nclose=np.argsort(abs(np.add(-par1,pararr)))
    return nclose[0]

def pad(seq, target_length, padding=np.nan):
    seq = list(seq)
    length = len(seq)               
    seq.extend([padding] * (target_length - length))
    return np.array(seq)


# Rotation
def rot_euler(v, xyz):
    ''' Rotate vector v (or array of vectors) by the euler angles xyz '''
    # https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    for theta, axis in zip(xyz, np.eye(3)):
        e = expm(np.cross(np.eye(3), axis*-theta))
        v = np.dot(e, np.array(v))
    return v.T

def rotation():
    """ Rotates xyz and corresponding velocites by angles alphai, bettai, and gammai respectively.
    -Provides different lines-of-sight of the same object
    -Also calculates rotation for clouds if applicable
    -Futher returns distance to the stars
    """

    xstar,ystar,zstar,vxstar,vystar,vzstar=np.loadtxt(params.filestar,usecols=(0,1,2,3,4,5),unpack=True)

    nstar=len(xstar)
    positionvector=[xstar,ystar,zstar]
    velocityvector=[vxstar,vystar,vzstar]

    zeinab=rot_euler(positionvector,np.multiply([params.alphai,params.bettai,params.gammai],pi/180.))
    xstar=zeinab[0:nstar,0]
    ystar=zeinab[0:nstar,1]
    zstar=zeinab[0:nstar,2]

    zeinabv=rot_euler(velocityvector,np.multiply([params.alphai,params.bettai,params.gammai],pi/180.))
    vxstar=zeinabv[0:nstar,0]
    vystar=zeinabv[0:nstar,1]
    vzstar=zeinabv[0:nstar,2]

    distancestar=np.add(params.distance,-zstar)

    pc2pixstar=206264.806247/distancestar/params.res
    newx=xstar*pc2pixstar #convert x[pc] into pixel position
    newy=ystar*pc2pixstar
    newz=zstar*pc2pixstar
    print('READ STARs: ',nstar)

    #######################################################################################################################
    # 

    newxcloud,newycloud,newzcloud,newhcloud=0,0,0,0

    if (params.Columndensities == 'sph'):
        myso_logo('cloud')
        xcloud,ycloud,zcloud,vxcloud,vycloud,vzcloud,hpar=np.loadtxt(params.filecloud,usecols=(0,1,2,3,4,5,7),unpack=True)
        ncloud=len(xcloud)
        
        positioncvector=[xcloud,ycloud,zcloud]
        zeinabc=rot_euler(positioncvector,np.multiply([params.alphai,params.bettai,params.gammai],pi/180.))
        xcloud=zeinabc[0:ncloud,0]
        ycloud=zeinabc[0:ncloud,1]
        zcloud=zeinabc[0:ncloud,2]

        distancecloud=np.add(params.distance,-zcloud)
        pc2pixcloud=206264.806247/distancecloud/params.res
        newxcloud=xcloud*pc2pixcloud #convert x[pc] into pixel position
        newycloud=ycloud*pc2pixcloud
        newzcloud=zcloud*pc2pixcloud
        newhcloud=np.multiply(hpar,pc2pixcloud)
    
    return nstar,newx,newy,newz,vxstar,vystar,vzstar,distancestar,newxcloud,newycloud,newzcloud,newhcloud


# Dust and Gas
def draine_dust_model():
    DrainearrLam, DrainearrK=[0.0],[0.0]
    if (params.EXTmodel == 'Dmodel'):
        DrainearrLamu,drainealbedo,drainecos,draineC,DrainearrKu,drainecos2=np.loadtxt(directories.Drainemodel,usecols=(0,1,2,3,4,5),unpack=True)
        DrainearrLamu = DrainearrLamu*1.0E+4
        DrainearrKu = DrainearrKu/constants.DraineKappaV
        DrainearrLam, DrainearrK = zip(*sorted(zip(DrainearrLamu,DrainearrKu)))
    return np.array(DrainearrLam), np.array(DrainearrK)

def Kernel_fn(dr, h):
    pi=np.pi
    if (dr/h < 1.0):
        w = (1/(pi*h**3.)) * (1-1.5*(dr/h)**2. + 0.75*(dr/h)**3.)
    elif (dr/h < 2.0):
        w = (1/(pi*h**3.)) * (0.25 * (2.0-(dr/h))**3.)
    else: w = 0.0

    return w
  
def COLUMN_DENSITY(xstar,ystar,zstar,xarr,yarr,zarr,marr,harr):
    nc=len(marr)
    den2dtot=0.0
    for ii in range(nc):
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

# Filters
def filters():
    lambdaF,weight=np.loadtxt(params.filterfile,unpack=True)
    wavelengthvega,fluxvega=np.loadtxt(directories.sed_vega,unpack=True)

    wavelength_weightvega=np.interp(wavelengthvega,lambdaF,weight)
    par2vega=0.0
    for i in range(1, (len(wavelengthvega))):
      par2vega += fluxvega[i]*wavelengthvega[i]*wavelength_weightvega[i]*(wavelengthvega[i]-wavelengthvega[i-1])

    return lambdaF,weight,par2vega

# Spectroscopy
def spec_wave_out():
    Lspecarr=[]
    sceneimFL = []
    myso_logo('spec')
    Lspecarr.append(params.lminspec)

    while (Lspecarr[-1] <= params.lmaxspec): Lspecarr.append(Lspecarr[-1]+Lspecarr[-1]/params.Rspec)

    nspec=len(Lspecarr)
    sceneL=np.zeros(nspec)

    sceneimFL=np.full((nspec,int(params.ypix),int(params.xpix)),0.0)
    lun=open(directories.outputspecL,"w")
    lun.write('# Lambda[A] \n')
    for ll in range(nspec): lun.write("%f\n" %(Lspecarr[ll]))  
    lun.close()
    return sceneimFL,Lspecarr,nspec

def spectroscopy(sceneimFL,nspec,wavelength,lambdaF,weight,flux,Lspecarr,AVstar,Teffstar,DrainearrLam,DrainearrK,loglstar,distancestar,newx,newy,vzstar):
        for ll in range(2,nspec-3): 
            #linterp,lambda,weight,wavelength,wavelength_weight
            wavelength_weight=np.interp(wavelength,lambdaF,weight)
            # print(np.max(wavelength_weight))
            # break

            bc3=BCcals(wavelength,flux*wavelength_weight,[Lspecarr[ll-1],Lspecarr[ll],Lspecarr[ll+1],Lspecarr[ll+2]],[0.0,1.0,1.0,0.0],AVstar,params.Rv,Teffstar,params.EXTmodel,DrainearrLam,DrainearrK)
            mag3=constants.Mbolsun-2.5*loglstar-bc3+5.0*nlog10(distancestar/10.0)
            fluxstar3=float(10.**(mag3/(-2.5)))
            if (params.PSFtype == 'gaussian'):
                airy3=makeGaussian((params.xpix,params.ypix),params.fwhm/params.res,center=(newx+params.xpix/2.,newy+params.ypix/2.))

            airy3=airy3/(np.sum(airy3))

            if (params.velocitydis == 'yes'):
                shiftv=vzstar*Lspecarr[ll]/3.0E+5 #shift is in A
                lambdashift=Lspecarr[ll]-shiftv        #if vz>0 ==> source comes toward observer ==> lambdashift<lamda0 (blue-shift)
                llchannel=FINDCLOSE(lambdashift,Lspecarr)
                if (params.Adaptiveoptics == 'yes'): sceneimFL[llchannel,:,:] += fluxstar3*(params.SR*airy3+ (1.0-params.SR)*halo)
                if (params.Adaptiveoptics == 'no' ): sceneimFL[llchannel,:,:] += fluxstar3*airy3

            elif (params.velocitydis == 'no'):
                if (params.Adaptiveoptics == 'yes'): sceneimFL[ll,:,:] += fluxstar3*(params.SR*airy3+ (1.0-params.SR)*halo)
                if (params.Adaptiveoptics == 'no' ): sceneimFL[ll,:,:] += fluxstar3*airy3
        return sceneimFL 


# Bolometric Correction
def BCcal(wavelength,flux,lambdaF,weight,AVstar,Rv,Teff,par2vega,EXTmodel,DrainearrLam,DrainearrK): # sed, filter
    #  BC= Mbol_sun -2.5 log10 [4!PI 10pc^2 Sigma Teff^4 / Lsun] + 2.5 .... Girardi + 2002
 
    n_wave=len(wavelength)
    # I put the zero weight at the edges just to solve the problem of edges in linear interpolation
    weight[0]=0.0
    weight[len(weight)-1]=0.0
    wavelength_weight=np.interp(wavelength,lambdaF,weight)
    par1=0.0
    par2=0.0
    alamarr=np.zeros(len(wavelength)) #fltarr((size(wavelength))[-1])

    if (EXTmodel == 'Dmodel'):
     kappased=np.interp(wavelength,DrainearrLam,DrainearrK)
     alamarr=AVstar*kappased

  
    for i in range(1, len(wavelength)):
      ltemp=wavelength[i]*1.0e-4 
      if (EXTmodel == 'Fmodel' and wavelength[i] >= 1250. and wavelength[i] <= 33333.):
          alamarr[i]=extinctions(ltemp,AVstar,Rv)
      par1 += wavelength[i]*flux[i]*(10.0**(-0.4*alamarr[i]))*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])
    #     par2 += wavelength[i]*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])  ;!!!! par2 will be calculating from Vega flux
    BCfilter = constants.Mbolsun + constants.bolconstant - 10.*log10(Teff)+2.5*log10(par1/par2vega)
    return BCfilter

def BCcals(wavelength,flux,lambdaF,weight,AVstar,Rv,Teff,EXTmodel,DrainearrLam,DrainearrK): # sed, filter
    #  BC= Mbol_sun -2.5 log10 [4!PI 10pc^2 Sigma Teff^4 / Lsun] + 2.5 .... Girardi + 2002
 
    n_wave=len(wavelength)
    # I put the zero weight at the edges just to solve the problem of edges in linear interpolation
    weight[0]=0.0
    weight[len(weight)-1]=0.0
    wavelength_weight=np.interp(wavelength,lambdaF,weight)
    par1=0.0
    par2=0.0
    alamarr=np.zeros(len(wavelength)) #fltarr((size(wavelength))[-1])

    if (EXTmodel == 'Dmodel'):
     kappased=np.interp(wavelength,DrainearrLam,DrainearrK)
     alamarr=AVstar*kappased

  
    for i in range(1, len(wavelength)):
        ltemp=wavelength[i]*1.0e-4 
        if (EXTmodel == 'Fmodel' and wavelength[i] >= 1250. and wavelength[i] <= 33333.):
            alamarr[i]=extinctions(ltemp,AVstar,Rv)

        if par1==0: par1=1e-9
        if par2==0: par2=1e-9
        par1 += wavelength[i]*flux[i]*(10.0**(-0.4*alamarr[i]))*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])
        par2 += wavelength[i]*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])  #!!!! par2 will be calculating from Vega flux
    BCfilter = constants.Mbolsun + constants.bolconstant - 10.*log10(Teff)+2.5*log10(par1/par2)
    return BCfilter


# PSF
def makeGaussian(size, fwhm, center):
    """ Make a square gaussian kernel.
    size is the length-array of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    x = np.arange(0, size[0], 1, float)
    y = np.arange(0, size[1], 1, float)
    y = y[:,np.newaxis]
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    return np.exp(-4*log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


# Isochrones
def get_iso_params(logagestar,massstar,kzstar,ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso):
    # Reading the isochrones
    # Finding the best match isochrones based on logage and mass
    # Returns teff, logg, logl of best match


    flag = False #Flag for if the mass or age of the star differs by more than 10% from the iso

    teffiso=10.**logteff
    niso=len(miniiso)

    nage=FINDCLOSE(logagestar,logageiso)
    selectedage=logageiso[nage]

    if abs(logageiso[nage]-logagestar)/logagestar > 0.1:
        flag=True
    # logagestar[ii]=logageiso[nage]


    nmetalicity=FINDCLOSE(kzstar,ziso)
    selectedz=ziso[nmetalicity]

    marrtemp=np.full(niso,999.99)

    for kk in range(niso):  
        if ((ziso[kk] == selectedz) and (logageiso[kk] == selectedage)):  marrtemp[kk]=mactiso[kk]
    ns=FINDCLOSE(massstar,marrtemp)

    if abs(marrtemp[ns]-massstar)/massstar > 0.1:
        flag=True

    return teffiso[ns], loggiso[ns], logliso[ns], flag


# Output
def noise_for_image(fluxstar):
    fwhm,res,SNR = params.fwhm,params.res,params.SNR

    faintestflux=min(i for i in fluxstar if i > 0)
    if (params.SNR != 0.0): 
        noise=faintestflux*4.*log(2.0)/(pi*(fwhm/res)*(fwhm/res))/SNR 
    else: 
        noise=params.noise2add

    noise2addim=noise*np.random.rand(int(params.ypix),int(params.xpix))

    return faintestflux, noise2addim, noise


def log_output(nfovstars,noise,faintestflux,massstar,logagestar,kzstar,Teffstar,loggstar,loglstar,AVstar,mag,xpos,ypos,readsed,outputstarinfo):
    lunstarinfo=open(outputstarinfo,"w")

    lunstarinfo.write("#     mass[mo] ,   logage[yr]  ,     Z      , log[Teff[k]] ,    logg    ,   logL/Lo    ,  Av_star   ,      mag     ,    Xpix     ,   Ypix       ,       assignedSED \n")
    lunstarinfo.write("#       1              2             3            4              5              6           7               8            9            10                   11 \n")

    for ii in range(len(massstar)):
        lunstarinfo.write("%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %60s \n" %(massstar[ii], logagestar[ii],kzstar[ii],log10(Teffstar[ii]),loggstar[ii],loglstar[ii],AVstar[ii],mag[ii],xpos[ii],ypos[ii],readsed[ii]))

    lunstarinfo.write('#nstars: %f \n' %(len(massstar)))

    lunstarinfo.write('#faintestflux: %f \n' %(faintestflux)) #calc
    lunstarinfo.write('#noise:  %f \n' %(noise)) #calc

    lunstarinfo.write('#SNR:  %f \n' %(params.SNR)) #params
    lunstarinfo.write('#Filter:  %s \n' %(params.filterfile)) #params
    lunstarinfo.write('#distance:  %f \n' %(params.distance)) #params
    lunstarinfo.write('#res:  %f \n' %(params.res)) #params
    lunstarinfo.write('#FoV:  %d %d \n' %(params.fovx,params.fovy)) #params
    lunstarinfo.write('#fwhm:  %f \n' %(params.fwhm)) #params
    lunstarinfo.write('#OBtreatment:  %s \n' %(params.OBtreatment)) #params
    lunstarinfo.write('#Adaptiveoptics:  %s \n' %(params.Adaptiveoptics)) #params
    lunstarinfo.write('#SR:  %f \n' %(params.SR)) #params
    lunstarinfo.write('#Seeing:  %f \n' %(params.seeing)) #params
    lunstarinfo.write('#Spectroscopy:  %s \n' %(params.spectroscopy)) #params
    if params.spectroscopy=='yes':
        lunstarinfo.write('#Spectroscopic Resolution:  %f \n' %(params.Rspec)) #params
        lunstarinfo.write('#Lambda_min:  %f \n' %(params.lminspec)) #params
        lunstarinfo.write('#Lambda_max:  %f \n' %(params.lmaxspec)) #params
    lunstarinfo.write('#filestar:  %s \n' %(params.filestar)) #params
    lunstarinfo.write('#Columndensities:  %s \n' %(params.Columndensities)) #params
    lunstarinfo.write('#filecloud:  %s \n' %(params.filecloud)) #params
    lunstarinfo.write('#EXTmodel:  %s \n' %(params.EXTmodel)) #params
    lunstarinfo.write('#Rv:  %f \n' %(params.Rv)) #params
    lunstarinfo.close()

def create_image(sceneim,noise2addim,sceneimFL):

    if params.filetype=='fits':
        hdu = fits.PrimaryHDU(sceneim)
        hdu.writeto(directories.outputim+'.fits')

        if params.getnoise:
            hdunoise = fits.PrimaryHDU(noise2addim)
            hdunoise.writeto(directories.outputimnoise+'.fits')

    if params.filetype=='png':

        hdu = fits.PrimaryHDU(sceneim)
        hdu.writeto(directories.outputim+'.fits')

        if params.getnoise:
            hdunoise = fits.PrimaryHDU(noise2addim)
            hdunoise.writeto(directories.outputimnoise+'.fits')

        data = fits.open(directories.outputim+'.fits')[0].data.T

        if params.stretched:
            dd = 2*np.arctan(data*params.stretchfactor*10**params.stretchpower)/pi
            plt.imsave(directories.outputim+'.png',dd,cmap=params.colormap)
        else:
            plt.imsave(directories.outputim+'.png',data,cmap=params.colormap)

    else:
        print("!! No image output created. Please specify either 'fits' or 'png'")

    if (params.spectroscopy == 'yes'):
        hdu = fits.PrimaryHDU(sceneimFL)
        hdu.writeto(directories.outputspecFL)


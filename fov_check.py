import params_clean as params
import numpy as np
import functions

def rotation_abridged():
    """ Rotates xyz and corresponding velocites by angles alphai, bettai, and gammai respectively.
    -Provides different lines-of-sight of the same object
    -Also calculates rotation for clouds if applicable
    -Futher returns distance to the stars
    """

    xstar,ystar,zstar=np.loadtxt(params.filestar,usecols=(0,1,2),unpack=True)

    nstar=len(xstar)
    positionvector=[xstar,ystar,zstar]

    zeinab=functions.rot_euler(positionvector,np.multiply([params.alphai,params.bettai,params.gammai],np.pi/180.))
    xstar=zeinab[0:nstar,0]
    ystar=zeinab[0:nstar,1]
    zstar=zeinab[0:nstar,2]

    distancestar=np.add(params.distance,-zstar)

    pc2pixstar=206264.806247/distancestar/params.res
    newx=xstar*pc2pixstar #convert x[pc] into pixel position
    newy=ystar*pc2pixstar
    newz=zstar*pc2pixstar

    print('READ STARs: ',nstar)
    
    return nstar,newx,newy,newz


print('FoV: ',params.fovx,'" x ',params.fovy,'" ::',int(params.xpix),'[pix] x',int(params.ypix),'[pix]')

nstar,newx_,newy_,newz_ = rotation_abridged()


nfovstars=0
jj = 0
for ii in range(nstar):
    if ((abs(newx_[ii]) < (params.xpix/2)-1) and (abs(newy_[ii]) < (params.ypix/2)-1)): # Check if star is within chosen FOV
        nfovstars += 1
        jj+=1
    else:
        jj=jj

if nfovstars==0:
    raise Exception("No stars in FOV")
else:
    print('Number of valid stars in the FoV: ',nfovstars)
    print('Number of stars outside of FOV: ',nstar-nfovstars)
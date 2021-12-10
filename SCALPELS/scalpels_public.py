# -*- coding: utf-8 -*-
#
#   scalpels_public.py - Tools for decomposing radial-velocity measurements 
#   into shape-driven and shift-driven components.
#
#   Copyright (C) 2021  Prof Andrew Collier Cameron, University of St Andrews
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#


import numpy as np
import scipy
import scipy.optimize as opt
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import glob
import pickle


def acf2d(ccf2d):
    # X dimension of data
    nroll = np.shape(ccf2d)[1]
    half = np.round(nroll/2,0)
    # Normalise (divide by np.median of each row)
    normval = np.percentile(ccf2d,75,axis=1)
    ccfNorm = np.transpose(np.transpose(ccf2d) / normval)
    tr = []
    # Simple autocorrelation: circular shift and co-multiply with unshifted profile.
    denom =np.sum(ccfNorm*ccfNorm,axis=1)
    for iroll in np.arange(nroll):
        shifted = np.roll(ccfNorm,(iroll-half).astype(int),axis=1)
        tr.append(np.sum(ccfNorm*shifted,axis=1)/denom)
    
    # Normalise, transpose and return as numpy array
    tr = tr/np.mean(tr,axis=0)
    acf = np.array(tr).T    
    return acf

def median_abs_dev(array2d):
    dev = array2d - np.median(array2d,axis=0)
    absdev = np.abs(dev)
    mad = np.median(absdev,axis=0)
    return(mad)

def goodflag(array2d,madthresh):
    dev = array2d - np.median(array2d,axis=0)
    absdev = np.abs(dev)
    mad = np.median(absdev,axis=0)
    good2d = absdev < madthresh*mad
    good = np.prod(good2d,axis=1)
    return(good.astype('bool'))
    
def rearrange_uloocv(uloocv,alploo,rvorth,rverr,jit=0.0):
    # Initialise arrays and lists
    karray = np.arange(len(uloocv.T))
    good = np.ones_like(karray).astype('bool')
    res = rvorth
    invar = 1/(rverr**2+jit**2)
    chs = np.sum(res**2*invar)
    chslist = []
    kdxlist=[]

    # Loop over columns of U.
    for j in range(len(karray)):
        newchs = []
        idx=[]
        rvklist=[]
        klist=[]
        # Compute next chs obtained for  each of the as-yet unused columns of U
        #print(karray[good],np.shape(alploo),np.shape(uloocv))
        for k in karray[good]:
            rvk = alploo[:,k]*uloocv[:,k]
            klist.append(k)
            rvklist.append(rvk)
            newchs.append(np.sum((res-rvk)**2*invar))  
        # Pick the one that gives the best improvement.
        kdx = np.argmin(newchs) # This column gives best improvement
        kdxlist.append(klist[kdx])
        res = res - rvklist[kdx]
        chs = newchs[kdx]
        chslist.append(chs)
        # Mask out the column we just accumulated into chisq.
        good[klist[kdx]] = False
        
    # Return the index array
    return np.array(kdxlist).astype('int'), np.array(chslist)

def loocv(acfarray,rvorth):
    nobs,nvel = np.shape(acfarray)
    k = np.min((nvel,nobs-1))
    #k = np.min((nvel,nobs))
    uloocv = []
    alploo = []
    acfmodel = []
    u,s,v = np.linalg.svd(acfarray,full_matrices=False)
    obslist = np.arange(nobs)
    for j in obslist:
        mask = obslist != j
        acfloo = acfarray[mask]
        uloo,sloo,vloo = np.linalg.svd(acfloo,full_matrices=False)
        sv = np.dot(np.diag(sloo[:k]),vloo[:k])
        vsv = np.dot(sv,vloo[:k].T)
        asv = np.dot(acfarray[j],vloo[:k].T)
        uj = np.linalg.solve(vsv,asv)
        uloocv.append(uj[:k])
        acfj = np.dot(uj[:k],sv[:k])
        acfmodel.append(acfj)
        rvloo = rvorth[mask]
        alpha = np.dot(rvloo,uloo[:,:k])
        alploo.append(alpha[:k])
        
    # Note that uloocv has one column fewer than u if nobs<ncol
    return np.array(uloocv),np.array(alploo),np.array(acfmodel)

def vscalpels_loocv(ccfarray,rvarray,rverr,npc=None,sort=True,jit=0.0):
    
    # Set maximum number of principal components to use if npc=None.
    if npc == None:
        nobs,nvel = np.shape(ccfarray)
        kmax = np.min((nvel,nobs-1))
        #kmax = np.min((nvel,nobs))
        npc = kmax
    else:
        kmax = npc
    
    
    # Compute ACF
    acfarray = acf2d(ccfarray)
    # Orthogonalise RVs
    rvorth = rvarray - np.mean(rvarray)
    
    # Compute leave-one-out cross-validation estimate of uloocv
    uloocv,alploo,acfmodel = loocv(acfarray,rvorth)
    
    # Use leave-one-out cross-validation to determine how 
    # many columns of U can be used without overfitting.
        
    if sort:
        # Sort columns of U in fastest-descending chisq order.
        idx, chs = rearrange_uloocv(uloocv[:,:npc],alploo[:,:npc],rvorth,rverr,jit)
        alploo = alploo.T[idx].T
        uloocv = uloocv.T[idx].T
        # Define kmax to be argmin(chs), where chisq of the cleaned RV is minimised.
        kmax = np.argmin(chs) # NB if argmin(chs)=1, there are 2 useful PCs.
        npc = kmax+1
         
    # Compute matrix of velocity projections on to LOOCV basis vectors
    print('Using ',npc,' principal components (kmax = %0.0i).'%kmax)
    vresp = alploo.T[:npc].T*uloocv.T[:npc].T
    
    rvshape = np.sum(vresp.T[:npc],axis=0)    
    rvclean = rvorth - rvshape
                
    return rvorth,np.array(rvshape),np.array(rvclean),np.array(uloocv),np.array(alploo)

def vscalpels_recover_loocv(bjd,ccfarray,rvarray,rverr,period,npc=None,sort=True,jit=0.0):
    
    # Set maximum number of principal components to use if npc=None.
    if npc == None:
        nobs,nvel = np.shape(ccfarray)
        kmax = np.min((nvel,nobs-1))
        #kmax = np.min((nvel,nobs))
        npc = kmax
    else:
        kmax = npc
    
    # Compute ACF
    acfarray = acf2d(ccfarray)
    # Orthogonalise RVs
    rvorth = rvarray - np.mean(rvarray)
    
    # Compute leave-one-out cross-validation estimate of uloocv
    uloocv,alploo,acfmodel = loocv(acfarray,rvorth)
    
    # Use leave-one-out cross-validation to determine how 
    # many columns of U can be used without overfitting.
        
    if sort:
        # Sort columns of U in fastest-descending chisq order.
        idx, chs = rearrange_uloocv(uloocv[:,:npc],alploo[:,:npc],rvorth,rverr,jit)
        alploo = alploo.T[idx].T
        uloocv = uloocv.T[idx].T
        # Define kmax to be argmin(chs), where chisq of the cleaned RV is minimised.
        if npc == 0:
            kmax = 0
        else:
            kmax = np.argmin(chs) # NB if argmin(chs)=1, there are 2 useful PCs.
            npc = kmax+1
    
    # Compute matrix of velocity projections on to LOOCV basis vectors
    #vresp = alploo*uloocv
    print('Using ',npc,' principal components (kmax = %0.0i).'%npc)
    
    # Define array F of orbit basis functions 
    fftrn = []
    nplanets = len(period)
    omega = 2.*np.pi / period
    nrow,ncol = np.shape(ccfarray)                
    for planet in range(nplanets):
        # Compute orbital phases, sines and cosines
    
        phi = omega[planet] * (bjd - np.mean(bjd))
        cphi = np.cos(phi)
        sphi = np.sin(phi)
        #print (np.shape(cphi))
        fftrn.append(cphi)
        fftrn.append(sphi)
    
    fftrn = np.asarray(fftrn)
    ffuncs = fftrn.T
    
    # Define pperp projection operator
    pperp = np.identity(nrow)-np.dot(uloocv[:,:npc],uloocv[:,:npc].T)    

    # Project ffuncs and rvorth into the orthogonal complement of the shape subspace
    fperp = np.dot(pperp,ffuncs)
    vperp = np.dot(pperp,rvorth)
    
    # Define inverse variance matrix for independent RV variances
    covar = np.diag(rverr*rverr)
    invar = np.linalg.pinv(covar)
    
    # Set up least-squares problem
    amat = np.dot(np.dot(fperp.T,invar),fperp)
    bvec = np.dot(np.dot(fperp.T,invar),vperp)
    
    # Solve for model parameter vector theta and its uncertainties
    theta = np.linalg.lstsq(amat,bvec,rcond=None)[0]
    theterr = np.sqrt(np.diag(np.linalg.inv(amat)))
    
    # Compute orbital velocity amplitudes and uncertainties.
    icol = 0
    amp = []
    amperr = []
    for planet in range(nplanets):
        Kx, dKx = theta[icol], theterr[icol]
        Ky, dKy = theta[icol+1], theterr[icol+1]
        K = np.sqrt(Kx*Kx+Ky*Ky)
        dK = np.sqrt((Kx*dKx)**2 + (Ky*dKy)**2)/K        
        amp.append(K)
        amperr.append(dK)
        icol+=1
        icol+=1    

    # Separate out the different RV components
    rvorbit = np.dot(fftrn.T,theta)
    rvclean = vperp
    rvshape = rvorth - rvclean
    rvresid = vperp - np.dot(fperp,theta) # = rvorth - rvshape - rvorbit
    
    return rvorth,np.array(rvshape),np.array(rvclean),np.array(rvorbit),np.array(rvresid),period,np.array(amp),np.array(amperr),np.array(uloocv),np.array(alploo)

def mask_outliers(ccfNorm,rvcList,rveList,uclip_threshold=7):
    acfarray = np.array(acf2d(ccfNorm))
    u,s,v = np.linalg.svd(acfarray,full_matrices=False)
    rvorth,rvshape,rvclean,uloocv,alploo = vscalpels_loocv(ccfNorm,rvcList,rveList,npc=None,sort=False)
    nobs,nvel = np.shape(acfarray)
    kmax = np.min((nvel,nobs-1))
    #kmax = np.min((nvel,nobs))
        
    colmed_u = np.median(u[:,:kmax],axis=0)
    absdev_u = np.abs(u[:,:kmax]-colmed_u)
    colmad_u = np.median(absdev_u,axis=0)
    
    badfrac=[]
    
    for badthresh in range(20):
        badmask = absdev_u>badthresh*colmad_u
        goodmask = 1-badmask
        rowmask = np.prod(goodmask,axis=1).astype('bool')
        badfrac.append(np.sum(1-rowmask)/np.size(rowmask))
    
    badfrac = np.array(badfrac)
    
    # Default threshold for clipping outliers is 7*MAD(U-Uhat)_k
    goodmask = 1-(absdev_u>uclip_threshold*colmad_u)
    rowmask = np.prod(goodmask,axis=1).astype('bool')  
    
    return(rowmask,badfrac)

def quality_control(ccfNorm,rvcList,rveList,madfac=3,uclip_threshold=7):
    # Start by weeding out bad rows
    rowmask,badfrac = mask_outliers(ccfNorm,rvcList,rveList,uclip_threshold=7)
    
    #Perform SVD on masked dataset
    acfarray = np.array(acf2d(ccfNorm[rowmask]))
    u,s,v = np.linalg.svd(acfarray,full_matrices=False)
    # Perform LOOCV on masked dataset
    rvorth,rvshape,rvclean,uloocv,alploo = vscalpels_loocv(ccfNorm[rowmask],rvcList[rowmask],rveList[rowmask],npc=None,sort=False)
    nobs,nvel = np.shape(ccfNorm[rowmask])
    kmax = np.min((nvel,nobs-1))
    #kmax = np.min((nvel,nobs))
    
    # Compute ratio MAD(U-Uhat)/MAD(Uhat) 
    h = uloocv      # 'hat'
    d = u[:,:kmax]-h         # 'difference'
    colmed_d = np.median(d,axis=0)
    absdev_d = np.abs(d-colmed_d)
    colmad_d = np.median(absdev_d,axis=0)
    colmed_h = np.median(h,axis=0)
    absdev_h = np.abs(h-colmed_h)
    colmad_h = np.median(absdev_h,axis=0)
    madratio = colmad_d/colmad_h

    # Determine the size of the useful U basis found by LOOCV.
    # The ratio MAD(U-Uhat)/MAD(Uhat) << 1 for good Uhat vectors.
    
    for kk in range(0,len(madratio)-1):
        if madratio[kk] > .25:
            if madratio[kk] > 1-madfac*median_abs_dev(madratio[kk+1:]):
                kopt = kk
                break 

    return kopt, madratio, badfrac, rowmask

def sinefit(tp, vp, ep, P):
    # Simple sine fitter. ep needs to be an array of same length as tp, vp
    from scipy.optimize import leastsq
    invar = 1/ep/ep
    om = 2*np.pi/P
    Kest = np.std(vp)
    # Orthogonalise
    that = np.sum(tp*invar)/np.sum(invar)
    t = tp - that
    vhat = np.sum(vp*invar)/np.sum(invar)
    v = vp - vhat
    cwt = np.cos(om*t)
    swt = np.sin(om*t)
    optimize_func = lambda x: (x[0] + x[1]*cwt + x[2]*swt-v)/ep
    dv0, C, S = leastsq(optimize_func,[0,Kest,Kest])[0]
    K = np.sqrt(S*S+C*C)
    dt = np.arctan2(S,C)/om 
    T0 = that + dt + P/4
    #print('dbg>',T0)
    pha = np.mod(tp-T0,P)/P 
    m = -np.sin(2*np.pi*pha)
    Kerr = np.sqrt(1/np.sum(m*m*invar))
    # r = (dv0 + C*np.cos(om*t) + S*np.sin(om*t)-v)/ep
    r = (dv0 + K*m - v)/ep
    chs = np.dot(r,r)
    
    return(T0,K,Kerr,chs)

def parabfit(x,y):
    r = (y[2]-y[1])/(x[2]-x[1])
    q = (y[2]-y[0])/(x[2]-x[0])
    aa = (r - q)/(x[1] - x[0])
    bb = q - aa * (x[0] + x[2])
    xmin = -bb/2./aa
    cc = y[1] - (aa*x[1] + bb)*x[1]
    ymin = (aa*xmin + bb)*xmin + cc
    return xmin,ymin

def find_extrema(x,y,thresh=.001):
    n = len(x)
    result = []
    for i in range(1,n-2):
        if (y[i]-y[i-1])*(y[i+1]-y[i]) < 0:
            x3 = x[i-1:i+2]
            y3 = y[i-1:i+2]
            if abs(parabfit(x3,y3)[1]) > thresh:
                result.append(parabfit(x3,y3))
    return np.array(result) 


def derivs(ccfIn,cdelt1=0.82,npc=10,use_svd=True):
    if use_svd:
        u,s,v = np.linalg.svd(ccfIn,full_matrices = False)
        ccf2d = np.dot(u[:,:npc] * s[:npc], v[:npc])
    else:
        ccf2d = ccfIn
        
    nvel = np.shape(ccf2d)[1]
    fm = np.roll(ccf2d,-1)
    f0 = np.array(ccf2d)
    fp = np.roll(ccf2d,1)
    # First derivative
    dfdv = (fm-fp) / 2 / cdelt1
    # Fix the ends
    dfdv[:,0] = dfdv[:,1]
    dfdv[:,nvel-1] = dfdv[:,nvel-2]
    # Second derivative
    dfdvm = np.roll(dfdv,-1)
    dfdvp = np.roll(dfdv,1)
    d2fdv2 = (dfdvm-dfdvp) / 2 / cdelt1
    # Fix the ends
    d2fdv2[:,0] = d2fdv2[:,1]
    d2fdv2[:,nvel-1] = d2fdv2[:,nvel-2]
    # third derivative
    d2fdv2m = np.roll(dfdv,-1)
    d2fdv2p = np.roll(dfdv,1)
    d3fdv3 = (d2fdv2m-d2fdv2p) / 2 / cdelt1
    # Fix the ends
    d3fdv3[:,0] = d3fdv3[:,1]
    d3fdv3[:,nvel-1] = d3fdv3[:,nvel-2]
    return dfdv,d2fdv2,d3fdv3

def vel2ccf(vseries,dfdv,d2fdv2,d3fdv3):
    ccfVel = np.transpose(vseries*np.transpose(dfdv)
                          +  vseries*vseries*np.transpose(d2fdv2)/2    
                          +  vseries*vseries*vseries*np.transpose(d3fdv3)/6
                         )
    return ccfVel

def bary2hel(ccfBary, bsrv, cdelt1=0.82, frac=1.000):
    # Normalise (divide by np.mean of each row)
    ccfNorm = np.transpose(np.transpose(ccfBary) / np.mean(ccfBary,axis=1))

    dfdv,d2fdv2,d3fdv3 = derivs(ccfNorm,cdelt1)
    # Orthogonalise bsrv and scale if necessary as in getvel_matrix
    dbsrv = (bsrv - np.mean(bsrv)) * frac
    reflex = vel2ccf(dbsrv,dfdv,d2fdv2,d3fdv3)

    # Subtract (some fraction of) solar reflex motion 
    ccfHel = ccfNorm - reflex
    
    return np.transpose(np.transpose(ccfHel) * np.mean(ccfBary,axis=1))


def getweights(ccf2d,npc=10):
    # Subtract column np.means to get residual map
    ccfRes = ccf2d-np.mean(ccf2d,axis=0)
    u,s,v = np.linalg.svd(ccfRes)
    ccfSmoo = np.dot(u[:,:npc]*s[:npc],v[:npc])
    ccfDiff = ccfRes-ccfSmoo
    stdvec = np.std(ccfDiff,axis=1)
    invar = 1/stdvec/stdvec
    ccfWeight = (np.zeros_like(ccf2d).T+invar).T
    return ccfWeight

def covtrgl(vel,rowvar,hdvbase):
    # Scale a triangular profile of width hdvbase by the row variances.
    matrix=[]
    elem = -1
    for v1 in vel:
        elem = elem + 1
        v2 = vel
        # Don't normalise the ridge area
        profile = rowvar * np.clip(1 - np.abs(v2-v1)/hdvbase,0,None)
        # Normalise the ridge area
        #profile = rowvar * np.clip(1 - np.abs(v2-v1)/hdvbase,0,None)/hdvbase/hdvbase
        matrix.append(profile)
    # Convert to numpy array and add transpose
    matrix = np.array(matrix) 
    return (matrix+np.transpose(matrix))/2

def covdiag(vel,rowvar):
    ccfVar = 1/rowvar
    matrix = np.diag(ccfVar) 
    return matrix

def eval(dv,ccfRes,dfdv,d2fdv2,invar):
    res = ccfRes - dv * (dfdv + dv*d2fdv2/2)
    chs = np.dot(res,np.dot(invar,res))
    return chs

def getvel_matrix(vel,ccf2d,err2d,hdvbase=0.82,cdelt1=0.82,k=10,sc=1.000,use_ccfShift=False,empirical=False):
    # General-purpose velocity measurement routine 
    # for an arbitrary covariance matrix.
    # sc is a scale factor which will differ from unity only for pixel-sampled solar ccfs.
    #ccf2d = forceAlign(ccf2d)
    if empirical:
        # Estimate errors from RMS of normalised CCF rows
        ccfWeight = getweights(ccf2d,npc=k)
    else:
        # Use row means of normalised CCF error array
        stdvec = np.mean(err2d,axis=1)
        invar = 1/stdvec/stdvec
        ccfWeight = (np.zeros_like(ccf2d).T+invar).T
    # Use default use_svd=true to get dfdv=Ctilde'
    dfdv,d2fdv2,d3fdv3 = derivs(ccf2d,cdelt1,npc=k)
    # Uncomment to use derivative of mean profile only
    #dfdv,d2fdv2,d3fdv3 = derivs(ccf2d,cdelt1,npc=1)
    # Subtract mean CCF profile
    ccfRes = ccf2d - np.mean(ccf2d,axis=0)
    uu,ss,vv = np.linalg.svd(ccfRes)
    # Compute reduced-rank covariance matrix
    covar = np.dot(vv[:k].T,np.dot(np.diag(ss[:k]*ss[:k]),vv[:k]))/len(ccf2d)
    
    # Initialise output arrays
    vseries = []
    eseries = []
    cseries = []
    sseries = []
    for row in range(len(ccf2d)):  
        # all the elements of ccfVar are identical.
        rowvar = 1/(ccfWeight[row])        
        # Triangular kernel, default base half-width 1.5 km/sec
        matrix = covtrgl(vel,rowvar,hdvbase)
        #rowvar = err2d[row]*err2d[row]
        #matrix = covdiag(vel,rowvar)
        # Full matrix model including reduced-rank profile-covariance model
        matrix2 = matrix+covar 
        # Invert the covariance matrices
        invar = np.linalg.inv(matrix)
        invar2 = np.linalg.inv(matrix2)
        if use_ccfShift:
            x = ccfShift[row]
            p = dfdv[row]
        else:
            x = ccfRes[row]
            p = dfdv[row]
        # Measure velocity using row covariances only
        vhat = -np.dot(x.T,np.dot(invar,p)) / np.dot(p.T,np.dot(invar,p))
        # Measure formal photon-noise uncertainty on velocity 
        var = 1 / np.dot(p.T,np.dot(invar,p))  
        # The old DRS adds 43 cm/s of drift noise to the photon-noise variances.
        # In practice this is unnecessary - any such systematics will contribute
        # to varsys.
        # Measure velocity uncertainty including large-scale time-varying systematics
        varsys = 1 / np.dot(p.T,np.dot(invar2,p)) 
        vseries.append(vhat)
        eseries.append(np.sqrt(var))
        sseries.append(np.sqrt(varsys))
    return np.array(vseries)/sc,np.array(eseries)/sc,np.array(sseries)/sc

def getvel_matrix_quad(vel,ccf2d,err2d,hdvbase=0.82,cdelt1=0.82,k=10,sc=1.000,use_ccfShift=False,empirical=False):
    # General-purpose velocity measurement routine 
    # for an arbitrary covariance matrix.
    # sc is a scale factor which will differ from unity only for pixel-sampled solar ccfs.
    #ccf2d = forceAlign(ccf2d)
    #import scipy
    
    if empirical:
        # Estimate errors from RMS of normalised CCF rows
        ccfWeight = getweights(ccf2d,npc=k)
    else:
        # Use row means of normalised CCF error array
        stdvec = np.mean(err2d,axis=1)
        invar = 1/stdvec/stdvec
        ccfWeight = (np.zeros_like(ccf2d).T+invar).T
    # Use default use_svd=true to get dfdv=Ctilde'
    dfdv,d2fdv2,d3fdv3 = derivs(ccf2d,cdelt1,npc=k)
    # Uncomment to use derivative of mean profile only
    #dfdv,d2fdv2,d3fdv3 = derivs(ccf2d,cdelt1,npc=1)
    # Subtract mean CCF profile
    ccfRes = ccf2d - np.mean(ccf2d,axis=0)
    uu,ss,vv = np.linalg.svd(ccfRes)
    # Compute reduced-rank covariance matrix
    covar = np.dot(vv[:k].T,np.dot(np.diag(ss[:k]*ss[:k]),vv[:k]))/len(ccf2d)
    
    # Initialise output arrays
    vseries = []
    eseries = []
    cseries = []
    sseries = []
    for row in range(len(ccf2d)):  
        # all the elements of ccfVar are identical.
        rowvar = 1/(ccfWeight[row])        
        # Triangular kernel, default base half-width 1.5 km/sec
        matrix = covtrgl(vel,rowvar,hdvbase)
        #rowvar = err2d[row]*err2d[row]
        #matrix = covdiag(vel,rowvar)
        # Full matrix model including reduced-rank profile-covariance model
        matrix2 = matrix+covar 
        # Invert the covariance matrices
        invar = np.linalg.inv(matrix)
        invar2 = np.linalg.inv(matrix2)
        if use_ccfShift:
            x = ccfShift[row]
            p1 = dfdv[row]
            p2 = d2fdv2[row]
        else:
            x = ccfRes[row]
            p1 = dfdv[row]
            p2 = d2fdv2[row]
        # Measure velocity using row covariances only
        solv = opt.minimize(eval,0,args=(x,p1,p2,invar))
        vhat = -solv.x[0]
        # Measure formal photon-noise uncertainty on velocity 
        var = solv.hess_inv[0][0]
        # The old DRS adds 43 cm/s of drift noise to the photon-noise variances.
        # In practice this is unnecessary - any such systematics will contribute
        # to varsys.
        # Measure velocity uncertainty including large-scale time-varying systematics
        solv2 = opt.minimize(eval,-vhat,args=(x,p1,p2,invar2))
        varsys = solv2.hess_inv[0][0]
        vseries.append(vhat)
        eseries.append(np.sqrt(var))
        sseries.append(np.sqrt(varsys))
    return np.array(vseries)/sc,np.array(eseries)/sc,np.array(sseries)/sc

def inject(ccfHel, bjd, T0, P, K, inject=True, cdelt1=0.82):
    # Similar to align, but starting from the heliocentric map.
    # Normalise (divide by np.mean of each row)
   
    ccfNorm = (ccfHel.T / np.mean(ccfHel,axis=1)).T

    # Compute first and second derivatives of CCF.
    dfdv,d2fdv2,d3fdv3 = derivs(ccfNorm,cdelt1)
    
    # Compute velocity timeseries and np.mean profile derivative
    pha = np.mod(bjd-T0,P)/P
    vInject = - K * np.sin (2 * np.pi * pha)

    # If inject is set to False, remove the signal instead of adding it.
    if inject:
        reflex = vel2ccf(vInject,dfdv,d2fdv2,d3fdv3)
    
    else:
        reflex = -vel2ccf(vInject,dfdv,d2fdv2,d3fdv3)

    # Subtract or add noise-free model of secular reflex motion and rescale 
    ccfNew = ccfNorm - reflex

    # Return corrected map with original row scaling
    return (ccfNew.T * np.mean(ccfHel,axis=1)).T




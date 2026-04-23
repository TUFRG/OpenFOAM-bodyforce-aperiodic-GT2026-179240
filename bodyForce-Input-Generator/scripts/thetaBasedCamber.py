#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 12:24:05 2025

@author: adekola
"""
import argparse
import numpy as np
#import matplotlib.pyplot as plt


#%% Used Functions
def surfnorm_python(X, Y, Z):
    """
    Calculates the 3D normal vectors for a surface defined by X, Y, Z coordinates.

    This function replicates the core functionality of MATLAB's `surfnorm`.
    It computes the normal vector at each point on the surface by taking
    the cross product of the tangent vectors in the two primary grid directions.

    Args:
        X (np.ndarray): A 2D array of x-coordinates.
        Y (np.ndarray): A 2D array of y-coordinates.
        Z (np.ndarray): A 2D array of z-coordinates (surface height).

    Returns:
        tuple: A tuple containing three 2D numpy arrays (Nx, Ny, Nz),
               representing the x, y, and z components of the unit normal
               vectors at each corresponding point in the X, Y, Z arrays.
    """
    dX_u, dX_v = np.gradient(X, edge_order=2)
    dY_u, dY_v = np.gradient(Y, edge_order=2)
    dZ_u, dZ_v = np.gradient(Z, edge_order=2)

    # Compute the components of the normal vector using the cross product formula:    
    Nx = dZ_u * dY_v - dY_u * dZ_v 
    Ny = dX_u * dZ_v - dZ_u * dX_v
    Nz = dY_u * dX_v - dX_u * dY_v

    # Calculate the magnitude (norm) of each normal vector
    norm = np.sqrt(Nx**2 + Ny**2 + Nz**2)
    Nx = np.where(norm != 0, Nx / norm, 0)
    Ny = np.where(norm != 0, Ny / norm, 0)
    Nz = np.where(norm != 0, Nz / norm, 0)

    return Nx, Ny, Nz


def computeCamber(As, Bs, Cs, Nbs, sSections, Ns, pOrA):
    A = []
    Br = []
    Bc = []
    Bth = []
    sNth = np.zeros([Nbs, sSections, Ns+2, 4])
    sNr = np.zeros([Nbs, sSections, Ns+2, 4])
    sNc = np.zeros([Nbs, sSections, Ns+2, 4])
    if pOrA == 0:  # Periodic
        # Compute normals ONCE for reference blade
        naS_ref, nbS_ref, ncS_ref = surfnorm_python(As[0], Bs[0], Cs[0])
        a_ref = As[0,:,:]
        b_ref = Bs[0,:,:]
        theta_ref = np.arctan2(b_ref, a_ref)
        r_ref = np.sqrt(a_ref**2 + b_ref**2)
        for p in range(Nbs):
            blade_angle = p * 2*np.pi / Nbs
            cos_p = np.cos(blade_angle)
            sin_p = np.sin(blade_angle)
            
            # Rotate reference normals
            naS = cos_p * naS_ref - sin_p * nbS_ref
            nbS = sin_p * naS_ref + cos_p * nbS_ref
            ncS = ncS_ref.copy()
            
            # Transform to cylindrical
            theta_pts = theta_ref + blade_angle
            if np.max(theta_pts)  - np.min(theta_pts ) > np.pi:
                theta_pts = np.unwrap(theta_pts.flatten()).reshape(theta_pts.shape)
                print('Yes')
            # theta_pts = np.arctan2(np.sin(theta_pts), np.cos(theta_pts))
            r_pts =r_ref
            
            nr = np.cos(theta_pts) * naS + np.sin(theta_pts) * nbS
            nth = -np.sin(theta_pts) * naS + np.cos(theta_pts) * nbS
            norm = np.sqrt(nr**2 + nth**2 + ncS**2)
            # bladeAngles[p,:] = np.degrees(np.arctan2(-nth, ncS))
            # Avoid division by zero and normalize
            norm = np.where(norm != 0, norm, 1.0)  # Replace 0 with 1 to avoid division issues
            nr = nr / norm
            nth = nth / norm
            ncS = ncS / norm
            m, n = naS_ref.shape
            # Store results
            for q in range(m):
                for rq in range(n):
                    sNr[p,rq,q,:] = [theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nr[q,rq]]
                    sNth[p,rq,q,:] = [theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nth[q,rq]]
                    sNc[p,rq,q,:] = [theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], ncS[q,rq]]
                    
            A.append(r_pts)   
            Br.append(nr) 
            Bth.append(nth)
            Bc.append(ncS) 
    else:  # Non-periodic (aperiodic)
        for p in range(Nbs):
            # Compute normals for EACH blade separately
            naS, nbS, ncS = surfnorm_python(As[p], Bs[p], Cs[p])
            
            # Transform to cylindrical
            a = As[p,:,:]
            b = Bs[p,:,:]
            theta_pts = np.arctan2(b, a)
            r_pts = np.sqrt(a**2 + b**2)
            
            nr = np.cos(theta_pts) * naS + np.sin(theta_pts) * nbS
            nth = -np.sin(theta_pts) * naS + np.cos(theta_pts) * nbS
            
            norm = np.sqrt(nr**2 + nth**2 + ncS**2)
    
            # Avoid division by zero and normalize
            norm = np.where(norm != 0, norm, 1.0)  # Replace 0 with 1 to avoid division issues
            nr = nr / norm
            nth = nth / norm
            ncS = ncS / norm
            
            m, n = naS.shape
            # Store results
            for q in range(m):
                for rq in range(n):
                    sNr[p,rq,q,:] = [theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nr[q,rq]]
                    sNth[p,rq,q,:] = [theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nth[q,rq]]
                    sNc[p,rq,q,:] = [theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], ncS[q,rq]]
            A.append(r_pts)   
            Br.append(nr) 
            Bth.append(nth)
            Bc.append(ncS)  
    return sNr, sNth, sNc, A, Br, Bth, Bc
#%%
parser = argparse.ArgumentParser(description="Getting the Geometry parameters.")
parser.add_argument(
    "Nbr", 
    type=int,
    help="number of rotor blades"
)
parser.add_argument(
    "Nbs", 
    type=int,
    help="number of stator blades"
)
parser.add_argument(
    "rSections", 
    type=int,
    help="number of rotor blade profiles"
)
parser.add_argument(
    "sSections", 
    type=int,
    help="number of stator blade profiles"
)
parser.add_argument(
    "periodicOrAperiodic", 
    type=int,
    help="periodic (0) or nonPeriodic (1)"
)

# Parse arguments
args = parser.parse_args()

# Use parsed arguments
Nbr = args.Nbr  # Total number of blades in the annulus for rotor
Nbs = args.Nbs  # Total number of blades in the annulus for stator
rSections = args.rSections  # Number of rotor blade profiles
sSections = args.sSections  # Number of stator blade profiles
periodicOrAperiodic = args.periodicOrAperiodic  # if periodic select 0 otherwise select 1


Nr=101
Ns=101
periodicORaperiodic = 1
if periodicORaperiodic == 0:
    dataPath = '../processedData/periodic/'
    outputDataPath = '../inputData/periodic/'
else:
    dataPath = '../processedData/nonPeriodic/'
    outputDataPath = '../inputData/nonPeriodic/'
    
statorCamber = np.empty([Nbs, sSections, Ns+2, 3]) #x,y,z
rotorCamber = np.empty([Nbr, rSections, Nr+2, 3])

for a in range(Nbr):
    for b in range(rSections):    
        rFile = np.loadtxt(dataPath + '/rotor/blade{}'.format(a) + '/camber{}.txt'.format(b), delimiter=',') 
        vecRotor = rFile[0, :] - rFile[1, :]
        newRLE = rFile[0, :] + 2 * vecRotor
        vecRotor = rFile[-1, :] - rFile[-2,:] 
        newRTE = rFile[-1, :] + 2 * vecRotor
        rotorCamber[a,b,:] =  np.concatenate(([newRLE], rFile, [newRTE]))
for c in range(Nbs):
    for d in range(sSections):
        sfile = np.loadtxt(dataPath + '/stator/blade{}'.format(c) + '/camber{}.txt'.format(d), delimiter=',') 
        # Extend leading edge
        vecStator = sfile[0, :] - sfile[1, :]
        newSLE = sfile[0, :] + 2 * vecStator    
        vecStator = sfile[-1, :] - sfile[-2, :]
        newSTE = sfile[-1, :] + 2 * vecStator
        statorCamber[c,d,:] = np.concatenate(([newSLE], sfile, [newSTE]))

#%% Determine the Camber Normals        
Ar = np.zeros((Nbr, Nr+2, rSections))
Br = np.zeros((Nbr, Nr+2, rSections))
Cr = np.zeros((Nbr, Nr+2, rSections))
for e in range(Nbr):
    for f in range(rSections):
        Ar[e,:,f] = rotorCamber[e,f,:,0]
        Br[e,:,f] = rotorCamber[e,f,:,1]
        Cr[e,:,f] = rotorCamber[e,f,:,2]       

As = np.zeros((Nbs, Ns+2, sSections))
Bs = np.zeros((Nbs, Ns+2, sSections))
Cs = np.zeros((Nbs, Ns+2, sSections))
for g in range(Nbs):
    for h in range(sSections):
        As[g,:,h] = statorCamber[g,h,:,0]
        Bs[g,:,h] = statorCamber[g,h,:,1]
        Cs[g,:,h] = statorCamber[g,h,:,2]
        
GG = Cr[0]
# GA = Cr[1]
#%%
Ra = []
Rr = []
Rc = []
Rth = []
rNth = np.zeros([Nbr, rSections, Nr+2,4])
rNr = np.zeros([Nbr, rSections, Nr+2, 4])
rNc = np.zeros([Nbr, rSections, Nr+2, 4])

rNr, rNth, rNc, Ra, Rr, Rth, Rc = computeCamber(Ar, Br, Cr, Nbr, rSections, Nr, periodicORaperiodic)

A = []
Br = []
Bc = []
Bth = []
sNth = np.zeros([Nbs, sSections, Ns+2, 4])
sNr = np.zeros([Nbs, sSections, Ns+2, 4])
sNc = np.zeros([Nbs, sSections, Ns+2, 4])

sNr, sNth, sNc, A, Br, Bth, Bc = computeCamber(As, Bs, Cs, Nbs, sSections, Ns, periodicORaperiodic)
#%%
'''
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
xx = 0
# Nr contour
im1 = ax1.contourf(Cs[xx], A[xx], Br[xx], levels=200, cmap='viridis')
ax1.set_xlabel('C')
ax1.set_ylabel('R')
ax1.set_title(r'Rotor $|\hat{n}_r|$')
ax1.set_aspect('equal')
plt.colorbar(im1, ax=ax1)
im1.set_clim(-1, 1)

# Nth contour
im2 = ax2.contourf(Cs[xx], A[xx], Bth[xx], levels=200, cmap='viridis')
ax2.set_xlabel('C')
ax2.set_ylabel('R')
ax2.set_title(r'Rotor $|\hat{n}_{\theta}|$')
ax2.set_aspect('equal')
plt.colorbar(im2, ax=ax2)
im2.set_clim(-1, 1)

# Nc contour
im3 = ax3.contourf(Cs[xx], A[xx], Bc[xx], levels=200, cmap='viridis')
ax3.set_xlabel('C')
ax3.set_ylabel('R')
ax3.set_title(r'Rotor $|\hat{n}_c|$')
ax3.set_aspect('equal')
plt.colorbar(im3, ax=ax3)
im3.set_clim(-1, 1)

plt.tight_layout()
plt.show()  

#%%
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
xx = 10
# Nr contour
im1 = ax1.contourf(Cr[xx], Ra[xx], Rr[xx], levels=200, cmap='viridis')
ax1.set_xlabel('C')
ax1.set_ylabel('R')
ax1.set_title(r'Rotor $|\hat{n}_r|$')
ax1.set_aspect('equal')
plt.colorbar(im1, ax=ax1)
im1.set_clim(-1, 1)

# Nth contour
im2 = ax2.contourf(Cr[xx], Ra[xx], Rth[xx], levels=200, cmap='viridis')
ax2.set_xlabel('C')
ax2.set_ylabel('R')
ax2.set_title(r'Rotor $|\hat{n}_{\theta}|$')
ax2.set_aspect('equal')
plt.colorbar(im2, ax=ax2)
im2.set_clim(-1, 1)

# Nc contour
im3 = ax3.contourf(Cr[xx], Ra[xx], Rc[xx], levels=200, cmap='viridis')
ax3.set_xlabel('C')
ax3.set_ylabel('R')
ax3.set_title(r'Rotor $|\hat{n}_c|$')
ax3.set_aspect('equal')
plt.colorbar(im3, ax=ax3)
im3.set_clim(-1, 1)

plt.tight_layout()
plt.show() 
'''
#%%
rotorNr = rNr.reshape(-1,4)   #theta, radius, axial, camberNormal
rotorNc = rNc.reshape(-1,4)
rotorNth = rNth.reshape(-1,4)
statorNr = sNr.reshape(-1,4)
statorNc = sNc.reshape(-1,4)
statorNth = sNth.reshape(-1,4)

np.savetxt(outputDataPath + '/rotorNr.txt', rotorNr, delimiter=',')
np.savetxt(outputDataPath + '/rotorNc.txt', rotorNc, delimiter=',')
np.savetxt(outputDataPath + '/rotorNth.txt', rotorNth, delimiter=',')

np.savetxt(outputDataPath + '/statorNr.txt', statorNr, delimiter=',')
np.savetxt(outputDataPath + '/statorNc.txt', statorNc, delimiter=',')
np.savetxt(outputDataPath + '/statorNth.txt', statorNth, delimiter=',')



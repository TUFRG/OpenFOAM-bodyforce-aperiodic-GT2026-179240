#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 12:24:05 2025

@author: adekola
"""
import numpy as np
import matplotlib.pyplot as plt


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
#%%
statorPath = '../processedData/'
filePath = '../inputData'
periodicORaperiodic = 0 #choose 0 if periodic and 1 if aperiodic 

sSections = 5
Nbs = 60 #Number of blades
Ns = 101
statorCamber = np.empty([Nbs, sSections, Ns+2, 3]) #x,y,z

for c in range(Nbs):
    for d in range(sSections):
        if periodicORaperiodic == 0:
            sfile = np.loadtxt(statorPath + '/blade{}'.format(c+1) + '/blade{}.txt'.format(d+1), delimiter=',') 
            # filePath = './inputData/periodic/'
        else:
            sfile = np.loadtxt(statorPath + 'nonPeriodic/blade{}'.format(c) + '/camber{}.txt'.format(d), delimiter=',') 
            filePath = './inputData/nonPeriodic/'
        # Extend leading edge
        vecStator = sfile[0, :] - sfile[1, :]
        newSLE = sfile[0, :] + 2 * vecStator    
        vecStator = sfile[-1, :] - sfile[-2, :]
        newSTE = sfile[-1, :] + 2 * vecStator
        statorCamber[c,d,:] = np.concatenate(([newSLE], sfile, [newSTE]))

#%% Determine the Camber Normals        
As = np.zeros((Nbs, Ns+2, sSections))
Bs = np.zeros((Nbs, Ns+2, sSections))
Cs = np.zeros((Nbs, Ns+2, sSections))
for g in range(Nbs):
    for h in range(sSections):
        As[g,:,h] = statorCamber[g,h,:,0]
        Bs[g,:,h] = statorCamber[g,h,:,1]
        Cs[g,:,h] = statorCamber[g,h,:,2]
        
GG = As[0]
# GA = Cr[1]

#%%
# A = []
# Br = []
# Bc = []
# Bth = []
# sNth = np.zeros([Nbs, sSections, Ns+2, 4])
# sNr = np.zeros([Nbs, sSections, Ns+2, 4])
# sNc = np.zeros([Nbs, sSections, Ns+2, 4])
# for p in range(Nbs):
#    naS, nbS, ncS = surfnorm_python(As[p], Bs[p], Cs[p])
#    m, n = naS.shape
#    # Initialize cylindrical normal arrays
#    nr = np.zeros((m, n))
#    nth = np.zeros((m, n))
#    r_pts = np.zeros((m, n))
#    theta_pts = np.zeros((m, n))
#    for q in range(m):
#        for rq in range(n):
#            a = As[p,q,rq]
#            b = Bs[p,q,rq]
#            theta = np.arctan2(b, a)
#            theta_pts[q,rq] = theta
#            r = np.sqrt(a**2 + b**2)
#            r_pts[q,rq] = r
#            # Transform normal components
#            nr[q,rq] = np.cos(theta) * naS[q,rq] + np.sin(theta) * nbS[q,rq]
#            nth[q,rq] = -np.sin(theta) * naS[q,rq] + np.cos(theta) * nbS[q,rq]  
#            sNr[p,rq,q,:] = np.array([theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nr[q,rq]])
#            sNth[p,rq,q,:] = np.array([theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nth[q,rq]])
#            sNc[p,rq,q,:] = np.array([theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], ncS[q,rq]])  
#    A.append(r_pts)   
#    Br.append(nr)   
#    Bth.append(nth)
#    Bc.append(ncS)
# naS_ref, nbS_ref, ncS_ref = surfnorm_python(As[0], Bs[0], Cs[0])
# for p in range(Nbs):
#     # Rotation angle for this blade
#     blade_angle = p * 2*np.pi / Nbs
#     cos_p = np.cos(blade_angle)
#     sin_p = np.sin(blade_angle)
    
#     # Rotate the ENTIRE normal field at once (vectorized)
#     naS = cos_p * naS_ref - sin_p * nbS_ref
#     nbS = sin_p * naS_ref + cos_p * nbS_ref
#     ncS = ncS_ref.copy()  # Axial doesn't change
    
#     # Get theta and r for all points on THIS blade (vectorized)
#     a = As[p,:,:]
#     b = Bs[p,:,:]
#     theta_pts = np.arctan2(b, a)
#     r_pts = np.sqrt(a**2 + b**2)
#     m, n = naS_ref.shape
#     # Transform to cylindrical (vectorized)
#     nr = np.cos(theta_pts) * naS + np.sin(theta_pts) * nbS
#     nth = -np.sin(theta_pts) * naS + np.cos(theta_pts) * nbS
    
#     # Store results
#     for q in range(m):
#         for rq in range(n):
#             sNr[p,rq,q,:] = np.array([theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nr[q,rq]])
#             sNth[p,rq,q,:] = np.array([theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], nth[q,rq]])
#             sNc[p,rq,q,:] = np.array([theta_pts[q,rq], r_pts[q,rq], Cs[p,q,rq], ncS[q,rq]])
# # for p in range(Nbs):
# #     # Rotation angle for this blade
# #     blade_angle = p * 2*np.pi / Nbs
# #     # Rotation matrix for this blade
# #     cos_p = np.cos(blade_angle)
# #     sin_p = np.sin(blade_angle)
# #     m, n = naS_ref.shape
# #     nr = np.zeros((m, n))
# #     nth = np.zeros((m, n))
# #     r_pts = np.zeros((m, n))
# #     for q in range(m):
# #         for rq in range(n):
# #             # Rotate the REFERENCE normal to this blade's position
# #             naS_p = cos_p * naS_ref[q,rq] - sin_p * nbS_ref[q,rq]
# #             nbS_p = sin_p * naS_ref[q,rq] + cos_p * nbS_ref[q,rq]
# #             ncS_p = ncS_ref[q,rq]
            
# #             # Get position for THIS blade
# #             a = As[p,q,rq]
# #             b = Bs[p,q,rq]
# #             theta = np.arctan2(b, a)
# #             r = np.sqrt(a**2 + b**2)
# #             r_pts[q,rq] = r
# #             # Transform to cylindrical at THIS position
# #             nr[q,rq] = np.cos(theta) * naS_p + np.sin(theta) * nbS_p
# #             nth[q,rq] = -np.sin(theta) * naS_p + np.cos(theta) * nbS_p
            
# #             sNr[p,rq,q,:] = np.array([theta, r, Cs[p,q,rq], nr[q,rq]])
# #             sNth[p,rq,q,:] = np.array([theta, r, Cs[p,q,rq], nth[q,rq]])
# #             sNc[p,rq,q,:] = np.array([theta, r, Cs[p,q,rq], ncS_p])
#     A.append(r_pts)   
#     Br.append(nr) 
#     Bth.append(nth)
#     Bc.append(ncS) 
#%%
A = []
Br = []
Bc = []
Bth = []
Normal = []
sNth = np.zeros([Nbs, sSections, Ns+2, 4])
sNr = np.zeros([Nbs, sSections, Ns+2, 4])
sNc = np.zeros([Nbs, sSections, Ns+2, 4])
bladeAngles = np.zeros([Nbs, Ns+2, sSections])

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
    norm1 = np.sqrt(sNr[p,rq,:,3]**2 + sNth[p,rq,:,3]**2 + sNc[p,rq,:,3]**2)
    A.append(r_pts)   
    Br.append(nr) 
    Bth.append(nth)
    Bc.append(ncS)  
    Normal.append(norm1)
#%%
# Get all values with their positions
violations = []
for i, sublist in enumerate(Normal):
    for j, value in enumerate(sublist):
        if value > 1:
            violations.append((i, j, value))

print(f"Found {len(violations)} value(s) > 1:")
for i, j, value in violations:
    print(f"  [{i}][{j}] = {value}")
#%%
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
xx = 5
# Nr contour
im1 = ax1.contourf(Cs[xx], A[xx], Br[xx], levels=200, cmap='jet')
ax1.set_xlabel('C')
ax1.set_ylabel('R')
ax1.set_title(r'Rotor $|\hat{n}_r|$')
ax1.set_aspect('equal')
plt.colorbar(im1, ax=ax1)
im1.set_clim(-1, 1)

# Nth contour
im2 = ax2.contourf(Cs[xx], A[xx], Bth[xx], levels=200, cmap='jet')
ax2.set_xlabel('C')
ax2.set_ylabel('R')
ax2.set_title(r'Rotor $|\hat{n}_{\theta}|$')
ax2.set_aspect('equal')
plt.colorbar(im2, ax=ax2)
im2.set_clim(-1, 1)

# Nc contour
im3 = ax3.contourf(Cs[xx], A[xx], Bc[xx], levels=200, cmap='jet')
ax3.set_xlabel('C')
ax3.set_ylabel('R')
ax3.set_title(r'Rotor $|\hat{n}_c|$')
ax3.set_aspect('equal')
plt.colorbar(im3, ax=ax3)
im3.set_clim(-1, 1)

plt.tight_layout()
# plt.show()  
# plt.savefig('/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/plot/camberSurfaceNormals.jpg', dpi=500)#, bbox_inches='tight')
#%%
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
xx = 0
Anorm = A[xx]/np.max(A[xx])
Cnorm = (Cs[xx]-np.min(Cs[xx]))/np.max(A[xx])

# Define common levels from -1 to 1
levels = np.linspace(-1, 1, 200, endpoint=True)

# Nr contour
im1 = ax1.contourf(Cnorm, Anorm, Br[xx], levels=levels, cmap='viridis', vmin=-1, vmax=1, extend='neither')
ax1.set_xlabel(r'c/$r_{\mathrm{tip}}$', fontsize=18)
ax1.set_ylabel(r'r/$r_{\mathrm{tip}}$', fontsize=18)
ax1.set_title(r'Rotor $|\hat{n}_r|$', fontsize=18)
ax1.set_aspect('equal')
# plt.colorbar(im1, ax=ax1)

# Nth contour
im2 = ax2.contourf(Cnorm, Anorm, Bth[xx], levels=levels, cmap='viridis', vmin=-1, vmax=1, extend='neither')
ax2.set_xlabel(r'c/$r_{\mathrm{tip}}$', fontsize=18)
# ax2.set_ylabel('R')
ax2.set_title(r'Rotor $|\hat{n}_{\theta}|$', fontsize=18)
ax2.set_aspect('equal')
# plt.colorbar(im2, ax=ax2)

# Nc contour
im3 = ax3.contourf(Cnorm, Anorm, Bc[xx], levels=levels, cmap='viridis', vmin=-1, vmax=1, extend='neither')
ax3.set_xlabel(r'c/$r_{\mathrm{tip}}$', fontsize=18)
# ax3.set_ylabel('R')
ax3.set_title(r'Rotor $|\hat{n}_c|$', fontsize=18)
ax3.set_aspect('equal')
# plt.colorbar(im3, ax=ax3)
fig.colorbar(im3, ax=[ax1, ax2, ax3], fraction=0.046, pad=0.04)
# plt.savefig('/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/plot/camberSurfaceNormals.jpg', dpi=500)
#%%
thetaB = np.linspace(0, 2*np.pi,Nbs, endpoint=False)
angle = np.zeros(Nbs)
for a in range(Nbs):
    for b in range(5):
        prev_blade = (a-1) % Nbs
        next_blade = (a+1) % Nbs
        pC = statorCamber[a,b,0,:]
        pN = statorCamber[next_blade,b,0,:]
        rC = np.sqrt(pC[0]**2 + pC[1]**2)
        rN = np.sqrt(pN[0]**2 + pN[1]**2)
        rAvg = 0.5*(rC+rN)
        qC = np.arctan2(pC[1], pC[0])
        qN = np.arctan2(pN[1], pN[0])
        # Handle angle wrapping
        dQ = qN - qC
        if dQ < 0:
            dQ += 2*np.pi
        # angle[a] = rAvg * np.degrees(dQ)
        angle[a] = np.degrees(dQ)
        # angle[a] = dQ
plt.plot(thetaB, angle,'k', label='varingAngle')
#%%
# After line 334, add stagger calculation
stagger = np.zeros(Nbs)
normalAngle = np.zeros(Nbs)
for a in range(Nbs):
    for b in range(5):  # Loop over sections, or choose one specific section
        # Leading edge (first point after extension, so index 1)
        pLE = statorCamber[a, b, 1, :]  # [x, y, z]
        # Trailing edge (second to last point before extension, so index -2)
        pTE = statorCamber[a, b, -2, :]  # [x, y, z]
        # Chord vector
        chord_vec = pTE - pLE
        # Axial component (z-direction)
        axial_comp = chord_vec[2]
        # Tangential component in cylindrical coords
        # Project chord onto tangential direction at midchord
        rLE = np.sqrt(pLE[0]**2 + pLE[1]**2)
        rTE = np.sqrt(pTE[0]**2 + pTE[1]**2)
        rMid = 0.5 * (rLE + rTE)
        thetaLE = np.arctan2(pLE[1], pLE[0])
        thetaTE = np.arctan2(pTE[1], pTE[0])
        dTheta = thetaTE - thetaLE
        if dTheta < -np.pi:
            dTheta += 2*np.pi
        elif dTheta > np.pi:
            dTheta -= 2*np.pi
        tangential_comp = rMid * dTheta
        # tangential_comp =  0.9*dTheta
        # Stagger angle (positive = leaning in direction of rotation)
        # stagger[a] = np.degrees(
        stagger[a] = np.degrees(np.arctan2(tangential_comp, axial_comp))
        normalAngle[a] = np.degrees(np.arctan2(-sNth[a,b,-1,3], sNc[a,b,-1,3]))#+90

# plt.figure()
# plt.plot(thetaB, stagger, 'b', label='Stagger Angle')
# # plt.plot(thetaB, normalAngle, 'r', label='Normal Angle')
# plt.xlabel('Circumferential Position (rad)')
# plt.ylabel('Stagger Angle (deg)')
# plt.legend()
# plt.grid(True)
#%%

#%%
statorNr = sNr.reshape(-1,4)
statorNc = sNc.reshape(-1,4)
statorNth = sNth.reshape(-1,4)

np.savetxt(filePath + '/statorNr.txt', statorNr, delimiter=',')
np.savetxt(filePath + '/statorNc.txt', statorNc, delimiter=',')
np.savetxt(filePath + '/statorNth.txt', statorNth, delimiter=',')

#%%






            
            
            
            
            
            
            

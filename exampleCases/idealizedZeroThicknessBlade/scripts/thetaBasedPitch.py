#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 14:21:30 2026

@author: adekola
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def cart2pol(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.atan2(y,x)
    z = z
    return theta, r, z

def cart2pol2(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.atan2(y,x)
    theta = np.unwrap(theta)
    z = z
    return theta, r, z

periodicORaperiodic = 0
dataPath = '../processedData/'
filePath = '../inputData/'
sSections = 5
saveData = 1
Ns = 101
Nbs = 60

# Load blade geometry
blade = np.empty([Nbs, sSections, Ns, 3])
for a in range(Nbs):
    for b in range(sSections):
        blade[a, b, :] = np.loadtxt(
            dataPath + '/blade{}/blade{}.txt'.format(a+1, b+1), 
            delimiter=','
        )

# Convert to cylindrical coordinates
bladeCyl = np.empty([Nbs, sSections, Ns, 3])
for c in range(Nbs):
    for d in range(sSections):
        for e in range(Ns):
            bladeCyl[c, d, e, :] = np.array(
                cart2pol(blade[c, d, e, 0], blade[c, d, e, 1], blade[c, d, e, 2])
            ).T
#%%
bladePitch = np.zeros([Nbs, sSections, Ns, 4])

for f in range(Nbs):
    nextBlade = (f + 1) % Nbs
    for g in range(sSections):
        qC = bladeCyl[f, g, :, 0] #current blade
        rC = bladeCyl[f, g, :, 1]
        zC = bladeCyl[f, g, :, 2]
        qN = bladeCyl[nextBlade, g, :, 0] #next blade
        rN = bladeCyl[nextBlade, g, :, 1]
        zN = bladeCyl[nextBlade, g, :, 2]        
        
        # Use actual arc length as the interpolation coordinate
        # Calculate cumulative arc length for current blade
        dxC = np.diff(blade[f, g, :, 0])
        dyC = np.diff(blade[f, g, :, 1])
        dzC = np.diff(blade[f, g, :, 2])
        dsC = np.sqrt(dxC**2 + dyC**2 + dzC**2) #Arc length
        sC = np.concatenate([[0], np.cumsum(dsC)])
        sC = sC / sC[-1]  # Normalize to [0, 1]
        
        # Calculate arc length for next blade
        dxN = np.diff(blade[nextBlade, g, :, 0])
        dyN = np.diff(blade[nextBlade, g, :, 1])
        dzN = np.diff(blade[nextBlade, g, :, 2])
        dsN = np.sqrt(dxN**2 + dyN**2 + dzN**2) #Arc length
        sN = np.concatenate([[0], np.cumsum(dsN)])
        sN = sN / sN[-1]  # Normalize to [0, 1]        
        
        
        
        # Interpolate next blade onto current blade's parametrization
        interpqN= interp1d(sN, qN, kind='cubic', fill_value='extrapolate')
        interprN = interp1d(sN, rN, kind='cubic', fill_value='extrapolate')
        interpzN = interp1d(sN, zN, kind='cubic', fill_value='extrapolate')
        
        qInterp = interpqN(sC)
        rInterp = interprN(sC)
        zInterp = interpzN(sC)
        
        # Average r and z for passage midpoint
        rMid = 0.5 * (rC + rInterp)
        zMid = 0.5 * (zC + zInterp)
        
        # Compute pitch (angular spacing)
        pitch = qInterp - qC
        # Handle periodicity for last blade wrapping to first
        if f == Nbs - 1:
            pitch = pitch + 2 * np.pi
        # Ensure pitch is positive and reasonable
        pitch = np.where(pitch < 0, pitch + 2*np.pi, pitch)
        pitch = np.where(pitch > 1.5*np.pi, pitch - 2*np.pi, pitch)
        pitch = np.where(pitch < 0, pitch + 2*np.pi, pitch)
        
        # Midpoint theta (passage center)
        qMid = qC + 0.5 * pitch
        
        # Store results: [theta_mid, r_mid, z_mid, pitch]
        bladePitch[f, g, :, 0] = qMid
        bladePitch[f, g, :, 1] = rMid
        bladePitch[f, g, :, 2] = zMid
        bladePitch[f, g, :, 3] = pitch
        # Compute pitch (angular spacing between consecutive blades)


# Reshape and save
bladePitchShaped = bladePitch.reshape(-1, 4)

# Optional: wrap theta_mid back to [0, 2π] for output
# bladePitchShaped[:, 0] = np.mod(bladePitchShaped[:, 0], 2*np.pi)

np.savetxt(filePath + "/beta.txt", bladePitchShaped, delimiter=",")
#%%
# plt.scatter(bladePitchShaped[:,1]*np.cos(bladePitchShaped[:,0]), 
#            bladePitchShaped[:,1]*np.sin(bladePitchShaped[:,0]),
#            c=np.degrees(bladePitchShaped[:,3]), 
#            cmap='viridis', s=10)
# plt.colorbar(label='blade pitch s (deg)')
# plt.axis('equal')
# # plt.xlabel('Axial distance (m)')
# # plt.ylabel('Radial distance (m)')
# # plt.xlabel(r'c/$r_{\mathrm{tip}}$')
# # plt.ylabel(r'r/$r_{\mathrm{tip}}$')
# plt.title('Rotor blade pitch variation ')




























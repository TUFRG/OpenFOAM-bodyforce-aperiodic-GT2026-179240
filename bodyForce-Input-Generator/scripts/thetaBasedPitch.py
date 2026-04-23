#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 14:21:30 2026

@author: adekola
"""
import argparse
import numpy as np
#import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def cart2pol(x,y,z):
    """
    Convert Cartesian to cylindrical WITHOUT arbitrary offsets
    """
    r = np.sqrt(x**2 + y**2)
    theta = np.atan2(y, x)
    # Normalize to [0, 2π)
    # theta = theta % (2 * np.pi)
    return theta, r, z

def cart2pol2(x,y,z):
    import numpy as np
    theta = np.arctan2(y,x) % (2 * np.pi)
    # theta = np.unwrap(theta)
    r = np.sqrt(x**2 + y**2)
    return theta, r, z

def reOrderBlades(bladeCyl, blade):
    """
    Reorder blades in ascending angular order based on their leading edge position
    """    
    # Get reference angle for each blade (use LE of mid-span)
    midSpan = bladeCyl.shape[1] // 2
    bladeAngles = bladeCyl[:, midSpan, 0, 0]  # theta at LE
    # Sort indices
    sorted_indices = np.argsort(bladeAngles)
    # Reorder
    bladeCylSorted = bladeCyl[sorted_indices]
    bladeSorted = blade[sorted_indices]
    
    print(f"\nBlade reordering:")
    print(f"Original order -> Sorted order:")
    for i, orig_idx in enumerate(sorted_indices):
        print(f"  Position {i}: Blade {orig_idx} (angle: {np.degrees(bladeAngles[orig_idx]):.1f}°)")
    return bladeCylSorted, bladeSorted, sorted_indices

def computePitch(bladeCyl, blade, Nbs, sSections, Ns):
    """
    Robust pitch computation with proper handling of periodic boundaries
    """
    bladePitch = np.zeros([Nbs, sSections, Ns, 4])    
    for f in range(Nbs):
        nextBlade = (f + 1) % Nbs
        for g in range(sSections):
            # Current blade
            qC = bladeCyl[f, g, :, 0].copy()
            rC = bladeCyl[f, g, :, 1]
            zC = bladeCyl[f, g, :, 2]
            # Next blade
            qN = bladeCyl[nextBlade, g, :, 0].copy()
            rN = bladeCyl[nextBlade, g, :, 1]
            zN = bladeCyl[nextBlade, g, :, 2]
            # Handle wraparound: if next blade has smaller angle, add 2π
            # if qN[0] < qC[0]:  # Next blade wrapped around
            #     qN = qN + 2 * np.pi
            if np.mean(qN) < np.mean(qC):  # Next blade wrapped around
                qN = qN + 2 * np.pi                
            # Unwrap to handle internal discontinuities
            qC = np.unwrap(qC)
            qN = np.unwrap(qN)
            
            # offset = np.mean(qC)

            # # Shift qC to be near 0, unwrap, shift back
            # qC = np.unwrap(qC - offset) + offset
            
            # # Shift qN to be close to qC (within π), then unwrap
            # qN = qN - offset
            # # Bring qN within (-π, π] of qC's mean
            # qN = (qN + np.pi) % (2 * np.pi) - np.pi
            # qN = np.unwrap(qN) + offset
            
            # Ensure next blade is ahead of current blade angularly
            if np.mean(qN) < np.mean(qC):
                qN = qN + 2 * np.pi
            # Compute arc length for interpolation
            dxC = np.diff(blade[f, g, :, 0])
            dyC = np.diff(blade[f, g, :, 1])
            dzC = np.diff(blade[f, g, :, 2])
            dsC = np.sqrt(dxC**2 + dyC**2 + dzC**2)
            sC = np.concatenate([[0], np.cumsum(dsC)])
            sC = sC / sC[-1] if sC[-1] > 0 else sC
            
            dxN = np.diff(blade[nextBlade, g, :, 0])
            dyN = np.diff(blade[nextBlade, g, :, 1])
            dzN = np.diff(blade[nextBlade, g, :, 2])
            dsN = np.sqrt(dxN**2 + dyN**2 + dzN**2)
            sN = np.concatenate([[0], np.cumsum(dsN)])
            sN = sN / sN[-1] if sN[-1] > 0 else sN
            
            # Interpolate next blade
            interpqN = interp1d(sN, qN, kind='cubic', fill_value='extrapolate')
            interprN = interp1d(sN, rN, kind='cubic', fill_value='extrapolate')
            interpzN = interp1d(sN, zN, kind='cubic', fill_value='extrapolate')
            
            qInterp = interpqN(sC)
            rInterp = interprN(sC)
            zInterp = interpzN(sC)
            # Compute pitch
            rawPitch = qInterp - qC
            # Ensure pitch is positive 
            pitch = np.where(rawPitch < 0, rawPitch + 2*np.pi, rawPitch)
 
            # Mid-passage values
            rMid = 0.5 * (rC + rInterp)
            zMid = 0.5 * (zC + zInterp)
            qMid = (qC + qInterp) / 2
            
            # Wrap back to [0, 2π)
            qMid = qMid % (2*np.pi)
                        
            bladePitch[f, g, :, 0] = qMid
            bladePitch[f, g, :, 1] = rMid
            bladePitch[f, g, :, 2] = zMid
            bladePitch[f, g, :, 3] = pitch
    
    bladePitchShaped = bladePitch.reshape(-1, 4)
    
    # Diagnostic: check pitch statistics
    print(f"\nPitch statistics:")
    print(f"  Mean pitch: {np.degrees(np.mean(bladePitchShaped[:,3])):.4f}°")
    print(f"  Min pitch: {np.degrees(np.min(bladePitchShaped[:,3])):.4f}°")
    print(f"  Max pitch: {np.degrees(np.max(bladePitchShaped[:,3])):.4f}°")
    print(f"  Std pitch: {np.degrees(np.std(bladePitchShaped[:,3])):.4f}°")
    print(f"  Range: {np.degrees(np.max(bladePitchShaped[:,3]) - np.min(bladePitchShaped[:,3])):.4f}°")
    return bladePitchShaped

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

args = parser.parse_args()

# Use parsed arguments
Nbr = args.Nbr  # Total number of blades in the annulus for rotor
Nbs = args.Nbs  # Total number of blades in the annulus for stator
rSections = args.rSections  # Number of rotor blade profiles
sSections = args.sSections  # Number of stator blade profiles
periodicOrAperiodic = args.periodicOrAperiodic  # if periodic select 0 otherwise select 1

Ns = 101 #number o
Nr = 101            
if periodicORaperiodic == 0:
    filePath = '../inputData/periodic/'
    dataPath = '../processedData/periodic/'
else:
    filePath = '../inputData/nonPeriodic/'
    dataPath = '../processedData/nonPeriodic/'
print("Loading blade geometry data...")

statorCamber = np.zeros([Nbs, sSections, Ns, 3]) #x,y,z
rotorCamber = np.zeros([Nbr, rSections, Nr, 3]) #x,y,z

for a in range(Nbs):
    for b in range(sSections):
        statorCamber[a,b,:] = np.loadtxt(dataPath + 'stator/blade{}/camber{}.txt'.format(a,b), delimiter=',')
                
for aa in range(Nbr):
    for bb in range(rSections):
        rotorCamber[aa,bb,:] = np.loadtxt(dataPath + 'rotor/blade{}/camber{}.txt'.format(aa,bb), delimiter=',')
        
statorCamberCyl = np.zeros([Nbs, sSections, Ns, 3]) #theta, r, z
for c in range(Nbs):
    for d in range(sSections):
        statorCamberCyl[c,d,:,:] = np.array(cart2pol(statorCamber[c,d,:,0], statorCamber[c,d,:,1], statorCamber[c,d,:,2])).T     
        if max(statorCamberCyl[c,d,:,0]) - min(statorCamberCyl[c,d,:,0]) > np.pi:
            statorCamberCyl[c,d,:,:] = np.array(cart2pol2(statorCamber[c,d,:,0], statorCamber[c,d,:,1], statorCamber[c,d,:,2])).T

rotorCamberCyl  = np.zeros([Nbr, rSections, Nr, 3])
for cc in range(Nbr):
    for dd in range(rSections):
        rotorCamberCyl[cc,dd,:,:] = np.array(cart2pol(rotorCamber[cc,dd,:,0], rotorCamber[cc,dd,:,1], rotorCamber[cc,dd,:,2])).T    
        if max(rotorCamberCyl[cc,dd,:,0]) - min(rotorCamberCyl[cc,dd,:,0]) > np.pi:
            rotorCamberCyl[cc,dd,:,:] = np.array(cart2pol2(rotorCamber[cc,dd,:,0], rotorCamber[cc,dd,:,1], rotorCamber[cc,dd,:,2])).T
   
#%%
rotorCamberCylSorted, rotorCamberSorted, _ = reOrderBlades(rotorCamberCyl, rotorCamber)
statorCamberCylSorted, statorCamberSorted, _ = reOrderBlades(statorCamberCyl, statorCamber)

rotorPitch = computePitch(rotorCamberCylSorted, rotorCamberSorted, Nbr, rSections, Nr)
statorPitch = computePitch(statorCamberCylSorted, statorCamberSorted, Nbs, sSections, Ns)

np.savetxt(filePath + "/rotorBeta.txt", rotorPitch, delimiter=",")
np.savetxt(filePath + "/statorBeta.txt", statorPitch, delimiter=",")

























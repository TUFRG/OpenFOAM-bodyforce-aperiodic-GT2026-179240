#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 04:42:58 2025

@author: adekola
"""
import os
import argparse
import numpy as np 
from scipy.interpolate import interp1d, CubicSpline
#import matplotlib.pyplot as plt
import getCamber as gC
#%%
def rMatrix(points, R):
    r11,r12,r13,r21,r22,r23,r31,r32,r33 = R
    R = np.array([
        [ r11, r12, r13],
        [ r21, r22, r23],
        [r31, r32, r33]])
    return points @ R 

def vectorRot3D (x,y,z,theta):
    import numpy as np
    xRot = x*np.cos(theta) - y*np.sin(theta)
    yRot = x*np.sin(theta) + y*np.cos(theta)
    return xRot, yRot, z

def cart2pol(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    z = z
    return theta, r, z

def cart2pol2(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    theta = np.unwrap(theta)
    z = z
    return theta, r, z

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
    "Nr", 
    type=int,
    help="number of point on rotor blade profile"
)
parser.add_argument(
    "Ns", 
    type=int,
    help="number of points stator blade profile"
)
parser.add_argument(
    "periodicOrAperiodic", 
    type=int,
    help="periodic (0) or nonPeriodic (1)"
)
parser.add_argument(
    "--RmatrixCoeff",
    type=float,
    nargs=9,
    default=[1,0,0,0,1,0,0,0,1],
    help="9 coefficients for rotation matrix (r11 r12 r13 r21 r22 r23 r31 r32 r33)"
)

# Parse arguments
args = parser.parse_args()

# Use parsed arguments
Nbr = args.Nbr  # Total number of blades in the annulus for rotor
Nbs = args.Nbs  # Total number of blades in the annulus for stator
rSections = args.rSections  # Number of rotor blade profiles
sSections = args.sSections  # Number of stator blade profiles
Nr = args.Nr  # Number of points on rotor blade profiles
Ns = args.Ns  # Number of points on stator blade profiles
periodicOrAperiodic = args.periodicOrAperiodic  # if periodic select 0 otherwise select 1
RmatrixCoeff = args.RmatrixCoeff  # Define based on your machine axis

Ncr = 101
Ncs = 101
thetaR = np.linspace(0, 2*np.pi,Nbr, endpoint=False)
thetaS = np.linspace(0, 2*np.pi,Nbs, endpoint=False)


rBladeData = np.zeros([Nbr, rSections, Nr, 3])
sBladeData = np.zeros([Nbs, sSections, Ns, 3])

rBlade = np.zeros((rSections, Nr, 3))
sBlade = np.zeros((sSections, Ns, 3))
if periodicOrAperiodic == 0:
    for a in range(rSections):
        rBlade[a,:] = np.loadtxt('../rawData/periodic/bladeSurface/rotor/blade{}.txt'.format(a), delimiter=',')
    for b in range(sSections):
        sBlade[b,:] = np.loadtxt('../rawData/periodic/bladeSurface/stator/blade{}.txt'.format(b), delimiter=',')
else:
    for a in range(rSections):
       rBlade[a,:] = np.loadtxt('../rawData/periodic/bladeSurface/rotor/blade{}.txt'.format(a), delimiter=',')   
    
for c in range(Nbr):
    for d in range(rSections):
        if periodicOrAperiodic == 0:
            abcBlade = rMatrix(rBlade[d,:], RmatrixCoeff)
            rBladeData[c,d,:] = np.array(vectorRot3D(abcBlade[:,0], abcBlade[:,1], abcBlade[:,2],thetaR[c])).T
            inputPath = '../processedData/periodic/rotor/blade{}'.format(c)
            if not os.path.exists(inputPath):
                os.makedirs(inputPath)
                print(inputPath)
            np.savetxt(inputPath + '/blade{}.txt'.format(d),rBladeData[c,d,:] , delimiter=',')
            # np.savetxt(inputPath + '/camber{}.txt'.format(b),camberData[a,b,:] , delimiter=',')
        else:
            abcBlade = rMatrix(rBlade[d,:], RmatrixCoeff)
            rBladeData[c,d,:] = np.array(vectorRot3D(abcBlade[:,0], abcBlade[:,1], abcBlade[:,2],thetaR[c])).T
            inputPath = '../processedData/nonPeriodic/rotor/blade{}'.format(c)
            if not os.path.exists(inputPath):
                os.makedirs(inputPath)
            np.savetxt(inputPath + '/blade{}.txt'.format(d),rBladeData[c,d,:] , delimiter=',')
            

for c in range(Nbs):
    for d in range(sSections):
        if periodicOrAperiodic == 0:
            abcBlade = rMatrix(sBlade[d,:], RmatrixCoeff)
            sBladeData[c,d,:] = np.array(vectorRot3D(abcBlade[:,0], abcBlade[:,1], abcBlade[:,2],thetaS[c])).T
            inputPath = '../processedData/periodic/stator/blade{}'.format(c)
            if not os.path.exists(inputPath):
                os.makedirs(inputPath)
            np.savetxt(inputPath + '/blade{}.txt'.format(d),sBladeData[c,d,:] , delimiter=',')
        else:
            nonAxiBlade = np.loadtxt('../rawData/nonPeriodic/bladeSurface/stator/blade{}/blade{}.txt'.format(c,d), delimiter=',')
            inputPath = '../processedData/nonPeriodic/stator/blade{}'.format(c)
            if not os.path.exists(inputPath):
                os.makedirs(inputPath)
            abcBlade = rMatrix( nonAxiBlade, RmatrixCoeff)
            sBladeData[c,d,:] = abcBlade
            np.savetxt(inputPath + '/blade{}.txt'.format(d),sBladeData[c,d,:] , delimiter=',')

#%% Get camber
rCamberData = np.zeros([Nbr, rSections, Ncr, 3])
sCamberData = np.zeros([Nbs, sSections, Ncs, 3])

for e in range(Nbr):
    for f in range(rSections):
        rCamber, _ = gC.getCamber(rBladeData[0,f,:], Ncr-1)
        if periodicOrAperiodic == 0:
            rCamberData[e,f,:] = np.array(vectorRot3D(rCamber[:,0], rCamber[:,1], rCamber[:,2],thetaR[e])).T
            inputPath = '../processedData/periodic/rotor/blade{}'.format(e)
            np.savetxt(inputPath + '/camber{}.txt'.format(f),rCamberData[e,f,:] , delimiter=',')
        else:
            rCamberData[e,f,:] = np.array(vectorRot3D(rCamber[:,0], rCamber[:,1], rCamber[:,2],thetaR[e])).T
            inputPath = '../processedData/nonPeriodic/rotor/blade{}'.format(e)
            np.savetxt(inputPath + '/camber{}.txt'.format(f), rCamberData[e,f,:] , delimiter=',')

for e in range(Nbs):
    for f in range(sSections):
        if periodicOrAperiodic == 0:
            sCamber, _ = gC.getCamber(sBladeData[0,f,:], Ncs-1)
            sCamberData[e,f,:] = np.array(vectorRot3D(sCamber[:,0], sCamber[:,1], sCamber[:,2],thetaS[e])).T
            inputPath = '../processedData/periodic/stator/blade{}'.format(e)
            np.savetxt(inputPath + '/camber{}.txt'.format(f),sCamberData[e,f,:] , delimiter=',')
        else:
            sCamber, _ = gC.getCamber(sBladeData[e,f,:], Ncs-1)
            inputPath = '../processedData/nonPeriodic/stator/blade{}'.format(e)
            np.savetxt(inputPath + '/camber{}.txt'.format(f), sCamber , delimiter=',')


#%% Adjust nonPeriodic cases to fit the meridional extent
#load the streamline data
def computeArcLength(curve):
    """
    Compute meridional arc length of a 3D curve (using r–z plane)
    """
    r = np.sqrt(curve[:, 0]**2 + curve[:, 1]**2)
    z = curve[:, 2]
    meridional = np.column_stack([r, z])
    diffs = np.diff(meridional, axis=0)
    seg_lengths = np.linalg.norm(diffs, axis=1)
    fracSum = np.concatenate(([0.0], np.cumsum(seg_lengths)))
    frac = fracSum/fracSum[-1]
    return frac, fracSum

def scaleMeridionalExtent(curveA, curveB):
    '''
    This will basically make CurveA and CurveB follow the same meridional Path and meridional extent.
    That is the only way to guarantee that the same grid is valid for both cases.
    '''
    # Compute fractional meridional arc length for both curves
    fracA, _ = computeArcLength(curveA)
    fracB, _ = computeArcLength(curveB)
    # Extract theta from A
    thetaA = np.arctan2(curveA[:, 1], curveA[:, 0])
    # Extract r, z from B
    rB = np.sqrt(curveB[:, 0]**2 + curveB[:, 1]**2)
    zB = curveB[:, 2]
    # Interpolate B's r,z onto A's meridional fractions
    # This gives us B's geometry at the same meridional stations as A
    # rB_interp = np.interp(fracA, fracB, rB)
    # zB_interp = np.interp(fracA, fracB, zB)
    rB_interp = CubicSpline(fracB, rB)(fracA)
    zB_interp = CubicSpline(fracB, zB)(fracA)
    # Reconstruct 3D curve using B's r,z and A's theta
    x = rB_interp * np.cos(thetaA)
    y = rB_interp * np.sin(thetaA)
    z = zB_interp
    
    return np.column_stack((x, y, z))

if periodicOrAperiodic != 0:
    baselineBlade = np.zeros((sSections, Ncr, 3))
    for g in range(sSections):
        baselineBlade[g] = np.loadtxt('../processedData/periodic/stator/blade0/camber{}.txt'.format(g), delimiter=',')
    for h in range(Nbs):
        bladeData = np.zeros((sSections, Ncr, 3))
        for i in range(sSections):
            bladeData[i] = np.loadtxt('../processedData/nonPeriodic/stator/blade{}/camber{}.txt'.format(h,i), delimiter=',')
        newBladeData = np.zeros((sSections, Ncr, 3))
        for j in range(sSections):
            mapped = scaleMeridionalExtent(bladeData[j], baselineBlade[j])
            newBladeData[j] = mapped
            np.savetxt('../processedData/nonPeriodic/stator/blade{}/camber{}.txt'.format(h,j), mapped , delimiter=',')










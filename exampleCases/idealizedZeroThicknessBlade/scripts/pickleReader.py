#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 23:39:26 2025

@author: adekola
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp1d, CubicSpline

#%%
filePath = '/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/pickleFiles.curved/'

def lossCoeff(Pin, Pout, Pstat):
    loss = (Pin - Pout)/(Pout - Pstat)
    return loss 

def vectorRotP(x,y,P,theta):
    xRot = ((x-P[0])*np.cos(theta) - (P[1]-y)*np.sin(theta)) + P[0]
    yRot = P[1] -((x-P[0])*np.sin(theta) + (P[1]-y)*np.cos(theta))
    return xRot, yRot

def massAvgQty(Pt, Uz, Nt, Nr, x, y, theta, rho=1):
    # Rest of your code remains the same
    radius = np.sqrt(x**2 + y**2)
    hub = min(radius)*1.001
    cas = max(radius)*0.999
    spanFrac = np.linspace(0,1,Nr)
    radi = hub + spanFrac*(cas - hub)
    # theta = np.linspace(0, 2*np.pi, Nt)
    index = 0
    data = np.zeros([Nr*Nt, 2])
    X = np.zeros((Nt, Nr))
    Y = np.zeros((Nt, Nr))
    for b in range(Nr):
        for c in range(Nt):   
            data[index] = (radi[b] * np.sin(theta[c])),(radi[b] * np.cos(theta[c]))
            X[c,b] = (radi[b] * np.sin(theta[c]))
            Y[c,b] = (radi[b] * np.cos(theta[c])) 
            index +=1
    Area = np.zeros(Nr*Nt) #A=0.5*​(router2​−rinner2​)Δθ
    dr = np.diff(radi)[0]  # Radial spacing
    dtheta = np.diff(theta)[0]  # Angular spacing
    count = 0
    for c in range(Nr):
        for d in range(Nt):
            r = radi[c]
            # Define cell boundaries
            rInner = r - 0.5 * dr
            rOuter = r + 0.5 * dr
            # Handle boundary cases
            if c == 0:  # Inner boundary
                rInner = hub
            elif c == Nr - 1:  # Outer boundary  
                rOuter = cas       
            # Ensure positive radii
            rInner = max(rInner, 0)
            rOuter = max(rOuter, rInner + 1e-10)      
            # Correct annular sector area
            Area[count] = 0.5 * (rOuter**2 - rInner**2) * dtheta
            count += 1
    print(Area.shape)
    newPt = griddata((x,y), Pt, (X, Y), method='cubic')
    newUz = griddata((x,y), Uz, (X, Y), method='cubic')
    massAvgPt = np.zeros(Nr)
    for i in range(Nr):
        AvgPt = np.zeros(Nt)
        mDot= np.zeros(Nt)
        for j in range(Nt):
            # print(f"j={j}, i={i}")
            # print(f"newPt[j,i] shape: {np.shape(newPt[j,i])}, value: {newPt[j,i]}")
            # print(f"Area[j+Nt*i] shape: {np.shape(Area[j+Nt*i])}, value: {Area[j+Nt*i]}")
            # print(f"newUz[j,i] shape: {np.shape(newUz[j,i])}, value: {newUz[j,i]}")
            # print(f"newRho[j,i] shape: {np.shape(newRho[j,i])}, value: {newRho[j,i]}")
            AvgPt[j] = newPt[j,i]*Area[j+Nt*i]*newUz[j,i]
            mDot[j] = Area[j+Nt*i]*newUz[j,i]
        massAvgPt[i] = sum(AvgPt) / sum(mDot)
    mass_flux = np.zeros(Nt)
    for j in range(Nt):
        sum_mdot = 0
        sum_area = 0
        for i in range(Nr):
            idx = j + Nt * i
            mdot = rho * newUz[j, i] * Area[idx]
            sum_area += Area[idx]
            sum_mdot += mdot
        if sum_area > 0:
            mass_flux[j] = sum_mdot / sum_area  # Area-weighted average mass flux (ρ * Uz_avg)
        else:
            mass_flux[j] = 0    
    return massAvgPt, newPt, mass_flux, spanFrac

def densifyCurve(points, n_points_new, distribution='both'):
    if n_points_new <= len(points):
        return points[:n_points_new]
    points = np.atleast_2d(points)
    if points.shape[0] == 1:
        points = points.T
    # Step 1: Create parametric representation using cumulative distance
    diffs = np.diff(points, axis=0)
    segment_lengths = np.linalg.norm(diffs, axis=1)
    cumulative_length = np.concatenate([[0], np.cumsum(segment_lengths)])
    # Step 2: Fit cubic splines for each dimension
    n_dims = points.shape[1]
    splines = []
    for dim in range(n_dims):
        # splines.append(CubicSpline(cumulative_length, points[:, dim]))
        splines.append(interp1d(cumulative_length, points[:, dim]))
    # Step 3: Create parameter values with specified distribution
    total_length = cumulative_length[-1]
    if distribution == 'both':
        # Cosine clustering at both ends (like airfoil distribution)
        theta = np.linspace(0, np.pi, n_points_new)
        # Map from [0, π] to [0, 1], then to [0, total_length]
        xi = (1 - np.cos(theta)) / 2  # Maps to [0, 1]
        new_params = xi * total_length
    elif distribution == 'TE':
        # Clustering at start (trailing edge)
        theta = np.linspace(0, np.pi/2, n_points_new)
        xi = np.sin(theta)  # Maps to [0, 1] with clustering at start
        new_params = xi * total_length
    elif distribution == 'LE':
        # Clustering at end (leading edge)
        theta = np.linspace(0, np.pi/2, n_points_new)
        xi = 1 - np.cos(theta)  # Maps to [0, 1] with clustering at end
        new_params = xi * total_length
    elif distribution == 'uniform':
        # Uniform distribution (original behavior)
        new_params = np.linspace(0, total_length, n_points_new)
        
    else:
        raise ValueError(f"Unknown distribution: {distribution}")
    # Step 4: Evaluate the splines at all parameter values
    new_coords = []
    for dim in range(n_dims):
        new_coords.append(splines[dim](new_params))
    return np.column_stack(new_coords)

Nt = 90#360
Nr = 5
theta = np.linspace(0, 2*np.pi, Nt)
thetaB = np.linspace(0, 2*np.pi,Nt, endpoint=True)
# meanTheta = 360/Nt
# A = 0.5 #amplitude  ratio of max to min is 1.25
# perturbation = A * meanTheta * np.sin(thetaB)
# varyingPitch = thetaB + perturbation

#%%
with open(filePath + '/bodyForceBaseline.pkl', 'rb') as f:
    bFBaseline = pickle.load(f)  

bFBaselinePtIn = bFBaseline['surfaceMassAvgInlet']['totalP']
bFBaselineAlphaFlowIn = bFBaseline['surfaceMassAvgInlet']['alphaFlow']
bFBaselineUIn = bFBaseline['surfaceMassAvgInlet']['Uz']
xbFBaselineIn = bFBaseline['surfaceMassAvgInlet']['y']
ybFBaselineIn = bFBaseline['surfaceMassAvgInlet']['x']
bFBaselineUyIn = bFBaseline['surfaceMassAvgInlet']['Utangent']
bFBaselineUzIn = bFBaseline['surfaceMassAvgInlet']['Uaxial']
bFBaselinePtBarIn, bfBaselineSurfIn, _,spanFrac = massAvgQty(bFBaselinePtIn, bFBaselineUIn, Nt, Nr, xbFBaselineIn, ybFBaselineIn, theta)
bFBaselineAlphaFlowBarIn, bfBaselineAplhaFlowSurfIn, _,_ = massAvgQty(bFBaselineAlphaFlowIn, bFBaselineUIn, Nt, Nr, xbFBaselineIn, ybFBaselineIn, theta)
# bFBaselineUyBarIn, bFBaselineUySurfIn, _,_ = massAvgQty(bFBaselineUyIn, bFBaselineUIn, Nt, Nr, xbFBaselineIn, ybFBaselineIn, theta)
# bFBaselineUzBarIn, bFBaselineUzSurfIn,bFBaselineMassFluxIn, _ = massAvgQty(bFBaselineUzIn, bFBaselineUIn, Nt, Nr, xbFBaselineIn, ybFBaselineIn, theta)

bFBaselinePtRotIn = bFBaseline['surfaceMassAvgRotIn']['totalP']
bFBaselineAlphaFlowRotIn = bFBaseline['surfaceMassAvgRotIn']['alphaFlow']
bFBaselineURotIn = bFBaseline['surfaceMassAvgRotIn']['Uz']
xbFBaselineRotIn = bFBaseline['surfaceMassAvgRotIn']['y']
ybFBaselineRotIn = bFBaseline['surfaceMassAvgRotIn']['x']
bFBaselineUyRotIn = bFBaseline['surfaceMassAvgRotIn']['Utangent']
bFBaselineUzRotIn = bFBaseline['surfaceMassAvgRotIn']['Uaxial']
bFBaselineAlphaFlowBarRotIn, bfBaselineAplhaFlowSurfRotIn, _,_ = massAvgQty(bFBaselineAlphaFlowRotIn, bFBaselineURotIn, Nt, Nr, xbFBaselineRotIn, ybFBaselineRotIn, theta)
bFBaselinePtBarRotIn, bfBaselineSurfRotIn, _,_ = massAvgQty(bFBaselinePtRotIn, bFBaselineURotIn, Nt, Nr, xbFBaselineRotIn, ybFBaselineRotIn, theta)
# bFBaselineUyBarRotIn, bFBaselineUySurfRotIn, _,_ = massAvgQty(bFBaselineUyRotIn, bFBaselineURotIn, Nt, Nr, xbFBaselineRotIn, ybFBaselineRotIn, theta)
# bFBaselineUzBarRotIn, bFBaselineUzSurfRotIn,bFBaselineMassFluxRotIn, _ = massAvgQty(bFBaselineUzRotIn, bFBaselineURotIn, Nt, Nr, xbFBaselineRotIn, ybFBaselineRotIn, theta)

bFBaselinePtRotOut = bFBaseline['surfaceMassAvgRotOut']['p']
bFBaselineAlphaFlowRotOut = bFBaseline['surfaceMassAvgRotOut']['alphaFlow']
bFBaselineURotOut = bFBaseline['surfaceMassAvgRotOut']['Uz']
xbFBaselineRotOut = bFBaseline['surfaceMassAvgRotOut']['y']
ybFBaselineRotOut = bFBaseline['surfaceMassAvgRotOut']['x']
bFBaselineUyRotOut = bFBaseline['surfaceMassAvgRotOut']['Utangent']
bFBaselineUzRotOut = bFBaseline['surfaceMassAvgRotOut']['Uaxial']
bFBaselinePtBarRotOut, bfBaselineSurfRotOut, _,_ = massAvgQty(bFBaselinePtRotOut, bFBaselineURotOut, Nt, Nr, xbFBaselineRotOut, ybFBaselineRotOut, theta)
bFBaselineAlphaFlowBarRotOut, bfBaselineAplhaFlowSurfRotOut, _,_ = massAvgQty(bFBaselineAlphaFlowRotOut, bFBaselineURotOut, Nt, Nr, xbFBaselineRotOut, ybFBaselineRotOut, theta)
# bFBaselineUyBarRotOut, bFBaselineUySurfRotOut, _,_ = massAvgQty(bFBaselineUyRotOut, bFBaselineURotOut, Nt, Nr, xbFBaselineRotOut, ybFBaselineRotOut, theta)
# bFBaselineUzBarRotOut, bFBaselineUzSurfRotOut,bFBaselineMassFluxRotOut, _ = massAvgQty(bFBaselineUzRotOut, bFBaselineURotOut, Nt, Nr, xbFBaselineRotOut, ybFBaselineRotOut, theta)

bFBaselinePtOut = bFBaseline['surfaceMassAvgOutlet']['p']
bFBaselineAlphaFlowOut = bFBaseline['surfaceMassAvgOutlet']['alphaFlow']
bFBaselineUOut = bFBaseline['surfaceMassAvgOutlet']['Uz']
xbFBaselineOut = bFBaseline['surfaceMassAvgOutlet']['y']
ybFBaselineOut = bFBaseline['surfaceMassAvgOutlet']['x']
bFBaselineUyOut = bFBaseline['surfaceMassAvgOutlet']['Utangent']
bFBaselineUzOut = bFBaseline['surfaceMassAvgOutlet']['Uaxial']
bFBaselinePtBarOut, bfBaselineSurfOut, _,spanFrac = massAvgQty(bFBaselinePtOut, bFBaselineUOut, Nt, Nr, xbFBaselineOut, ybFBaselineOut, theta)
bFBaselineAlphaFlowBarOut, bfBaselineAplhaFlowSurfOut, _,_ = massAvgQty(bFBaselineAlphaFlowOut, bFBaselineUOut, Nt, Nr, xbFBaselineOut, ybFBaselineOut, theta)
# bFBaselineUyBarOut, bFBaselineUySurfOut, _,_ = massAvgQty(bFBaselineUyOut, bFBaselineUOut, Nt, Nr, xbFBaselineOut, ybFBaselineOut, theta)
# bFBaselineUzBarOut, bFBaselineUzSurfOut,bFBaselineMassFluxOut, _ = massAvgQty(bFBaselineUzOut, bFBaselineUOut, Nt, Nr, xbFBaselineOut, ybFBaselineOut, theta)

bFBaselineLossCoeffOut = lossCoeff(bfBaselineSurfIn, bfBaselineSurfRotOut, 0)
bFBaselineLossCoeffOut = lossCoeff(bfBaselineSurfIn, bfBaselineSurfOut, 0)

with open(filePath + '/bodyForceVaryingPitch.pkl', 'rb') as f:
    bFVaryingPitch = pickle.load(f)  

bFVaryingPitchPtIn = bFVaryingPitch['surfaceMassAvgInlet']['totalP']
bFVaryingPitchAlphaFlowIn = bFVaryingPitch['surfaceMassAvgInlet']['alphaFlow']
bFVaryingPitchUIn = bFVaryingPitch['surfaceMassAvgInlet']['Uz']
xbFVaryingPitchIn = bFVaryingPitch['surfaceMassAvgInlet']['y']
ybFVaryingPitchIn = bFVaryingPitch['surfaceMassAvgInlet']['x']
bFVaryingPitchUyIn = bFVaryingPitch['surfaceMassAvgInlet']['Utangent']
bFVaryingPitchUzIn = bFVaryingPitch['surfaceMassAvgInlet']['Uaxial']
bFVaryingPitchPtBarIn, bfVaryingPitchSurfIn, _,spanFrac = massAvgQty(bFVaryingPitchPtIn, bFVaryingPitchUIn, Nt, Nr, xbFVaryingPitchIn, ybFVaryingPitchIn, theta)
bFVaryingPitchAlphaFlowBarIn, bfVaryingPitchAplhaFlowSurfIn, _,_ = massAvgQty(bFVaryingPitchAlphaFlowIn, bFVaryingPitchUIn, Nt, Nr, xbFVaryingPitchIn, ybFVaryingPitchIn, theta)
# bFVaryingPitchUyBarIn, bFVaryingPitchUySurfIn, _,_ = massAvgQty(bFVaryingPitchUyIn, bFVaryingPitchUIn, Nt, Nr, xbFVaryingPitchIn, ybFVaryingPitchIn, theta)
# bFVaryingPitchUzBarIn, bFVaryingPitchUzSurfIn,bFVaryingPitchMassFluxIn, _ = massAvgQty(bFVaryingPitchUzIn, bFVaryingPitchUIn, Nt, Nr, xbFVaryingPitchIn, ybFVaryingPitchIn, theta)

bFVaryingPitchPtRotIn = bFVaryingPitch['surfaceMassAvgRotIn']['totalP']
bFVaryingPitchAlphaFlowRotIn = bFVaryingPitch['surfaceMassAvgRotIn']['alphaFlow']
bFVaryingPitchURotIn = bFVaryingPitch['surfaceMassAvgRotIn']['Uz']
xbFVaryingPitchRotIn = bFVaryingPitch['surfaceMassAvgRotIn']['y']
ybFVaryingPitchRotIn = bFVaryingPitch['surfaceMassAvgRotIn']['x']
bFVaryingPitchUyRotIn = bFVaryingPitch['surfaceMassAvgRotIn']['Utangent']
bFVaryingPitchUzRotIn = bFVaryingPitch['surfaceMassAvgRotIn']['Uaxial']
bFVaryingPitchAlphaFlowBarRotIn, bfVaryingPitchAplhaFlowSurfRotIn, _,_ = massAvgQty(bFVaryingPitchAlphaFlowRotIn, bFVaryingPitchURotIn, Nt, Nr, xbFVaryingPitchRotIn, ybFVaryingPitchRotIn, theta)
bFVaryingPitchPtBarRotIn, bfVaryingPitchSurfRotIn, _,_ = massAvgQty(bFVaryingPitchPtRotIn, bFVaryingPitchURotIn, Nt, Nr, xbFVaryingPitchRotIn, ybFVaryingPitchRotIn, theta)
# bFVaryingPitchUyBarRotIn, bFVaryingPitchUySurfRotIn, _,_ = massAvgQty(bFVaryingPitchUyRotIn, bFVaryingPitchURotIn, Nt, Nr, xbFVaryingPitchRotIn, ybFVaryingPitchRotIn, theta)
# bFVaryingPitchUzBarRotIn, bFVaryingPitchUzSurfRotIn,bFVaryingPitchMassFluxRotIn, _ = massAvgQty(bFVaryingPitchUzRotIn, bFVaryingPitchURotIn, Nt, Nr, xbFVaryingPitchRotIn, ybFVaryingPitchRotIn, theta)

bFVaryingPitchPtRotOut = bFVaryingPitch['surfaceMassAvgRotOut']['totalP']
bFVaryingPitchAlphaFlowRotOut = bFVaryingPitch['surfaceMassAvgRotOut']['alphaFlow']
bFVaryingPitchURotOut = bFVaryingPitch['surfaceMassAvgRotOut']['Uz']
xbFVaryingPitchRotOut = bFVaryingPitch['surfaceMassAvgRotOut']['y']
ybFVaryingPitchRotOut = bFVaryingPitch['surfaceMassAvgRotOut']['x']
bFVaryingPitchUyRotOut = bFVaryingPitch['surfaceMassAvgRotOut']['Utangent']
bFVaryingPitchUzRotOut = bFVaryingPitch['surfaceMassAvgRotOut']['Uaxial']
bFVaryingPitchPtBarRotOut, bfVaryingPitchSurfRotOut, _,_ = massAvgQty(bFVaryingPitchPtRotOut, bFVaryingPitchURotOut, Nt, Nr, xbFVaryingPitchRotOut, ybFVaryingPitchRotOut, theta)
bFVaryingPitchAlphaFlowBarRotOut, bfVaryingPitchAplhaFlowSurfRotOut, _,_ = massAvgQty(bFVaryingPitchAlphaFlowRotOut, bFVaryingPitchURotOut, Nt, Nr, xbFVaryingPitchRotOut, ybFVaryingPitchRotOut, theta)
# bFVaryingPitchUyBarRotOut, bFVaryingPitchUySurfRotOut, _,_ = massAvgQty(bFVaryingPitchUyRotOut, bFVaryingPitchURotOut, Nt, Nr, xbFVaryingPitchRotOut, ybFVaryingPitchRotOut, theta)
# bFVaryingPitchUzBarRotOut, bFVaryingPitchUzSurfRotOut,bFVaryingPitchMassFluxRotOut, _ = massAvgQty(bFVaryingPitchUzRotOut, bFVaryingPitchURotOut, Nt, Nr, xbFVaryingPitchRotOut, ybFVaryingPitchRotOut, theta)

bFVaryingPitchPtOut = bFVaryingPitch['surfaceMassAvgOutlet']['totalP']
bFVaryingPitchAlphaFlowOut = bFVaryingPitch['surfaceMassAvgOutlet']['alphaFlow']
bFVaryingPitchUOut = bFVaryingPitch['surfaceMassAvgOutlet']['Uz']
xbFVaryingPitchOut = bFVaryingPitch['surfaceMassAvgOutlet']['y']
ybFVaryingPitchOut = bFVaryingPitch['surfaceMassAvgOutlet']['x']
bFVaryingPitchUyOut = bFVaryingPitch['surfaceMassAvgOutlet']['Utangent']
bFVaryingPitchUzOut = bFVaryingPitch['surfaceMassAvgOutlet']['Uaxial']
bFVaryingPitchPtBarOut, bfVaryingPitchSurfOut, _,spanFrac = massAvgQty(bFVaryingPitchPtOut, bFVaryingPitchUOut, Nt, Nr, xbFVaryingPitchOut, ybFVaryingPitchOut, theta)
bFVaryingPitchAlphaFlowBarOut, bfVaryingPitchAplhaFlowSurfOut, _,_ = massAvgQty(bFVaryingPitchAlphaFlowOut, bFVaryingPitchUOut, Nt, Nr, xbFVaryingPitchOut, ybFVaryingPitchOut, theta)
# bFVaryingPitchUyBarOut, bFVaryingPitchUySurfOut, _,_ = massAvgQty(bFVaryingPitchUyOut, bFVaryingPitchUOut, Nt, Nr, xbFVaryingPitchOut, ybFVaryingPitchOut, theta)
# bFVaryingPitchUzBarOut, bFVaryingPitchUzSurfOut,bFVaryingPitchMassFluxOut, _ = massAvgQty(bFVaryingPitchUzOut, bFVaryingPitchUOut, Nt, Nr, xbFVaryingPitchOut, ybFVaryingPitchOut, theta)

bFVaryingPitchLossCoeffOut = lossCoeff(bfVaryingPitchSurfIn, bfVaryingPitchSurfRotOut, 0)
bFVaryingPitchLossCoeffOut = lossCoeff(bfVaryingPitchSurfIn, bfVaryingPitchSurfOut, 0)


with open(filePath + '/bodyForceStagger.pkl', 'rb') as f:
    bFVaryingStagger = pickle.load(f)  

bFVaryingStaggerPtIn = bFVaryingStagger['surfaceMassAvgInlet']['totalP']
bFVaryingStaggerAlphaFlowIn = bFVaryingStagger['surfaceMassAvgInlet']['alphaFlow']
bFVaryingStaggerUIn = bFVaryingStagger['surfaceMassAvgInlet']['Uz']
xbFVaryingStaggerIn = bFVaryingStagger['surfaceMassAvgInlet']['y']
ybFVaryingStaggerIn = bFVaryingStagger['surfaceMassAvgInlet']['x']
bFVaryingStaggerUyIn = bFVaryingStagger['surfaceMassAvgInlet']['Utangent']
bFVaryingStaggerUzIn = bFVaryingStagger['surfaceMassAvgInlet']['Uaxial']
bFVaryingStaggerPtBarIn, bfVaryingStaggerSurfIn, _,spanFrac = massAvgQty(bFVaryingStaggerPtIn, bFVaryingStaggerUIn, Nt, Nr, xbFVaryingStaggerIn, ybFVaryingStaggerIn, theta)
bFVaryingStaggerAlphaFlowBarIn, bfVaryingStaggerAplhaFlowSurfIn, _,_ = massAvgQty(bFVaryingStaggerAlphaFlowIn, bFVaryingStaggerUIn, Nt, Nr, xbFVaryingStaggerIn, ybFVaryingStaggerIn, theta)
# bFVaryingStaggerUyBarIn, bFVaryingStaggerUySurfIn, _,_ = massAvgQty(bFVaryingStaggerUyIn, bFVaryingStaggerUIn, Nt, Nr, xbFVaryingStaggerIn, ybFVaryingStaggerIn, theta)
# bFVaryingStaggerUzBarIn, bFVaryingStaggerUzSurfIn,bFVaryingStaggerMassFluxIn, _ = massAvgQty(bFVaryingStaggerUzIn, bFVaryingStaggerUIn, Nt, Nr, xbFVaryingStaggerIn, ybFVaryingStaggerIn, theta)

bFVaryingStaggerPtRotIn = bFVaryingStagger['surfaceMassAvgRotIn']['totalP']
bFVaryingStaggerAlphaFlowRotIn = bFVaryingStagger['surfaceMassAvgRotIn']['alphaFlow']
bFVaryingStaggerURotIn = bFVaryingStagger['surfaceMassAvgRotIn']['Uz']
xbFVaryingStaggerRotIn = bFVaryingStagger['surfaceMassAvgRotIn']['y']
ybFVaryingStaggerRotIn = bFVaryingStagger['surfaceMassAvgRotIn']['x']
bFVaryingStaggerUyRotIn = bFVaryingStagger['surfaceMassAvgRotIn']['Utangent']
bFVaryingStaggerUzRotIn = bFVaryingStagger['surfaceMassAvgRotIn']['Uaxial']
bFVaryingStaggerAlphaFlowBarRotIn, bfVaryingStaggerAplhaFlowSurfRotIn, _,_ = massAvgQty(bFVaryingStaggerAlphaFlowRotIn, bFVaryingStaggerURotIn, Nt, Nr, xbFVaryingStaggerRotIn, ybFVaryingStaggerRotIn, theta)
bFVaryingStaggerPtBarRotIn, bfVaryingStaggerSurfRotIn, _,_ = massAvgQty(bFVaryingStaggerPtRotIn, bFVaryingStaggerURotIn, Nt, Nr, xbFVaryingStaggerRotIn, ybFVaryingStaggerRotIn, theta)
# bFVaryingStaggerUyBarRotIn, bFVaryingStaggerUySurfRotIn, _,_ = massAvgQty(bFVaryingStaggerUyRotIn, bFVaryingStaggerURotIn, Nt, Nr, xbFVaryingStaggerRotIn, ybFVaryingStaggerRotIn, theta)
# bFVaryingStaggerUzBarRotIn, bFVaryingStaggerUzSurfRotIn,bFVaryingStaggerMassFluxRotIn, _ = massAvgQty(bFVaryingStaggerUzRotIn, bFVaryingStaggerURotIn, Nt, Nr, xbFVaryingStaggerRotIn, ybFVaryingStaggerRotIn, theta)

bFVaryingStaggerPtRotOut = bFVaryingStagger['surfaceMassAvgRotOut']['totalP']
bFVaryingStaggerAlphaFlowRotOut = bFVaryingStagger['surfaceMassAvgRotOut']['alphaFlow']
bFVaryingStaggerURotOut = bFVaryingStagger['surfaceMassAvgRotOut']['Uz']
xbFVaryingStaggerRotOut = bFVaryingStagger['surfaceMassAvgRotOut']['y']
ybFVaryingStaggerRotOut = bFVaryingStagger['surfaceMassAvgRotOut']['x']
bFVaryingStaggerUyRotOut = bFVaryingStagger['surfaceMassAvgRotOut']['Utangent']
bFVaryingStaggerUzRotOut = bFVaryingStagger['surfaceMassAvgRotOut']['Uaxial']
bFVaryingStaggerPtBarRotOut, bfVaryingStaggerSurfRotOut, _,_ = massAvgQty(bFVaryingStaggerPtRotOut, bFVaryingStaggerURotOut, Nt, Nr, xbFVaryingStaggerRotOut, ybFVaryingStaggerRotOut, theta)
bFVaryingStaggerAlphaFlowBarRotOut, bfVaryingStaggerAplhaFlowSurfRotOut, _,_ = massAvgQty(bFVaryingStaggerAlphaFlowRotOut, bFVaryingStaggerURotOut, Nt, Nr, xbFVaryingStaggerRotOut, ybFVaryingStaggerRotOut, theta)
# bFVaryingStaggerUyBarRotOut, bFVaryingStaggerUySurfRotOut, _,_ = massAvgQty(bFVaryingStaggerUyRotOut, bFVaryingStaggerURotOut, Nt, Nr, xbFVaryingStaggerRotOut, ybFVaryingStaggerRotOut, theta)
# bFVaryingStaggerUzBarRotOut, bFVaryingStaggerUzSurfRotOut,bFVaryingStaggerMassFluxRotOut, _ = massAvgQty(bFVaryingStaggerUzRotOut, bFVaryingStaggerURotOut, Nt, Nr, xbFVaryingStaggerRotOut, ybFVaryingStaggerRotOut, theta)

bFVaryingStaggerPtOut = bFVaryingStagger['surfaceMassAvgOutlet']['totalP']
bFVaryingStaggerAlphaFlowOut = bFVaryingStagger['surfaceMassAvgOutlet']['alphaFlow']
bFVaryingStaggerUOut = bFVaryingStagger['surfaceMassAvgOutlet']['Uz']
xbFVaryingStaggerOut = bFVaryingStagger['surfaceMassAvgOutlet']['y']
ybFVaryingStaggerOut = bFVaryingStagger['surfaceMassAvgOutlet']['x']
bFVaryingStaggerUyOut = bFVaryingStagger['surfaceMassAvgOutlet']['Utangent']
bFVaryingStaggerUzOut = bFVaryingStagger['surfaceMassAvgOutlet']['Uaxial']
bFVaryingStaggerPtBarOut, bfVaryingStaggerSurfOut, _,spanFrac = massAvgQty(bFVaryingStaggerPtOut, bFVaryingStaggerUOut, Nt, Nr, xbFVaryingStaggerOut, ybFVaryingStaggerOut, theta)
bFVaryingStaggerAlphaFlowBarOut, bfVaryingStaggerAplhaFlowSurfOut, _,_ = massAvgQty(bFVaryingStaggerAlphaFlowOut, bFVaryingStaggerUOut, Nt, Nr, xbFVaryingStaggerOut, ybFVaryingStaggerOut, theta)
# bFVaryingStaggerUyBarOut, bFVaryingStaggerUySurfOut, _,_ = massAvgQty(bFVaryingStaggerUyOut, bFVaryingStaggerUOut, Nt, Nr, xbFVaryingStaggerOut, ybFVaryingStaggerOut, theta)
# bFVaryingStaggerUzBarOut, bFVaryingStaggerUzSurfOut,bFVaryingStaggerMassFluxOut, _ = massAvgQty(bFVaryingStaggerUzOut, bFVaryingStaggerUOut, Nt, Nr, xbFVaryingStaggerOut, ybFVaryingStaggerOut, theta)

bFVaryingStaggerLossCoeffOut = lossCoeff(bfVaryingStaggerSurfIn, bfVaryingStaggerSurfRotOut, 0)
bFVaryingStaggerLossCoeffOut = lossCoeff(bfVaryingStaggerSurfIn, bfVaryingStaggerSurfOut, 0)

with open(filePath + '/bodyForceStaggerOldMethod.pkl', 'rb') as f:
    bFVaryingStaggerOld = pickle.load(f) 
    
bFVaryingStaggerOldPtRotOut = bFVaryingStaggerOld['surfaceMassAvgRotOut']['totalP']
bFVaryingStaggerOldAlphaFlowRotOut = bFVaryingStaggerOld['surfaceMassAvgRotOut']['alphaFlow']
bFVaryingStaggerOldURotOut = bFVaryingStaggerOld['surfaceMassAvgRotOut']['Uz']
xbFVaryingStaggerOldRotOut = bFVaryingStaggerOld['surfaceMassAvgRotOut']['y']
ybFVaryingStaggerOldRotOut = bFVaryingStaggerOld['surfaceMassAvgRotOut']['x']
bFVaryingStaggerOldUyRotOut = bFVaryingStaggerOld['surfaceMassAvgRotOut']['Utangent']
bFVaryingStaggerOldUzRotOut = bFVaryingStaggerOld['surfaceMassAvgRotOut']['Uaxial']
bFVaryingStaggerOldPtBarRotOut, bfVaryingStaggerOldSurfRotOut, _,_ = massAvgQty(bFVaryingStaggerOldPtRotOut, bFVaryingStaggerOldURotOut, Nt, Nr, xbFVaryingStaggerOldRotOut, ybFVaryingStaggerOldRotOut, theta)
bFVaryingStaggerOldAlphaFlowBarRotOut, bfVaryingStaggerOldAplhaFlowSurfRotOut, _,_ = massAvgQty(bFVaryingStaggerOldAlphaFlowRotOut, bFVaryingStaggerOldURotOut, Nt, Nr, xbFVaryingStaggerOldRotOut, ybFVaryingStaggerOldRotOut, theta)


#%%
# plt.plot(theta, bfBaselineAplhaFlowSurfRotOut[:,2], 'k', label='outlet')
# plt.plot(theta, bfBaselineAplhaFlowSurfIn[:,2], 'k--', label='inlet')
# plt.plot(theta, bfBaselineAplhaFlowSurfRotIn[:,2], 'k-.', label='inlet')
# plt.plot(theta, bfBaselineAplhaFlowSurfRotOut[:,2], 'k--.', label='inlet')


# plt.plot(theta, bfVaryingPitchAplhaFlowSurfIn[:,2], 'k--.', label='inlet')
# plt.plot(theta, bfVaryingPitchAplhaFlowSurfOut[:,2], 'k--', label='outlet')
# plt.plot(theta, bfVaryingPitchAplhaFlowSurfRotOut[:,2], 'k--.', label='TE')
# plt.plot(theta, bfVaryingPitchAplhaFlowSurfRotIn[:,2], 'k', label=' LE')

plt.plot(theta, bfVaryingStaggerAplhaFlowSurfOut[:,2], 'k--', label='outlet')
# plt.plot(theta, bfVaryingStaggerAplhaFlowSurfIn[:,2], 'k--', label='nonUniform')
plt.plot(theta, bfVaryingStaggerAplhaFlowSurfRotOut[:,2], 'k--.', label='TE')
# plt.plot(theta, bfVaryingStaggerAplhaFlowSurfRotIn[:,2], 'k', label=' LE')

plt.grid()
plt.legend(loc=0)
# plt.axis('equal')
plt.xlabel('circumferential direction')
plt.ylabel('Flow angle (deg)')
# plt.ylabel('Loss Coeff (Y)')
# plt.ylabel('stagger angle (deg)')
# plt.ylim(0,11)
# plt.savefig('/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/plot/varyingStaggerBFFlowAngle.png', dpi=500, bbox_inches='tight')
# plt.plot(theta, bfVaryingStaggerAlphaFlowSurf[:,0], 'b')
#%% Bladed

# with open(filePath + '/bladedStraightSlipWalls.pkl', 'rb') as f:
#     bLStraight = pickle.load(f)  

# bLStraightPt = bLStraight['surfaceMassAvgRotOut']['totalP']
# bLStraightAlphaFlow = bLStraight['surfaceMassAvgRotOut']['alphaFlow']
# bLStraightU = bLStraight['surfaceMassAvgRotOut']['Uz']
# xbLStraight = bLStraight['surfaceMassAvgRotOut']['x']
# ybLStraight = bLStraight['surfaceMassAvgRotOut']['y']

# bLStraightPtBar, blStraightSurf, spanFrac = massAvgQty(bLStraightPt, bLStraightU, Nt, Nr, xbLStraight, ybLStraight)
# bLStraightAlphaFlowBar, blStraightAplhaFlowSurf, _ = massAvgQty(bLStraightAlphaFlow, bLStraightU, Nt, Nr, xbLStraight, ybLStraight)
#%%
with open(filePath + '/bladedBaselineSlip.pkl', 'rb') as f:
# with open(filePath + '/bladedStraightSlipWalls.pkl', 'rb') as f:    
    bLBaselineSlip = pickle.load(f)  
    
bLBaselineSlipPtIn = bLBaselineSlip['surfaceMassAvgInlet']['totalP']
bLBaselineSlipAlphaFlowIn = bLBaselineSlip['surfaceMassAvgInlet']['alphaFlow']
bLBaselineSlipUIn = bLBaselineSlip['surfaceMassAvgInlet']['Uz']
xbLBaselineSlipIn = bLBaselineSlip['surfaceMassAvgInlet']['y']
ybLBaselineSlipIn = bLBaselineSlip['surfaceMassAvgInlet']['x']
bLBaselineSlipUyIn = bLBaselineSlip['surfaceMassAvgInlet']['Utangent']
bLBaselineSlipUzIn = bLBaselineSlip['surfaceMassAvgInlet']['Uaxial']

bLBaselineSlipPtBarIn, bLBaselineSlipSurfIn, _, stotalPanFrac = massAvgQty(bLBaselineSlipPtIn, bLBaselineSlipUIn, Nt, Nr, xbLBaselineSlipIn, ybLBaselineSlipIn, theta)
bLBaselineSlipAlphaFlowBarIn, bLBaselineSlipAlphaFlowSurfIn, _,_ = massAvgQty(bLBaselineSlipAlphaFlowIn, bLBaselineSlipUIn, Nt, Nr, xbLBaselineSlipIn, ybLBaselineSlipIn, theta)
# bLBaselineSlipUyBarIn, bLBaselineSlipUySurfIn, _,_ = massAvgQty(bLBaselineSlipUyIn, bLBaselineSlipUIn, Nt, Nr, xbLBaselineSlipIn, ybLBaselineSlipIn, theta)
# bLBaselineSlipUzBarIn, bLBaselineSlipUzSurfIn, bLBaselineSlipMassFluxIn, _ = massAvgQty(bLBaselineSlipUzIn, bLBaselineSlipUIn, Nt, Nr, xbLBaselineSlipIn, ybLBaselineSlipIn, theta)    

bLBaselineSlipPtRotIn = bLBaselineSlip['surfaceMassAvgRotIn']['totalP']
bLBaselineSlipAlphaFlowRotIn = bLBaselineSlip['surfaceMassAvgRotIn']['alphaFlow']
bLBaselineSlipURotIn = bLBaselineSlip['surfaceMassAvgRotIn']['Uz']
xbLBaselineSlipRotIn = bLBaselineSlip['surfaceMassAvgRotIn']['y']
ybLBaselineSlipRotIn = bLBaselineSlip['surfaceMassAvgRotIn']['x']
bLBaselineSlipUyRotIn = bLBaselineSlip['surfaceMassAvgRotIn']['Utangent']
bLBaselineSlipUzRotIn = bLBaselineSlip['surfaceMassAvgRotIn']['Uaxial']

bLBaselineSlipPtBarRotIn, bLBaselineSlipSurfRotIn, _, stotalPanFrac = massAvgQty(bLBaselineSlipPtRotIn, bLBaselineSlipURotIn, Nt, Nr, xbLBaselineSlipRotIn, ybLBaselineSlipRotIn, theta)
bLBaselineSlipAlphaFlowBarRotIn, bLBaselineSlipAlphaFlowSurfRotIn, _,_ = massAvgQty(bLBaselineSlipAlphaFlowRotIn, bLBaselineSlipURotIn, Nt, Nr, xbLBaselineSlipRotIn, ybLBaselineSlipRotIn, theta)
# bLBaselineSlipUyBarRotIn, bLBaselineSlipUySurfRotIn, _,_ = massAvgQty(bLBaselineSlipUyRotIn, bLBaselineSlipURotIn, Nt, Nr, xbLBaselineSlipRotIn, ybLBaselineSlipRotIn, theta)
# bLBaselineSlipUzBarRotIn, bLBaselineSlipUzSurfRotIn, bLBaselineSlipMassFluxRotIn, _ = massAvgQty(bLBaselineSlipUzRotIn, bLBaselineSlipURotIn, Nt, Nr, xbLBaselineSlipRotIn, ybLBaselineSlipRotIn, theta)

bLBaselineSlipPtRotOut = bLBaselineSlip['surfaceMassAvgRotOut']['p']
bLBaselineSlipAlphaFlowRotOut = bLBaselineSlip['surfaceMassAvgRotOut']['alphaFlow']
bLBaselineSlipURotOut = bLBaselineSlip['surfaceMassAvgRotOut']['Uz']
xbLBaselineSlipRotOut = bLBaselineSlip['surfaceMassAvgRotOut']['y']
ybLBaselineSlipRotOut = bLBaselineSlip['surfaceMassAvgRotOut']['x']
bLBaselineSlipUyRotOut = bLBaselineSlip['surfaceMassAvgRotOut']['Utangent']
bLBaselineSlipUzRotOut = bLBaselineSlip['surfaceMassAvgRotOut']['Uaxial']

bLBaselineSlipPtBarRotOut, bLBaselineSlipSurfRotOut, _, stotalPanFrac = massAvgQty(bLBaselineSlipPtRotOut, bLBaselineSlipURotOut, Nt, Nr, xbLBaselineSlipRotOut, ybLBaselineSlipRotOut, theta)
bLBaselineSlipAlphaFlowBarRotOut, bLBaselineSlipAlphaFlowSurfRotOut, _,_ = massAvgQty(bLBaselineSlipAlphaFlowRotOut, bLBaselineSlipURotOut, Nt, Nr, xbLBaselineSlipRotOut, ybLBaselineSlipRotOut, theta)
# bLBaselineSlipUyBarRotOut, bLBaselineSlipUySurfRotOut, _,_ = massAvgQty(bLBaselineSlipUyRotOut, bLBaselineSlipURotOut, Nt, Nr, xbLBaselineSlipRotOut, ybLBaselineSlipRotOut, theta)
# bLBaselineSlipUzBarRotOut, bLBaselineSlipUzSurfRotOut, bLBaselineSlipMassFluxRotOut, _ = massAvgQty(bLBaselineSlipUzRotOut, bLBaselineSlipURotOut, Nt, Nr, xbLBaselineSlipRotOut, ybLBaselineSlipRotOut, theta)


bLBaselineSlipPtOut = bLBaselineSlip['surfaceMassAvgOutlet']['p']
bLBaselineSlipAlphaFlowOut = bLBaselineSlip['surfaceMassAvgOutlet']['alphaFlow']
bLBaselineSlipUOut = bLBaselineSlip['surfaceMassAvgOutlet']['Uz']
xbLBaselineSlipOut = bLBaselineSlip['surfaceMassAvgOutlet']['y']
ybLBaselineSlipOut = bLBaselineSlip['surfaceMassAvgOutlet']['x']
bLBaselineSlipUyOut = bLBaselineSlip['surfaceMassAvgOutlet']['Utangent']
bLBaselineSlipUzOut = bLBaselineSlip['surfaceMassAvgOutlet']['Uaxial']

bLBaselineSlipPtBarOut, bLBaselineSlipSurfOut, _, _ = massAvgQty(bLBaselineSlipPtOut, bLBaselineSlipUOut, Nt, Nr, xbLBaselineSlipOut, ybLBaselineSlipOut, theta)
bLBaselineSlipAlphaFlowBarOut, bLBaselineSlipAlphaFlowSurfOut, _,_ = massAvgQty(bLBaselineSlipAlphaFlowOut, bLBaselineSlipUOut, Nt, Nr, xbLBaselineSlipOut, ybLBaselineSlipOut, theta)
# bLBaselineSlipUyBarOut, bLBaselineSlipUySurfOut, _,_ = massAvgQty(bLBaselineSlipUyOut, bLBaselineSlipUOut, Nt, Nr, xbLBaselineSlipOut, ybLBaselineSlipOut, theta)
# bLBaselineSlipUzBarOut, bLBaselineSlipUzSurfOut, bLBaselineSlipMassFluxOut, _ = massAvgQty(bLBaselineSlipUzOut, bLBaselineSlipUOut, Nt, Nr, xbLBaselineSlipOut, ybLBaselineSlipOut, theta)

#%%
# with open(filePath + '/bladedBaselineNoSlip.pkl', 'rb') as f:
# # with open(filePath + '/bladedStraightNoSlipWalls.pkl', 'rb') as f:    
#     bLBaselineNoSlip = pickle.load(f)  
    
# bLBaselineNoSlipPtIn = bLBaselineNoSlip['surfaceMassAvgInlet']['totalP']
# bLBaselineNoSlipAlphaFlowIn = bLBaselineNoSlip['surfaceMassAvgInlet']['alphaFlow']
# bLBaselineNoSlipUIn = bLBaselineNoSlip['surfaceMassAvgInlet']['Uz']
# xbLBaselineNoSlipIn = bLBaselineNoSlip['surfaceMassAvgInlet']['y']
# ybLBaselineNoSlipIn = bLBaselineNoSlip['surfaceMassAvgInlet']['x']
# bLBaselineNoSlipUyIn = bLBaselineNoSlip['surfaceMassAvgInlet']['Utangent']
# bLBaselineNoSlipUzIn = bLBaselineNoSlip['surfaceMassAvgInlet']['Uaxial']

# bLBaselineNoSlipPtBarIn, bLBaselineNoSlipSurfIn, _, stotalPanFrac = massAvgQty(bLBaselineNoSlipPtIn, bLBaselineNoSlipUIn, Nt, Nr, xbLBaselineNoSlipIn, ybLBaselineNoSlipIn, theta)
# bLBaselineNoSlipAlphaFlowBarIn, bLBaselineNoSlipAlphaFlowSurfIn, _,_ = massAvgQty(bLBaselineNoSlipAlphaFlowIn, bLBaselineNoSlipUIn, Nt, Nr, xbLBaselineNoSlipIn, ybLBaselineNoSlipIn, theta)
# bLBaselineNoSlipUyBarIn, bLBaselineNoSlipUySurfIn, _,_ = massAvgQty(bLBaselineNoSlipUyIn, bLBaselineNoSlipUIn, Nt, Nr, xbLBaselineNoSlipIn, ybLBaselineNoSlipIn, theta)
# bLBaselineNoSlipUzBarIn, bLBaselineNoSlipUzSurfIn, bLBaselineNoSlipMassFluxIn, _ = massAvgQty(bLBaselineNoSlipUzIn, bLBaselineNoSlipUIn, Nt, Nr, xbLBaselineNoSlipIn, ybLBaselineNoSlipIn, theta)    

# bLBaselineNoSlipPtRotIn = bLBaselineNoSlip['surfaceMassAvgRotIn']['totalP']
# bLBaselineNoSlipAlphaFlowRotIn = bLBaselineNoSlip['surfaceMassAvgRotIn']['alphaFlow']
# bLBaselineNoSlipURotIn = bLBaselineNoSlip['surfaceMassAvgRotIn']['Uz']
# xbLBaselineNoSlipRotIn = bLBaselineNoSlip['surfaceMassAvgRotIn']['y']
# ybLBaselineNoSlipRotIn = bLBaselineNoSlip['surfaceMassAvgRotIn']['x']
# bLBaselineNoSlipUyRotIn = bLBaselineNoSlip['surfaceMassAvgRotIn']['Utangent']
# bLBaselineNoSlipUzRotIn = bLBaselineNoSlip['surfaceMassAvgRotIn']['Uaxial']

# bLBaselineNoSlipPtBarRotIn, bLBaselineNoSlipSurfRotIn, _, stotalPanFrac = massAvgQty(bLBaselineNoSlipPtRotIn, bLBaselineNoSlipURotIn, Nt, Nr, xbLBaselineNoSlipRotIn, ybLBaselineNoSlipRotIn, theta)
# bLBaselineNoSlipAlphaFlowBarRotIn, bLBaselineNoSlipAlphaFlowSurfRotIn, _,_ = massAvgQty(bLBaselineNoSlipAlphaFlowRotIn, bLBaselineNoSlipURotIn, Nt, Nr, xbLBaselineNoSlipRotIn, ybLBaselineNoSlipRotIn, theta)
# bLBaselineNoSlipUyBarRotIn, bLBaselineNoSlipUySurfRotIn, _,_ = massAvgQty(bLBaselineNoSlipUyRotIn, bLBaselineNoSlipURotIn, Nt, Nr, xbLBaselineNoSlipRotIn, ybLBaselineNoSlipRotIn, theta)
# bLBaselineNoSlipUzBarRotIn, bLBaselineNoSlipUzSurfRotIn, bLBaselineNoSlipMassFluxRotIn, _ = massAvgQty(bLBaselineNoSlipUzRotIn, bLBaselineNoSlipURotIn, Nt, Nr, xbLBaselineNoSlipRotIn, ybLBaselineNoSlipRotIn, theta)

# bLBaselineNoSlipPtRotOut = bLBaselineNoSlip['surfaceMassAvgRotOut']['totalP']
# bLBaselineNoSlipAlphaFlowRotOut = bLBaselineNoSlip['surfaceMassAvgRotOut']['alphaFlow']
# bLBaselineNoSlipURotOut = bLBaselineNoSlip['surfaceMassAvgRotOut']['Uz']
# xbLBaselineNoSlipRotOut = bLBaselineNoSlip['surfaceMassAvgRotOut']['y']
# ybLBaselineNoSlipRotOut = bLBaselineNoSlip['surfaceMassAvgRotOut']['x']
# bLBaselineNoSlipUyRotOut = bLBaselineNoSlip['surfaceMassAvgRotOut']['Utangent']
# bLBaselineNoSlipUzRotOut = bLBaselineNoSlip['surfaceMassAvgRotOut']['Uaxial']

# bLBaselineNoSlipPtBarRotOut, bLBaselineNoSlipSurfRotOut, _, stotalPanFrac = massAvgQty(bLBaselineNoSlipPtRotOut, bLBaselineNoSlipURotOut, Nt, Nr, xbLBaselineNoSlipRotOut, ybLBaselineNoSlipRotOut, theta)
# bLBaselineNoSlipAlphaFlowBarRotOut, bLBaselineNoSlipAlphaFlowSurfRotOut, _,_ = massAvgQty(bLBaselineNoSlipAlphaFlowRotOut, bLBaselineNoSlipURotOut, Nt, Nr, xbLBaselineNoSlipRotOut, ybLBaselineNoSlipRotOut, theta)
# bLBaselineNoSlipUyBarRotOut, bLBaselineNoSlipUySurfRotOut, _,_ = massAvgQty(bLBaselineNoSlipUyRotOut, bLBaselineNoSlipURotOut, Nt, Nr, xbLBaselineNoSlipRotOut, ybLBaselineNoSlipRotOut, theta)
# bLBaselineNoSlipUzBarRotOut, bLBaselineNoSlipUzSurfRotOut, bLBaselineNoSlipMassFluxRotOut, _ = massAvgQty(bLBaselineNoSlipUzRotOut, bLBaselineNoSlipURotOut, Nt, Nr, xbLBaselineNoSlipRotOut, ybLBaselineNoSlipRotOut, theta)


# bLBaselineNoSlipPtOut = bLBaselineNoSlip['surfaceMassAvgOutlet']['totalP']
# bLBaselineNoSlipAlphaFlowOut = bLBaselineNoSlip['surfaceMassAvgOutlet']['alphaFlow']
# bLBaselineNoSlipUOut = bLBaselineNoSlip['surfaceMassAvgOutlet']['Uz']
# xbLBaselineNoSlipOut = bLBaselineNoSlip['surfaceMassAvgOutlet']['y']
# ybLBaselineNoSlipOut = bLBaselineNoSlip['surfaceMassAvgOutlet']['x']
# bLBaselineNoSlipUyOut = bLBaselineNoSlip['surfaceMassAvgOutlet']['Utangent']
# bLBaselineNoSlipUzOut = bLBaselineNoSlip['surfaceMassAvgOutlet']['Uaxial']

# bLBaselineNoSlipPtBarOut, bLBaselineNoSlipSurfOut, _, _ = massAvgQty(bLBaselineNoSlipPtOut, bLBaselineNoSlipUOut, Nt, Nr, xbLBaselineNoSlipOut, ybLBaselineNoSlipOut, theta)
# bLBaselineNoSlipAlphaFlowBarOut, bLBaselineNoSlipAlphaFlowSurfOut, _,_ = massAvgQty(bLBaselineNoSlipAlphaFlowOut, bLBaselineNoSlipUOut, Nt, Nr, xbLBaselineNoSlipOut, ybLBaselineNoSlipOut, theta)
# bLBaselineNoSlipUyBarOut, bLBaselineNoSlipUySurfOut, _,_ = massAvgQty(bLBaselineNoSlipUyOut, bLBaselineNoSlipUOut, Nt, Nr, xbLBaselineNoSlipOut, ybLBaselineNoSlipOut, theta)
# bLBaselineNoSlipUzBarOut, bLBaselineNoSlipUzSurfOut, bLBaselineNoSlipMassFluxOut, _ = massAvgQty(bLBaselineNoSlipUzOut, bLBaselineNoSlipUOut, Nt, Nr, xbLBaselineNoSlipOut, ybLBaselineNoSlipOut, theta)


with open(filePath + '/bladedVaryingPitchSlip.pkl', 'rb') as f:
# with open(filePath + '/bladedStraightSlipWalls.pkl', 'rb') as f:    
    bLVaryingPitchSlip = pickle.load(f)  
    
bLVaryingPitchSlipPtIn = bLVaryingPitchSlip['surfaceMassAvgInlet']['totalP']
bLVaryingPitchSlipAlphaFlowIn = bLVaryingPitchSlip['surfaceMassAvgInlet']['alphaFlow']
bLVaryingPitchSlipUIn = bLVaryingPitchSlip['surfaceMassAvgInlet']['Uz']
xbLVaryingPitchSlipIn = bLVaryingPitchSlip['surfaceMassAvgInlet']['y']
ybLVaryingPitchSlipIn = bLVaryingPitchSlip['surfaceMassAvgInlet']['x']
bLVaryingPitchSlipUyIn = bLVaryingPitchSlip['surfaceMassAvgInlet']['Utangent']
bLVaryingPitchSlipUzIn = bLVaryingPitchSlip['surfaceMassAvgInlet']['Uaxial']

bLVaryingPitchSlipPtBarIn, blVaryingPitchSlipSurfIn, _, stotalPanFrac = massAvgQty(bLVaryingPitchSlipPtIn, bLVaryingPitchSlipUIn, Nt, Nr, xbLVaryingPitchSlipIn, ybLVaryingPitchSlipIn, theta)
bLVaryingPitchSlipAlphaFlowBarIn, blVaryingPitchSlipAlphaFlowSurfIn, _,_ = massAvgQty(bLVaryingPitchSlipAlphaFlowIn, bLVaryingPitchSlipUIn, Nt, Nr, xbLVaryingPitchSlipIn, ybLVaryingPitchSlipIn, theta)
# bLVaryingPitchSlipUyBarIn, blVaryingPitchSlipUySurfIn, _,_ = massAvgQty(bLVaryingPitchSlipUyIn, bLVaryingPitchSlipUIn, Nt, Nr, xbLVaryingPitchSlipIn, ybLVaryingPitchSlipIn, theta)
# bLVaryingPitchSlipUzBarIn, blVaryingPitchSlipUzSurfIn, bLVaryingPitchSlipMassFluxIn, _ = massAvgQty(bLVaryingPitchSlipUzIn, bLVaryingPitchSlipUIn, Nt, Nr, xbLVaryingPitchSlipIn, ybLVaryingPitchSlipIn, theta)    

bLVaryingPitchSlipPtRotIn = bLVaryingPitchSlip['surfaceMassAvgRotIn']['totalP']
bLVaryingPitchSlipAlphaFlowRotIn = bLVaryingPitchSlip['surfaceMassAvgRotIn']['alphaFlow']
bLVaryingPitchSlipURotIn = bLVaryingPitchSlip['surfaceMassAvgRotIn']['Uz']
xbLVaryingPitchSlipRotIn = bLVaryingPitchSlip['surfaceMassAvgRotIn']['y']
ybLVaryingPitchSlipRotIn = bLVaryingPitchSlip['surfaceMassAvgRotIn']['x']
bLVaryingPitchSlipUyRotIn = bLVaryingPitchSlip['surfaceMassAvgRotIn']['Utangent']
bLVaryingPitchSlipUzRotIn = bLVaryingPitchSlip['surfaceMassAvgRotIn']['Uaxial']

bLVaryingPitchSlipPtBarRotIn, blVaryingPitchSlipSurfRotIn, _, sPanFrac = massAvgQty(bLVaryingPitchSlipPtRotIn, bLVaryingPitchSlipURotIn, Nt, Nr, xbLVaryingPitchSlipRotIn, ybLVaryingPitchSlipRotIn, theta)
bLVaryingPitchSlipAlphaFlowBarRotIn, blVaryingPitchSlipAlphaFlowSurfRotIn, _,_ = massAvgQty(bLVaryingPitchSlipAlphaFlowRotIn, bLVaryingPitchSlipURotIn, Nt, Nr, xbLVaryingPitchSlipRotIn, ybLVaryingPitchSlipRotIn, theta)
# bLVaryingPitchSlipUyBarRotIn, blVaryingPitchSlipUySurfRotIn, _,_ = massAvgQty(bLVaryingPitchSlipUyRotIn, bLVaryingPitchSlipURotIn, Nt, Nr, xbLVaryingPitchSlipRotIn, ybLVaryingPitchSlipRotIn, theta)
# bLVaryingPitchSlipUzBarRotIn, blVaryingPitchSlipUzSurfRotIn, bLVaryingPitchSlipMassFluxRotIn, _ = massAvgQty(bLVaryingPitchSlipUzRotIn, bLVaryingPitchSlipURotIn, Nt, Nr, xbLVaryingPitchSlipRotIn, ybLVaryingPitchSlipRotIn, theta)

bLVaryingPitchSlipPtRotOut = bLVaryingPitchSlip['surfaceMassAvgRotOut']['totalP']
bLVaryingPitchSlipAlphaFlowRotOut = bLVaryingPitchSlip['surfaceMassAvgRotOut']['alphaFlow']
bLVaryingPitchSlipURotOut = bLVaryingPitchSlip['surfaceMassAvgRotOut']['Uz']
xbLVaryingPitchSlipRotOut = bLVaryingPitchSlip['surfaceMassAvgRotOut']['y']
ybLVaryingPitchSlipRotOut = bLVaryingPitchSlip['surfaceMassAvgRotOut']['x']
bLVaryingPitchSlipUyRotOut = bLVaryingPitchSlip['surfaceMassAvgRotOut']['Utangent']
bLVaryingPitchSlipUzRotOut = bLVaryingPitchSlip['surfaceMassAvgRotOut']['Uaxial']

bLVaryingPitchSlipPtBarRotOut, bLVaryingPitchSlipSurfRotOut, _, sPanFrac = massAvgQty(bLVaryingPitchSlipPtRotOut, bLVaryingPitchSlipURotOut, Nt, Nr, xbLVaryingPitchSlipRotOut, ybLVaryingPitchSlipRotOut, theta)
bLVaryingPitchSlipAlphaFlowBarRotOut, bLVaryingPitchSlipAlphaFlowSurfRotOut, _,_ = massAvgQty(bLVaryingPitchSlipAlphaFlowRotOut, bLVaryingPitchSlipURotOut, Nt, Nr, xbLVaryingPitchSlipRotOut, ybLVaryingPitchSlipRotOut, theta)
# bLVaryingPitchSlipUyBarRotOut, bLVaryingPitchSlipUySurfRotOut, _,_ = massAvgQty(bLVaryingPitchSlipUyRotOut, bLVaryingPitchSlipURotOut, Nt, Nr, xbLVaryingPitchSlipRotOut, ybLVaryingPitchSlipRotOut, theta)
# bLVaryingPitchSlipUzBarRotOut, bLVaryingPitchSlipUzSurfRotOut, bLVaryingPitchSlipMassFluxRotOut, _ = massAvgQty(bLVaryingPitchSlipUzRotOut, bLVaryingPitchSlipURotOut, Nt, Nr, xbLVaryingPitchSlipRotOut, ybLVaryingPitchSlipRotOut, theta)


bLVaryingPitchSlipPtOut = bLVaryingPitchSlip['surfaceMassAvgOutlet']['totalP']
bLVaryingPitchSlipAlphaFlowOut = bLVaryingPitchSlip['surfaceMassAvgOutlet']['alphaFlow']
bLVaryingPitchSlipUOut = bLVaryingPitchSlip['surfaceMassAvgOutlet']['Uz']
xbLVaryingPitchSlipOut = bLVaryingPitchSlip['surfaceMassAvgOutlet']['y']
ybLVaryingPitchSlipOut = bLVaryingPitchSlip['surfaceMassAvgOutlet']['x']
bLVaryingPitchSlipUyOut = bLVaryingPitchSlip['surfaceMassAvgOutlet']['Utangent']
bLVaryingPitchSlipUzOut = bLVaryingPitchSlip['surfaceMassAvgOutlet']['Uaxial']

bLVaryingPitchSlipPtBarOut, bLVaryingPitchSlipSurfOut, _, _ = massAvgQty(bLVaryingPitchSlipPtOut, bLVaryingPitchSlipUOut, Nt, Nr, xbLVaryingPitchSlipOut, ybLVaryingPitchSlipOut, theta)
bLVaryingPitchSlipAlphaFlowBarOut, bLVaryingPitchSlipAlphaFlowSurfOut, _,_ = massAvgQty(bLVaryingPitchSlipAlphaFlowOut, bLVaryingPitchSlipUOut, Nt, Nr, xbLVaryingPitchSlipOut, ybLVaryingPitchSlipOut, theta)
# bLVaryingPitchSlipUyBarOut, bLVaryingPitchSlipUySurfOut, _,_ = massAvgQty(bLVaryingPitchSlipUyOut, bLVaryingPitchSlipUOut, Nt, Nr, xbLVaryingPitchSlipOut, ybLVaryingPitchSlipOut, theta)
# bLVaryingPitchSlipUzBarOut, bLVaryingPitchSlipUzSurfOut, bLVaryingPitchSlipMassFluxOut, _ = massAvgQty(bLVaryingPitchSlipUzOut, bLVaryingPitchSlipUOut, Nt, Nr, xbLVaryingPitchSlipOut, ybLVaryingPitchSlipOut, theta)


# with open(filePath + '/bladedVaryingPitchNoSlip.pkl', 'rb') as f:
# # with open(filePath + '/bladedStraightNoSlipWalls.pkl', 'rb') as f:    
#     bLVaryingPitchNoSlip = pickle.load(f)  
    
# bLVaryingPitchNoSlipPtIn = bLVaryingPitchNoSlip['surfaceMassAvgInlet']['totalP']
# bLVaryingPitchNoSlipAlphaFlowIn = bLVaryingPitchNoSlip['surfaceMassAvgInlet']['alphaFlow']
# bLVaryingPitchNoSlipUIn = bLVaryingPitchNoSlip['surfaceMassAvgInlet']['Uz']
# xbLVaryingPitchNoSlipIn = bLVaryingPitchNoSlip['surfaceMassAvgInlet']['y']
# ybLVaryingPitchNoSlipIn = bLVaryingPitchNoSlip['surfaceMassAvgInlet']['x']
# bLVaryingPitchNoSlipUyIn = bLVaryingPitchNoSlip['surfaceMassAvgInlet']['Utangent']
# bLVaryingPitchNoSlipUzIn = bLVaryingPitchNoSlip['surfaceMassAvgInlet']['Uaxial']

# bLVaryingPitchNoSlipPtBarIn, bLVaryingPitchNoSlipSurfIn, _, sPanFrac = massAvgQty(bLVaryingPitchNoSlipPtIn, bLVaryingPitchNoSlipUIn, Nt, Nr, xbLVaryingPitchNoSlipIn, ybLVaryingPitchNoSlipIn, theta)
# bLVaryingPitchNoSlipAlphaFlowBarIn, bLVaryingPitchNoSlipAlphaFlowSurfIn, _,_ = massAvgQty(bLVaryingPitchNoSlipAlphaFlowIn, bLVaryingPitchNoSlipUIn, Nt, Nr, xbLVaryingPitchNoSlipIn, ybLVaryingPitchNoSlipIn, theta)
# bLVaryingPitchNoSlipUyBarIn, bLVaryingPitchNoSlipUySurfIn, _,_ = massAvgQty(bLVaryingPitchNoSlipUyIn, bLVaryingPitchNoSlipUIn, Nt, Nr, xbLVaryingPitchNoSlipIn, ybLVaryingPitchNoSlipIn, theta)
# bLVaryingPitchNoSlipUzBarIn, bLVaryingPitchNoSlipUzSurfIn, bLVaryingPitchNoSlipMassFluxIn, _ = massAvgQty(bLVaryingPitchNoSlipUzIn, bLVaryingPitchNoSlipUIn, Nt, Nr, xbLVaryingPitchNoSlipIn, ybLVaryingPitchNoSlipIn, theta)    

# bLVaryingPitchNoSlipPtRotIn = bLVaryingPitchNoSlip['surfaceMassAvgRotIn']['totalP']
# bLVaryingPitchNoSlipAlphaFlowRotIn = bLVaryingPitchNoSlip['surfaceMassAvgRotIn']['alphaFlow']
# bLVaryingPitchNoSlipURotIn = bLVaryingPitchNoSlip['surfaceMassAvgRotIn']['Uz']
# xbLVaryingPitchNoSlipRotIn = bLVaryingPitchNoSlip['surfaceMassAvgRotIn']['y']
# ybLVaryingPitchNoSlipRotIn = bLVaryingPitchNoSlip['surfaceMassAvgRotIn']['x']
# bLVaryingPitchNoSlipUyRotIn = bLVaryingPitchNoSlip['surfaceMassAvgRotIn']['Utangent']
# bLVaryingPitchNoSlipUzRotIn = bLVaryingPitchNoSlip['surfaceMassAvgRotIn']['Uaxial']

# bLVaryingPitchNoSlipPtBarRotIn, bLVaryingPitchNoSlipSurfRotIn, _, sPanFrac = massAvgQty(bLVaryingPitchNoSlipPtRotIn, bLVaryingPitchNoSlipURotIn, Nt, Nr, xbLVaryingPitchNoSlipRotIn, ybLVaryingPitchNoSlipRotIn, theta)
# bLVaryingPitchNoSlipAlphaFlowBarRotIn, bLVaryingPitchNoSlipAlphaFlowSurfRotIn, _,_ = massAvgQty(bLVaryingPitchNoSlipAlphaFlowRotIn, bLVaryingPitchNoSlipURotIn, Nt, Nr, xbLVaryingPitchNoSlipRotIn, ybLVaryingPitchNoSlipRotIn, theta)
# bLVaryingPitchNoSlipUyBarRotIn, bLVaryingPitchNoSlipUySurfRotIn, _,_ = massAvgQty(bLVaryingPitchNoSlipUyRotIn, bLVaryingPitchNoSlipURotIn, Nt, Nr, xbLVaryingPitchNoSlipRotIn, ybLVaryingPitchNoSlipRotIn, theta)
# bLVaryingPitchNoSlipUzBarRotIn, bLVaryingPitchNoSlipUzSurfRotIn, bLVaryingPitchNoSlipMassFluxRotIn, _ = massAvgQty(bLVaryingPitchNoSlipUzRotIn, bLVaryingPitchNoSlipURotIn, Nt, Nr, xbLVaryingPitchNoSlipRotIn, ybLVaryingPitchNoSlipRotIn, theta)

# bLVaryingPitchNoSlipPtRotOut = bLVaryingPitchNoSlip['surfaceMassAvgRotOut']['totalP']
# bLVaryingPitchNoSlipAlphaFlowRotOut = bLVaryingPitchNoSlip['surfaceMassAvgRotOut']['alphaFlow']
# bLVaryingPitchNoSlipURotOut = bLVaryingPitchNoSlip['surfaceMassAvgRotOut']['Uz']
# xbLVaryingPitchNoSlipRotOut = bLVaryingPitchNoSlip['surfaceMassAvgRotOut']['y']
# ybLVaryingPitchNoSlipRotOut = bLVaryingPitchNoSlip['surfaceMassAvgRotOut']['x']
# bLVaryingPitchNoSlipUyRotOut = bLVaryingPitchNoSlip['surfaceMassAvgRotOut']['Utangent']
# bLVaryingPitchNoSlipUzRotOut = bLVaryingPitchNoSlip['surfaceMassAvgRotOut']['Uaxial']

# bLVaryingPitchNoSlipPtBarRotOut, bLVaryingPitchSurfNoSlipRotOut, _, sPanFrac = massAvgQty(bLVaryingPitchNoSlipPtRotOut, bLVaryingPitchNoSlipURotOut, Nt, Nr, xbLVaryingPitchNoSlipRotOut, ybLVaryingPitchNoSlipRotOut, theta)
# bLVaryingPitchNoSlipAlphaFlowBarRotOut, bLVaryingPitchNoSlipAlphaFlowSurfRotOut, _,_ = massAvgQty(bLVaryingPitchNoSlipAlphaFlowRotOut, bLVaryingPitchNoSlipURotOut, Nt, Nr, xbLVaryingPitchNoSlipRotOut, ybLVaryingPitchNoSlipRotOut, theta)
# bLVaryingPitchNoSlipUyBarRotOut, bLVaryingPitchNoSlipUySurfRotOut, _,_ = massAvgQty(bLVaryingPitchNoSlipUyRotOut, bLVaryingPitchNoSlipURotOut, Nt, Nr, xbLVaryingPitchNoSlipRotOut, ybLVaryingPitchNoSlipRotOut, theta)
# bLVaryingPitchNoSlipUzBarRotOut, bLVaryingPitchNoSlipUzSurfRotOut, bLVaryingPitchNoSlipMassFluxRotOut, _ = massAvgQty(bLVaryingPitchNoSlipUzRotOut, bLVaryingPitchNoSlipURotOut, Nt, Nr, xbLVaryingPitchNoSlipRotOut, ybLVaryingPitchNoSlipRotOut, theta)


# bLVaryingPitchNoSlipPtOut = bLVaryingPitchNoSlip['surfaceMassAvgOutlet']['totalP']
# bLVaryingPitchNoSlipAlphaFlowOut = bLVaryingPitchNoSlip['surfaceMassAvgOutlet']['alphaFlow']
# bLVaryingPitchNoSlipUOut = bLVaryingPitchNoSlip['surfaceMassAvgOutlet']['Uz']
# xbLVaryingPitchNoSlipOut = bLVaryingPitchNoSlip['surfaceMassAvgOutlet']['y']
# ybLVaryingPitchNoSlipOut = bLVaryingPitchNoSlip['surfaceMassAvgOutlet']['x']
# bLVaryingPitchNoSlipUyOut = bLVaryingPitchNoSlip['surfaceMassAvgOutlet']['Utangent']
# bLVaryingPitchNoSlipUzOut = bLVaryingPitchNoSlip['surfaceMassAvgOutlet']['Uaxial']

# bLVaryingPitchNoSlipPtBarOut, bLVaryingPitchNoSlipSurfOut, _, _ = massAvgQty(bLVaryingPitchNoSlipPtOut, bLVaryingPitchNoSlipUOut, Nt, Nr, xbLVaryingPitchNoSlipOut, ybLVaryingPitchNoSlipOut, theta)
# bLVaryingPitchNoSlipAlphaFlowBarOut, bLVaryingPitchNoSlipAlphaFlowSurfOut, _,_ = massAvgQty(bLVaryingPitchNoSlipAlphaFlowOut, bLVaryingPitchNoSlipUOut, Nt, Nr, xbLVaryingPitchNoSlipOut, ybLVaryingPitchNoSlipOut, theta)
# bLVaryingPitchNoSlipUyBarOut, bLVaryingPitchNoSlipUySurfOut, _,_ = massAvgQty(bLVaryingPitchNoSlipUyOut, bLVaryingPitchNoSlipUOut, Nt, Nr, xbLVaryingPitchNoSlipOut, ybLVaryingPitchNoSlipOut, theta)
# bLVaryingPitchNoSlipUzBarOut, bLVaryingPitchNoSlipUzSurfOut, bLVaryingPitchNoSlipMassFluxOut, _ = massAvgQty(bLVaryingPitchNoSlipUzOut, bLVaryingPitchNoSlipUOut, Nt, Nr, xbLVaryingPitchNoSlipOut, ybLVaryingPitchNoSlipOut, theta)

with open(filePath + '/bladedVaryingStaggerSlip.pkl', 'rb') as f:
# with open(filePath + '/bladedStraightSlipWalls.pkl', 'rb') as f:    
    bLVaryingStaggerSlip = pickle.load(f)  
    
bLVaryingStaggerSlipPtIn = bLVaryingStaggerSlip['surfaceMassAvgInlet']['totalP']
bLVaryingStaggerSlipAlphaFlowIn = bLVaryingStaggerSlip['surfaceMassAvgInlet']['alphaFlow']
bLVaryingStaggerSlipUIn = bLVaryingStaggerSlip['surfaceMassAvgInlet']['Uz']
xbLVaryingStaggerSlipIn = bLVaryingStaggerSlip['surfaceMassAvgInlet']['y']
ybLVaryingStaggerSlipIn = bLVaryingStaggerSlip['surfaceMassAvgInlet']['x']
bLVaryingStaggerSlipUyIn = bLVaryingStaggerSlip['surfaceMassAvgInlet']['Utangent']
bLVaryingStaggerSlipUzIn = bLVaryingStaggerSlip['surfaceMassAvgInlet']['Uaxial']

bLVaryingStaggerSlipPtBarIn, bLVaryingStaggerSlipSurfIn, _, _ = massAvgQty(bLVaryingStaggerSlipPtIn, bLVaryingStaggerSlipUIn, Nt, Nr, xbLVaryingStaggerSlipIn, ybLVaryingStaggerSlipIn, theta)
bLVaryingStaggerSlipAlphaFlowBarIn, bLVaryingStaggerSlipAlphaFlowSurfIn, _,_ = massAvgQty(bLVaryingStaggerSlipAlphaFlowIn, bLVaryingStaggerSlipUIn, Nt, Nr, xbLVaryingStaggerSlipIn, ybLVaryingStaggerSlipIn, theta)
# bLVaryingStaggerSlipUyBarIn, bLVaryingStaggerSlipUySurfIn, _,_ = massAvgQty(bLVaryingStaggerSlipUyIn, bLVaryingStaggerSlipUIn, Nt, Nr, xbLVaryingStaggerSlipIn, ybLVaryingStaggerSlipIn, theta)
# bLVaryingStaggerSlipUzBarIn, bLVaryingStaggerSlipUzSurfIn, bLVaryingStaggerSlipMassFluxIn, _ = massAvgQty(bLVaryingStaggerSlipUzIn, bLVaryingStaggerSlipUIn, Nt, Nr, xbLVaryingStaggerSlipIn, ybLVaryingStaggerSlipIn, theta)    

bLVaryingStaggerSlipPtRotIn = bLVaryingStaggerSlip['surfaceMassAvgRotIn']['totalP']
bLVaryingStaggerSlipAlphaFlowRotIn = bLVaryingStaggerSlip['surfaceMassAvgRotIn']['alphaFlow']
bLVaryingStaggerSlipURotIn = bLVaryingStaggerSlip['surfaceMassAvgRotIn']['Uz']
xbLVaryingStaggerSlipRotIn = bLVaryingStaggerSlip['surfaceMassAvgRotIn']['y']
ybLVaryingStaggerSlipRotIn = bLVaryingStaggerSlip['surfaceMassAvgRotIn']['x']
bLVaryingStaggerSlipUyRotIn = bLVaryingStaggerSlip['surfaceMassAvgRotIn']['Utangent']
bLVaryingStaggerSlipUzRotIn = bLVaryingStaggerSlip['surfaceMassAvgRotIn']['Uaxial']

bLVaryingStaggerSlipPtBarRotIn, bLVaryingStaggerSlipSurfRotIn, _, _ = massAvgQty(bLVaryingStaggerSlipPtRotIn, bLVaryingStaggerSlipURotIn, Nt, Nr, xbLVaryingStaggerSlipRotIn, ybLVaryingStaggerSlipRotIn, theta)
bLVaryingStaggerSlipAlphaFlowBarRotIn, bLVaryingStaggerSlipAlphaFlowSurfRotIn, _,_ = massAvgQty(bLVaryingStaggerSlipAlphaFlowRotIn, bLVaryingStaggerSlipURotIn, Nt, Nr, xbLVaryingStaggerSlipRotIn, ybLVaryingStaggerSlipRotIn, theta)
# bLVaryingStaggerSlipUyBarRotIn, bLVaryingStaggerSlipUySurfRotIn, _,_ = massAvgQty(bLVaryingStaggerSlipUyRotIn, bLVaryingStaggerSlipURotIn, Nt, Nr, xbLVaryingStaggerSlipRotIn, ybLVaryingStaggerSlipRotIn, theta)
# bLVaryingStaggerSlipUzBarRotIn, bLVaryingStaggerSlipUzSurfRotIn, bLVaryingStaggerSlipMassFluxRotIn, _ = massAvgQty(bLVaryingStaggerSlipUzRotIn, bLVaryingStaggerSlipURotIn, Nt, Nr, xbLVaryingStaggerSlipRotIn, ybLVaryingStaggerSlipRotIn, theta)

bLVaryingStaggerSlipPtRotOut = bLVaryingStaggerSlip['surfaceMassAvgRotOut']['totalP']
bLVaryingStaggerSlipAlphaFlowRotOut = bLVaryingStaggerSlip['surfaceMassAvgRotOut']['alphaFlow']
bLVaryingStaggerSlipURotOut = bLVaryingStaggerSlip['surfaceMassAvgRotOut']['Uz']
xbLVaryingStaggerSlipRotOut = bLVaryingStaggerSlip['surfaceMassAvgRotOut']['y']
ybLVaryingStaggerSlipRotOut = bLVaryingStaggerSlip['surfaceMassAvgRotOut']['x']
bLVaryingStaggerSlipUyRotOut = bLVaryingStaggerSlip['surfaceMassAvgRotOut']['Utangent']
bLVaryingStaggerSlipUzRotOut = bLVaryingStaggerSlip['surfaceMassAvgRotOut']['Uaxial']

bLVaryingStaggerSlipPtBarRotOut, bLVaryingStaggerSlipSurfRotOut, _, _ = massAvgQty(bLVaryingStaggerSlipPtRotOut, bLVaryingStaggerSlipURotOut, Nt, Nr, xbLVaryingStaggerSlipRotOut, ybLVaryingStaggerSlipRotOut, theta)
bLVaryingStaggerSlipAlphaFlowBarRotOut, bLVaryingStaggerSlipAlphaFlowSurfRotOut, _,_ = massAvgQty(bLVaryingStaggerSlipAlphaFlowRotOut, bLVaryingStaggerSlipURotOut, Nt, Nr, xbLVaryingStaggerSlipRotOut, ybLVaryingStaggerSlipRotOut, theta)
# bLVaryingStaggerSlipUyBarRotOut, bLVaryingStaggerSlipUySurfRotOut, _,_ = massAvgQty(bLVaryingStaggerSlipUyRotOut, bLVaryingStaggerSlipURotOut, Nt, Nr, xbLVaryingStaggerSlipRotOut, ybLVaryingStaggerSlipRotOut, theta)
# bLVaryingStaggerSlipUzBarRotOut, bLVaryingStaggerSlipUzSurfRotOut, bLVaryingStaggerSlipMassFluxRotOut, _ = massAvgQty(bLVaryingStaggerSlipUzRotOut, bLVaryingStaggerSlipURotOut, Nt, Nr, xbLVaryingStaggerSlipRotOut, ybLVaryingStaggerSlipRotOut, theta)


bLVaryingStaggerSlipPtOut = bLVaryingStaggerSlip['surfaceMassAvgOutlet']['totalP']
bLVaryingStaggerSlipAlphaFlowOut = bLVaryingStaggerSlip['surfaceMassAvgOutlet']['alphaFlow']
bLVaryingStaggerSlipUOut = bLVaryingStaggerSlip['surfaceMassAvgOutlet']['Uz']
xbLVaryingStaggerSlipOut = bLVaryingStaggerSlip['surfaceMassAvgOutlet']['y']
ybLVaryingStaggerSlipOut = bLVaryingStaggerSlip['surfaceMassAvgOutlet']['x']
bLVaryingStaggerSlipUyOut = bLVaryingStaggerSlip['surfaceMassAvgOutlet']['Utangent']
bLVaryingStaggerSlipUzOut = bLVaryingStaggerSlip['surfaceMassAvgOutlet']['Uaxial']

bLVaryingStaggerSlipPtBarOut, bLVaryingStaggerSlipSurfOut, _, _ = massAvgQty(bLVaryingStaggerSlipPtOut, bLVaryingStaggerSlipUOut, Nt, Nr, xbLVaryingStaggerSlipOut, ybLVaryingStaggerSlipOut, theta)
bLVaryingStaggerSlipAlphaFlowBarOut, bLVaryingStaggerSlipAlphaFlowSurfOut, _,_ = massAvgQty(bLVaryingStaggerSlipAlphaFlowOut, bLVaryingStaggerSlipUOut, Nt, Nr, xbLVaryingStaggerSlipOut, ybLVaryingStaggerSlipOut, theta)
# bLVaryingStaggerSlipUyBarOut, bLVaryingStaggerSlipUySurfOut, _,_ = massAvgQty(bLVaryingStaggerSlipUyOut, bLVaryingStaggerSlipUOut, Nt, Nr, xbLVaryingStaggerSlipOut, ybLVaryingStaggerSlipOut, theta)
# bLVaryingStaggerSlipUzBarOut, bLVaryingStaggerSlipUzSurfOut, bLVaryingStaggerSlipMassFluxOut, _ = massAvgQty(bLVaryingStaggerSlipUzOut, bLVaryingStaggerSlipUOut, Nt, Nr, xbLVaryingStaggerSlipOut, ybLVaryingStaggerSlipOut, theta)


# with open(filePath + '/bladedVaryingStaggerNoSlip.pkl', 'rb') as f:
# # with open(filePath + '/bladedStraightNoSlipWalls.pkl', 'rb') as f:    
#     bLVaryingStaggerNoSlip = pickle.load(f)  
    
# bLVaryingStaggerNoSlipPtIn = bLVaryingStaggerNoSlip['surfaceMassAvgInlet']['totalP']
# bLVaryingStaggerNoSlipAlphaFlowIn = bLVaryingStaggerNoSlip['surfaceMassAvgInlet']['alphaFlow']
# bLVaryingStaggerNoSlipUIn = bLVaryingStaggerNoSlip['surfaceMassAvgInlet']['Uz']
# xbLVaryingStaggerNoSlipIn = bLVaryingStaggerNoSlip['surfaceMassAvgInlet']['y']
# ybLVaryingStaggerNoSlipIn = bLVaryingStaggerNoSlip['surfaceMassAvgInlet']['x']
# bLVaryingStaggerNoSlipUyIn = bLVaryingStaggerNoSlip['surfaceMassAvgInlet']['Utangent']
# bLVaryingStaggerNoSlipUzIn = bLVaryingStaggerNoSlip['surfaceMassAvgInlet']['Uaxial']

# bLVaryingStaggerNoSlipPtBarIn, bLVaryingStaggerNoSlipSurfIn, _, _ = massAvgQty(bLVaryingStaggerNoSlipPtIn, bLVaryingStaggerNoSlipUIn, Nt, Nr, xbLVaryingStaggerNoSlipIn, ybLVaryingStaggerNoSlipIn, theta)
# bLVaryingStaggerNoSlipAlphaFlowBarIn, bLVaryingStaggerNoSlipAlphaFlowSurfIn, _,_ = massAvgQty(bLVaryingStaggerNoSlipAlphaFlowIn, bLVaryingStaggerNoSlipUIn, Nt, Nr, xbLVaryingStaggerNoSlipIn, ybLVaryingStaggerNoSlipIn, theta)
# bLVaryingStaggerNoSlipUyBarIn, bLVaryingStaggerNoSlipUySurfIn, _,_ = massAvgQty(bLVaryingStaggerNoSlipUyIn, bLVaryingStaggerNoSlipUIn, Nt, Nr, xbLVaryingStaggerNoSlipIn, ybLVaryingStaggerNoSlipIn, theta)
# bLVaryingStaggerNoSlipUzBarIn, bLVaryingStaggerNoSlipUzSurfIn, bLVaryingStaggerNoSlipMassFluxIn, _ = massAvgQty(bLVaryingStaggerNoSlipUzIn, bLVaryingStaggerNoSlipUIn, Nt, Nr, xbLVaryingStaggerNoSlipIn, ybLVaryingStaggerNoSlipIn, theta)    

# bLVaryingStaggerNoSlipPtRotIn = bLVaryingStaggerNoSlip['surfaceMassAvgRotIn']['totalP']
# bLVaryingStaggerNoSlipAlphaFlowRotIn = bLVaryingStaggerNoSlip['surfaceMassAvgRotIn']['alphaFlow']
# bLVaryingStaggerNoSlipURotIn = bLVaryingStaggerNoSlip['surfaceMassAvgRotIn']['Uz']
# xbLVaryingStaggerNoSlipRotIn = bLVaryingStaggerNoSlip['surfaceMassAvgRotIn']['y']
# ybLVaryingStaggerNoSlipRotIn = bLVaryingStaggerNoSlip['surfaceMassAvgRotIn']['x']
# bLVaryingStaggerNoSlipUyRotIn = bLVaryingStaggerNoSlip['surfaceMassAvgRotIn']['Utangent']
# bLVaryingStaggerNoSlipUzRotIn = bLVaryingStaggerNoSlip['surfaceMassAvgRotIn']['Uaxial']

# bLVaryingStaggerNoSlipPtBarRotIn, bLVaryingStaggerNoSlipSurfRotIn, _, _ = massAvgQty(bLVaryingStaggerNoSlipPtRotIn, bLVaryingStaggerNoSlipURotIn, Nt, Nr, xbLVaryingStaggerNoSlipRotIn, ybLVaryingStaggerNoSlipRotIn, theta)
# bLVaryingStaggerNoSlipAlphaFlowBarRotIn, bLVaryingStaggerNoSlipAlphaFlowSurfRotIn, _,_ = massAvgQty(bLVaryingStaggerNoSlipAlphaFlowRotIn, bLVaryingStaggerNoSlipURotIn, Nt, Nr, xbLVaryingStaggerNoSlipRotIn, ybLVaryingStaggerNoSlipRotIn, theta)
# bLVaryingStaggerNoSlipUyBarRotIn, bLVaryingStaggerNoSlipUySurfRotIn, _,_ = massAvgQty(bLVaryingStaggerNoSlipUyRotIn, bLVaryingStaggerNoSlipURotIn, Nt, Nr, xbLVaryingStaggerNoSlipRotIn, ybLVaryingStaggerNoSlipRotIn, theta)
# bLVaryingStaggerNoSlipUzBarRotIn, bLVaryingStaggerNoSlipUzSurfRotIn, bLVaryingStaggerNoSlipMassFluxRotIn, _ = massAvgQty(bLVaryingStaggerNoSlipUzRotIn, bLVaryingStaggerNoSlipURotIn, Nt, Nr, xbLVaryingStaggerNoSlipRotIn, ybLVaryingStaggerNoSlipRotIn, theta)

# bLVaryingStaggerNoSlipPtRotOut = bLVaryingStaggerNoSlip['surfaceMassAvgRotOut']['totalP']
# bLVaryingStaggerNoSlipAlphaFlowRotOut = bLVaryingStaggerNoSlip['surfaceMassAvgRotOut']['alphaFlow']
# bLVaryingStaggerNoSlipURotOut = bLVaryingStaggerNoSlip['surfaceMassAvgRotOut']['Uz']
# xbLVaryingStaggerNoSlipRotOut = bLVaryingStaggerNoSlip['surfaceMassAvgRotOut']['y']
# ybLVaryingStaggerNoSlipRotOut = bLVaryingStaggerNoSlip['surfaceMassAvgRotOut']['x']
# bLVaryingStaggerNoSlipUyRotOut = bLVaryingStaggerNoSlip['surfaceMassAvgRotOut']['Utangent']
# bLVaryingStaggerNoSlipUzRotOut = bLVaryingStaggerNoSlip['surfaceMassAvgRotOut']['Uaxial']

# bLVaryingStaggerNoSlipPtBarRotOut, bLVaryingStaggerSurfRotOut, _, _ = massAvgQty(bLVaryingStaggerNoSlipPtRotOut, bLVaryingStaggerNoSlipURotOut, Nt, Nr, xbLVaryingStaggerNoSlipRotOut, ybLVaryingStaggerNoSlipRotOut, theta)
# bLVaryingStaggerNoSlipAlphaFlowBarRotOut, bLVaryingStaggerNoSlipAlphaFlowSurfRotOut, _,_ = massAvgQty(bLVaryingStaggerNoSlipAlphaFlowRotOut, bLVaryingStaggerNoSlipURotOut, Nt, Nr, xbLVaryingStaggerNoSlipRotOut, ybLVaryingStaggerNoSlipRotOut, theta)
# bLVaryingStaggerNoSlipUyBarRotOut, bLVaryingStaggerNoSlipUySurfRotOut, _,_ = massAvgQty(bLVaryingStaggerNoSlipUyRotOut, bLVaryingStaggerNoSlipURotOut, Nt, Nr, xbLVaryingStaggerNoSlipRotOut, ybLVaryingStaggerNoSlipRotOut, theta)
# bLVaryingStaggerNoSlipUzBarRotOut, bLVaryingStaggerNoSlipUzSurfRotOut, bLVaryingStaggerNoSlipMassFluxRotOut, _ = massAvgQty(bLVaryingStaggerNoSlipUzRotOut, bLVaryingStaggerNoSlipURotOut, Nt, Nr, xbLVaryingStaggerNoSlipRotOut, ybLVaryingStaggerNoSlipRotOut, theta)


# bLVaryingStaggerNoSlipPtOut = bLVaryingStaggerNoSlip['surfaceMassAvgOutlet']['totalP']
# bLVaryingStaggerNoSlipAlphaFlowOut = bLVaryingStaggerNoSlip['surfaceMassAvgOutlet']['alphaFlow']
# bLVaryingStaggerNoSlipUOut = bLVaryingStaggerNoSlip['surfaceMassAvgOutlet']['Uz']
# xbLVaryingStaggerNoSlipOut = bLVaryingStaggerNoSlip['surfaceMassAvgOutlet']['y']
# ybLVaryingStaggerNoSlipOut = bLVaryingStaggerNoSlip['surfaceMassAvgOutlet']['x']
# bLVaryingStaggerNoSlipUyOut = bLVaryingStaggerNoSlip['surfaceMassAvgOutlet']['Utangent']
# bLVaryingStaggerNoSlipUzOut = bLVaryingStaggerNoSlip['surfaceMassAvgOutlet']['Uaxial']

# bLVaryingStaggerNoSlipPtBarOut, bLVaryingStaggerNoSlipSurfOut, _, _ = massAvgQty(bLVaryingStaggerNoSlipPtOut, bLVaryingStaggerNoSlipUOut, Nt, Nr, xbLVaryingStaggerNoSlipOut, ybLVaryingStaggerNoSlipOut, theta)
# bLVaryingStaggerNoSlipAlphaFlowBarOut, bLVaryingStaggerNoSlipAlphaFlowSurfOut, _,_ = massAvgQty(bLVaryingStaggerNoSlipAlphaFlowOut, bLVaryingStaggerNoSlipUOut, Nt, Nr, xbLVaryingStaggerNoSlipOut, ybLVaryingStaggerNoSlipOut, theta)
# bLVaryingStaggerNoSlipUyBarOut, bLVaryingStaggerNoSlipUySurfOut, _,_ = massAvgQty(bLVaryingStaggerNoSlipUyOut, bLVaryingStaggerNoSlipUOut, Nt, Nr, xbLVaryingStaggerNoSlipOut, ybLVaryingStaggerNoSlipOut, theta)
# bLVaryingStaggerNoSlipUzBarOut, bLVaryingStaggerNoSlipUzSurfOut, bLVaryingStaggerNoSlipMassFluxOut, _ = massAvgQty(bLVaryingStaggerNoSlipUzOut, bLVaryingStaggerNoSlipUOut, Nt, Nr, xbLVaryingStaggerNoSlipOut, ybLVaryingStaggerNoSlipOut, theta)
#%%
# plt.plot(theta, blBaselineAplhaFlowSurfIn[:,2], 'k', label='blade LE')
plt.plot(theta, bLBaselineSlipAlphaFlowSurfRotOut[:,2], 'k--', label='blade -TE')
# plt.plot(theta, bLBaselineSlipAlphaFlowSurfOut[:,2], 'k--.', label='blade -Outlet')
plt.plot(theta, bfBaselineAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')

# plt.plot(theta, bLBaselineSlipSurfRotOut[:,2], 'k--', label='blade -TE')
# plt.plot(theta, bLBaselineSlipSurfOut[:,2], 'k--.', label='blade -Outlet')
# plt.plot(theta, bfBaselineSurfRotOut[:,2], 'k', label='bodyForce')
# plt.plot(theta, bLVaryingPitchSlipAlphaFlowSurfRotOut[:,2], 'k--', label='blade')
# plt.plot(theta, bfVaryingPitchAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')

# plt.plot(theta, bLVaryingStaggerSlipAlphaFlowSurfRotOut[:,2], 'k--', label='blade')
# plt.plot(theta, bfVaryingStaggerAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')
# plt.plot(theta, bLBaselineNoSlipAlphaFlowSurfOut[:,2], 'k--.', label='outlet')

# plt.plot(theta, bLBaselineMassFluxIn, 'k', label='blade LE')
# plt.plot(theta, bLBaselineMassFlux, 'k--', label='blade TE')
# plt.plot(theta, bLBaselineMassFluxOut, 'k--.', label='outlet')

# plt.plot(theta, blBaselineUzSurfIn[:,4], 'k', label='blade LE')
# plt.plot(theta, blBaselineUzSurf[:,4], 'k--', label='blade TE')
# plt.plot(theta, blBaselineUzSurfOut[:,4], 'k--.', label='outlet')

# plt.plot(theta, (blBaselineSurfIn[:,2] - bLBaselinePtBarIn[2])/(0.5*bLBaselineUzBarIn[2]**2), 'k', label='blade LE')
# plt.plot(theta, (blBaselineSurf[:,2] - bLBaselinePtBarIn[2])/(0.5*bLBaselineUzBarIn[2]**2), 'k--', label='blade TE')
# plt.plot(theta, (blBaselineSurfOut[:,2] - bLBaselinePtBarIn[2])/(0.5*bLBaselineUzBarIn[2]**2), 'k--.', label='outlet')

# plt.plot(theta, blBaselineLossCoeff[:,2], 'k--')#, label='bladed'  )
# plt.plot(theta,blVaryingPitchSlipAlphaFlowSurfRotOut[:,2], 'r')
# plt.plot(theta,blVaryingPitchNoSlipAlphaFlowSurfRotOut[:,2], 'b')

# plt.plot(theta,blVaryingStaggerSlipAlphaFlowSurfRotOut[:,2], 'r')
# plt.plot(theta,blVaryingStaggerNoSlipAlphaFlowSurfRotOut[:,2], 'b')
# plt.plot(theta, blBaselineUySurf[:,5], 'g')
plt.grid()
# plt.axis('equal')
# plt.plot(theta, blBaselineAplhaFlowSurf[:,0], 'k')
# plt.plot(theta, blVaryingPitchAlphaFlowSurf[:,9], 'r')
plt.legend()
plt.xlabel('circumferential direction')
# plt.ylabel('stagger angle (deg)')
plt.ylabel('flow angle (deg)')
# plt.ylabel('mass flux (Uaxial) (m/s)')
# plt.ylabel('pressure coefficient')
# plt.ylabel('static pressure (guage)')
# plt.savefig('/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/plot/baselineMassFlux.jpg', dpi=500, bbox_inches='tight')

#%%
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

# Light smoothing to remove measurement/numerical noise only
# V_tan_smooth = gaussian_filter1d(bLVaryingPitchSlipAlphaFlowSurfRotOut[:,2], sigma=5, mode='wrap')
baselineSmooth = gaussian_filter1d(bLBaselineSlipAlphaFlowSurfRotOut[:,2], sigma=25, mode='wrap')
varyingPitchSmooth = gaussian_filter1d(bLVaryingPitchSlipAlphaFlowSurfRotOut[:,2], sigma=15, mode='wrap')
varyingStaggerSmooth = gaussian_filter1d(bLVaryingStaggerSlipAlphaFlowSurfRotOut[:,2], sigma=35, mode='wrap')


# Plot
plt.figure(figsize=(10, 5))
plt.plot(theta, bLBaselineSlipAlphaFlowSurfRotOut[:,2], 'ko-', alpha=0.4, linewidth=1, markersize=3, label='blade Raw')
plt.plot(theta, bLVaryingPitchSlipAlphaFlowSurfRotOut[:,2], 'ko-', alpha=0.4, linewidth=1, markersize=3, label='blade Raw')
plt.plot(theta, bfBaselineAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')
plt.plot(theta, baselineSmooth, 'k--', linewidth=2, label='Smoothed blade(σ=5)')
plt.plot(theta, varyingPitchSmooth, 'k--', linewidth=2, label='Smoothed blade(σ=5)')
plt.plot(theta, bfVaryingPitchAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')
plt.xlabel('Circumferential Direction ')
plt.ylabel('flow angle (deg)')
plt.grid(True, alpha=0.7)
plt.legend()
# plt.ylim(8,12)
# plt.title('Tangential Velocity at Blade Inlet - Varying Pitch')
#%%
fig, axes1 = plt.subplots(3, 1, figsize=(10, 8))
for ax in axes1:
    ax.set_ylim(2,15)
    ax.grid(True)
axes1[0].plot(theta, bfBaselineAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')
axes1[0].plot(theta, bLBaselineSlipAlphaFlowSurfRotOut[:,2], 'ko-', alpha=0.4, linewidth=1, markersize=3, label='blade Raw')
axes1[0].plot(theta, baselineSmooth, 'k--', linewidth=2, label='Smoothed blade(σ=5)')
axes1[0].set_ylabel('flow angle')
axes1[0].legend()
axes1[1].plot(theta, bfVaryingPitchAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')
axes1[1].plot(theta, bLVaryingPitchSlipAlphaFlowSurfRotOut[:,2], 'ko-', alpha=0.4, linewidth=1, markersize=3, label='blade Raw')
axes1[1].plot(theta, varyingPitchSmooth, 'k--', linewidth=2, label='Smoothed blade(σ=5)')
axes1[1].set_ylabel('flow angle')
axes1[1].legend(loc=2)
axes1[2].plot(theta, bfVaryingStaggerAplhaFlowSurfRotOut[:,2], 'k', label='bodyForce')
axes1[2].plot(theta, bLVaryingStaggerSlipAlphaFlowSurfRotOut[:,2], 'ko-', alpha=0.4, linewidth=1, markersize=3, label='blade Raw')
axes1[2].plot(theta, varyingStaggerSmooth, 'k--', linewidth=2, label='Smoothed blade(σ=5)')
axes1[2].set_ylabel('flow angle')
axes1[2].legend(loc=1)
#%%
baselineMetalAngle = np.loadtxt('../processedData/baseLine.txt')
baseExitAngle = densifyCurve(baselineMetalAngle, Nt, 'uniform').T[0]
staggerMetalAngle = np.loadtxt('../processedData/varyingStagger.txt')
XX = np.linspace(0,2*np.pi, 60)
func = CubicSpline(XX, staggerMetalAngle)
staggerExitAngle = func(theta)

#%%

baseLineDeviationBL = baseExitAngle - bLBaselineSlipAlphaFlowSurfRotOut[:,2] 
baseLineDeviationBF = baseExitAngle - bfBaselineAplhaFlowSurfRotOut[:,2]
varyingPitchDeviationBL = baseExitAngle - bLVaryingPitchSlipAlphaFlowSurfRotOut[:,2]
varyingPitchDeviationBF = baseExitAngle - bfVaryingPitchAplhaFlowSurfRotOut[:,2]
varyingStaggerDeviationBL = staggerExitAngle - bLVaryingStaggerSlipAlphaFlowSurfRotOut[:,2] 
varyingStaggerDeviationBF = staggerExitAngle -  bfVaryingStaggerAplhaFlowSurfRotOut[:,2]

baselineSmoothD = gaussian_filter1d(baseLineDeviationBL, sigma=25, mode='wrap')
varyingPitchSmoothD = gaussian_filter1d(varyingPitchDeviationBL, sigma=15, mode='wrap')
varyingStaggerSmoothD = gaussian_filter1d(varyingStaggerDeviationBL, sigma=35, mode='wrap')

plt.plot(theta, baselineSmoothD, 'k--', linewidth=2, label='Smoothed blade(σ=5)')
plt.plot(theta, varyingPitchSmoothD, 'g', linewidth=2, label='Smoothed blade(σ=5)')
plt.plot(theta, varyingStaggerSmoothD, 'k--', linewidth=2, label='Smoothed blade(σ=5)')
plt.plot(theta, baseLineDeviationBF, 'k')
plt.plot(theta, varyingPitchDeviationBF, 'k')
plt.plot(theta, varyingStaggerDeviationBF, 'k')
#%%
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Your existing smoothing
# baselineSmooth = gaussian_filter1d(baseLineDeviationBL, sigma=10, mode='wrap')
# varyingPitchSmooth = gaussian_filter1d(varyingPitchDeviationBL, sigma=5, mode='wrap')
# varyingStaggerSmooth = gaussian_filter1d(varyingStaggerDeviationBL, sigma=5, mode='wrap')

# baselineSmooth1 = gaussian_filter1d(bLBaselineSlipAlphaFlowSurfRotOut[:,2], sigma=10, mode='wrap')
# varyingPitchSmooth1 = gaussian_filter1d(bLVaryingPitchSlipAlphaFlowSurfRotOut[:,2], sigma=5, mode='wrap')
# varyingStaggerSmooth1 = gaussian_filter1d(bLVaryingStaggerSlipAlphaFlowSurfRotOut[:,2] , sigma=5, mode='wrap')

# # Create figure
# fig, ax = plt.subplots(figsize=(10, 6))

# # Case 1: Baseline (circles)
# ax.plot(theta, baselineSmooth, 'k--o', linewidth=2, markersize=8, 
#         markevery=60, fillstyle='none')
# ax.plot(theta, baseLineDeviationBF, 'k-o', linewidth=1.5, markersize=8, 
#         markevery=60)

# # Case 2: Varying Pitch (squares)
# ax.plot(theta, varyingPitchSmooth, 'k--s', linewidth=2, markersize=8, 
#         markevery=60, fillstyle='none')
# ax.plot(theta, varyingPitchDeviationBF, 'k-s', linewidth=1.5, markersize=8, 
#         markevery=60)

# # Case 3: Varying Stagger (triangles)
# ax.plot(theta, varyingStaggerSmooth, 'k--^', linewidth=2, markersize=8, 
#         markevery=60, fillstyle='none')
# ax.plot(theta, varyingStaggerDeviationBF, 'k-^', linewidth=1.5, markersize=8, 
#         markevery=60)

# # Custom legend
# # Symbols for cases
# symbol_handles = [
#     Line2D([0], [0], color='k', marker='o', linestyle='none', markersize=6, 
#            label='Baseline', fillstyle='none'),
#     Line2D([0], [0], color='k', marker='s', linestyle='none', markersize=6, 
#            label='Varying Pitch', fillstyle='none'),
#     Line2D([0], [0], color='k', marker='^', linestyle='none', markersize=6, 
#            label='Varying Stagger', fillstyle='none')
# ]

# # Line styles for methods
# line_handles = [
#     Line2D([0], [0], color='k', linestyle='--', linewidth=2, label='RANS'),
#     Line2D([0], [0], color='k', linestyle='-', linewidth=2.5, label='Body Force')
# ]

# # Combine legends
# first_legend = ax.legend(handles=symbol_handles, loc='lower left', 
#                          title='Cases', fontsize=18)
# ax.add_artist(first_legend)  # Keep first legend when adding second
# ax.legend(handles=line_handles, loc='upper right', 
#           title='Method', fontsize=18)

# ax.set_xlabel('Theta (radians)', fontsize=20)
# ax.set_ylabel('Deviation (deg)', fontsize=20)
# ax.grid(True, alpha=0.5)
# plt.tight_layout()
# # plt.savefig('/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/plot/All.png', dpi=500, bbox_inches='tight')
# # plt.savefig('/home/adekola/Documents/New_PhD/laTexDocuments/bodyForcePaper/All.png', dpi=500, bbox_inches='tight')
# plt.show()
#%%
plt.plot(bFBaselineAlphaFlowBarRotOut, spanFrac, 'k' , label='bodyforce')
plt.plot(bLBaselineSlipAlphaFlowBarRotOut, spanFrac, 'k--', label='bladed')
# # plt.plot(bLBaselineNoSlipAlphaFlowBarRotOut, spanFrac, 'k--', label='bladed baseline slipWalls')
# # plt.plot(bFBaselineUyBar, spanFrac, 'k-.' , label='bodyforce baseline')
# # plt.plot(bLBaselineUyBarOut, spanFrac,  'k', label='bladed baseline')
# # plt.plot(theta, bfVaryingPitchAlphaFlowSurf[:,2], 'k--', label='bodyfoce VaryingPitch')
# # plt.plot(theta, blVaryingPitchAlphaFlowSurfOut[:,2], 'k', label='bladed varyingPitch')
# # plt.plot(theta, bfVaryingStaggerAlphaFlowSurf[:,2], 'k--', label='bodyfoce VaryingStagger')
# # plt.plot(theta, blVaryingStaggerAlphaFlowSurfOut[:,2], 'k', label='bladed varyingStagger')
plt.grid()
# plt.legend(loc=0)
# plt.ylabel(r'$r/r_{tip}$', fontsize=10)
plt.ylabel('span fraction', fontsize=20)
plt.xlabel('flow angle (deg)',  fontsize=20)
# plt.savefig('/home/adekola/Documents/New_PhD/laTexDocuments/bodyForcePaper/spanwiseZero.png', dpi=500, bbox_inches='tight')
# plt.xlabel('circumferential direction')
# plt.ylabel('flow angle (deg)')
# plt.xlabel('tangential velocity (m/s)')
plt.xlim(0,10)
plt.savefig('/home/adekola/Documents/New_PhD/papers/ConferenceSubmission/TE2026/BF-OpenFOAM/plots/spanwise.png', dpi=500, bbox_inches='tight')
# plt.savefig('/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/plot/bodyForceBladedVaryingStagger.jpg', dpi=500, bbox_inches='tight')

#%%%
baselineSmooth1 = gaussian_filter1d(baseLineDeviationBL, sigma=25, mode='wrap')
varyingPitchSmooth1 = gaussian_filter1d(varyingPitchDeviationBL, sigma=15, mode='wrap')
varyingStaggerSmooth1 = gaussian_filter1d(varyingStaggerDeviationBL, sigma=35, mode='wrap')
pitchData = np.loadtxt('../pitch.txt', delimiter=',')
staggerVal = np.loadtxt('../stagger.txt', delimiter=',')
deltaBaseline = baseLineDeviationBF - baselineSmooth1
deltaPitch = varyingPitchDeviationBF - varyingPitchSmooth1
deltaStagger = varyingStaggerDeviationBF - varyingStaggerSmooth1
#%%
fig, axes1 = plt.subplots(3, 1, figsize=(10, 15), sharex=True)
for ax in axes1:
    # ax.set_ylim(2,17)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=17)

axes1[0].plot(np.degrees(pitchData[:,0]), pitchData[:,1], 'b.', label='blade position at TE (pitch)')
axes1[0].plot(np.degrees(staggerVal[:,0]),staggerVal[:,1], 'r.', label='blade position at TE (stagger)')
axes1[0].plot
axes1[0].tick_params(labelbottom=False)
axes1[1].plot(np.degrees(theta), baselineSmooth,  'k', linewidth=2, label='baseline')
axes1[1].plot(np.degrees(theta), varyingPitchSmooth, 'b', linewidth=2, label='varyingPitch')
axes1[1].tick_params(labelbottom=False)
axes1[1].plot(np.degrees(theta), varyingStaggerSmooth, 'r', linewidth=2, label='varyingStagger')
axes1[2].plot(np.degrees(theta), np.round(deltaBaseline, 1), 'k', )
axes1[2].plot(np.degrees(theta), deltaPitch, 'b', )
axes1[2].plot(np.degrees(theta), deltaStagger, 'r', )

# axes1[2].plot(theta, baseLineDeviationBF, 'k', )
# axes1[2].plot(theta, baselineSmooth, 'k--', )
# axes1[2].plot(theta, varyingPitchSmooth, 'b--')
# axes1[2].plot(theta, varyingPitchDeviationBF, 'b')
# axes1[2].plot(theta, varyingStaggerSmooth, 'r--')
# axes1[2].plot(theta, varyingStaggerDeviationBF, 'r')
axes1[0].set_ylabel('angles (deg)', fontsize=20)
# axes1[0].legend(fontsize=25)
axes1[1].set_ylabel('flow angles (deg)', fontsize=20)
# axes1[2].legend(fontsize=25)
axes1[2].set_ylabel('$\Delta$ deviation (deg)', fontsize=20)
# axes1[2].legend(fontsize=25)
axes1[2].set_xlabel('$\Theta$ (deg)' , fontsize=20)
plt.subplots_adjust(hspace=0.0)
# plt.tick_params(axis='both', labelsize=14)
# plt.savefig('/home/adekola/Documents/New_PhD/laTexDocuments/All.png', dpi=500, bbox_inches='tight')
# plt.savefig('/home/adekola/Documents/New_PhD/laTexDocuments/bodyForcePaper/All.png', dpi=500, bbox_inches='tight')
#%%
plt.plot(theta, baseLineDeviationBF - baselineSmoothD)
plt.plot(theta, varyingPitchDeviationBF-varyingPitchSmooth1)
plt.plot(theta, varyingStaggerDeviationBF-varyingStaggerSmooth1)
#%%
fig, axes1 = plt.subplots(2, 1, figsize=(20, 15), sharex=True)
for ax in axes1:
    # ax.set_ylim(2,17)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=20)
    
axes1[0].plot(np.degrees(pitchData[:,0]), pitchData[:,1], 'k.', label='blade position at TE (pitch)')   
axes1[1].plot(np.degrees(theta), bfVaryingPitchAplhaFlowSurfRotOut[:,2], 'k', linewidth=2, label='varyingPitch')
axes1[1].plot(np.degrees(theta), bfBaselineAplhaFlowSurfRotOut[:,2], 'k--', linewidth=2, label='varyingPitch')
plt.subplots_adjust(hspace=0.05)
axes1[0].set_ylabel('angular pitch (deg)', fontsize=30)
axes1[1].set_ylabel('flow angle (deg)', fontsize=30)
axes1[1].set_xlabel('$\Theta$ (deg)' , fontsize=30) 
axes1[0].tick_params(labelsize=30)
axes1[1].tick_params(labelsize=30)
plt.savefig('/home/adekola/Documents/New_PhD/papers/ConferenceSubmission/TE2026/BF-OpenFOAM/plots/pitchVariation.png', dpi=500, bbox_inches='tight')
    
#%%
fig, axes1 = plt.subplots(2, 1, figsize=(20, 15), sharex=True)
for ax in axes1:
    # ax.set_ylim(2,17)
    ax.grid(True)
    ax.tick_params(axis='both', labelsize=20)
    
axes1[0].plot(np.degrees(staggerVal[:,0]), staggerVal[:,1], 'k.', label='blade position at TE (pitch)')   
axes1[1].plot(np.degrees(theta), bfVaryingStaggerAplhaFlowSurfRotOut[:,2], 'k', linewidth=2, label='newFormulation')
axes1[1].plot(np.degrees(theta), bfVaryingStaggerOldAplhaFlowSurfRotOut[:,2], 'r--.', linewidth=2, label='oldFormulation')
plt.subplots_adjust(hspace=0.05)
axes1[0].set_ylabel('angular pitch (deg)', fontsize=30)
axes1[1].set_ylabel('flow angle (deg)', fontsize=30)
axes1[1].set_xlabel('$\Theta$ (deg)' , fontsize=30) 
    
    
    
    
    
    
    
    
    
    


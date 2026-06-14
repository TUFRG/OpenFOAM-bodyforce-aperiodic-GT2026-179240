#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 10:27:41 2025

@author: adekola
"""
import os
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve, brentq
#%%
def pol2cart(theta, r, Z):
    import numpy as np
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = Z
    return x, y, z  
def vectorRotZ3D (x,y,z,theta):
    import numpy as np
    xRot = x*np.cos(theta) - y*np.sin(theta)
    yRot = x*np.sin(theta) + y*np.cos(theta)
    return xRot, yRot, z
def meridRot(z, r, theta):
    c = np.cos(theta)
    s = np.sin(theta)
    xRot = r * c   # x in Cartesian
    yRot = r * s   # y in Cartesian
    zRot = z       # unchanged (rotation about z-axis)
    return xRot, yRot, zRot

def yPlusCalc(U, L, yPlus):
    rho = 1.225
    mu = 0.000014607
    Re = rho * U * L / mu
    #using Schlichting eqn https://www.cfd-online.com/Wiki/Skin_friction_coefficient
    Cf = (2*np.log(Re) - 0.65)**-2.3
    delThick = 0.37 * L * Re**-0.2
    tauWall = 0.5 * Cf * rho * U**2
    tauU = np.sqrt(tauWall / rho)
    y1 = yPlus * mu / (rho * tauU)
    return delThick, y1

def makeStl(Xvalues, Yvalues, Zvalues, filename):
    rows = Zvalues.shape[0]  ## Use Z because it is the axis of rotation
    columns = Zvalues.shape[1]
    X = Xvalues
    Y = Yvalues
    Z = Zvalues
    
    numFacets = 0
    
    file = open(filename, 'w')
    file.write('solid \n')
    
    def unitVector(file, p1, p2, p3):
        # VECTORS TANGENT TO FACET
        vector1 = p3 - p2
        vector2 = p3 - p1
        
        normalVec = np.cross(vector1, vector2)
        # magnitude = (normalVec[0]**2 + normalVec[1]**2 + normalVec[2]**2)**0.5
        magnitude = np.linalg.norm(normalVec)
        unitVec = normalVec / magnitude
        file.write(f'facet normal {unitVec[0]} {unitVec[1]} {unitVec[2]} \n'
                   f'outer loop \n'
                   f'vertex {p1[0]} {p1[1]} {p1[2]} \n'
                   f'vertex {p2[0]} {p2[1]} {p2[2]} \n'
                   f'vertex {p3[0]} {p3[1]} {p3[2]} \n'
                   f'endloop \n'
                   f'endfacet \n')
        return
    for i in range(rows - 1):
        for j in range(columns - 1):
            #FACET A VERTICES            
            p1 = np.asarray([X[i, j], Y[i, j], Z[i, j]])
            p2 = np.asarray([X[i, j+1], Y[i, j+1], Z[i, j+1]])
            p3 = np.asarray([X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]])
            
            unitVector(file, p1, p2, p3)
            
            p1 = np.asarray([X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]])
            p2 = np.asarray([X[i+1, j], Y[i+1, j], Z[i+1, j]])
            p3 = np.asarray([X[i, j], Y[i, j], Z[i, j]]) 
            
            unitVector(file, p1, p2, p3)
            
            numFacets += 2
    
    file.write('endsolid')
    file.close()
    print(f'Number of Facets = {numFacets} \n')
    return file

def makeStlNamed(Xvalues, Yvalues, Zvalues, file, solidName):
    """
    Write STL data to an already-open file with named solid
    """
    rows = Zvalues.shape[0]
    columns = Zvalues.shape[1]
    X = Xvalues
    Y = Yvalues
    Z = Zvalues
    numFacets = 0
    
    file.write(f'solid {solidName}\n')
    
    def unitVector(file, p1, p2, p3):
        vector1 = p3 - p2
        vector2 = p3 - p1
        normalVec = np.cross(vector1, vector2)
        magnitude = np.linalg.norm(normalVec)
        unitVec = normalVec / magnitude
        
        file.write(f'facet normal {unitVec[0]} {unitVec[1]} {unitVec[2]} \n'
                   f'outer loop \n'
                   f'vertex {p1[0]} {p1[1]} {p1[2]} \n'
                   f'vertex {p2[0]} {p2[1]} {p2[2]} \n'
                   f'vertex {p3[0]} {p3[1]} {p3[2]} \n'
                   f'endloop \n'
                   f'endfacet \n')
        return
    
    for i in range(rows - 1):
        for j in range(columns - 1):
            p1 = np.asarray([X[i, j], Y[i, j], Z[i, j]])
            p2 = np.asarray([X[i, j+1], Y[i, j+1], Z[i, j+1]])
            p3 = np.asarray([X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]])
            unitVector(file, p1, p2, p3)
            
            p1 = np.asarray([X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]])
            p2 = np.asarray([X[i+1, j], Y[i+1, j], Z[i+1, j]])
            p3 = np.asarray([X[i, j], Y[i, j], Z[i, j]])
            unitVector(file, p1, p2, p3)
            
            numFacets += 2
    
    file.write(f'endsolid {solidName}\n')
    print(f'{solidName}: {numFacets} facets')
    
    return numFacets

def commonRatio(y1, delta, n):
    # Define equation 1: y1*r^n - delta*r + (delta - y1) = 0
    def equation(r):
        if abs(r - 1.0) < 1e-10:
            # L'Hôpital's rule for r→1
            return y1 * n - delta
        return y1 * r**n - delta * r + (delta - y1)
    # Find bounds for r
    # Minimum r: uniform spacing would give delta = y1 * n
    # So r_min is slightly above 1
    r_min = 1.001
    # Maximum r: check if equation has solution in reasonable range
    # For large r: y1*r^n dominates, so r^n ≈ delta*r/y1, r^(n-1) ≈ delta/y1
    r_max_estimate = (delta / y1) ** (1.0 / (n - 1)) * 1.5
    r_max = min(r_max_estimate, 10.0)  # Cap at 10 for safety
    # Check if solution exists in range
    try:
        r = brentq(equation, r_min, r_max)
    except ValueError:
        # Try fsolve with better initial guess
        r_guess = (delta / (y1 * n)) ** (1.0 / (n - 1)) + 1
        result = fsolve(equation, r_guess, full_output=True)
        r = result[0][0]
        if result[2] != 1:  # Check if fsolve converged
            raise ValueError(f"Could not find solution. Check if delta={delta}, y1={y1}, n={n} are compatible")
    # Calculate yn
    yn = y1 * r**(n - 1)
    return r, yn


# def circularArc(zLE, sLE, zTE, sTE, betaLE, betaTE):
#     """
#     Find circular arc passing through (z_LE, s_LE) and (z_TE, s_TE)
#     with tangent angles beta_LE and beta_TE (in degrees).
#     Returns: center (z_c, s_c), radius R
#     """
#     # Convert angles to radians
#     betaLE = np.deg2rad(betaLE)
#     betaTE = np.deg2rad(betaTE)
#     # Normal angles (perpendicular to tangent)
#     alphaLE = betaLE + np.pi/2  # 90° if beta_LE = 0°
#     alphaTE = betaTE + np.pi/2  # 100° if beta_TE = 10°
#     # Normal directions (unit vectors pointing toward center)
#     nLE = np.array([np.cos(alphaLE), np.sin(alphaLE)])
#     nTE = np.array([np.cos(alphaTE), np.sin(alphaTE)])
#     # Solve for center: LE + t1*n_LE = TE + t2*n_TE
#     # This is a 2x2 linear system
#     A = np.column_stack([nLE, -nTE])
#     b = np.array([zTE - zLE, sTE - sLE])
#     # Check if normals are parallel (no solution)
#     if np.abs(np.linalg.det(A)) < 1e-10:
#         raise ValueError("Normals are parallel - no circular arc solution")
#     t = np.linalg.solve(A, b)
#     t1 = t[0]
#     # Center location
#     zC = zLE + t1 * np.cos(alphaLE)
#     sC = sLE + t1 * np.sin(alphaLE)
#     # Radius
#     R = np.sqrt((zLE - zC)**2 + (sLE - sC)**2)
#     # Verify (optional)
#     Rcheck = np.sqrt((zTE - zC)**2 + (sTE - sC)**2)
#     assert np.abs(R - Rcheck) < 1e-6, print(f"Warning: Radius mismatch R={R:.6f}, Rcheck={Rcheck:.6f}")#"Radius mismatch - check calculation"

#     return (zC, sC), R

def ArcPoints(zC, sC, R, zLE, sLE, zTE, sTE, res=101):
    """Generate points along the circular arc"""
    # Angles at LE and TE relative to center
    thetaLE = np.arctan2(sLE - sC, zLE - zC)
    thetaTE = np.arctan2(sTE - sC, zTE - zC)
    # Ensure we go the right direction around the arc
    if thetaTE < thetaLE:
        thetaTE += 2*np.pi
    # Generate points
    theta = np.linspace(thetaLE, thetaTE, res)
    z = zC + R * np.cos(theta)
    s = sC + R * np.sin(theta)
    return z, s

# def plotArcAndCircle(zLE, sLE, zTE, sTE, betaLE, betaTE, res=101):
#     (zC, sC), R = circularArc(zLE, sLE, zTE, sTE, betaLE, betaTE)
    
#     # Arc
#     zArc, sArc = ArcPoints(zC, sC, R, zLE, sLE, zTE, sTE, res)
    
#     # Full circle
#     theta = np.linspace(0, 2*np.pi, 360)
#     zCircle = zC + R * np.cos(theta)
#     sCircle = sC + R * np.sin(theta)
    
#     plt.plot(zArc, sArc, 'b-', linewidth=2, label='Arc')
#     plt.plot(zCircle, sCircle, 'b--', linewidth=1, label='Circle')
#     plt.plot([zLE, zTE], [sLE, sTE], 'ko', markersize=5, label='LE/TE')
#     plt.plot(zC, sC, 'r+', markersize=10, label='Center')
#     plt.axis('equal')
#     plt.legend()
#     plt.grid(True)
#     plt.show()

# def plotArcAndCircle(zLE, sLE, zTE, sTE, betaLE, betaTE, res=101, margin=0.1):
#     (zC, sC), R = circularArc(zLE, sLE, zTE, sTE, betaLE, betaTE)
    
#     # Arc
#     zArc, sArc = ArcPoints(zC, sC, R, zLE, sLE, zTE, sTE, res)
    
#     # Portion of circle — centered around the arc's angular midpoint
#     thetaLE = np.arctan2(sLE - sC, zLE - zC)
#     thetaTE = np.arctan2(sTE - sC, zTE - zC)
#     if thetaTE < thetaLE:
#         thetaTE += 2*np.pi
#     thetaMid = (thetaLE + thetaTE) / 2
#     arcSpan = (thetaTE - thetaLE)
    
#     # Plot a portion wider than the arc itself
#     theta = np.linspace(thetaMid - arcSpan, thetaMid + arcSpan, 360)
#     zCircle = zC + R * np.cos(theta)
#     sCircle = sC + R * np.sin(theta)
    
#     # plt.plot([0,0], [-0.01,0.1], 'r--')
#     # plt.plot([-0.01,0.2], [0,0], 'r--')
#     plt.plot(zArc, sArc, 'b-', linewidth=2, label='Arc')
#     plt.plot(zArc, sArc+0.01, 'b-', linewidth=2, label='Arc')
#     plt.plot(zArc, sArc+0.015, 'b-', linewidth=2, label='Arc')
#     plt.plot(zArc, sArc+0.018, 'b-', linewidth=2, label='Arc')
#     plt.plot(zArc, sArc+0.06, 'b-', linewidth=2, label='Arc') 
#     plt.plot(zArc, sArc+0.1, 'b-', linewidth=2, label='Arc') 
#     # plt.plot(zCircle, sCircle, 'b--', linewidth=1, label='Circle')
#     # plt.plot([zLE, zTE], [sLE, sTE], 'ko', markersize=5, label='LE/TE')
#     plt.xticks([])
#     plt.yticks([])
#     plt.axis('equal')
#     # plt.legend()
#     # plt.grid(True)
#     plt.show()
# Generate arc
def circularArc(zLE, sLE, zTE, sTE, betaLE, betaTE):
    betaLE = np.deg2rad(betaLE)
    betaTE = np.deg2rad(betaTE)
    
    alphaLE = betaLE + np.pi/2
    alphaTE = betaTE + np.pi/2
    
    nLE = np.array([np.cos(alphaLE), np.sin(alphaLE)])
    nTE = np.array([np.cos(alphaTE), np.sin(alphaTE)])
    
    A = np.column_stack([nLE, -nTE])
    b = np.array([zTE - zLE, sTE - sLE])
    
    if np.abs(np.linalg.det(A)) < 1e-10:
        raise ValueError("Normals are parallel - no circular arc solution")
    
    t = np.linalg.solve(A, b)
    t1, t2 = t[0], t[1]
    
    # Compute center from both ends and average for robustness
    zC = 0.5 * ((zLE + t1 * np.cos(alphaLE)) + (zTE + t2 * np.cos(alphaTE)))
    sC = 0.5 * ((sLE + t1 * np.sin(alphaLE)) + (sTE + t2 * np.sin(alphaTE)))
    
    # Average radius from both endpoints
    R = 0.5 * (np.sqrt((zLE - zC)**2 + (sLE - sC)**2) + 
               np.sqrt((zTE - zC)**2 + (sTE - sC)**2))
    
    Rcheck = np.sqrt((zTE - zC)**2 + (sTE - sC)**2)
    if np.abs(np.sqrt((zLE - zC)**2 + (sLE - sC)**2) - Rcheck) > 1e-4:
        print(f"Warning: Residual radius mismatch R_LE={np.sqrt((zLE-zC)**2+(sLE-sC)**2):.6f}, R_TE={Rcheck:.6f}")
    
    return (zC, sC), R


# Call it
# plotBladeStack(zLE=0, sLE=0, zTE=0.1, sTE=0.01, betaLE=0, betaTE=10)
# Example usage for your case:
# zLE = 0.0  # axialBladeIn
# sLE = 0.0  # starting arc length (rHub * 0)
# zTE = 0.1  # axialBladeOut
# sTE = 0.1 * np.tan(np.deg2rad(5))  # from your stagger formula

# (zc, sc), R = circularArc(zLE, sLE, zTE, sTE, 
#                                            betaLE=0, betaTE=10)

# zarc, s_arc = ArcPoints(z_c, s_c, R, z_LE, s_LE, z_TE, s_TE)
# print(f"Center: ({z_c:.6f}, {s_c:.6f})")
# print(f"Radius: {R:.6f}")
# print(f"Camber/Chord: {s_c/np.sqrt((z_TE-z_LE)**2 + (s_TE-s_LE)**2):.6f}")

# plt.plot(z_arc, s_arc, 'k')
# plt.xlabel('axial (z)')
# plt.ylabel('arc length (s)')
# plt.grid()

# plotArcAndCircle(zLE, sLE, zTE, sTE,betaLE=0, betaTE=10)
# plotBladeStack(zLE=0, sLE=0, zTE=0.1, sTE=0.01, betaLE=0, betaTE=10)
#%% input Parameters

bladedOrbodyForce = 1 #if bladed set as 0 if body force set as 1
if bladedOrbodyForce == 0:
    dataPath1 = '../caseSetup/bladedCase/'
else:
    dataPath1 = '../caseSetup/bodyForceCase/'
velocity = 119.73
target_yPlus = 2
xL = 0.1 #characteristic length
res = 101
Nt = 360
Nb = 60 #number of blades
Np = 5 #number of profiles
# theta = np.linspace(-0.25*np.pi, 0.25*np.pi,Nt)
# thetaB = np.linspace(-np.pi, np.pi,Nb)

stagger = 5
axialIn = -1
axialOut = 1
axialBladeIn = 0
axialBladeOut = 0.1
rCas = 1
hubToTipRatio = 0.9
rHub = rCas*hubToTipRatio
ptsX = 10
ptsY = 6
ptsZ = 60
#if no varying pitch or stagger, put 0, if varying pitch but no varying stagger put 1, if no varying pitch but varying stagger put 2, if varying pitch and varying stagger put 3
nonAxiPitchOrStagger = 1 #if no varying pitch or stagger, 
#%% 

thetaB = np.linspace(0, 2*np.pi,Nb, endpoint=False)
meanTheta = 2*np.pi/Nb
# A = 0.25/2.25 #amplitude  ratio of max to min is 1.25
A = 1/2
# perturbation = A * meanTheta * np.sin(thetaB)
# varyingPitch = meanTheta + perturbation
pitch_spacing = meanTheta * (1 + A * np.sin(thetaB))

# Convert to actual blade positions (cumulative sum)
varyingPitch = np.zeros(Nb)
for i in range(1, Nb):
    varyingPitch[i] = varyingPitch[i-1] + pitch_spacing[i-1]

AA =  1/2
meanStagger = stagger 
perturbedStagger = AA * meanStagger * np.sin(thetaB)
perturbedStaggerH = AA * meanStagger * np.sin(thetaB)
varyingStagger = stagger + perturbedStagger


inletHz = np.linspace(axialIn,axialBladeIn,res)
inletCz = np.linspace(axialIn,axialBladeIn,res)
# zC, sC = circularArc(z_LE, s_LE, z_TE, s_TE, betaLE=0, betaTE=10)
# bladeHz, sC = ArcPoints(zC, sC, R, axialBladeIn, 0, axialBladeOut, s_TE)
# bladeHz = np.linspace(axialBladeIn, axialBladeOut, res)
# bladeCz = np.linspace(axialBladeIn, axialBladeOut, res)
outletHz = np.linspace(axialBladeOut, axialOut,res)
outletCz = np.linspace(axialBladeOut, axialOut,res)

inletH = np.column_stack((inletHz, np.full(res, rHub)))
inletC = np.column_stack((inletCz, np.full(res, rCas)))

outletH = np.column_stack((outletHz, np.full(res, rHub)))
outletC = np.column_stack((outletCz, np.full(res, rCas)))

hubSurfaceIn = np.zeros((Nb, res, 3)) #x,y,z
casSurfaceIn = np.zeros((Nb, res, 3))
hubSurfaceBlade = np.zeros((Nb, res, 3))
casSurfaceBlade = np.zeros((Nb, res, 3))
hubSurfaceOut = np.zeros((Nb, res, 3))
casSurfaceOut = np.zeros((Nb, res, 3))
pPeriodicSurface = np.zeros((Nb, 2, res*3-2, 3)) #hub&casing, number of points on each, the coordinates
nPeriodicSurface = np.zeros((Nb, 2, res*3-2, 3))
staggerHub = np.zeros((Nb,res,3)) # This contains all the staggered blades in z,r, theta
staggerCas = np.zeros((Nb,res,3))
allBladeSurface = np.zeros((Nb, Np, res, 3))
radius = np.linspace(rHub, rCas,Np)
exitAngle = np.zeros(Nb)
inletAngle = np.zeros(Nb)
pitchAngle = np.zeros(Nb)
thetaRad = np.zeros(Nb)
pos = np.zeros(Nb)
posIn = np.zeros(Nb)
tempAngle = np.zeros((Nb,5))
if nonAxiPitchOrStagger == 1 or nonAxiPitchOrStagger == 3:
    # For varying pitch, use varyingPitch
    midAngle = np.diff(np.append(varyingPitch, varyingPitch[0] + 2*np.pi))
else:
    # For constant pitch, use thetaB
    midAngle = np.diff(np.append(thetaB, thetaB[0] + 2*np.pi))
if nonAxiPitchOrStagger == 0:
    thetaHIn = thetaB
    thetaCIn = thetaB
    for a in range(Nb):
        sLE = 0
        sTE = (axialBladeOut - axialBladeIn) * np.tan(np.deg2rad(stagger))
        betaTE = 2*stagger
        (zC, sC), R = circularArc(axialBladeIn, sLE, axialBladeOut, sTE, 0, betaTE)
        bladeHz, bladeHs = ArcPoints(zC, sC, R, axialBladeIn, sLE, axialBladeOut, sTE, res)
        staggerHub[a,:,0] = bladeHz #axial direction 
        staggerHub[a,:,1] = np.full(res, rHub) #radius
        staggerHub[a,:,2] = bladeHs #angular shift in theta using s=rtheta
        staggerCas[a,:,0] = bladeHz
        staggerCas[a,:,1] = np.full(res, rCas)
        staggerCas[a,:,2] = bladeHs    
        angleH = staggerHub[a,:,2]/rHub #remove radius to take back to pure angles and not arc length
        angleC = staggerCas[a,:,2]/rCas
        thetHB = thetaB[a] + angleH
        thetCB = thetaB[a] + angleC
        thetaHOut = thetaB[a] + angleH[-1]
        thetaCOut = thetaB[a] + angleC[-1]
        hubSurfaceIn[a,:,:] = np.array(meridRot(inletH[:,0], inletH[:,1], thetaHIn[a])).T
        casSurfaceIn[a,:,:] = np.array(meridRot(inletC[:,0], inletC[:,1], thetaCIn[a])).T
        hubSurfaceBlade[a,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHub[a,:,1], thetHB)).T
        casSurfaceBlade[a,:,:] = np.array(meridRot(staggerCas[a,:,0], staggerCas[a,:,1], thetCB)).T    
        hubSurfaceOut[a,:,:] = np.array(meridRot(outletH[:,0], outletH[:,1], thetaHOut)).T
        casSurfaceOut[a,:,:] = np.array(meridRot(outletC[:,0], outletC[:,1], thetaCOut)).T 
        dz   = staggerHub[a,-1,0] - staggerHub[a,-2,0]   # axial increment (bladeHz)
        drt  = staggerHub[a,-1,2] - staggerHub[a,-2,2]   # increment in (r*theta)
        dzIn   = staggerHub[a,1,0] - staggerHub[a,0,0]   # axial increment (bladeHz)
        drtIn  = staggerHub[a,1,2] - staggerHub[a,0,2]   # increment in (r*theta) 
        inletAngle[a] = np.rad2deg(np.arctan2(drtIn, dzIn))
        exitAngle[a] = np.rad2deg(np.arctan2(drt, dz))
        pos[a] =  thetCB[-1]
        for b in range(Np):
            staggerHubradius = np.full(res, radius[b])
            angleAtRadius = staggerHub[a,:,2]/radius[b]  # Correct angle for this radius
            thetAtRadius = thetaB[a] + angleAtRadius  # or varyingPitch[a] for cases 1,3
            allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetAtRadius)).T
            tempAngle[a,b] = thetAtRadius[-1]
            # allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetHB)).T
        for a in range(Nb):
            prev_blade = (a-1) % Nb
            next_blade = (a+1) % Nb
            pPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[a,:,:], hubSurfaceBlade[a,:,:][1:-1], hubSurfaceOut[a,:,:]))
            pPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[a,:,:], casSurfaceBlade[a,:,:][1:-1], casSurfaceOut[a,:,:]))    
            nPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[prev_blade,:,:], hubSurfaceBlade[prev_blade,:,:][1:-1], hubSurfaceOut[prev_blade,:,:]))
            nPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[prev_blade,:,:], casSurfaceBlade[prev_blade,:,:][1:-1], casSurfaceOut[prev_blade,:,:]))
            pC = allBladeSurface[a,b,0,:]
            pN = allBladeSurface[next_blade,b,0,:]
            rC = np.sqrt(pC[0]**2 + pC[1]**2)
            rN = np.sqrt(pN[0]**2 + pN[1]**2)
            rAvg = 0.5*(rC+rN)
            qC = np.arctan2(pC[1], pC[0])
            qN = np.arctan2(pN[1], pN[0])
            # Handle angle wrapping
            dQ = qN - qC
            if dQ < 0:
                dQ += 2*np.pi
            # pitchAngle[a] = rAvg * np.degrees(dQ)
            pitchAngle[a] = np.degrees(dQ)

            
elif nonAxiPitchOrStagger == 1:
    thetaHIn = varyingPitch
    thetaCIn = varyingPitch
    for a in range(Nb):
        sLE = 0
        sTE = (axialBladeOut - axialBladeIn) * np.tan(np.deg2rad(stagger))
        betaTE = 2*stagger
        (zC, sC), R = circularArc(axialBladeIn, sLE, axialBladeOut, sTE, 0, betaTE)
        bladeHz, bladeHs = ArcPoints(zC, sC, R, axialBladeIn, sLE, axialBladeOut, sTE, res)
        staggerHub[a,:,0] = bladeHz
        staggerHub[a,:,1] = np.full(res, rHub)
        staggerHub[a,:,2] = bladeHs
        staggerCas[a,:,0] = bladeHz
        staggerCas[a,:,1] = np.full(res, rCas)
        staggerCas[a,:,2] = bladeHs
        angleH = staggerHub[a,:,2]/rHub
        angleC = staggerCas[a,:,2]/rCas
        thetHB = varyingPitch[a] + angleH
        # thetaRad[a] = thetHB[-1]
        thetCB = varyingPitch[a] + angleC
        thetaHOut = varyingPitch[a] + angleH[-1]
        thetaCOut = varyingPitch[a] + angleC[-1]
        hubSurfaceIn[a,:,:] = np.array(meridRot(inletH[:,0], inletH[:,1], thetaHIn[a])).T
        casSurfaceIn[a,:,:] = np.array(meridRot(inletC[:,0], inletC[:,1], thetaCIn[a])).T
        hubSurfaceBlade[a,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHub[a,:,1], thetHB)).T
        casSurfaceBlade[a,:,:] = np.array(meridRot(staggerCas[a,:,0], staggerCas[a,:,1], thetCB)).T    
        hubSurfaceOut[a,:,:] = np.array(meridRot(outletH[:,0], outletH[:,1], thetaHOut)).T
        casSurfaceOut[a,:,:] = np.array(meridRot(outletC[:,0], outletC[:,1], thetaCOut)).T 
        dz   = staggerCas[a,-1,0] - staggerCas[a,-2,0]   # axial increment (bladeHz)
        drt  = staggerCas[a,-1,2] - staggerCas[a,-2,2]   # increment in (r*theta)
        dzIn   = staggerHub[a,1,0] - staggerHub[a,0,0]   # axial increment (bladeHz)
        drtIn  = staggerHub[a,1,2] - staggerHub[a,0,2]   # increment in (r*theta) 
        inletAngle[a] = np.rad2deg(np.arctan2(drtIn, dzIn))
        exitAngle[a] = np.rad2deg(np.arctan2(drt, dz))
        for b in range(Np):
            staggerHubradius = np.full(res, radius[b])
            angleAtRadius = staggerHub[a,:,2]/radius[b]  # Correct angle for this radius
            thetAtRadius = varyingPitch[a] + angleAtRadius #for cases 1,3
            allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetAtRadius)).T
            # allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetHB)).T
        for a in range(Nb):
            prev_blade = (a-1) % Nb
            next_blade = (a+1) % Nb
            pPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[a,:,:], hubSurfaceBlade[a,:,:][1:-1], hubSurfaceOut[a,:,:]))
            pPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[a,:,:], casSurfaceBlade[a,:,:][1:-1], casSurfaceOut[a,:,:]))    
            nPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[prev_blade,:,:], hubSurfaceBlade[prev_blade,:,:][1:-1], hubSurfaceOut[prev_blade,:,:]))
            nPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[prev_blade,:,:], casSurfaceBlade[prev_blade,:,:][1:-1], casSurfaceOut[prev_blade,:,:]))
            pC = allBladeSurface[a,b,0,:]
            pN = allBladeSurface[next_blade,b,0,:]
            rC = np.sqrt(pC[0]**2 + pC[1]**2)
            rN = np.sqrt(pN[0]**2 + pN[1]**2)
            rAvg = 0.5*(rC+rN)
            qC = np.arctan2(pC[1], pC[0])
            qN = np.arctan2(pN[1], pN[0])
            # Handle angle wrapping
            dQ = qN - qC
            if dQ < 0:
                dQ += 2*np.pi
            # pitchAngle[a] = rAvg * np.degrees(dQ)
            pitchAngle[a] = np.degrees(dQ)
            
elif nonAxiPitchOrStagger == 2:
    thetaHIn = thetaB
    thetaCIn = thetaB
    # midAngle = 0.5*np.diff(thetaB)
    for a in range(Nb):
        sLE = 0
        sTE = (axialBladeOut - axialBladeIn) * np.tan(np.deg2rad(varyingStagger[a]))
        betaTE = 2*varyingStagger[a]
        (zC, sC), R = circularArc(axialBladeIn, sLE, axialBladeOut, sTE, 0, betaTE)
        bladeHz, bladeHs = ArcPoints(zC, sC, R, axialBladeIn, sLE, axialBladeOut, sTE, res)
        staggerHub[a,:,0] = bladeHz
        staggerHub[a,:,1] = np.full(res, rHub)
        staggerHub[a,:,2] = bladeHs
        staggerCas[a,:,0] = bladeHz
        staggerCas[a,:,1] = np.full(res, rCas)
        staggerCas[a,:,2] = bladeHs
        angleH = staggerHub[a,:,2]/rHub
        angleC = staggerCas[a,:,2]/rCas
        thetHB = thetaB[a] + angleH
        thetCB = thetaB[a] + angleC
        thetaHOut = thetaB[a] + angleH[-1]
        thetaCOut = thetaB[a] + angleC[-1]
        hubSurfaceIn[a,:,:] = np.array(meridRot(inletH[:,0], inletH[:,1], thetaHIn[a])).T
        casSurfaceIn[a,:,:] = np.array(meridRot(inletC[:,0], inletC[:,1], thetaCIn[a])).T
        hubSurfaceBlade[a,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHub[a,:,1], thetHB)).T
        casSurfaceBlade[a,:,:] = np.array(meridRot(staggerCas[a,:,0], staggerCas[a,:,1], thetCB)).T    
        hubSurfaceOut[a,:,:] = np.array(meridRot(outletH[:,0], outletH[:,1], thetaHOut)).T
        casSurfaceOut[a,:,:] = np.array(meridRot(outletC[:,0], outletC[:,1], thetaCOut)).T 
        dz   = staggerCas[a,-1,0] - staggerCas[a,-2,0]   # axial increment (bladeHz)
        drt  = staggerCas[a,-1,2] - staggerCas[a,-2,2]   # increment in (r*theta)
        dzIn   = staggerHub[a,1,0] - staggerHub[a,0,0]   # axial increment (bladeHz)
        drtIn  = staggerHub[a,1,2] - staggerHub[a,0,2]   # increment in (r*theta) 
        inletAngle[a] = np.rad2deg(np.arctan2(drtIn, dzIn))
        exitAngle[a] = np.rad2deg(np.arctan2(drt, dz))
        pos[a] = thetCB[-1]
        posIn[a] = thetCB[0]
        for b in range(Np):
            staggerHubradius = np.full(res, radius[b])
            angleAtRadius = staggerHub[a,:,2]/radius[b]  # Correct angle for this radius
            thetAtRadius = thetaB[a] + angleAtRadius  # or varyingPitch[a] for cases 1,3
            allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetAtRadius)).T
            # allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetaB[a])).T
            # allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetHB)).T
        for a in range(Nb):
            prev_blade = (a-1) % Nb
            next_blade = (a+1) % Nb
            pPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[a,:,:], hubSurfaceBlade[a,:,:][1:-1], hubSurfaceOut[a,:,:]))
            pPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[a,:,:], casSurfaceBlade[a,:,:][1:-1], casSurfaceOut[a,:,:]))    
            nPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[prev_blade,:,:], hubSurfaceBlade[prev_blade,:,:][1:-1], hubSurfaceOut[prev_blade,:,:]))
            nPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[prev_blade,:,:], casSurfaceBlade[prev_blade,:,:][1:-1], casSurfaceOut[prev_blade,:,:]))
            pC = allBladeSurface[a,b,0,:]
            pN = allBladeSurface[next_blade,b,0,:]
            rC = np.sqrt(pC[0]**2 + pC[1]**2)
            rN = np.sqrt(pN[0]**2 + pN[1]**2)
            rAvg = 0.5*(rC+rN)
            qC = np.arctan2(pC[1], pC[0])
            qN = np.arctan2(pN[1], pN[0])
            # Handle angle wrapping
            dQ = qN - qC
            if dQ < 0:
                dQ += 2*np.pi
            # pitchAngle[a] = rAvg * np.degrees(dQ)
            pitchAngle[a] = np.degrees(dQ)
        
# else:
#     thetaHIn = varyingPitch
#     thetaCIn = varyingPitch
#     for a in range(Nb):
#         staggerHub[a,:,0] = bladeHz
#         staggerHub[a,:,1] = np.full(res, rHub)
#         staggerHub[a,:,2] = (bladeHz - axialBladeIn) * np.tan(np.deg2rad(varyingStagger[a]))
#         staggerCas[a,:,0] = bladeCz
#         staggerCas[a,:,1] = np.full(res, rCas)
#         staggerCas[a,:,2] = (bladeCz - axialBladeIn) * np.tan(np.deg2rad(varyingStagger[a]))    
#         angleH = staggerHub[a,:,2]/rHub
#         # angleC = staggerCas[a,:,2]/rCas
#         thetHB = varyingPitch[a] + angleH
#         thetCB = varyingPitch[a] + angleH#angleC
#         thetaHOut = varyingPitch[a] + angleH[-1]
#         thetaCOut = varyingPitch[a] + angleH[-1]#angleC[-1]
#         hubSurfaceIn[a,:,:] = np.array(meridRot(inletH[:,0], inletH[:,1], thetaHIn[a])).T
#         casSurfaceIn[a,:,:] = np.array(meridRot(inletC[:,0], inletC[:,1], thetaCIn[a])).T
#         hubSurfaceBlade[a,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHub[a,:,1], thetHB)).T
#         casSurfaceBlade[a,:,:] = np.array(meridRot(staggerCas[a,:,0], staggerCas[a,:,1], thetCB)).T    
#         hubSurfaceOut[a,:,:] = np.array(meridRot(outletH[:,0], outletH[:,1], thetaHOut)).T
#         casSurfaceOut[a,:,:] = np.array(meridRot(outletC[:,0], outletC[:,1], thetaCOut)).T 
#         dz   = staggerHub[a,-1,0] - staggerHub[a,-2,0]   # axial increment (bladeHz)
#         drt  = staggerHub[a,-1,2] - staggerHub[a,-2,2]   # increment in (r*theta)
#         exitAngle[a] = np.rad2deg(np.arctan2(drt, dz))
#         for b in range(Np):
#             staggerHubradius = np.full(res, radius[b])
#             angleAtRadius = staggerHub[a,:,2]/rHub#radius[b]  # Correct angle for this radius
#             thetAtRadius = varyingPitch[a] + angleAtRadius #for cases 1,3
#             allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetAtRadius)).T
#             # allBladeSurface[a,b,:,:] = np.array(meridRot(staggerHub[a,:,0], staggerHubradius, thetHB)).T
#         for a in range(Nb):
#             prev_blade = (a-1) % Nb
#             next_blade = (a+1) % Nb
#             pPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[a,:,:], hubSurfaceBlade[a,:,:][1:-1], hubSurfaceOut[a,:,:]))
#             pPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[a,:,:], casSurfaceBlade[a,:,:][1:-1], casSurfaceOut[a,:,:]))    
#             nPeriodicSurface[a,0,:,:] = np.vstack((hubSurfaceIn[prev_blade,:,:], hubSurfaceBlade[prev_blade,:,:][1:-1], hubSurfaceOut[prev_blade,:,:]))
#             nPeriodicSurface[a,1,:,:] = np.vstack((casSurfaceIn[prev_blade,:,:], casSurfaceBlade[prev_blade,:,:][1:-1], casSurfaceOut[prev_blade,:,:]))
#             pC = allBladeSurface[a,b,0,:]
#             pN = allBladeSurface[next_blade,b,0,:]
#             rC = np.sqrt(pC[0]**2 + pC[1]**2)
#             rN = np.sqrt(pN[0]**2 + pN[1]**2)
#             rAvg = 0.5*(rC+rN)
#             qC = np.arctan2(pC[1], pC[0])
#             qN = np.arctan2(pN[1], pN[0])
#             # Handle angle wrapping
#             dQ = qN - qC
#             if dQ < 0:
#                 dQ += 2*np.pi
#             # pitchAngle[a] = rAvg * np.degrees(dQ)
#             pitchAngle[a] = np.degrees(dQ)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(pPeriodicSurface[:,0,:,0], pPeriodicSurface[:,0,:,1], pPeriodicSurface[:,0,:,2], 'k.')
#%% save the blade surface/camber 
np.savetxt('../processedData/varyingStagger.txt', exitAngle, delimiter=',')
for c in range(Nb):
    for d in range(Np):
        bladeData = '../processedData/blade{}'.format(c+1)
        np.savetxt(bladeData + '/blade{}.txt'.format(d+1), allBladeSurface[c,d,:,:], delimiter=',')
        
#%% Defining the blocks 
delThick,y1 = yPlusCalc(velocity, xL, target_yPlus)
domHeight = rCas - rHub
Nbc = 0.75*ptsX
# for i in range(Nb):
#     dataPath = dataPath1 + '/case{}/constant/geometry/'.format(i)

#     XhB = np.zeros((Nb, res))
#     YhB = np.zeros((Nb, res))
#     ZhB = np.zeros((Nb, res))   
    
#     XcB = np.zeros((Nb, res))
#     YcB = np.zeros((Nb, res))
#     ZcB = np.zeros((Nb, res))
#     for j in range(Nb):
#         for f in range(res):
#             XhB[j,f] = hubSurfaceBlade[j,f][0]
#             YhB[j,f] = hubSurfaceBlade[j,f][1]
#             ZhB[j,f] = hubSurfaceBlade[j,f][2]
#             XcB[j,f] = casSurfaceBlade[j,f][0]
#             YcB[j,f] = casSurfaceBlade[j,f][1]
#             ZcB[j,f] = casSurfaceBlade[j,f][2]
#     filenamesPer = ['hub', 'casing']
#     Xvalues = [XhB, XcB]
#     Yvalues = [YhB, YcB]
#     Zvalues = [ZhB, ZcB] 
#     for m in range(2):
#         solidName = dataPath + f'{filenamesPer[m]}.stl'
#         X = Xvalues[m]
#         Y = Yvalues[m]
#         Z = Zvalues[m]
#         makeStl(X, Y, Z, solidName)
        #%%
for e in range(Nb):
# for e in range(1):  
#     e=1
    dataPath = dataPath1 + '/case{}/constant/geometry/'.format(e+1)
    XhB = np.zeros((Nb, res))
    YhB = np.zeros((Nb, res))
    ZhB = np.zeros((Nb, res))   

    XcB = np.zeros((Nb, res))
    YcB = np.zeros((Nb, res))
    ZcB = np.zeros((Nb, res))
    for j in range(Nb):
        for f in range(res):
            XhB[j,f] = hubSurfaceBlade[j,f][0]
            YhB[j,f] = hubSurfaceBlade[j,f][1]
            ZhB[j,f] = hubSurfaceBlade[j,f][2]
            XcB[j,f] = casSurfaceBlade[j,f][0]
            YcB[j,f] = casSurfaceBlade[j,f][1]
            ZcB[j,f] = casSurfaceBlade[j,f][2]
    filenamesPer = ['hub', 'casing']
    Xvalues = [XhB, XcB]
    Yvalues = [YhB, YcB]
    Zvalues = [ZhB, ZcB] 
    for m in range(2):
        solidName = dataPath + f'{filenamesPer[m]}.stl'
        X = Xvalues[m]
        Y = Yvalues[m]
        Z = Zvalues[m]
        makeStl(X, Y, Z, solidName)
    Xp = np.zeros((2, res*3-2))
    Yp= np.zeros((2, res*3-2))
    Zp= np.zeros((2, res*3-2)) 

    Xn = np.zeros((2, res*3-2))
    Yn= np.zeros((2, res*3-2))
    Zn= np.zeros((2, res*3-2))
    
    for j in range(2):
        for k in range(res*3-2):
            Xp[j,k] = pPeriodicSurface[e,j,k][0]
            Yp[j,k] = pPeriodicSurface[e,j,k][1]
            Zp[j,k] = pPeriodicSurface[e,j,k][2]
            Xn[j,k] = nPeriodicSurface[e,j,k][0]
            Yn[j,k] = nPeriodicSurface[e,j,k][1]
            Zn[j,k] = nPeriodicSurface[e,j,k][2]
    filenamesPer = ['pPer', 'nPer']
    Xvalues = [Xp, Xn,]
    Yvalues = [Yp, Yn,]
    Zvalues = [Zp, Zn,]  
    for m in range(2):
        path = dataPath #+ '/constant/geometry/'
        if not os.path.exists(path):
            os.makedirs(path)
        solidName = path + f'{filenamesPer[m]}.stl'
        X = Xvalues[m]
        Y = Yvalues[m]
        Z = Zvalues[m]
        makeStl(X, Y, Z, solidName)
    #radial Grading 
    blFrac = 5*delThick/domHeight
    r, yn = commonRatio(y1, 5*delThick, Nbc)
    gradH = yn/y1
    gradC = y1/yn
    newNbc = 1 + np.log(yn/y1)/np.log(r)
            
    Hcore = domHeight - 2*5*delThick
    Ncore = Hcore/yn
    totalCellCount = np.round(Ncore + 2*newNbc)
    
    z1P = blFrac
    z2P = blFrac
    z3P = blFrac
    
    z1N = newNbc/totalCellCount
    z2N = newNbc/totalCellCount
    z3N = newNbc/totalCellCount
    
    z1upG = gradH
    z2upG = gradH
    z3upG = gradH
    
    z1dwG = gradC
    z2dwG = gradC
    z3dwG = gradC
    
    #Axial Grading 
    #for block1
    blkLen1 = axialBladeIn - axialIn
    NbcX1 = 0.5*(ptsZ)
    blkFrac1 = (0.25*blkLen1)/blkLen1
    r1, yn1 = commonRatio(y1, 0.25*blkLen1, NbcX1)
    gradHBlk1 = yn1/y1
    gradCBlk1 = y1/yn1
    newNbcX1 = 1 + np.log(yn1/y1)/np.log(r1)
    HcoreBlk1 = blkLen1 - 0.25*blkLen1
    NcoreBlk1 = HcoreBlk1/yn1
    blk1CellCount = NcoreBlk1 + newNbcX1
    #for block2
    blkLen2 = axialBladeOut - axialBladeIn
    NbcX2 = 0.3*(ptsZ)
    blkFrac2 = 0.25*blkLen2/blkLen2
    r2, yn2 = commonRatio(y1, 0.25*blkLen2, NbcX2)
    gradHBlk2 = yn2/y1
    gradCBlk2 = y1/yn2
    newNbcX2 = 1 + np.log(yn2/y1)/np.log(r2)
    HcoreBlk2 = blkLen2 - 2*0.25*blkLen2
    NcoreBlk2 = HcoreBlk2/yn2
    blk2CellCount = NcoreBlk2 + 2*newNbcX2
    #for block3
    blkLen3 = axialOut - axialBladeOut
    NbcX3 = 0.5*(ptsZ)
    blkFrac3 = 0.25*blkLen3/blkLen3
    r3, yn3 = commonRatio(y1,0.25*blkLen3 , NbcX3)
    gradHBlk3 = yn3/y1
    gradCBlk3 = y1/yn3
    newNbcX3 = 1 + np.log(yn3/y1)/np.log(r3)
    HcoreBlk3 = blkLen3 - 0.25*blkLen3
    NcoreBlk3 = HcoreBlk3/yn3
    blk3CellCount = NcoreBlk3 + newNbcX3
    
    # axialBlkCount = blk1CellCount + blk2CellCount + blk3CellCount
    
    ptX1 = np.round(blk1CellCount)
    ptX2 = np.round(blk2CellCount)
    ptX3 = np.round(blk3CellCount)
    
    axialBlkCount = ptX1 + ptX2 + ptX3
    print('This is the total cell Count {}'.format(axialBlkCount*totalCellCount))
    
    x1P = blkFrac1
    x2P = blkFrac2
    x3P = blkFrac3
    
    x1N = newNbcX1/blk1CellCount
    x2N = newNbcX2/blk2CellCount
    x3N = newNbcX3/blk3CellCount
    
    x1upG = gradCBlk1
    x2upG = gradCBlk2
    x3upG = gradCBlk3
    
    x1dwG = gradHBlk1
    x2dwG = gradHBlk2
    x3dwG = gradHBlk3
    
    #START HERE, complete the circumferential grading 
    #Circumferential Grading 
    Nbq = 0.3*ptsY
    #use atan, i don't need the angles to go beyond 90 degrees, I am working locally
    # midAngle = abs(0.5*(np.atan(pPeriodicSurface[e][0][0][1]/ pPeriodicSurface[e][0][0][0]) - np.atan(nPeriodicSurface[e][0][0][1]/ nPeriodicSurface[e][0][0][0])))
    circumLengthC = rCas*midAngle[e-1] #s=rTheta
    
    
    blFracC = 2*delThick/circumLengthC
    r, ynC = commonRatio(y1, 2*delThick, Nbq)
    gradCp = ynC/y1
    gradCn = y1/ynC
            
    Hcore = circumLengthC - 2*2*delThick
    Ncore = Hcore/ynC
    totalCellCountC = np.round(Ncore + 2*Nbq)
    
    y1P = blFracC
    y2P = blFracC
    y3P = blFracC
    
    y1N = Nbq/totalCellCountC
    y2N = Nbq/totalCellCountC
    y3N = Nbq/totalCellCountC
    
    y1p = gradCp
    y2p = gradCp
    y3p = gradCp
    
    y1n = gradCn
    y2n = gradCn
    y3n = gradCn
    
    xInH = axialIn #- 1.5*radius
    rInH = rHub
    xInC = axialIn #- 1.5*radius
    rInC = rCas
    xLEH = axialBladeIn
    rLEH = rHub
    xLEC = axialBladeIn
    rLEC = rCas
    xTEH = axialBladeOut
    rTEH = rHub
    xTEC = axialBladeOut
    rTEC = rCas
    xOutH = axialOut #+ 1.5*radius
    rOutH = rHub
    xOutC = axialOut #+ 1.5*radius
    rOutC = rCas
    
    
    x1r0p = pPeriodicSurface[e][0][0]
    x1r1p = pPeriodicSurface[e][1][0]
    x2r0p = pPeriodicSurface[e][0][res-1]
    x2r1p = pPeriodicSurface[e][1][res-1]
    x3r0p = pPeriodicSurface[e][0][2*res-2]
    x3r1p = pPeriodicSurface[e][1][2*res-2]
    x4r0p = pPeriodicSurface[e][0][-1]
    x4r1p = pPeriodicSurface[e][1][-1]
    
    x1r0n = nPeriodicSurface[e][0][0]
    x1r1n = nPeriodicSurface[e][1][0]
    x2r0n = nPeriodicSurface[e][0][res-1]
    x2r1n = nPeriodicSurface[e][1][res-1]
    x3r0n = nPeriodicSurface[e][0][2*res-2]
    x3r1n = nPeriodicSurface[e][1][2*res-2]
    x4r0n = nPeriodicSurface[e][0][-1]
    x4r1n = nPeriodicSurface[e][1][-1]
    
    # midAngle1 = 0.5*(np.atan2(x1r0n[1], x1r0n[0]) + np.atan2(x1r0p[1], x1r0p[0]))
    angle_p = np.atan2(x1r0p[1], x1r0p[0])
    angle_n = np.atan2(x1r0n[1], x1r0n[0])
    
    # Unwrap: ensure angles are within 2π of each other
    angle_diff = angle_p - angle_n
    if angle_diff > np.pi:
        angle_n += 2*np.pi
    elif angle_diff < -np.pi:
        angle_p += 2*np.pi
    
    midAngle1 = 0.5*(angle_p + angle_n)
    # x1r0m = [rHub*np.cos(midAngle[e-1]), rHub*np.sin(midAngle[e-1]), axialIn] 
    # x1r1m = [rCas*np.cos(midAngle[e-1]), rCas*np.sin(midAngle[e-1]), axialIn]
    # x2r0m = [rHub*np.cos(midAngle[e-1]), rHub*np.sin(midAngle[e-1]), axialBladeIn] 
    # x2r1m = [rCas*np.cos(midAngle[e-1]), rCas*np.sin(midAngle[e-1]), axialBladeIn]
    # x3r0m = [rHub*np.cos(midAngle[e-1]), rHub*np.sin(midAngle[e-1]), axialBladeOut] 
    # x3r1m = [rCas*np.cos(midAngle[e-1]), rCas*np.sin(midAngle[e-1]), axialBladeOut]
    # x4r0m = [rHub*np.cos(midAngle[e-1]), rHub*np.sin(midAngle[e-1]), axialOut]
    # x4r1m = [rCas*np.cos(midAngle[e-1]), rCas*np.sin(midAngle[e-1]), axialOut]
    
    x1r0m = [rHub*np.cos(midAngle1), rHub*np.sin(midAngle1), axialIn] 
    x1r1m = [rCas*np.cos(midAngle1), rCas*np.sin(midAngle1), axialIn]
    x2r0m = [rHub*np.cos(midAngle1), rHub*np.sin(midAngle1), axialBladeIn] 
    x2r1m = [rCas*np.cos(midAngle1), rCas*np.sin(midAngle1), axialBladeIn]
    x3r0m = [rHub*np.cos(midAngle1), rHub*np.sin(midAngle1), axialBladeOut] 
    x3r1m = [rCas*np.cos(midAngle1), rCas*np.sin(midAngle1), axialBladeOut]
    x4r0m = [rHub*np.cos(midAngle1), rHub*np.sin(midAngle1), axialOut]
    x4r1m = [rCas*np.cos(midAngle1), rCas*np.sin(midAngle1), axialOut]   
    print(e)
    f = open(dataPath1 + '/case{}/system/pointData'.format(e+1), 'w')
    
    f.write('x1r0m ({} {} {}); \n'.format(x1r0m[0], x1r0m[1], x1r0m[2]))
    f.write('x1r1m ({} {} {}); \n'.format(x1r1m[0], x1r1m[1], x1r1m[2]))
    f.write('x2r0m ({} {} {}); \n'.format(x2r0m[0], x2r0m[1], x2r0m[2]))
    f.write('x2r1m ({} {} {}); \n'.format(x2r1m[0], x2r1m[1], x2r1m[2]))
    f.write('x3r0m ({} {} {}); \n'.format(x3r0m[0], x3r0m[1], x3r0m[2]))
    f.write('x3r1m ({} {} {}); \n'.format(x3r1m[0], x3r1m[1], x3r1m[2]))
    f.write('x1r0p ({} {} {}); \n'.format(x1r0p[0], x1r0p[1], x1r0p[2]))
    f.write('x1r0n ({} {} {}); \n'.format(x1r0n[0], x1r0n[1], x1r0n[2]))
    f.write('x1r1p ({} {} {}); \n'.format(x1r1p[0], x1r1p[1], x1r1p[2]))
    f.write('x1r1n ({} {} {}); \n'.format(x1r1n[0], x1r1n[1], x1r1n[2]))
    f.write('x2r0p ({} {} {}); \n'.format(x2r0p[0], x2r0p[1], x2r0p[2]))
    f.write('x2r0n ({} {} {}); \n'.format(x2r0n[0], x2r0n[1], x2r0n[2]))
    f.write('x2r1p ({} {} {}); \n'.format(x2r1p[0], x2r1p[1], x2r1p[2]))
    f.write('x2r1n ({} {} {}); \n'.format(x2r1n[0], x2r1n[1], x2r1n[2]))
    f.write('x3r0p ({} {} {}); \n'.format(x3r0p[0], x3r0p[1], x3r0p[2]))
    f.write('x3r0n ({} {} {}); \n'.format(x3r0n[0], x3r0n[1], x3r0n[2]))
    f.write('x3r1p ({} {} {}); \n'.format(x3r1p[0], x3r1p[1], x3r1p[2]))
    f.write('x3r1n ({} {} {}); \n'.format(x3r1n[0], x3r1n[1], x3r1n[2]))
    f.write('x4r0p ({} {} {}); \n'.format(x4r0p[0], x4r0p[1], x4r0p[2]))
    f.write('x4r0n ({} {} {}); \n'.format(x4r0n[0], x4r0n[1], x4r0n[2]))
    f.write('x4r1p ({} {} {}); \n'.format(x4r1p[0], x4r1p[1], x4r1p[2]))
    f.write('x4r1n ({} {} {}); \n'.format(x4r1n[0], x4r1n[1], x4r1n[2]))
    f.write('x4r0m ({} {} {}); \n'.format(x4r0m[0], x4r0m[1], x4r0m[2]))
    f.write('x4r1m ({} {} {}); \n'.format(x4r1m[0], x4r1m[1], x4r1m[2]))
    
    f.write('z1P {}; \n'.format(z1P))
    f.write('z2P {}; \n'.format(z2P))
    f.write('z3P {}; \n'.format(z3P))
    f.write('z1N {}; \n'.format(z1N))
    f.write('z2N {}; \n'.format(z2N))
    f.write('z3N {}; \n'.format(z3N))
    f.write('z1upG {}; \n'.format(z1upG))
    f.write('z2upG {}; \n'.format(z2upG))
    f.write('z3upG {}; \n'.format(z3upG))
    f.write('z1dwG {}; \n'.format(z1dwG))
    f.write('z2dwG {}; \n'.format(z2dwG))
    f.write('z3dwG {}; \n'.format(z3dwG))
    f.write('y1P {}; \n'.format(y1P))
    f.write('y2P {}; \n'.format(y2P))
    f.write('y3P {}; \n'.format(y3P))
    f.write('y1N {}; \n'.format(y1N))
    f.write('y2N {}; \n'.format(y2N))
    f.write('y3N {}; \n'.format(y3N))
    f.write('y1p {}; \n'.format(y1p))
    f.write('y2p {}; \n'.format(y2p))
    f.write('y3p {}; \n'.format(y3p))
    f.write('y1n {}; \n'.format(y1n))
    f.write('y2n {}; \n'.format(y2n))
    f.write('y3n {}; \n'.format(y3n))
    f.write('x1P {}; \n'.format(x1P))
    f.write('x3P {}; \n'.format(x2P))
    f.write('x3P {}; \n'.format(x3P))
    f.write('x1N {}; \n'.format(x1N))
    f.write('x2N {}; \n'.format(x2N))
    f.write('x3N {}; \n'.format(x3N))
    f.write('x1upG {}; \n'.format(x1upG))
    f.write('x2upG {}; \n'.format(x2upG))
    f.write('x3upG {}; \n'.format(x3upG))
    f.write('x1dwG {}; \n'.format(x1dwG))
    f.write('x2dwG {}; \n'.format(x2dwG))
    f.write('x3dwG {}; \n'.format(x3dwG))
    f.write('ptX1 {}; \n'.format(ptX1))
    f.write('ptX2 {}; \n'.format(ptX2))
    f.write('ptX3 {}; \n'.format(ptX3))
    f.write('ptsZ {}; \n'.format(totalCellCount))
    f.write('ptsY {}; \n'.format(totalCellCountC))
    f.close()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 13:14:44 2025

@author: adekola
"""
import argparse
import numpy as np
from scipy.interpolate import CubicSpline, interp1d
#%% Load input data
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
args = parser.parse_args()
# Use parsed arguments
Nbr = args.Nbr  # Total number of blades in the annulus for rotor
Nbs = args.Nbs  # Total number of blades in the annulus for stator
rSections = args.rSections  # Number of rotor blade profiles
sSections = args.sSections  # Number of stator blade profiles
Nr = args.Nr  # Number of points on rotor blade profiles
Ns = args.Ns  # Number of points on stator blade profiles
periodicOrAperiodic = args.periodicOrAperiodic  # if periodic select 0 otherwise select 1

Np = 10 #Number of interpolated profiles between each blade profiles
N = 101 #number of points on thickness
distro = 'cosine'
#%%
if periodicOrAperiodic == 0:
    dataPath = '../processedData/periodic/'
else:
    dataPath = '../processedData/nonPeriodic/'
    
    
rBlade = np.zeros([Nbr, rSections, Nr, 3])
sBlade = np.zeros([Nbs, sSections, Ns, 3])
for a in range(Nbr):
    for b in range(rSections):
        if periodicOrAperiodic == 0:
            rBlade[a,b,:] = np.loadtxt(dataPath + '/rotor/blade{}/blade{}.txt'.format(a,b), delimiter=',')
            outputDatapath = '../inputData/periodic/'
        else:
            rBlade[a,b,:] = np.loadtxt(dataPath + '/rotor/blade{}/blade{}.txt'.format(a,b), delimiter=',')          
            outputDatapath = '../inputData/nonPeriodic/'
            
for c in range(Nbs):
    for d in range(sSections):
        if periodicOrAperiodic == 0:
            sBlade[c,d,:] = np.loadtxt(dataPath + '/stator/blade{}/blade{}.txt'.format(c,d), delimiter=',')
            outputDatapath = '../inputData/periodic/'
        else:
            sBlade[c,d,:] = np.loadtxt(dataPath + '/stator/blade{}/blade{}.txt'.format(c,d), delimiter=',')           
            outputDatapath = '../inputData/nonPeriodic/'

#%% Used functions
def vectorRotZ3D (x,y,z,theta):
    xRot = x*np.cos(theta) - y*np.sin(theta)
    yRot = x*np.sin(theta) + y*np.cos(theta)
    return xRot, yRot, z
def cart2pol(x,y,z):
    import numpy as np
    theta = np.arctan2(y,x)
    rho = np.sqrt(np.square(x)+np.square(y))
    Z = z
    return(theta, rho, Z)
def cart2pol2(x,y,z):
    import numpy as np
    theta = np.arctan2(y,x)
    # theta = (theta+2*np.pi) % 2*np.pi
    theta = np.unwrap(theta)
    rho = np.sqrt(np.square(x)+np.square(y))
    Z = z
    return(theta, rho, Z)
def separate_blade_surfaces(blade_cylindrical):
    """Separate blade into pressure and suction surfaces"""
    min_axial_i = np.argmin(blade_cylindrical[:, 2])
    max_axial_i = np.argmax(blade_cylindrical[:, 2])
    if min_axial_i == 0:
        s1 = blade_cylindrical[:max_axial_i+1]
        s2 = blade_cylindrical[max_axial_i:]
        s2 = np.flip(s2, 0)
    elif (min_axial_i < max_axial_i) and (min_axial_i != 0):
        s1 = blade_cylindrical[min_axial_i:max_axial_i+1]
        s2 = np.concatenate((blade_cylindrical[max_axial_i:-1], 
                           blade_cylindrical[0:min_axial_i+1]))
        s2 = np.flip(s2, 0)
    else:
        s1 = blade_cylindrical[max_axial_i:min_axial_i+1]
        s2 = np.concatenate((blade_cylindrical[min_axial_i:-1], 
                             blade_cylindrical[0:max_axial_i+1]))
        s1 = np.flip(s1, 0)
    
    # Check if axial coordinates are monotonic
    if not (np.all(np.diff(s1[:, 2]) >= 0) or np.all(np.diff(s1[:, 2]) <= 0)):
        print("WARNING: Surface 1 is not monotonic in axial direction!")
    if not (np.all(np.diff(s2[:, 2]) >= 0) or np.all(np.diff(s2[:, 2]) <= 0)):
        print("WARNING: Surface 2 is not monotonic in axial direction!")
    return s1, s2

def calculate_thickness_distribution(bladeData, num_profiles, newLETE, n):
    """Calculate thickness distribution including theta coordinates"""
    thicknessData = np.zeros((num_profiles,n, 4))
    for blade_num, blade_cylindrical in enumerate(bladeData):
        if blade_cylindrical.size == 0:
            continue         
        s1, s2 = separate_blade_surfaces(blade_cylindrical)
        s1_theta_spline = interp1d(s1[:, 2], s1[:, 0],fill_value='extrapolate')  # theta vs z  
        s2_theta_spline = interp1d(s2[:, 2], s2[:, 0], fill_value='extrapolate')
        mn = blade_cylindrical[np.argmin(blade_cylindrical[:, 2]), 2]
        mx = blade_cylindrical[np.argmax(blade_cylindrical[:, 2]), 2]
        if distro == 'cosine':
            axial_distro = (mx - mn) * (0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))) + mn
        else:
            axial_distro = np.linspace(mn, mx, n)    
        s1_theta = s1_theta_spline(axial_distro) 
        s2_theta = s2_theta_spline(axial_distro) 
        thickness = np.abs(s1_theta - s2_theta)

        new_mn = newLETE[blade_num, 2]  # column 2 is new LE axial
        new_mx = newLETE[blade_num, 5]  # column 5 is new TE axial
        if distro == 'cosine':
            new_axial_distro = (new_mx - new_mn) * (0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))) + new_mn
        else:
            new_axial_distro = np.linspace(new_mn, new_mx, n)        
        r_interp = interp1d(blade_cylindrical[:,2], blade_cylindrical[:,1], fill_value='extrapolate') #handles floating point problems
        theta_interp = interp1d(blade_cylindrical[:,2], blade_cylindrical[:,0], fill_value='extrapolate')       
        new_mean_r = r_interp(new_axial_distro)
        new_mean_theta = theta_interp(new_axial_distro)
        # Store results
        thicknessData[blade_num, :, 0] = new_mean_theta      # theta
        thicknessData[blade_num, :, 1] = new_mean_r            # r (radial)
        thicknessData[blade_num, :, 2] = new_axial_distro        # x (axial)
        thicknessData[blade_num, :, 3] = thickness         # thickness   
    return thicknessData

def structuredInterpolation(thickness_data, num_sections, num_points, Np):
    n_blades = thickness_data.shape[0]
    # Calculate interpolated sections: ((original_sections-1)*Np)+1
    numProfilesInterp = ((num_sections-1)*Np)+1
    interpolated_data = np.zeros((n_blades, numProfilesInterp, num_points, 4)) 
    for blade in range(n_blades):
        thicknessCompressed = thickness_data[blade]  # [sections, points, 4] (theta, r, z, thickness)
        thicknessCompressedInterp = np.zeros((numProfilesInterp, num_points, 4))  # (z, r, theta, thickness, gradX, gradR)
        # Place existing values every Np rows in new array
        for r_idx in range(0, numProfilesInterp):
            if (r_idx % Np) == 0:
                thicknessCompressedInterp[r_idx, :, 0:4] = thicknessCompressed[int(r_idx/Np), :, 0:4]
        for p in range(0, num_points):
            for q in range(0, num_sections-1):
                rStart = thicknessCompressed[q, p, 1]  # r coordinate
                rStop = thicknessCompressed[q+1, p, 1]  # r coordinate  
                newR = np.linspace(rStart, rStop, 11)
                thicknessCompressedInterp[q*Np:q*Np+11, p, 1] = newR
        # Interpolate for new axial, theta, and thickness values
        for p in range(0, num_points):
            oldR = thicknessCompressed[:, p, 1]  # Original r coordinates
            newR = thicknessCompressedInterp[:, p, 1]  # New r coordinates
            cols = np.r_[2, 0, 3]  # new columns to fill in: z(0), theta(2), t_hat(3)
            old_cols = np.r_[2, 0, 3]  # old columns to pull data from: z, theta, thickness
            for s in range(0, len(cols)):
                old_col = old_cols[s]
                oldVar = thicknessCompressed[:, p, old_col]
                # csFunc = CubicSpline(oldR, oldVar) 
                csFunc = interp1d(oldR, oldVar) 
                newVar = csFunc(newR)
                col = cols[s]
                thicknessCompressedInterp[:, p, col] = newVar
        interpolated_data[blade] = thicknessCompressedInterp[:, :, :]
    return interpolated_data

# Function to reshape and save interpolated data
def saveThicknessData(rotor_interp, stator_interp, output_path):
    # Reshape to 2D for saving
    n_blades_r, n_sections_r, n_points_r, _ = rotor_interp.shape
    n_blades_s, n_sections_s, n_points_s, _ = stator_interp.shape
    
    t_hat_r = rotor_interp.reshape(n_blades_r * n_sections_r * n_points_r, 4)
    t_hat_s = stator_interp.reshape(n_blades_s * n_sections_s * n_points_s, 4)
    
    # Save only the first 4 columns [theta, r, z, thickness]
    np.savetxt(output_path + '/Rt_hat.txt', t_hat_r, delimiter=',', 
               header='theta,r,z,thickness')
    np.savetxt(output_path + '/St_hat.txt', t_hat_s, delimiter=',',
               header='theta,r,z,thickness')    
    return t_hat_r, t_hat_s
            
#%% Convert rotor blade data from Cartesian to Cylindrical
rBladeCyl = np.zeros([Nbr, rSections, Nr, 3])
rLETE = np.zeros([Nbr,rSections,6])
sLETE = np.zeros([Nbs,sSections,6]) 
for e in range(Nbr):
    for ee in range(rSections):
        rBladeCyl[e,ee,:,:] = np.array(cart2pol(rBlade[e,ee,:,0], rBlade[e,ee,:,1], rBlade[e,ee,:,2])).T 

        if max(rBladeCyl[e,ee,:,0]) - min(rBladeCyl[e,ee,:,0]) > np.pi:
            rBladeCyl[e,ee,:,:] = np.array(cart2pol2(rBlade[e,ee,:,0], rBlade[e,ee,:,1], rBlade[e,ee,:,2])).T 
        minIdx = np.argmin(rBladeCyl[e,ee,:,2]) 
        maxIdx = np.argmax(rBladeCyl[e,ee,:,2]) 
        rLETE[e,ee] = np.hstack((rBladeCyl[e,ee,minIdx], rBladeCyl[e,ee,maxIdx]))


# Convert stator blade data from Cartesian to Cylindrical
sBladeCyl = np.zeros([Nbs, sSections, Ns, 3])
for g in range(Nbs):
    for h in range(sSections):
        sBladeCyl[g,h,:,:] = np.array(cart2pol(sBlade[g,h,:,0], sBlade[g,h,:,1], sBlade[g,h,:,2])).T

        if max(sBladeCyl[g,h,:,0]) - min(sBladeCyl[g,h,:,0]) > np.pi:
            sBladeCyl[g,h,:,:] = np.array(cart2pol2(sBlade[g,h,:,0], sBlade[g,h,:,1], sBlade[g,h,:,2])).T
        minIdx = np.argmin(sBladeCyl[g,h,:,2]) 
        maxIdx = np.argmax(sBladeCyl[g,h,:,2])         
        sLETE[g,h] = np.hstack((sBladeCyl[g,h,minIdx], sBladeCyl[g,h,maxIdx]))
#%% Determine blade thickness based on Profiles 
rotorThickness = np.zeros([Nbr, rSections, N, 4])
statorThickness = np.zeros([Nbs, sSections, N, 4])

for j in range(Nbs):
    statorThickness[j,:,:,:] = calculate_thickness_distribution(sBladeCyl[j], sSections, sLETE[j], N)
    
for h in range(Nbr):
    rotorThickness[h,:,:,:] = calculate_thickness_distribution(rBladeCyl[h], rSections, rLETE[h], N)

#%% Now do the radial interpolation to ensure that the blockage is representative of the blade
rotorInterp = structuredInterpolation(rotorThickness, rSections, N, Np)
statorInterp = structuredInterpolation(statorThickness, sSections, N, Np)

rotor, stator = saveThicknessData(rotorInterp, statorInterp, outputDatapath)
















#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2024-09-18

@author: Amelia George

Modified on 2025-12-17

@author: Adekola Adeyemi
"""

import numpy as np
import os
from scipy import interpolate
import matplotlib.pyplot as plt
import geometryFunctions as mf
import LE_TE_Solver as lts
import camberFunctions as cf

def cart2pol(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)#%(2*np.pi)
    z = z
    return theta, r, z
def cart2pol2(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    theta = np.unwrap(theta)
    z = z
    return theta, r, z

#%%

# for a in range (numProfiles):
def getCamber(blade, n, ):
    newLETE = np.zeros(4)
    distro = 'cosine' #point distribution on generated camber line
    nIter = 10 #number of iterations in accurate camber loop
    tolIter = 0.001 #tolerance of accurate camber line
    tolLETE = 0.001 #tolerance of movement of LE and TE
    NLETE = 10 #number of loops of overall camber generation
    # blade = np.loadtxt(filePath + folderSelect + 'blade{}.txt'.format(a), delimiter=',')

    # Also get full blade in cylindrical coordinates for interpolation later
    bladeCylindrical = np.array(cart2pol(blade[:,0], blade[:,1], blade[:,2])).T
    if max(bladeCylindrical[:,0]) - min(bladeCylindrical[:,0]) > np.pi:
        bladeCylindrical = np.array(cart2pol2(blade[:,0], blade[:,1], blade[:,2])).T
    theta_mean = np.arctan2(np.mean(blade[:,1]), np.mean(blade[:,0]))
    # Rotate blade so its mean angle aligns with the x-axis
    cos_t, sin_t = np.cos(-theta_mean), np.sin(-theta_mean)
    x_rot = blade[:,0]*cos_t - blade[:,1]*sin_t
    y_rot = blade[:,0]*sin_t + blade[:,1]*cos_t
    # Now x_rot ≈ r (radial), y_rot ≈ r*dθ (tangential), blade[:,2] = axial
    # Use (z, y_rot) as the 2D working plane — no branch cut possible
    blade2D = np.column_stack((blade[:,2], y_rot))  
    #blade2D = np.column_stack((bladeCylindrical[:,2], bladeCylindrical[:,1]*bladeCylindrical[:,0]))

    # Find initial LE and TE indices based on furthest apart two points
    [dis, i, j] = mf.longest_distance_on_curve(blade2D[:,0], blade2D[:,1])

    #assign LE and TE values based on indices and axial location
    if blade2D[i,0] < blade2D [j,0]:
        le = i
        te = j
    else:
        le = j
        te = i
    LE = blade2D[le]
    TE = blade2D[te]
    
    #define largest thickness on blade as normalization value for certain calculations below
    maxThickness = max(blade2D[:,1]) - min(blade2D[:,1])

    # get approximate camber
    US, LS, camber = cf.getApproxCamber(blade2D, le, te, distro, n)
    PS = LS
    SS = US

    tol = 10
    count = 0 
    while tol >= tolLETE:
        if count > NLETE:
            break
        print('\nOverall camber iteration: ' + (str(count)))
        PS, SS, LEnew, TEnew, camberlineNewLETE = lts.LE_TE_Solver(PS, SS, LE, TE, camber, distro, n)
        camber = cf.getAccurateCamber(blade2D, camberlineNewLETE, tolIter, nIter, maxThickness, dis)
        LEdiff = (mf.dist2D(LE[0],LE[1],LEnew[0],LEnew[1]))
        TEdiff = (mf.dist2D(TE[0],TE[1],TEnew[0],TEnew[1]))
        tol = ((LEdiff + TEdiff)/2)/dis
        print('Overall camber tol: ' + str(tol))
        #LETE_stored += ((LEnew,TEnew))
        LE = LEnew
        TE = TEnew
        count += 1
    newLETE[0] = LEnew[0]
    newLETE[1] = LEnew[1]
    newLETE[2] = TEnew[0]
    newLETE[3] = TEnew[1]
    PS, SS, LEnew, TEnew, camber = lts.LE_TE_Solver(PS, SS, LE, TE, camber, distro, n)
    
    # if plotGraphs == 1:
    #     plt.plot(blade2D[:,0], blade2D[:,1], 'k', linewidth = 0.1)
    #     plt.plot(camber[:,0], camber[:,1], 'g', linewidth = 0.1)
    #     plt.show()

    nNew = len(camber)
    #convert camber lines back into cylindrical coordinates (columns: r-theta-axial)
    camberCylindrical = np.empty([nNew, 3])
    #axial coordiates are 1-1, no conversion needed
    camberCylindrical[:,2] = camber[:,0]
    #create interpolation object with axial and radial points of the original cylindrical coordinates of the blade data
    f = interpolate.interp1d(bladeCylindrical[:,2],bladeCylindrical[:,1])
    #fill new radial coordinates using interpolation object and new axial coordinates
    camberCylindrical[:,1] = f(camberCylindrical[:,2])
    y_rot_camber = camber[:,1]
    r_camber = camberCylindrical[:,1]
    x_rot_camber = np.sqrt(np.maximum(r_camber**2 - y_rot_camber**2, 0))
    
    # theta in local frame
    theta_local = np.arctan2(y_rot_camber, x_rot_camber)
    
    # Rotate back to global frame
    theta_global = theta_local + theta_mean
    camberCylindrical[:,0] = theta_global
    #theta coordinate obtained by dividing radial coordinate from r-theta coordinate
    #camberCylindrical[:,0] = (camber[:,1]) / (camberCylindrical[:,1])

    #convert camber lines to a-b-c coordinates for MATLAB camber surface normals generation
    camberCartesian = np.empty([nNew, 3])
    #axial coordiates are 1-1, no conversion needed
    camberCartesian[:,2] = camberCylindrical[:,2]
    #standard cylindrical transformation to get a and b coordinates from r and theta
    camberCartesian[:,0] = camberCylindrical[:,1] * np.cos(camberCylindrical[:,0])
    camberCartesian[:,1] = camberCylindrical[:,1] * np.sin(camberCylindrical[:,0])
    # camberCart[a,:,:] = camberCylindrical
    # camberCart[a,:,:] = camberCartesian
    return camberCartesian, newLETE
    





























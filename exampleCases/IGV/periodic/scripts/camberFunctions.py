#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2024-09-18

@author: Amelia George
"""

#imports embedded in each function

#%%
def getApproxCamber(blade, le, te, distro, n):

    from scipy.interpolate import CubicSpline
    import numpy as np
    import geometryFunctions as mf

    LE = blade[le]
    TE = blade[te]

    theta = np.arctan2((TE[1] - LE[1]), (TE[0] - LE[0])) ##Angle of rotation
    RotBlade = np.zeros(blade.shape)
    for j in range(blade[:,0].size):
        RotBlade[j,:] = mf.vectorRotP(blade[j,0], blade[j,1], TE, theta)[0:2]

    if (le < te) and (le != 0):
        LS_rot = RotBlade[le:te+1]
        US_rot = np.concatenate((RotBlade[te:-1], RotBlade[0:le+1]))
        US_rot = np.flip(US_rot,0)
    elif (le < te) and (le == 0):
        LS_rot = RotBlade[le:te+1]
        US_rot = RotBlade[te:-1]
        US_rot = np.flip(US_rot,0)
    elif (te < le) and (te != 0):
        US_rot = RotBlade[te:le+1]
        LS_rot = np.concatenate((RotBlade[le:-1], RotBlade[0:te+1]))
        US_rot = np.flip(US_rot,0)
    elif (te < le) and (te == 0):
        US_rot = RotBlade[te:le+1]
        LS_rot = RotBlade[le:-1]
        US_rot = np.flip(US_rot,0)

    mn = min(RotBlade[:,0])
    mx = max(RotBlade[:,0])
    if distro == 'cosine':
        query_points = (mx - mn)*(0.5*(1 - np.cos(np.linspace(0, np.pi, n)))) + mn
    elif distro == 'linear':
        query_points = np.linspace(mn, mx, n)
    
    US_spline = CubicSpline(US_rot[:,0], US_rot[:,1])
    US_newY = US_spline(query_points)
    LS_spline = CubicSpline(LS_rot[:,0], LS_rot[:,1])
    LS_newY = LS_spline(query_points)   

    camberRotated = np.zeros((n,2))
    camberRotated[:,0] = query_points
    camberRotated[:,1] = 0.5*(US_newY + LS_newY)

    approxCamber = np.zeros((n,2))
    for j in range(0,n):
        approxCamber[j,:] = mf.vectorRotP(camberRotated[j,0], camberRotated[j,1], TE, (-theta))[0:2]

    US = np.zeros((np.shape(US_rot)))
    for j in range(0,len(US)):
        US[j,:] = mf.vectorRotP(US_rot[j,0], US_rot[j,1], TE, (-theta))[0:2]

    LS = np.zeros((np.shape(LS_rot)))
    for j in range(0,len(LS)):
        LS[j,:] = mf.vectorRotP(LS_rot[j,0], LS_rot[j,1], TE, (-theta))[0:2]

    return US, LS, approxCamber

#%%
def getAccurateCamber(blade, approxCamber, tolIter, N, maxThickness, dis):
    import numpy as np
    import geometryFunctions as mf

    count = 0
    tol = 10
    n = len(approxCamber)
    oneThird = int(n/3)
    
    cam_prev = approxCamber.T.copy()
    cam_new = approxCamber.T.copy()
    
    a = np.flip(np.arange(1,oneThird*2))
    b = np.arange(oneThird,n)
    while tol >= tolIter:# or count < N:
        print('Accurate camber iteration: ' + str(count))
        if count > N:
            break
        leftCutoffRequired = 0
        rightCutoffRequired = 0
        for i in range(a.size):
            k = int(a[i])
            p1 = cam_new[:, k+2]
            p2 = cam_new[:, k+1]
            pm = cam_new[:,k]
    
            slope = -1/mf.Slope(p2[0], p2[1], p1[0], p1[1]) #perpendicular slope
            interC = pm[1] - slope*pm[0]
            x1, y1, x2, y2 = mf.pointPlusDistSlope(maxThickness, slope, pm)
    
            testLine = np.zeros((10,2))
            testLine[:,1] = np.linspace(y2, y1, 10)
            testLine[:,0] = (testLine[:,1] - interC) / slope
            numP, P = mf.TwoLinesIntersect(blade,testLine)
            if numP == 2:
                PM = ((P[0][0] + P[1][0]) / 2) , ((P[0][1] + P[1][1]) / 2)
            if (numP != 2) or (np.abs(mf.angleBetween(p2, p1, PM)) < 3): #to prevent the new point making too large an angle with the exitsting camber
                LE1 = cam_new[:,k+1]
                LE2 = cam_new[:,k+2]
                testLine = mf.extrapolateLinearLine(LE1, LE2, maxThickness)
                numP, P = mf.TwoLinesIntersect(blade,testLine)
                finalLE = P
                left_cutoff = k+1
                leftCutoffRequired = 1
                break
            cam_new[:,k] = PM
        
        
        for i in range(b.size):
            k = int(b[i])
            p1 = cam_new[:, k-2]
            p2 = cam_new[:, k-1]
            pm = cam_new[:,k]
    
            slope = -1/mf.Slope(p2[0], p2[1], p1[0], p1[1]) #perpendicular slope
            interC = pm[1] - slope*pm[0]
            x1, y1, x2, y2 = mf.pointPlusDistSlope(maxThickness, slope, pm)
    
            testLine = np.zeros((10,2))
            testLine[:,1] = np.linspace(y2, y1, 10)
            testLine[:,0] = (testLine[:,1] - interC) / slope
    
            numP, P = mf.TwoLinesIntersect(blade,testLine)
            if numP == 2:
                PM = ((P[0][0] + P[1][0]) / 2) , ((P[0][1] + P[1][1]) / 2)
            if (numP != 2) or (np.abs(mf.angleBetween(p2, p1, PM)) < 3): #to prevent the new point making too large an angle with the exitsting camber
                TE1 = cam_new[:,k-1]
                TE2 = cam_new[:,k-2]
                testLine = mf.extrapolateLinearLine(TE1, TE2, maxThickness)
                numP, P = mf.TwoLinesIntersect(blade,testLine)
                #plt.plot(blade[:,0], blade[:,1])
                #plt.plot(testLine[:,0], testLine[:,1])
                finalTE = P
                right_cutoff = k-1
                rightCutoffRequired = 1
                break
            cam_new[:,k] = PM
    
        if leftCutoffRequired == 0:
            left_cutoff = 1
            LE1 = cam_new[:,1]
            LE2 = cam_new[:,2]
            testLine = mf.extrapolateLinearLine(LE1, LE2, maxThickness)
            numP, P = mf.TwoLinesIntersect(blade,testLine)
            finalLE = P
        if rightCutoffRequired == 0:
            right_cutoff = n-2
            TE1 = cam_new[:,98]
            TE2 = cam_new[:,97]
            testLine = mf.extrapolateLinearLine(TE1, TE2, maxThickness)
            numP, P = mf.TwoLinesIntersect(blade,testLine)
            finalTE = P
    
        cam_new[0,0:left_cutoff+1] = np.linspace(finalLE[0], cam_new[0,left_cutoff],left_cutoff+1)
        cam_new[1,0:left_cutoff+1] = np.linspace(finalLE[1], cam_new[1,left_cutoff],left_cutoff+1)
        cam_new[0,right_cutoff:] = np.linspace(cam_new[0,right_cutoff], finalTE[0],n-right_cutoff)
        cam_new[1,right_cutoff:] = np.linspace(cam_new[1,right_cutoff], finalTE[1],n-right_cutoff)
        
        cam_diff = cam_prev - cam_new
        tol = (np.max(cam_diff[0]**2 + cam_diff[1]**2))/dis
        print('Accurate camber tol: ' + str(tol))
        cam_prev = cam_new.copy()
        
        count += 1  
    
    camber = cam_new.T
    camber[0] = approxCamber[0]
    camber[-1] = approxCamber[-1]

    return camber





























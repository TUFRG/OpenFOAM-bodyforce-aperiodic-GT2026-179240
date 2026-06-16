#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2024-09-18

@author: Justin Smart
"""

import numpy as np
import geometryFunctions as mf
from scipy.interpolate import CubicSpline

#%%
def LE_TE_Solver(PS, SS, LE, TE, camberLineEstimate, distro, n):
    blade = np.vstack((PS[0:len(PS)+1],SS[::-1][1:len(SS)+1]))
    
    # Rotate the SS and PS such that the LE and TE are on the x-axis
    thetaPrime = np.arctan2((TE[1] - LE[1]), (TE[0] - LE[0]))
    SS_rotated = mf.vectorRot(SS[:, 0], SS[:, 1], -thetaPrime)
    SS_rotated = np.column_stack(SS_rotated)
    PS_rotated = mf.vectorRot(PS[:, 0], PS[:, 1], -thetaPrime)
    PS_rotated = np.column_stack(PS_rotated)
    SS_sorted = SS_rotated[SS_rotated[:, 0].argsort()]
    PS_sorted = PS_rotated[PS_rotated[:, 0].argsort()]
    
    # Rotate the blade and camber such that the LE and TE are on the x-axis
    blade_rotated = mf.vectorRot(blade[:, 0], blade[:, 1], -thetaPrime)
    blade_rotated = np.column_stack(blade_rotated)
    camber_rotated = mf.vectorRot(camberLineEstimate[:, 0], camberLineEstimate[:, 1], -thetaPrime)
    camber_rotated = np.column_stack(camber_rotated)
    minCamberX = np.min(camber_rotated[:,0])
    maxCamberX = np.max(camber_rotated[:,0])
    try:
        tempSpline = CubicSpline(camber_rotated[:, 0], camber_rotated[:, 1])
    except:
        camber_rotated = camber_rotated[camber_rotated[:, 0].argsort()]
        tempSpline = CubicSpline(camber_rotated[:, 0], camber_rotated[:, 1])
    camber_rotated = np.zeros((500,2))
    camber_rotated[:,0] = np.linspace(minCamberX, maxCamberX, 500)
    camber_rotated[:,1] = tempSpline(camber_rotated[:,0])
    camber_rotated = camber_rotated[camber_rotated[:, 0].argsort()]
    
    LE_rotated = mf.vectorRot(LE[0], LE[1], -thetaPrime)
    TE_rotated = mf.vectorRot(TE[0], TE[1], -thetaPrime)

    
    # Assigning the sorted arrays back to the variables
    SS = SS_sorted
    PS = PS_sorted
    camber_line = camber_rotated
    blade = blade_rotated
    LE = LE_rotated
    TE = TE_rotated

    # Make spline for camber and take 2nd derr.
    cs_camber_line = CubicSpline(camber_line[:, 0], camber_line[:, 1])
    d2ydx2_camber = cs_camber_line.derivative(nu=2)(camber_line[:, 0])
    numCamberPoints = len(camber_line)
    
    # Find all boundary points where the 2nd derivative changes sign
    boundary_points = np.where(np.diff(np.sign(d2ydx2_camber)))[0]
    if any(np.in1d(boundary_points, np.arange(50,450))) == True:
        boundary_points = (np.where( abs(np.diff(d2ydx2_camber)) > 0.01)[0])

    # Seperate boundary points into left and right sides
    left_boundary_points = boundary_points[boundary_points[:] < (numCamberPoints/2)]
    right_boundary_points = boundary_points[boundary_points[:] > (numCamberPoints/2)]

    # Remove boundary points except in the first and last 4%
    left_boundary_points = left_boundary_points[left_boundary_points[:] < (numCamberPoints/25)]
    right_boundary_points = right_boundary_points[right_boundary_points[:] > (24*numCamberPoints/25)] 
    
    # Find x and y values of the left and right boundary points
    left_boundary_points_x = camber_line[left_boundary_points, 0]
    left_boundary_points_y = camber_line[left_boundary_points, 1]
    right_boundary_points_x = camber_line[right_boundary_points, 0]
    right_boundary_points_y = camber_line[right_boundary_points, 1]
    
    # stack
    left_boundary_points = np.column_stack((left_boundary_points_x, left_boundary_points_y))
    right_boundary_points = np.column_stack((right_boundary_points_x, right_boundary_points_y))
    
    if left_boundary_points.size > 0 and right_boundary_points.size > 0:
        camber_line = camber_line[
            (camber_line[:, 0] > left_boundary_points[:, 0].max()) & 
            (camber_line[:, 0] < right_boundary_points[:, 0].min())    
        ]    
    elif left_boundary_points.size > 0 and right_boundary_points.size == 0:
        camber_line = camber_line[
            (camber_line[:, 0] > left_boundary_points[:, 0].max())
        ]        
    elif left_boundary_points.size == 0 and right_boundary_points.size > 0:
        camber_line = camber_line[
            (camber_line[:, 0] < right_boundary_points[:, 0].min())
        ]  
    
    # Split SS, PS, and camber into front half and back half
    midway = (np.max(blade[:,0]) + np.min(blade[:,0]))/2
    SSleft = SS[SS[:,0] < midway]
    SSright = SS[SS[:,0] > midway]
    PSleft = PS[PS[:,0] < midway]
    PSright = PS[PS[:,0] > midway]
    camberLeft = camber_line[camber_line[:,0] < midway]
    camberRight = camber_line[camber_line[:,0] > midway]
    
    def circle_from_points(p1, p2, p3):
        num_points=10000
        temp = p2[0]**2 + p2[1]**2
        bc = (p1[0]**2 + p1[1]**2 - temp) / 2
        cd = (temp - p3[0]**2 - p3[1]**2) / 2
        cx = (bc * (p2[1] - p3[1]) - cd * (p1[1] - p2[1])) / det
        cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det
        radius = np.sqrt((cx - p1[0])**2 + (cy - p1[1])**2)
        
        angles = np.linspace(0, 2 * np.pi, num_points)
        circle_points = np.zeros((num_points, 2))
        circle_points[:, 0] = cx + radius * np.cos(angles)
        circle_points[:, 1] = cy + radius * np.sin(angles)
        return circle_points

    # BEGIN ROTATION - LEFT
    # Rotate SS and PS so that the leftmost 2 points of the camber are lying on the x axis
    leftmostCamberPoint1 = camber_line[0]
    leftmostCamberPoint2 = camber_line[1]
    
    theta = np.arctan2((leftmostCamberPoint2[1] - leftmostCamberPoint1[1]), (leftmostCamberPoint2[0] - leftmostCamberPoint1[0]))
    rotatedSSLeft = mf.vectorRot(SSleft[:, 0], SSleft[:, 1], -theta)
    rotatedSSLeft = np.column_stack(rotatedSSLeft)
    rotatedPSLeft = mf.vectorRot(PSleft[:, 0], PSleft[:, 1], -theta)
    rotatedPSLeft = np.column_stack(rotatedPSLeft)
    rotatedCamberLeft = mf.vectorRot(camberLeft[:, 0], camberLeft[:, 1], -theta)
    rotatedCamberLeft = np.column_stack(rotatedCamberLeft)
    rotatedLE = mf.vectorRot(LE[0], LE[1], -theta)
    
    # Split out leftmost points that exceed beyond the cut-off camber
    leftmostSS = rotatedSSLeft[rotatedSSLeft[:,0] < np.min(rotatedCamberLeft[:,0])]
    leftmostPS = rotatedPSLeft[rotatedPSLeft[:,0] < np.min(rotatedCamberLeft[:,0])]
    left = np.concatenate((leftmostSS ,leftmostPS))
    left = left[left[:, 1].argsort()]
    
    # Remove the leftmost points from the SS and PS dataset
    rotatedSSLeft = rotatedSSLeft[rotatedSSLeft[:,0] > np.min(rotatedCamberLeft[:,0])]
    rotatedPSLeft = rotatedPSLeft[rotatedPSLeft[:,0] > np.min(rotatedCamberLeft[:,0])]
   
    # Keep increasing point interval until they are no longer colinear (by checking determinant)
    i=1
    while True:
        p1 = rotatedCamberLeft[0]
        p2 = rotatedCamberLeft[i]
        p3 = rotatedCamberLeft[2*i]
        det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])
        if abs(det) < 1.0e-12:
            i += 1
            continue
        else:
            break
    
    circle_points = circle_from_points(p1, p2, p3)
    circle_points = np.array(circle_points)
    try:
        numP, intersection_point_LE = mf.TwoLinesIntersect(circle_points, left)
        tempLE = np.column_stack((intersection_point_LE[0], intersection_point_LE[1]))
        # Get new leftmost SS and PS from splitting at new LE
        leftmostSS = left[left[:,1] > intersection_point_LE[1]]
        leftmostPS = left[left[:,1] < intersection_point_LE[1]]
        
        #Combine to get new SS and PS left
        rotatedSSLeft = np.concatenate((tempLE, leftmostSS, rotatedSSLeft))
        rotatedSSLeft = rotatedSSLeft[rotatedSSLeft[:, 0].argsort()]
        rotatedPSLeft = np.concatenate((tempLE, leftmostPS, rotatedPSLeft))
        rotatedPSLeft = rotatedPSLeft[rotatedPSLeft[:, 0].argsort()]
        
        # Get new camber left
        arc_leftmost = circle_points[(circle_points[:, 0] < p1[0]) &
                                     (circle_points[:, 0] > intersection_point_LE[0]) & 
                                     (circle_points[:, 1] > left[:, 1].min()) & 
                                     (circle_points[:, 1] < left[:, 1].max())]
        rotatedCamberLeft = np.concatenate((tempLE, arc_leftmost, rotatedCamberLeft))
        rotatedCamberLeft = rotatedCamberLeft[rotatedCamberLeft[:, 0].argsort()]
    except:
        intersection_point_LE = rotatedLE
    
    # Create spline for camber to get back to original point distribution
    cs_camber_line = CubicSpline(rotatedCamberLeft[:, 0], rotatedCamberLeft[:, 1])
    mn = min(rotatedCamberLeft[:,0])
    mx = max(rotatedCamberLeft[:,0])
    numPoints = int(n/2) + 1
    if distro == 'cosine':
        query_points = (mx - mn)*(1 - np.cos(np.linspace(0, (np.pi)/2, numPoints))) + mn
    elif distro == 'linear':
        query_points = np.linspace(mn, mx, numPoints)
    y_camber_line = cs_camber_line(query_points)
    rotatedCamberLeft = np.column_stack((query_points, y_camber_line))
    
    #Rotate back - new LE, SS, PS, camber (left)
    newLE = mf.vectorRot(intersection_point_LE[0], intersection_point_LE[1], theta)
    
    camberLeft = mf.vectorRot(rotatedCamberLeft[:, 0], rotatedCamberLeft[:, 1], theta)
    camberLeft = np.column_stack(camberLeft)
    SSLeft = mf.vectorRot(rotatedSSLeft[:, 0], rotatedSSLeft[:, 1], theta)
    SSLeft = np.column_stack(SSLeft)
    PSLeft = mf.vectorRot(rotatedPSLeft[:, 0], rotatedPSLeft[:, 1], theta)
    PSLeft = np.column_stack(PSLeft)
    # END ROTATION - LEFT
    
    # BEGIN ROTATION - RIGHT
    # Rotate SS and PS so that the rightmost 2 points of the camber are lying on the x axis
    rightmostCamberPoint1 = camber_line[-2]
    rightmostCamberPoint2 = camber_line[-1]
    
    theta = np.arctan2((rightmostCamberPoint2[1] - rightmostCamberPoint1[1]), (rightmostCamberPoint2[0] - rightmostCamberPoint1[0]))
    rotatedSSRight = mf.vectorRot(SSright[:, 0], SSright[:, 1], -theta)
    rotatedSSRight = np.column_stack(rotatedSSRight)
    rotatedPSRight = mf.vectorRot(PSright[:, 0], PSright[:, 1], -theta)
    rotatedPSRight = np.column_stack(rotatedPSRight)
    rotatedCamberRight = mf.vectorRot(camberRight[:, 0], camberRight[:, 1], -theta)
    rotatedCamberRight = np.column_stack(rotatedCamberRight)
    rotatedTE = mf.vectorRot(TE[0], TE[1], -theta)
    
    # Split out rightmost points that exceed beyond the cut-off camber
    rightmostSS = rotatedSSRight[rotatedSSRight[:,0] > np.max(rotatedCamberRight[:,0])]
    rightmostPS = rotatedPSRight[rotatedPSRight[:,0] > np.max(rotatedCamberRight[:,0])]
    right = np.concatenate((rightmostSS ,rightmostPS))
    right = right[right[:, 1].argsort()]
    
    # Remove the rightmost points from the SS and PS dataset
    rotatedSSRight = rotatedSSRight[rotatedSSRight[:,0] < np.max(rotatedCamberRight[:,0])]
    rotatedPSRight = rotatedPSRight[rotatedPSRight[:,0] < np.max(rotatedCamberRight[:,0])]
    
    # Keep increasing point interval until they are no longer colinear (by checking determinant)
    i=2
    while True:
        p3 = rotatedCamberRight[-1]
        p2 = rotatedCamberRight[-i]
        p1 = rotatedCamberRight[-2*i+1]
        det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])
        if abs(det) < 1.0e-12:
            i += 1
            continue
        else:
            break
    
    circle_points = circle_from_points(p1, p2, p3)
    circle_points = np.array(circle_points)
    try:
        numP, intersection_point_TE = mf.TwoLinesIntersect(circle_points, right)
        tempTE = np.column_stack((intersection_point_TE[0], intersection_point_TE[1]))
        # Get new rightmost SS and PS from splitting at new TE
        rightmostSS = right[right[:,1] > intersection_point_TE[1]]
        rightmostPS = right[right[:,1] < intersection_point_TE[1]]
        #Combine to get new SS and PS right
        rotatedSSRight = np.concatenate((rightmostSS, rotatedSSRight, tempTE))
        rotatedSSRight = rotatedSSRight[rotatedSSRight[:, 0].argsort()]
        rotatedPSRight = np.concatenate((rightmostPS, rotatedPSRight, tempTE))
        rotatedPSRight = rotatedPSRight[rotatedPSRight[:, 0].argsort()]
        # Get new camber right
        arc_rightmost = circle_points[(circle_points[:, 0] > p3[0]) &
                                     (circle_points[:, 0] < intersection_point_TE[0]) & 
                                     (circle_points[:, 1] > right[:, 1].min()) & 
                                     (circle_points[:, 1] < right[:, 1].max())]
        rotatedCamberRight = np.concatenate((arc_rightmost, rotatedCamberRight, tempTE))
        rotatedCamberRight = rotatedCamberRight[rotatedCamberRight[:, 0].argsort()]
    except:
        intersection_point_TE = rotatedTE
    
    # Create spline for camber to get back to original point distribution
    cs_camber_line = CubicSpline(rotatedCamberRight[:, 0], rotatedCamberRight[:, 1])
    mn = min(rotatedCamberRight[:,0])
    mx = max(rotatedCamberRight[:,0])
    numPoints = n - int(n/2) + 1
    if distro == 'cosine':
        query_points = (mx - mn)*(-1*(np.cos(np.linspace((np.pi)/2, np.pi, numPoints)))) + mn
    elif distro == 'linear':
        query_points = np.linspace(mn, mx, numPoints)
    y_camber_line = cs_camber_line(query_points)
    rotatedCamberRight = np.column_stack((query_points, y_camber_line))
    
    #Rotate back - new LE, SS, PS, camber (right)
    newTE = mf.vectorRot(intersection_point_TE[0], intersection_point_TE[1], theta)
    camberRight = mf.vectorRot(rotatedCamberRight[:, 0], rotatedCamberRight[:, 1], theta)
    camberRight = np.column_stack(camberRight)
    SSRight = mf.vectorRot(rotatedSSRight[:, 0], rotatedSSRight[:, 1], theta)
    SSRight = np.column_stack(SSRight)
    PSRight = mf.vectorRot(rotatedPSRight[:, 0], rotatedPSRight[:, 1], theta)
    PSRight = np.column_stack(PSRight)
    # END ROTATION - RIGHT
    
    # Connect left and right sides
    SS_new = np.concatenate((SSLeft, SSRight))
    PS_new = np.concatenate((PSLeft, PSRight))
    camber_line = np.concatenate((camberLeft, camberRight[1:]))
    
    # Final rotation back to original blade coordinates
    SS_new = mf.vectorRot(SS_new[:, 0], SS_new[:, 1], thetaPrime)
    SS_new = np.column_stack(SS_new)
    PS_new = mf.vectorRot(PS_new[:, 0], PS_new[:, 1], thetaPrime)
    PS_new = np.column_stack(PS_new)
    camber_line = mf.vectorRot(camber_line[:, 0], camber_line[:, 1], thetaPrime)
    camber_line = np.column_stack(camber_line)
    newLE = mf.vectorRot(newLE[0], newLE[1], thetaPrime)
    newTE = mf.vectorRot(newTE[0], newTE[1], thetaPrime)
    
    return PS_new, SS_new, newLE, newTE, camber_line

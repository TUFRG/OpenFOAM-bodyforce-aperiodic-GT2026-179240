#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2024-09-18

@author: Amelia George
"""

import numpy as np 

def Slope (x1,y1,x2,y2):
    grad = (y2 - y1)/ (x2 - x1)
    return grad

def vectorRot(x,y,theta):
    import numpy as np
    xRot = x*np.cos(theta) - y*np.sin(theta)
    yRot = x*np.sin(theta) + y*np.cos(theta)
    return xRot, yRot

def vectorRotP(x,y,P,theta):
    import numpy as np
    xRot = ((x-P[0])*np.cos(theta) - (P[1]-y)*np.sin(theta)) + P[0]
    yRot = P[1] -((x-P[0])*np.sin(theta) + (P[1]-y)*np.cos(theta))
    return xRot, yRot, 0

def rotate_points(points, angle):
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    return np.dot(points, rotation_matrix.T)

def extrapolateLinearLine(point1, point2, distance): #designed to extrapolate from point1, away from point2
    slope = Slope(point1[0],point1[1],point2[0],point2[1])
    interC = point1[1] - slope*point1[0]
    x1, y1, x2, y2 = pointPlusDistSlope(distance, slope, point1)
    d1 = dist2D(x1,y1,point2[0],point2[1])
    d2 = dist2D(x2,y2,point2[0],point2[1])
    testLine = np.zeros((10,2))
    if d1 > d2:
        testLine[:,0] = np.linspace(x1, point1[0], 10)
    elif d2 > d1:
        testLine[:,0] = np.linspace(x2, point1[0], 10)
    testLine[:,1] = slope*(testLine[:,0]) + interC
    return testLine

def TwoLinesIntersect(Line1, Line2): #Get the intersection point if the lines already intersects
    import numpy as np
    from shapely.geometry import LineString
    L1 = LineString(np.column_stack((Line1[:,0],Line1[:,1])))
    L2 = LineString(np.column_stack((Line2[:,0],Line2[:,1])))
    InterPt = L1.intersection(L2)
    if InterPt.is_empty:
        numP = 0
        P = 0
    elif InterPt.geom_type == 'Point':
        numP = 1
        P = InterPt.x, InterPt.y
    elif InterPt.geom_type == 'MultiPoint':
        numP = 0
        P = []
        for point in InterPt.geoms:
            numP += 1
            P.append((point.x, point.y))
    return numP, P

def pointPlusDistSlope(d, m, point):
    import numpy as np
    x0 = point[0]
    y0 = point[1]
    dx = d / (np.sqrt(d + m**2))
    dy = m * dx
    x1 = x0 + dx
    y1 = y0 + dy
    x2 = x0 - dx
    y2 = y0 - dy
    return x1, y1, x2, y2

def dist2D(x1,y1,x2,y2):
    import numpy as np
    distance = np.sqrt(np.square(x2 - x1) + np.square(y2 - y1))
    return distance

def longest_distance_on_curve(x, y):
    import numpy as np
    """
    Find the longest distance between two points on a closed curve defined by x and y coordinates.
    Return the maximum distance and the indices of the two points.
    """
    # Combine x and y coordinates into a single array of points
    points = np.column_stack([x, y])
    
    # Compute pairwise distances between points
    distances = pairwise_distances(points)
    
    # Find the indices of the maximum distance
    idx_i, idx_j = np.unravel_index(np.argmax(distances), distances.shape)
    
    # Find the maximum distance
    max_distance = distances[idx_i, idx_j]
    
    return max_distance, idx_i, idx_j

def pairwise_distances(points):
    import numpy as np
    """
    Calculate pairwise distances between points.
    """
    num_points = len(points)
    distances = np.zeros((num_points, num_points))
    for i in range(num_points):
        for j in range(i+1, num_points):
            distances[i, j] = distances[j, i] = np.linalg.norm(points[i] - points[j])
    return distances

def angleBetween(pM, pA, pB):
    import numpy as np
    vec1 = (( (pA[0] - pM[0]), (pA[1] - pM[1]) ))
    vec2 = (( (pB[0] - pM[0]), (pB[1] - pM[1]) ))
    dot = ((vec1[0]*vec2[0]) + (vec1[1]*vec2[1]))
    det = ((vec1[0]*vec2[1]) - (vec1[1]*vec2[0]))
    return np.arctan2(det,dot)

def cart2pol(x,y,z):
    import numpy as np
    theta = np.arctan(y/x)
    rho = np.sqrt(np.square(x)+np.square(y))
    Z = z
    return(theta, rho, Z)

#%%

'''

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False
def draw_curve(p1, p2):
    import numpy as np
    a = (p2[1] - p1[1]) / (np.cosh(p2[0]) - np.cosh(p1[0]))
    b = p1[1] - a * np.cosh(p1[0])
    x = np.linspace(p1[0], p2[0], 100)
    y = a * np.cosh(x) + b
    return x, y





def dist3D(x1,x2,y1,y2,z1,z2):
    import numpy as np
    distance = np.sqrt(np.square(x2 - x1) + np.square(y2 - y1) + np.square(z2- z1))
    return distance



def AngRadCenPts(x1,y1,x2,y2,slope,m1,m2):
    import numpy as np
    slope_perp = -1/slope
    Cinterp1 = y2 - slope_perp*x2
    P1 = [x1, x1+m1]
    P2 = [x2+m2, (slope_perp*(x2+m2)+Cinterp1)]
    Line1 = line(P1,[x1,y1])
    Line2 = line(P2,[x2,y2])
    InterPt = intersection(Line1, Line2)
    radius1 = np.sqrt((InterPt[0] - x1)**2 + (InterPt[1] - y1)**2)
    radius2 = np.sqrt((InterPt[0] - x2)**2 + (InterPt[1] - y2)**2)
    return (radius1, radius2, InterPt)
    
def ExtentionIntersection(x1,y1,Line2,theta,m1):
    import numpy as np
    from shapely.geometry import LineString
    slope = np.tan(theta*2*np.pi/180)
    slope_perp = -1/slope
    Cinterp = y1 - slope_perp*x1
    NewPt = [[x1-m1,slope_perp*(x1-m1)+Cinterp]]
    originalPt = [[x1,y1]]
    Line1 = np.concatenate((originalPt, NewPt))
    L1 = LineString(np.column_stack((Line1[:,0],Line1[:,1])))
    L2 = LineString(np.column_stack((Line2[:,0],Line2[:,1])))
    InterPt = L1.intersection(L2)
    [x,y] =InterPt.xy
    Center = [x[0],y[0],0]
    return (Center, Line1)

def simulEqn(x1,y1,x2,y2):
    import numpy as np
    A = np.array([[((2*x2)-(2*x1)), ((2*y2)-(2*y1))],[(2*x1), (2*y1)]])
    G = np.linalg.det(A)
    B = np.array([(x2**2+y2**2-x1**2-y1**2),(x1**2+y1**2)])
    result = np.linalg.inv(A).dot(B)
    return result

def RadAndAng(x1,y1,x2,y2):
    import numpy as np
    center = simulEqn(x1, y1, x2, y2)
    radius = np.sqrt((center[0] - x1)**2 + (center[1] - y1)**2)
    base = dist2D(x1, y1, x2, y2)
    cosAngle = ((radius**2 + radius**2) - base**2)/(2*radius*radius)
    Angle = np.arccos(cosAngle)
    AngleDeg = np.rad2deg(Angle)
    arcLen = Angle * radius
    return arcLen
    
def scale(rmin,rmax,tmin,tmax,m):
    NewM = (((rmin - m)/(rmax - rmin)) * (tmax - tmin)) + tmax
    return NewM



def MidPts (Line):
    x1 = Line[0][0]
    x2 = Line[len(Line)-1][0]
    y1 = Line[0][1]
    y2 = Line[len(Line)-1][1]
    x = (x1+x2)/2
    y = (y1+y2)/2
    MidPoint = [x,y,0]
    return MidPoint
    

def vectorRot3D (x,y,z,theta):
    import numpy as np
    yRot = y*np.cos(theta) - z*np.sin(theta)
    zRot = y*np.sin(theta) + z*np.cos(theta)
    return x, yRot, zRot
def TwoLinesIntersectExtsn (Line1, Line2, m): #Get the intersection point if the lines DO NOT intersects
    import numpy as np
    from shapely.geometry import LineString
    slope1 = (Line1[1][1] - Line1[0][1])/(Line1[1][0] - Line1[0][0])
    slope2 = (Line2[1][1] - Line2[0][1])/(Line2[1][0] - Line2[0][0])
    Intercept1 = Line1[1][1] - (slope1 * Line1[1][0])
    Intercept2 = Line2[1][1] - (slope2 * Line2[1][0])
    NLine1 = np.vstack(((Line1[1][0],Line1[1][1]),(Line1[1][0]+m, (slope1*(Line1[1][0]+m))+Intercept1)))
    NLine2 = np.vstack(((Line2[1][0],Line2[1][1]),(Line2[1][0]+m, (slope2*(Line2[1][0]+m))+Intercept2)))
    L1 = LineString(np.column_stack((NLine1[:,0],NLine1[:,1])))
    L2 = LineString(np.column_stack((NLine2[:,0],NLine2[:,1])))
    InterPt = L1.intersection(L2)
    [x,y] =InterPt.xy
    Center = [x[0],y[0],0]
    return Center    

def redistribute_point(x, y, num_points):
    import numpy as np
    from scipy.interpolate import interp1d
    from scipy.integrate import cumtrapz
    assert len(x) == len(y), "x and y arrays must have the same length"
    # Calculate cumulative distance along the curve
    distances = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    cumulative_distances = np.concatenate(([0], cumtrapz(distances)))
    cumulative_distances = np.append(cumulative_distances, cumulative_distances[-1])
    # Linearly interpolate the cumulative distances to get a function f(s) = t,
    # where s is the arc length parameter and t is the parameter along the curve
    interp_func = interp1d(cumulative_distances, np.arange(len(x)), kind='linear')
    
    # Equally distribute points along the curve
    new_distances = np.linspace(0, cumulative_distances[-1], num_points)
    new_indices = interp_func(new_distances).astype(int)
    
    # Extract redistributed points
    redistributed_x = x[new_indices]
    redistributed_y = y[new_indices]
    
    return redistributed_x, redistributed_y





def pol2cart(r, theta, z):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, z



def linearToX(p1,p2,x):
    slope = Slope(p2[0], p2[1], p1[0], p1[1])
    interC = p1[1] - slope*p1[0]
    y = slope*x + interC
    return y










'''















#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 19:15:33 2025

@author: adekola
"""
import argparse
import numpy as np
from scipy.interpolate import CubicSpline, interp1d
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve, brentq
#%%
def vectorRotZ3D (x,y,z,theta):
    import numpy as np
    xRot = x*np.cos(theta) - y*np.sin(theta)
    yRot = x*np.sin(theta) + y*np.cos(theta)
    return xRot, yRot, z
    
def TwoLinesIntersect (Line1, Line2): #Get the intersection point if the lines already intersects
    from shapely.geometry import LineString, Point, MultiPoint
    L1 = LineString(np.column_stack((Line1[:,0],Line1[:,1])))
    L2 = LineString(np.column_stack((Line2[:,0],Line2[:,1])))
    InterPt = L1.intersection(L2)
    if InterPt.is_empty:
        return 0
    else:
        if InterPt.geom_type == 'MultiPoint':
            return [[point.x, point.y] for point in InterPt.geoms]
        else:
            return [InterPt.x, InterPt.y]
        
def scale(rmin,rmax,tmin,tmax,m):
    NewM = (((rmin - m)/(rmax - rmin)) * (tmax - tmin)) + tmax
    return NewM

def Slope (x1,y1,x2,y2):
    grad = (y2 - y1)/ (x2 - x1)
    return grad

def cone_curve(r, T, Tr, rMid, interC, slope, Nb, num_points=100):
    """
    Parameters
    ----------
    r : This is the midPoint on the nose. Basically, the midpoint on the horizontal axis, the corresponding radial value is found, that is the value of r
    T : This is the radial value of the top of the diamond of the butterfly mesh
    Tr : So this the radial value of the mid level of the diamond. Look at it this way, the diamond has four points, the up/dwn left/right, Tr is radial value left/right. T is up, down is 0
    rMid : Because I am using equatioon of a cone. rMid is same as r assuming the curve formed a perfect cone. But it doesn't so rMid is the mid point assuming a straight line
    interC : intercept of the equation of a straight line
    slope : slope of the slant edge of the cone
    Nb : Number of blades so I can determine extent of rotation.
    num_points : TYPE, optional. he default is 100.
    Returns
    -------
    x , y , z

    """
    theta = np.linspace(0,  np.deg2rad(360/(0.5*Nb)), num_points)
    t = np.linspace(T-r, 0, num_points)
    th = np.linspace(Tr-rMid, 0, num_points)
    x = (r+t) * np.cos(theta)
    y = r * np.sin(theta)
    z = (rMid+th - interC)/slope
    return x, y, z

def plot_cone_curve(x, y, z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label='Curve on Cone Surface')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Curve on Cone Surface')
    ax.legend()
    plt.axis('equal')
    plt.show()
       
def scale_curve_endpoints(x, y, x1, y1, x2, y2):
    # Translation vector
    dx = x1 - x[0]
    dy = y1 - y[0]
    # Translate to start at (x1, y1)
    x_trans = x + dx
    y_trans = y + dy
    # Scale to reach (x2, y2)
    x_range_original = x[-1] - x[0]
    y_range_original = y[-1] - y[0]
    x_range_target = x2 - x1
    y_range_target = y2 - y1
    if abs(x_range_original) > 1e-10:
        x_scale = x_range_target / x_range_original
        x_scaled = x1 + (x - x[0]) * x_scale
    else:
        x_scaled = x_trans
    if abs(y_range_original) > 1e-10:
        y_scale = y_range_target / y_range_original
        y_scaled = y1 + (y - y[0]) * y_scale
    else:
        y_scaled = y_trans
    return x_scaled, y_scaled

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

def extndCurve(curve, length):
    #This fuction takes a curve an extend it based on the slope of the last two points at both ends
    #It expects a 2D curve z,r
    #Also, this expects the curve to go increase in r
    hubSlope = Slope(curve[0,1], curve[0,0], curve[1,1], curve[1,0])
    hubIntX = curve[0,0] - hubSlope*curve[0,1]
    hubExtn = [hubSlope*(curve[0,1]-length)+hubIntX, curve[0,1]-length]
    casSlope = Slope(curve[-1,1], curve[-1,0], curve[-2,1], curve[-2,0])
    casIntX = curve[-1,0] - casSlope*curve[-1,1]
    casExtn = [casSlope*(curve[-1,1]+length)+casIntX, curve[-1,1]+length]    
    newCurve = np.vstack((hubExtn, curve, casExtn))
    return newCurve

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

def insertPoints_batch(curve_points, new_points, closed_loop=False):
    curve = np.array(curve_points)
    new_pts = np.array(new_points)
    # Ensure 2D array
    if new_pts.ndim == 1:
        new_pts = new_pts.reshape(1, -1)
    # Find insertion indices for all points
    insertion_data = []
    n_segments = len(curve) if closed_loop else len(curve) - 1
    for new_pt in new_pts:
        min_dist = np.inf
        best_index = 0
        for i in range(n_segments):
            p1 = curve[i]
            p2 = curve[(i + 1) % len(curve)]
            closest_pt, dist = point_to_segment_distance(new_pt, p1, p2)
            if dist < min_dist:
                min_dist = dist
                best_index = i + 1
        insertion_data.append((best_index, new_pt))
    # Sort by index (descending) so we insert from back to front
    # This prevents index shifting from affecting later insertions
    insertion_data.sort(key=lambda x: x[0], reverse=True)
    # Insert all points
    for idx, pt in insertion_data:
        curve = np.insert(curve, idx, pt, axis=0)
    return curve
def point_to_segment_distance(point, seg_start, seg_end):
    point = np.array(point)
    seg_start = np.array(seg_start)
    seg_end = np.array(seg_end)
    # Vector from seg_start to seg_end
    segment_vec = seg_end - seg_start
    segment_length_sq = np.dot(segment_vec, segment_vec)
    if segment_length_sq == 0:
        # Degenerate segment (point)
        return seg_start, np.linalg.norm(point - seg_start)
    # Project point onto the line containing the segment
    # t is the parameter: closest_point = seg_start + t * segment_vec
    t = np.dot(point - seg_start, segment_vec) / segment_length_sq
    # Clamp t to [0, 1] to stay on the segment
    t = np.clip(t, 0, 1)
    # Find the closest point on the segment
    closest_point = seg_start + t * segment_vec
    dist = np.linalg.norm(point - closest_point)
    return closest_point, dist

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
#%%
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
    "res", 
    type=int,
    help="stl resolution"
)
parser.add_argument(
    "domainExtn", 
    type=float,
    help="domain extension"
)
parser.add_argument(
    "velocity", 
    type=float,
    help="farstream velocity"
)
parser.add_argument(
    "target_yPlus",
    type=float,
    help="target_yplus"
)
parser.add_argument(
    "NrCells", 
    type=int,
    help="number of cells radial direction"
)
parser.add_argument(
    "NzCells", 
    type=int,
    help="number of cells in axial direction"
)
parser.add_argument(
    "NqCells", 
    type=int,
    help="number of cells in circumferential direction"
)
parser.add_argument(
    "hubZ",
    type=float,
    help="Axial location where rotor hub stops rotating"
)
parser.add_argument(
    "spinnerIdx", 
    type=int,
    help="index on hub where spinner nose ends"
)
parser.add_argument(
    "scales",
    type=float,
    help="convert to meters or SI units"
)
# Parse arguments
args = parser.parse_args()

# Use parsed arguments
Nbr = args.Nbr  # Total number of blades in the annulus for rotor
Nbs = args.Nbs  # Total number of blades in the annulus for stator
rSections = args.rSections  # Number of rotor blade profiles
sSections = args.sSections  # Number of stator blade profiles
res = args.res  
domainExtn = args.domainExtn  
velocity = args.velocity 
target_yPlus = args.target_yPlus
Nr = args.NrCells
Nz = args.NzCells 
Nq = args.NqCells 
hubZ = args.hubZ
spinnerIdx = args.spinnerIdx
scales = args.scales

filePath = '../rawData/'

Np = 101
Ncp = 101
length = 20
Nb = 4 #Number of blades
percent = 0.313 #Set the percentage upstream and downstream of the domain where the flow becomes inviscid. This percentage is based on the domainExtension

pitchAngle = 2*np.pi/Nb
theta = np.linspace(-0.5*pitchAngle, 0.5*pitchAngle,res)
rotorCamberData = np.zeros((rSections, Np, 3))
statorCamberData = np.zeros((sSections, Np, 3))
hubData = np.loadtxt(filePath + '/gasPath/hub.txt', delimiter=',', skiprows=0)[:,[2,1,0]]
hubData = np.delete(hubData, 1, 1)
casData = np.loadtxt(filePath + '/gasPath/casing.txt', delimiter=',', skiprows=1)[:,[2,1,0]]
casData = np.delete(casData,1, 1)
rLEData = np.zeros((rSections,2)) #z,r
rTEData = np.zeros((rSections,2))
sLEData = np.zeros((sSections,2))
sTEData = np.zeros((sSections,2))
for a in range(rSections):
        rotorCamberData[a,:,:] = np.loadtxt(filePath + '/camberData/rotor/camber{}.txt'.format(a), delimiter=',')
        rLEData[a,:] = np.column_stack((rotorCamberData[a,:,:][0,2], np.sqrt(rotorCamberData[a,:,:][0,1]**2+rotorCamberData[a,:,:][0,0]**2)))
        rTEData[a,:] = np.column_stack((rotorCamberData[a,:,:][-1,2], np.sqrt(rotorCamberData[a,:,:][-1,1]**2+rotorCamberData[a,:,:][-1,0]**2)))

for b in range(sSections):
        statorCamberData[b,:,:] = np.loadtxt(filePath + '/camberData/stator/camber{}.txt'.format(b), delimiter=',')       
        sLEData[b,:] = np.column_stack((statorCamberData[b,:,:][0,2], np.sqrt(statorCamberData[b,:,:][0,1]**2+statorCamberData[b,:,:][0,0]**2)))
        sTEData[b,:] = np.column_stack((statorCamberData[b,:,:][-1,2], np.sqrt(statorCamberData[b,:,:][-1,1]**2+statorCamberData[b,:,:][-1,0]**2)))

 
hubData *= scales
casData *= scales
rLEData *= scales
rTEData *= scales
sLEData *= scales
sTEData *= scales
hubZ *= scales
#%% Ensure that the LE/TE of rotor and stator lies within hub and casing 
#first Extend the curves
rLEDataExtn = extndCurve(rLEData, length)
rTEDataExtn = extndCurve(rTEData, length)
sLEDataExtn = extndCurve(sLEData, length)
sTEDataExtn = extndCurve(sTEData, length)
#find the intersection on the hub and casing 
rLEInterXHub = TwoLinesIntersect(hubData, rLEDataExtn)
rTEInterXHub = TwoLinesIntersect(hubData, rTEDataExtn)
sLEInterXHub = TwoLinesIntersect(hubData, sLEDataExtn)
sTEInterXHub = TwoLinesIntersect(hubData, sTEDataExtn)

rLEInterXCas = TwoLinesIntersect(casData, rLEDataExtn)
rTEInterXCas = TwoLinesIntersect(casData, rTEDataExtn)
sLEInterXCas = TwoLinesIntersect(casData, sLEDataExtn)
sTEInterXCas = TwoLinesIntersect(casData, sTEDataExtn)
#newly defined LE/TE curves
rLEHubIdx = np.argmin(abs(rLEData[:,1]-rLEInterXHub[1]))
rLECasIdx = np.argmin(abs(rLEData[:,1]-rLEInterXCas[1]))
newRLEData = np.vstack((rLEInterXHub, rLEData[rLEHubIdx+1:rLECasIdx], rLEInterXCas))
newRLEData = densifyCurve(newRLEData, res-2, 'uniform') # I am subtracting 5 because I will add the intersection to the offset hub and casing 

rTEHubIdx = np.argmin(abs(rTEData[:,1]-rTEInterXHub[1]))
rTECasIdx = np.argmin(abs(rTEData[:,1]-rTEInterXCas[1]))
newRTEData = np.vstack((rTEInterXHub, rTEData[rTEHubIdx+1:rTECasIdx], rTEInterXCas))
newRTEData = densifyCurve(newRTEData, res-2, 'uniform')

sLEHubIdx = np.argmin(abs(sLEData[:,1]-sLEInterXHub[1]))
sLECasIdx = np.argmin(abs(sLEData[:,1]-sLEInterXCas[1]))
newSLEData = np.vstack((sLEInterXHub, sLEData[sLEHubIdx+1:sLECasIdx], sLEInterXCas))
newSLEData = densifyCurve(newSLEData, res-2, 'uniform')

sTEHubIdx = np.argmin(abs(sTEData[:,1]-sTEInterXHub[1]))
sTECasIdx = np.argmin(abs(sTEData[:,1]-sTEInterXCas[1]))
newSTEData = np.vstack((sTEInterXHub, sTEData[sTEHubIdx+1:sTECasIdx], sTEInterXCas))
newSTEData = densifyCurve(newSTEData, res-2, 'uniform')
#%%Building the hubCut surface

modHub = hubData
modCas = np.vstack(([[hubData[0,0]-3*casData[0,1], casData[0,1]], casData]))

hubFunc = interp1d(modHub[:,0], modHub[:,1],)
casFunc = interp1d(modCas[:,0], modCas[:,1],)


hubR = hubFunc(hubZ)
distH = abs(rTEInterXHub[0] - hubZ)
casZ = rTEInterXCas[0] + distH
casR = casFunc(casZ)

stageExitH = hubData[-1]
hub = np.zeros((res-5,2))
# hub[:,0] = np.linspace(modHub[0,0], hubData[-1,0]+domainExtn*2*casData[0,1], res-5)
hub[:,0] = np.linspace(modHub[0,0], modHub[-1,0], res-5)
hub[:,1] = hubFunc(hub[:,0])


hub = insertPoints_batch(hub, [rLEInterXHub, rTEInterXHub,  sLEInterXHub, sTEInterXHub, stageExitH])
cas = np.zeros((res-5,2))
cas[:,0] = np.linspace(modCas[0,0] , modCas[-1,0], res-5)
cas[:,1] = casFunc(cas[:,0])

stageExitC = casData[-1]
cas = insertPoints_batch(cas, [rLEInterXCas, rTEInterXCas,  sLEInterXCas, sTEInterXCas, stageExitC])
#define the position for inviscid region 

inletR = np.linspace(hub[0,1], cas[0,1], rSections)
inletX = np.linspace(cas[0,0], cas[0,0], rSections)

outletR = np.linspace(hub[-1,1], cas[-1,1], rSections)
outletX = np.linspace(hub[-1,0], cas[-1,0], rSections)

offsetTE = np.zeros((len(newRTEData), 2)) # I am subtracting 5 because I will add the intersection to the offset hub and casing 
offsetTE[:,0], offsetTE[:,1] = scale_curve_endpoints(newRTEData[:,0], newRTEData[:,1], hubZ, hubR, casZ, casR)
#%% Generate the curves on the nose (The butterfly curve)
nose = hubData[0:spinnerIdx] #This is the definition of the spinner nose. You can adjust as you wish, this is where I decide to cut mine
# please note that this is an approximate solution, but it WORKS! The nose is not a perfect cone, so I made modification to make it work.
f=interp1d(nose[:,0], nose[:,1])
z = np.linspace(nose[0][0], nose[len(nose)-1][0], res)
r = f(z)
r[0] = 1e-10 #To avoid 0 radius
# A = f(-13)

Theta, R = np.meshgrid(theta, r)
X = R * np.cos(Theta)
Y = R * np.sin(Theta)
Z = np.tile(z, (res, 1)).T
slope = Slope(nose[0,0], nose[0,1], nose[-1,0], nose[-1,1]) #This computes the slope of the slant edge as though it is a perfect cone
interC = nose[-1,1] - (slope*nose[-1,0]) # This computes the intercept, so I can define the slant edge as a straight line 

rMid = 0.5*(nose[0,1] + nose[-1,1])
zMid = 0.5*(nose[0,0] + nose[-1,0])

func = interp1d(nose[:,0], nose[:,1])
noseMid = func(zMid)

#Now be very careful here. I chose 0.6 and 0.4 arbitrarily. I found this values to give me good cell quality, so what exactly are these numbers 
#The Tr point, is the radial value as I move up along the nose surface, basically if you think of the butterfly as a kite or diamond, this point defines the top of the diamond, you can move it 
#up or down 

#So for the Tz, again, I chose this value arbitrarily, it is the point along the axial that determines where the apex of the diamond lies, this is why I then interpolate the value to get the 
#value on the nose.
Tr = 0.6 * (nose[0,1] + nose[-1,1])
Tz = 0.4 * (nose[0,0] - nose[-1,0]) + nose[-1,0]
T = func(Tz)

Zt = (nose[-1,1]*0.6 - interC)/slope
Zmid = (nose[-1,1]*0.5 - interC)/slope

xC, yC, zC = cone_curve(noseMid, T, Tr, rMid, interC, slope, 16, res)

pVal = np.concatenate((xC, yC, zC)).reshape((-1,3), order='F')
nVal = np.concatenate((xC, (-1*yC), zC)).reshape((-1,3), order='F')


mValz = np.linspace(nose[0,0], nose[-1,0], res)
mValr = func(mValz)
mVal = np.column_stack((mValz, mValr))

extrusion = np.linspace(0, res*scales, res)  # How far to extrude in X-direction

Xp = np.tile(pVal[:,0], (res, 1)).T + extrusion
Yp = np.tile(pVal[:,1], (res, 1)).T
Zp = np.tile(pVal[:,2], (res, 1)).T

Xn = np.tile(nVal[:,0], (res, 1)).T + extrusion
Yn = np.tile(nVal[:,1], (res, 1)).T
Zn = np.tile(nVal[:,2], (res, 1)).T


Xm = np.tile(mValr, (res, 1)).T + extrusion  # Shift in X
Ym = np.tile(np.zeros(res), (res, 1)).T                 # Y stays same
Zm = np.tile(mValz, (res, 1)).T 

#Defines where inviscid region ends at the inlet and outlet
inviscidR = casFunc(nose[0,0]-domainExtn*2*casData[0,1]*percent)
rotorInletR = casFunc(nose[-1,0])

inviscidOutRc = casFunc(cas[-1,0]-domainExtn*2*casData[0,1]*percent)
inviscidOutRh = hubFunc(cas[-1,0]-domainExtn*2*casData[0,1]*percent)

mid = int(0.5*rTEData.shape[0])

domLen = rTEData[mid,0] - rLEData[mid,0] 
delThick, y1 = yPlusCalc(velocity, domLen, target_yPlus)

#Determine the offset Curve
offCas = np.column_stack((cas[:,0], cas[:,1]-10*delThick))

offHub = np.column_stack((hub[:,0], hub[:,1]+10*delThick))

offsetHub = np.zeros((res,res+5,3))
offsetCas = np.zeros((res,res+5,3))

#Determine the intersection of the cross passage curve on the offset 
rLEOffsetH = TwoLinesIntersect(offHub, newRLEData)
rTEOffsetH = TwoLinesIntersect(offHub, newRTEData)
hubCutOffsetH = TwoLinesIntersect(offHub, offsetTE)
sLEOffsetH = TwoLinesIntersect(offHub, newSLEData) 
sTEOffsetH = TwoLinesIntersect(offHub, newSTEData) 

rLEOffsetC = TwoLinesIntersect(offCas, newRLEData)
rTEOffsetC = TwoLinesIntersect(offCas, newRTEData)
hubCutOffsetC = TwoLinesIntersect(offCas, offsetTE)
sLEOffsetC = TwoLinesIntersect(offCas, newSLEData) 
sTEOffsetC = TwoLinesIntersect(offCas, newSTEData) 

newRLEData = insertPoints_batch(newRLEData, [rLEOffsetH, rLEOffsetC])
newRTEData = insertPoints_batch(newRTEData, [rTEOffsetH, rTEOffsetC])
offsetTE = insertPoints_batch(offsetTE, [hubCutOffsetH, hubCutOffsetC])
newSLEData = insertPoints_batch(newSLEData, [sLEOffsetH, sLEOffsetC])
newSTEData = insertPoints_batch(newSTEData, [sTEOffsetH, sTEOffsetC])

offCas = insertPoints_batch(offCas, [rLEOffsetC,rTEOffsetC,hubCutOffsetC,sLEOffsetC,sTEOffsetC])
offHub = insertPoints_batch(offHub, [rLEOffsetH,rTEOffsetH,hubCutOffsetH,sLEOffsetH,sTEOffsetH])

for c in range(res):
    offsetHub[c,:,:] = np.array(vectorRotZ3D(offHub[:,1], 0, offHub[:,0], theta[c])).T
    offsetCas[c,:,:] = np.array(vectorRotZ3D(offCas[:,1], 0, offCas[:,0], theta[c])).T
    
XhO = np.zeros([res, res+5])
YhO = np.zeros([res, res+5])
ZhO = np.zeros([res, res+5])
XcO = np.zeros([res, res+5])
YcO = np.zeros([res, res+5])
ZcO = np.zeros([res, res+5])
for d in range(res):
    for e in range(res+5):
        XhO[d,e] = offsetHub[d,e][0]
        YhO[d,e] = offsetHub[d,e][1]
        ZhO[d,e] = offsetHub[d,e][2]
        XcO[d,e] = offsetCas[d,e][0]
        YcO[d,e] = offsetCas[d,e][1]
        ZcO[d,e] = offsetCas[d,e][2]
        
modHub = insertPoints_batch(hub, [newRLEData[0], newRTEData[0], offsetTE[0], newSLEData[0], newSTEData[0]])
modCas = insertPoints_batch(cas, [newRLEData[-1], newRTEData[-1], offsetTE[-1], newSLEData[-1], newSTEData[-1]])
#%% Rotate all surfaces 
newNsection = len(newRTEData)
rotInlet = np.zeros((res, rSections, 3))
rotOutlet = np.zeros((res, rSections, 3))
rotLE = np.zeros((res, newNsection, 3))
rotTE = np.zeros((res, newNsection, 3))
staLE = np.zeros((res, newNsection, 3))
staTE = np.zeros((res, newNsection, 3))
hubCut = np.zeros((res, newNsection, 3))
rotHub = np.zeros((res, len(modHub), 3))
rotCas = np.zeros((res, len(modCas), 3))


for b in range(res):
    rotInlet[b,:,:] = np.array(vectorRotZ3D(inletR,0,inletX,theta[b])).T
    rotOutlet[b,:,:] = np.array(vectorRotZ3D(outletR,0,outletX,theta[b])).T
    rotLE[b,:,:] = np.array(vectorRotZ3D(newRLEData[:,1],0,newRLEData[:,0],theta[b])).T
    rotTE[b,:,:] = np.array(vectorRotZ3D(newRTEData[:,1],0,newRTEData[:,0],theta[b])).T
    staLE[b,:,:] = np.array(vectorRotZ3D(newSLEData[:,1],0,newSLEData[:,0],theta[b])).T
    staTE[b,:,:] = np.array(vectorRotZ3D(newSTEData[:,1],0,newSTEData[:,0],theta[b])).T
    hubCut[b,:,:] = np.array(vectorRotZ3D(offsetTE[:,1],0,offsetTE[:,0],theta[b])).T
    rotHub[b,:,:] = np.array(vectorRotZ3D(modHub[:,1],0,modHub[:,0],theta[b])).T
    rotCas[b,:,:] = np.array(vectorRotZ3D(modCas[:,1],0,modCas[:,0],theta[b])).T
    
#%%
Xrle = np.zeros([res, newNsection])
Yrle = np.zeros([res, newNsection])
Zrle = np.zeros([res, newNsection])
Xrte = np.zeros([res, newNsection])
Yrte = np.zeros([res, newNsection])
Zrte = np.zeros([res, newNsection])
Xsle = np.zeros([res, newNsection])
Ysle = np.zeros([res, newNsection])
Zsle = np.zeros([res, newNsection])
Xste = np.zeros([res, newNsection])
Yste = np.zeros([res, newNsection])
Zste = np.zeros([res, newNsection])
XhC = np.zeros([res, newNsection])
YhC = np.zeros([res, newNsection])
ZhC = np.zeros([res, newNsection])

Xi = np.zeros([res, rSections])
Yi = np.zeros([res, rSections])
Zi = np.zeros([res, rSections])
Xo = np.zeros([res, rSections])
Yo = np.zeros([res, rSections])
Zo = np.zeros([res, rSections])

Xh = np.zeros([res, len(modHub)])
Yh = np.zeros([res, len(modHub)])
Zh = np.zeros([res, len(modHub)])
Xc = np.zeros([res, len(modCas)])
Yc = np.zeros([res, len(modCas)])
Zc = np.zeros([res, len(modCas)])

for c in range(res):
    for d in range(newNsection):
        Xrle[c,d] = rotLE[c,d][0]
        Yrle[c,d] = rotLE[c,d][1]
        Zrle[c,d] = rotLE[c,d][2]
        Xrte[c,d] = rotTE[c,d][0]
        Yrte[c,d] = rotTE[c,d][1]
        Zrte[c,d] = rotTE[c,d][2]   
        Xsle[c,d] = staLE[c,d][0]
        Ysle[c,d] = staLE[c,d][1]
        Zsle[c,d] = staLE[c,d][2]
        Xste[c,d] = staTE[c,d][0]
        Yste[c,d] = staTE[c,d][1]
        Zste[c,d] = staTE[c,d][2] 
        XhC[c,d] = hubCut[c,d][0]
        YhC[c,d] = hubCut[c,d][1]
        ZhC[c,d] = hubCut[c,d][2]         
        
    for e in range(len(modHub)):
        Xh[c,e] = rotHub[c,e][0]
        Yh[c,e] = rotHub[c,e][1]
        Zh[c,e] = rotHub[c,e][2]
        Xc[c,e] = rotCas[c,e][0]
        Yc[c,e] = rotCas[c,e][1]
        Zc[c,e] = rotCas[c,e][2]
      
    for f in range(rSections):
        Xi[c,f] = rotInlet[c,f][0]
        Yi[c,f] = rotInlet[c,f][1]
        Zi[c,f] = rotInlet[c,f][2]
        Xo[c,f] = rotOutlet[c,f][0]
        Yo[c,f] = rotOutlet[c,f][1]
        Zo[c,f] = rotOutlet[c,f][2]

#%%  Generate the stls 
filenames = ['casing', 'hubCut', 'Hub', 'LES', 'LE', 'TES', 'TE', 'mPer', 'nVal', 'pVal', 'outlet', 'inlet', 'offsetHub', 'offsetCas']

Xvalues = [Xc, XhC, Xh, Xsle, Xrle, Xste, Xrte, Xm, Xn, Xp, Xo, Xi, XhO, XcO]
Yvalues = [Yc, YhC, Yh, Ysle, Yrle, Yste, Yrte, Ym, Yn, Yp, Yo, Yi, YhO, YcO]
Zvalues = [Zc, ZhC, Zh, Zsle, Zrle, Zste, Zrte, Zm, Zn, Zp, Zo, Zi, ZhO, ZcO]

for qq in range(len(Xvalues)):
    filename = '../fullWheelRotorStator/constant/geometry/{}.stl'.format(filenames[qq])
    rows = Zvalues[qq].shape[0]
    columns = Zvalues[qq].shape[1]
    X = Xvalues[qq]
    Y = Yvalues[qq]
    Z = Zvalues[qq]

    numFacets = 0

    file = open(filename, 'w')
    file.write('solid \n')

    def unitVector(file, p1, p2, p3, tolerance=1e-10):
        """
        Calculate unit normal vector and write facet.
        Returns True if facet was written, False if degenerate.
        """
        # Check if points are coincident (happens at zero radius)
        if np.allclose(p1, p2, atol=tolerance) or \
           np.allclose(p2, p3, atol=tolerance) or \
           np.allclose(p1, p3, atol=tolerance):
            return False
        
        # VECTORS TANGENT TO FACET
        vector1 = p3 - p2
        vector2 = p3 - p1

        normalVec = np.cross(vector1, vector2)
        magnitude = np.linalg.norm(normalVec)
        
        # Check for zero or near-zero magnitude (degenerate triangle)
        if magnitude < tolerance:
            return False
        
        unitVec = normalVec / magnitude
        
        file.write(f'facet normal {unitVec[0]} {unitVec[1]} {unitVec[2]} \n'
                   f'outer loop \n'
                   f'vertex {p1[0]} {p1[1]} {p1[2]} \n'
                   f'vertex {p2[0]} {p2[1]} {p2[2]} \n'
                   f'vertex {p3[0]} {p3[1]} {p3[2]} \n'
                   f'endloop \n'
                   f'endfacet \n')
        return True

    for i in range(rows - 1):
        for j in range(columns - 1):
            # FACET A VERTICES
            p1 = np.asarray([X[i, j], Y[i, j], Z[i, j]])
            p2 = np.asarray([X[i, j+1], Y[i, j+1], Z[i, j+1]])
            p3 = np.asarray([X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]])

            if unitVector(file, p1, p2, p3):
                numFacets += 1

            # FACET B VERTICES
            p1 = np.asarray([X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1]])
            p2 = np.asarray([X[i+1, j], Y[i+1, j], Z[i+1, j]])
            p3 = np.asarray([X[i, j], Y[i, j], Z[i, j]])

            if unitVector(file, p1, p2, p3):
                numFacets += 1

    file.write('endsolid')
    file.close()
    
    print(f"Written {numFacets} facets to {filename}")


#%% Determine grading properties 
k = 10 # This gives more room for the BL thickness, allowing for more cells there and a smoother transition 
frac = 0.3
Nbc = int(frac*Nr) #boundary layer cells

domHeight = rotorInletR - Xm[-1][0]

blFrac = k*delThick/domHeight
r, yn = commonRatio(y1, k*delThick, Nbc)
gradH = yn/y1
gradC = y1/yn
newNbc = 1 + np.log(yn/y1)/np.log(r)

Hcore = domHeight - 2*k*delThick
Ncore = Hcore/yn
totalCellCount = Ncore + 2*newNbc
xdP = blFrac

xdN = newNbc/totalCellCount
xdupG = gradH
xddwG = gradC

#Define Grading properties 
blRadial = np.round(totalCellCount)
radial = np.round(xdN*blRadial)
coreRadial = np.round(blRadial - 2*radial)

diamondCellSize = Tr/Nq
lowMidGrad = yn/ diamondCellSize
butterflyGrad = yn/diamondCellSize
#%%
# Nbl = frac*Nz
mid = int(0.5*rTEData.shape[0])
frac1 = 0.15
Nbl = frac1*Nz

#for block2
blkLen2 = rLEData[0,0] - nose[-1,0]
NbcX2 = 2*(Nbl)
blkFrac2 = k*delThick/blkLen2
r2, yn2 = commonRatio(y1, k*delThick, NbcX2)
gradHBlk2 = yn2/y1
gradCBlk2 = y1/yn2
newNbcX2 = 1 + np.log(yn2/y1)/np.log(r2)
HcoreBlk2 = blkLen2 - k*delThick
NcoreBlk2 = HcoreBlk2/yn2
blk2CellCount = NcoreBlk2 + newNbcX2

#for blockS
blkLenS =  nose[-1,0] - casData[0,0]
NbcXS = 5*(Nbl)
blkFracS = 1
rS, ynS = commonRatio(yn2, blkLenS , NbcXS)
gradHBlkS = ynS/yn2
gradCBlkS = yn2/ynS
newNbcXS = 1 + np.log(ynS/yn2)/np.log(rS)
HcoreBlkS = blkLenS - k*delThick
NcoreBlkS = HcoreBlkS/ynS
blkSCellCount = 1*newNbcXS
#for block1

blkLen1 = casData[0,0]  - cas[0,0] 
NbcX1 = 5*(Nbl)
blkFrac1 = k*delThick/blkLen1
r1, yn1 = commonRatio(ynS, blkLen1, NbcX1)
# yn1 = 20*y1
gradHBlk1 = yn1/ynS
gradCBlk1 = ynS/yn1
newNbcX1 = 1 + np.log(yn1/ynS)/np.log(r1)
HcoreBlk1 = blkLen1 - k*delThick
NcoreBlk1 = HcoreBlk1/yn1
blk1CellCount = newNbcX1

#First find the cell height of the wedge block
blkHieghtBaseBlk = noseMid*np.cos(0.5*pitchAngle)
coreHeight = rotorInletR - Xm[-1][0] 

yFirst = blkHieghtBaseBlk/(0.5*Nq)
G = yn/yFirst
blkHeightMid = (nose[-1,1]+k*delThick - T)
AA = np.log((yn - blkHeightMid*G)/(yn - blkHeightMid))/np.log(G)
newNbcMid = AA/(AA-1)


#for block3
blkLen3 = rTEData[mid,0] - rLEData[mid,0] 
NbcX3 = 2.5*(Nbl)
blkFrac3 = k*delThick/blkLen3
r3, yn3 = commonRatio(y1, k*delThick, NbcX3)
gradHBlk3 = yn3/y1
gradCBlk3 = y1/yn3
newNbcX3 = 1 + np.log(yn3/y1)/np.log(r3)
HcoreBlk3 = blkLen3 - 2*k*delThick
NcoreBlk3 = HcoreBlk3/yn3
blk3CellCount = NcoreBlk3 + 2*newNbcX3
#for block4
blkLen4 = offsetTE[0,0] - rTEData[0,0] 
NbcX4 = 2*(Nbl)
blkFrac4 = delThick/blkLen4
r4, yn4 = commonRatio(y1, delThick, NbcX4)
gradHBlk4 = yn4/y1
gradCBlk4 = y1/yn4
newNbcX4 = 1 + np.log(yn4/y1)/np.log(r4)

HcoreBlk4 = blkLen4 - 0
NcoreBlk4 = HcoreBlk4/yn4
blk4CellCount = NcoreBlk4 + newNbcX4
#for block5
blkLen5 = sLEData[mid,0] - rTEData[mid,0]  
NbcX5 = 2*(Nbl)
blkFrac5 = k*delThick/blkLen5
r5, yn5 = commonRatio(y1, k*delThick, NbcX5)
gradHBlk5 = yn5/y1
gradCBlk5 = y1/yn5
newNbcX5 = 1 + np.log(yn5/y1)/np.log(r5)
HcoreBlk5 = blkLen5 - 2*k*delThick
NcoreBlk5 = HcoreBlk5/yn5
blk5CellCount = NcoreBlk5 + 2*newNbcX5
#for block6
blkLen6 = sTEData[mid,0] - sLEData[mid,0] 
NbcX6 = 2.5*(Nbl)
blkFrac6 = k*delThick/blkLen6
r6, yn6 = commonRatio(y1, k*delThick, NbcX6)
gradHBlk6 = yn6/y1
gradCBlk6 = y1/yn6
newNbcX6 = 1 + np.log(yn6/y1)/np.log(r6)
HcoreBlk6 = blkLen6 - 2*k*delThick
NcoreBlk6 = HcoreBlk6/yn6
blk6CellCount = NcoreBlk6 + 2*newNbcX6
#for block7

blkLen7 =hub[-1,0] - sTEData[0,0]
NbcX7 = 3*(Nbl)
blkFrac7 = 2*k*delThick/blkLen7
r7, yn7 = commonRatio(y1, 2*k*delThick, NbcX7)
gradHBlk7 = yn7/y1
gradCBlk7 = y1/yn7
newNbcX7 = 1 + np.log(yn7/y1)/np.log(r7)
HcoreBlk7 = blkLen7 - 2*k*delThick
NcoreBlk7 = HcoreBlk7/yn7
blk7CellCount =  NcoreBlk7 + 1*newNbcX7

#for block8
blkLen8 = hub[-1,0] - (hub[-1,0] -domainExtn*2*casData[0,1]*percent)
NbcX8 = 0.12*(Nbl)
blkFrac8 = delThick/blkLen8
r8, yn8 = commonRatio(yn7, blkLen8, NbcX8)
gradHBlk8 = yn8/yn7
gradCBlk8 = yn7/yn8
newNbcX8 = 1 + np.log(yn8/yn7)/np.log(r8)
HcoreBlk8 = blkLen8 - blkLen8
NcoreBlk8 = HcoreBlk8/yn8
blk8CellCount = NcoreBlk8 + newNbcX8

#diamond block

ptX1 = np.round(blk1CellCount)
ptXS = np.round(blkSCellCount)
ptX2 = np.round(blk2CellCount)
ptX3 = np.round(blk3CellCount)
ptX4 = np.round(blk4CellCount)
ptX5 = np.round(blk5CellCount)
ptX6 = np.round(blk6CellCount)
ptX7 = np.round(blk7CellCount)
ptX8 = np.round(blk8CellCount)

xS1 = blkFrac1
xSP = blkFracS
x2P = blkFrac2
x3P = blkFrac3
x4P = blkFrac4
x5P = blkFrac5
x6P = blkFrac6
x7P = blkFrac7
x8P = blkFrac8

x1N = newNbcX1/blk1CellCount
xSN = newNbcXS/blkSCellCount
x2N = newNbcX2/blk2CellCount
x3N = newNbcX3/blk3CellCount
x4N = newNbcX4/blk4CellCount
x5N = newNbcX5/blk5CellCount
x6N = newNbcX6/blk6CellCount
x7N = newNbcX7/blk7CellCount
x8N = newNbcX8/blk8CellCount

x1upG = gradHBlk1
xSupG = gradHBlkS
x2upG = gradHBlk2
x3upG = gradHBlk3
x4upG = gradHBlk4
x5upG = gradHBlk5
x6upG = gradHBlk6
x7upG = gradHBlk7
x8upG = gradHBlk8

x1dwG = gradCBlk1
xSdwG = gradCBlkS
x2dwG = gradCBlk2
x3dwG = gradCBlk3
x4dwG = gradCBlk4
x5dwG = gradCBlk5
x6dwG = gradCBlk6
x7dwG = gradCBlk7
x8dwG = gradCBlk8
#%% Define the vertices 
pAng = 0.5*pitchAngle
nAng = -0.5*pitchAngle
xAngP = np.cos((pAng))
xAngN = np.cos((nAng))
yAngP = np.sin((pAng))
yAngN = np.sin((nAng))

x0y0m = [inletR[0], 0, inletX[0]]
x1y0m = [inletR[0], 0, hub[0,0]-k*delThick]
x2y0m = [inletR[0], 0, hub[0,0]]

x0y1m = [T, 0, inletX[0]]
x0y1p = [noseMid*xAngP, noseMid*yAngP, inletX[0]]
x0y1n = [noseMid*xAngN, noseMid*yAngN, inletX[0]]
x1y1m = [T, 0, Zt-k*delThick]
x1y1p = [noseMid*xAngP, noseMid*yAngP, Zmid-k*delThick]
x1y1n = [noseMid*xAngN, noseMid*yAngN, Zmid-k*delThick]
x2y1m = [T, 0, Zt]
x2y1p = [noseMid*xAngP, noseMid*yAngP, Zmid]
x2y1n = [noseMid*xAngN, noseMid*yAngN, Zmid]
x0y2m = [nose[-1,1]+k*delThick, 0, inletX[0]]
x0y2p = [(nose[-1,1]+k*delThick)*xAngP, (nose[-1,1]+k*delThick)*yAngP, inletX[0]]
x0y2n = [(nose[-1,1]+k*delThick)*xAngN, (nose[-1,1]+k*delThick)*yAngN, inletX[0]]

x2y2m = [Xm[-1][0], Ym[-1][0], nose[-1,0]]
x2y2p = [Xm[-1][0]*xAngP, Xm[-1][0]*yAngP, nose[-1,0]]
x2y2n = [Xm[-1][0]*xAngN, Xm[-1][0]*yAngN, nose[-1,0]]

x0y4m = [cas[:,1][0]-k*delThick, 0, cas[:,0][0]]
x0y4p = [(cas[:,1][0]-k*delThick)*xAngP, (cas[:,1][0]-k*delThick)*yAngP, cas[:,0][0]]
x0y4n = [(cas[:,1][0]-k*delThick)*xAngN, (cas[:,1][0]-k*delThick)*yAngN, cas[:,0][0]]

x2y3m = [Xm[-1][0]+k*delThick, Ym[-1][0], nose[-1,0]]
x2y3p = [(Xm[-1][0]+k*delThick)*xAngP, (Xm[-1][0]+k*delThick)*yAngP, nose[-1,0]]
x2y3n = [(Xm[-1][0]+k*delThick)*xAngN, (Xm[-1][0]+k*delThick)*yAngN, nose[-1,0]]

x0y5m = [cas[:,1][0], 0, cas[:,0][0]]
x0y5p = [cas[:,1][0]*xAngP, cas[:,1][0]*yAngP, cas[:,0][0]]
x0y5n = [cas[:,1][0]*xAngN, cas[:,1][0]*yAngN, cas[:,0][0]]

x2y4m = [rotorInletR-k*delThick, 0, nose[-1,0]]
x2y4p = [(rotorInletR-k*delThick)*xAngP, (rotorInletR-k*delThick)*yAngP,  nose[-1,0]]
x2y4n = [(rotorInletR-k*delThick)*xAngN, (rotorInletR-k*delThick)*yAngN,  nose[-1,0]]

x2y5m = [rotorInletR, 0,nose[-1,0]]
x2y5p = [rotorInletR*xAngP, rotorInletR*yAngP, nose[-1,0]]
x2y5n = [rotorInletR*xAngN, rotorInletR*yAngN, nose[-1,0]]


xSy0m = [hub[0,1], 0, casData[0,0]]

xSy1m = [T, 0, casData[0,0]]
xSy1p = [noseMid*xAngP, noseMid*yAngP, casData[0,0]]
xSy1n = [noseMid*xAngN, noseMid*yAngN, casData[0,0]]

xSy2m = [nose[-1,1]+k*delThick, 0,  casData[0,0]]
xSy2p = [(nose[-1,1]+k*delThick)*xAngP, (nose[-1,1]+k*delThick)*yAngP,  casData[0,0]]
xSy2n = [(nose[-1,1]+k*delThick)*xAngN, (nose[-1,1]+k*delThick)*yAngN,  casData[0,0]]

xSy4m = [inviscidR-k*delThick, 0,  casData[0,0]]
xSy4p = [(inviscidR-k*delThick)*xAngP, (inviscidR-k*delThick)*yAngP,  casData[0,0]]
xSy4n = [(inviscidR-k*delThick)*xAngN, (inviscidR-k*delThick)*yAngN,  casData[0,0]]

xSy5m = [inviscidR, 0,  casData[0,0]]
xSy5p = [inviscidR*xAngP, inviscidR*yAngP,  casData[0,0]]
xSy5n = [inviscidR*xAngN, inviscidR*yAngN,  casData[0,0]]

x3y2m = [hub[-1,1], 0, hub[-1,0]]
x3y2p = [hub[-1,1]*xAngP, hub[-1,1]*yAngP, hub[-1,0]]
x3y2n = [hub[-1,1]*xAngN, hub[-1,1]*yAngN, hub[-1,0]]

x3y3m = [hub[-1,1]+k*delThick, 0, hub[-1,0]]
x3y3p = [(hub[-1,1]+k*delThick)*xAngP, (hub[-1,1]+k*delThick)*yAngP, hub[-1,0]]
x3y3n = [(hub[-1,1]+k*delThick)*xAngN, (hub[-1,1]+k*delThick)*yAngN, hub[-1,0]]

x4y2m = [hub[-1,1], 0, hub[-1,0]]
x4y2p = [hub[-1,1]*xAngP, hub[-1,1]*yAngP, hub[-1,0]]
x4y2n = [hub[-1,1]*xAngN, hub[-1,1]*yAngN, hub[-1,0]]

x4y3m = [hub[-1,1]+k*delThick, 0, hub[-1,0]]
x4y3p = [(hub[-1,1]+k*delThick)*xAngP, (hub[-1,1]+k*delThick)*yAngP, hub[-1,0]]
x4y3n = [(hub[-1,1]+k*delThick)*xAngN, (hub[-1,1]+k*delThick)*yAngN, hub[-1,0]]

x3y4m = [cas[-1,1]-k*delThick, 0, cas[-1,0]]
x3y4p = [(cas[-1,1]-k*delThick)*xAngP, (cas[-1,1]-k*delThick)*yAngP, cas[-1,0]]
x3y4n = [(cas[-1,1]-k*delThick)*xAngN, (cas[-1,1]-k*delThick)*yAngN, cas[-1,0]]

x3y5m = [cas[-1,1], 0, cas[-1,0]]
x3y5p = [cas[-1,1]*xAngP, cas[-1,1]*yAngP, cas[-1,0]]
x3y5n = [cas[-1,1]*xAngN, cas[-1,1]*yAngN, cas[-1,0]]

x4y4m = [cas[-1,1]-k*delThick, 0, cas[-1,0]]
x4y4p = [(cas[-1,1]-k*delThick)*xAngP, (cas[-1,1]-k*delThick)*yAngP, cas[-1,0]]
x4y4n = [(cas[-1,1]-k*delThick)*xAngN, (cas[-1,1]-k*delThick)*yAngN, cas[-1,0]]

x4y5m = [cas[-1,1], 0, cas[-1,0]]
x4y5p = [cas[-1,1]*xAngP, cas[-1,1]*yAngP, cas[-1,0]]
x4y5n = [cas[-1,1]*xAngN, cas[-1,1]*yAngN, cas[-1,0]]

xr1y2m = [newRLEData[0,1], 0, newRLEData[0,0]]
xr1y2p = [newRLEData[0,1]*xAngP, newRLEData[0,1]*yAngP, newRLEData[0,0]]
xr1y2n = [newRLEData[0,1]*xAngN, newRLEData[0,1]*yAngN, newRLEData[0,0]]

xr1y3m = [rLEOffsetH[1], 0, rLEOffsetH[0]]
xr1y3p = [rLEOffsetH[1]*xAngP, rLEOffsetH[1]*yAngP, rLEOffsetH[0]]
xr1y3n = [rLEOffsetH[1]*xAngN, rLEOffsetH[1]*yAngN, rLEOffsetH[0]]

xr2y2m = [newRTEData[0,1], 0, newRTEData[0,0]]
xr2y2p = [newRTEData[0,1]*xAngP, newRTEData[0,1]*yAngP, newRTEData[0,0]]
xr2y2n = [newRTEData[0,1]*xAngN, newRTEData[0,1]*yAngN, newRTEData[0,0]]

xr2y3m = [rTEOffsetH[1], 0, rTEOffsetH[0]]
xr2y3p = [rTEOffsetH[1]*xAngP, rTEOffsetH[1]*yAngP, rTEOffsetH[0]]
xr2y3n = [rTEOffsetH[1]*xAngN, rTEOffsetH[1]*yAngN, rTEOffsetH[0]]

xr1y4m = [rLEOffsetC[1], 0, rLEOffsetC[0]]
xr1y4p = [rLEOffsetC[1]*xAngP, rLEOffsetC[1]*yAngP, rLEOffsetC[0]]
xr1y4n = [rLEOffsetC[1]*xAngN, rLEOffsetC[1]*yAngN, rLEOffsetC[0]]

xr1y5m = [newRLEData[-1,1], 0, newRLEData[-1,0]]
xr1y5p = [newRLEData[-1,1]*xAngP, newRLEData[-1,1]*yAngP, newRLEData[-1,0]]
xr1y5n = [newRLEData[-1,1]*xAngN, newRLEData[-1,1]*yAngN, newRLEData[-1,0]]

xr2y4m = [rTEOffsetC[1], 0, rTEOffsetC[0]]
xr2y4p = [rTEOffsetC[1]*xAngP, rTEOffsetC[1]*yAngP, rTEOffsetC[0]]
xr2y4n = [rTEOffsetC[1]*xAngN, rTEOffsetC[1]*yAngN, rTEOffsetC[0]]

xr2y5m = [newRTEData[-1,1], 0, newRTEData[-1,0]]
xr2y5p = [newRTEData[-1,1]*xAngP, newRTEData[-1,1]*yAngP, newRTEData[-1,0]]
xr2y5n = [newRTEData[-1,1]*xAngN, newRTEData[-1,1]*yAngN, newRTEData[-1,0]]

xs1y2m = [newSLEData[0,1], 0, newSLEData[0,0]]
xs1y2p = [newSLEData[0,1]*xAngP, newSLEData[0,1]*yAngP, newSLEData[0,0]]
xs1y2n = [newSLEData[0,1]*xAngN, newSLEData[0,1]*yAngN, newSLEData[0,0]]

xs1y3m = [sLEOffsetH[1], 0, sLEOffsetH[0]]
xs1y3p = [sLEOffsetH[1]*xAngP, sLEOffsetH[1]*yAngP, sLEOffsetH[0]]
xs1y3n = [sLEOffsetH[1]*xAngN, sLEOffsetH[1]*yAngN, sLEOffsetH[0]]

xs2y2m = [newSTEData[0,1], 0, newSTEData[0,0]]
xs2y2p = [newSTEData[0,1]*xAngP, newSTEData[0,1]*yAngP, newSTEData[0,0]]
xs2y2n = [newSTEData[0,1]*xAngN, newSTEData[0,1]*yAngN, newSTEData[0,0]]

xs2y3m = [sTEOffsetH[1], 0, sTEOffsetH[0]]
xs2y3p = [sTEOffsetH[1]*xAngP, sTEOffsetH[1]*yAngP, sTEOffsetH[0]]
xs2y3n = [sTEOffsetH[1]*xAngN, sTEOffsetH[1]*yAngN, sTEOffsetH[0]]

xs1y4m = [sLEOffsetC[1], 0, sLEOffsetC[0]]
xs1y4p = [sLEOffsetC[1]*xAngP, sLEOffsetC[1]*yAngP, sLEOffsetC[0]]
xs1y4n = [sLEOffsetC[1]*xAngN, sLEOffsetC[1]*yAngN, sLEOffsetC[0]]

xs1y5m = [newSLEData[-1,1], 0, newSLEData[-1,0]]
xs1y5p = [newSLEData[-1,1]*xAngP, newSLEData[-1,1]*yAngP, newSLEData[-1,0]]
xs1y5n = [newSLEData[-1,1]*xAngN, newSLEData[-1,1]*yAngN, newSLEData[-1,0]]

xs2y4m = [sTEOffsetC[1], 0, sTEOffsetC[0]]
xs2y4p = [sTEOffsetC[1]*xAngP, sTEOffsetC[1]*yAngP, sTEOffsetC[0]]
xs2y4n = [sTEOffsetC[1]*xAngN, sTEOffsetC[1]*yAngN, sTEOffsetC[0]]

xs2y5m = [newSTEData[-1,1], 0, newSTEData[-1,0]]
xs2y5p = [newSTEData[-1,1]*xAngP, newSTEData[-1,1]*yAngP, newSTEData[-1,0]]
xs2y5n = [newSTEData[-1,1]*xAngN, newSTEData[-1,1]*yAngN, newSTEData[-1,0]]

xCy2m = [offsetTE[0,1], 0, offsetTE[0,0]]
xCy2p = [offsetTE[0,1]*xAngP, offsetTE[0,1]*yAngP, offsetTE[0,0]]
xCy2n = [offsetTE[0,1]*xAngN, offsetTE[0,1]*yAngN, offsetTE[0,0]]

xCy3m = [hubCutOffsetH[1], 0, hubCutOffsetH[0]]
xCy3p = [hubCutOffsetH[1]*xAngP, hubCutOffsetH[1]*yAngP, hubCutOffsetH[0]]
xCy3n = [hubCutOffsetH[1]*xAngN, hubCutOffsetH[1]*yAngN, hubCutOffsetH[0]]

xCy4m = [hubCutOffsetC[1], 0, hubCutOffsetC[0]]
xCy4p = [hubCutOffsetC[1]*xAngP, hubCutOffsetC[1]*yAngP, hubCutOffsetC[0]]
xCy4n = [hubCutOffsetC[1]*xAngN, hubCutOffsetC[1]*yAngN, hubCutOffsetC[0]]

xCy5m = [offsetTE[-1,1], 0, offsetTE[-1,0]]
xCy5p = [offsetTE[-1,1]*xAngP, offsetTE[-1,1]*yAngP, offsetTE[-1,0]]
xCy5n = [offsetTE[-1,1]*xAngN, offsetTE[-1,1]*yAngN, offsetTE[-1,0]]

x0y2ArcP =  [(nose[-1,1]+k*delThick)*np.cos(np.deg2rad(pAng*0.5)), (nose[-1,1]+k*delThick)*np.sin(np.deg2rad(pAng*0.5)), inletX[0]]
x0y2ArcN =  [(nose[-1,1]+k*delThick)*np.cos(np.deg2rad(nAng*0.5)), (nose[-1,1]+k*delThick)*np.sin(np.deg2rad(nAng*0.5)), inletX[0]]


#%%

f = open('../fullWheelRotorStator/system/pointData', 'w')

f.write('x0y0m ({} {} {}); \n'.format(x0y0m[0], x0y0m[1], x0y0m[2]))
f.write('x1y0m ({} {} {}); \n'.format(x1y0m[0], x1y0m[1], x1y0m[2]))
f.write('x2y0m ({} {} {}); \n'.format(x2y0m[0], x2y0m[1], x2y0m[2]))
f.write('x0y1m ({} {} {}); \n'.format(x0y1m[0], x0y1m[1], x0y1m[2]))
f.write('x0y1p ({} {} {}); \n'.format(x0y1p[0], x0y1p[1], x0y1p[2]))
f.write('x0y1n ({} {} {}); \n'.format(x0y1n[0], x0y1n[1], x0y1n[2]))
f.write('x1y1m ({} {} {}); \n'.format(x1y1m[0], x1y1m[1], x1y1m[2]))
f.write('x1y1p ({} {} {}); \n'.format(x1y1p[0], x1y1p[1], x1y1p[2]))
f.write('x1y1n ({} {} {}); \n'.format(x1y1n[0], x1y1n[1], x1y1n[2]))
f.write('x2y1m ({} {} {}); \n'.format(x2y1m[0], x2y1m[1], x2y1m[2]))
f.write('x2y1p ({} {} {}); \n'.format(x2y1p[0], x2y1p[1], x2y1p[2]))
f.write('x2y1n ({} {} {}); \n'.format(x2y1n[0], x2y1n[1], x2y1n[2]))
f.write('x0y2m ({} {} {}); \n'.format(x0y2m[0], x0y2m[1], x0y2m[2]))
f.write('x0y2p ({} {} {}); \n'.format(x0y2p[0], x0y2p[1], x0y2p[2]))
f.write('x0y2n ({} {} {}); \n'.format(x0y2n[0], x0y2n[1], x0y2n[2]))
f.write('x2y2m ({} {} {}); \n'.format(x2y2m[0], x2y2m[1], x2y2m[2]))
f.write('x2y2p ({} {} {}); \n'.format(x2y2p[0], x2y2p[1], x2y2p[2]))
f.write('x2y2n ({} {} {}); \n'.format(x2y2n[0], x2y2n[1], x2y2n[2]))
f.write('x0y4m ({} {} {}); \n'.format(x0y4m[0], x0y4m[1], x0y4m[2]))
f.write('x0y4p ({} {} {}); \n'.format(x0y4p[0], x0y4p[1], x0y4p[2]))
f.write('x0y4n ({} {} {}); \n'.format(x0y4n[0], x0y4n[1], x0y4n[2]))
f.write('x2y3m ({} {} {}); \n'.format(x2y3m[0], x2y3m[1], x2y3m[2]))
f.write('x2y3p ({} {} {}); \n'.format(x2y3p[0], x2y3p[1], x2y3p[2]))
f.write('x2y3n ({} {} {}); \n'.format(x2y3n[0], x2y3n[1], x2y3n[2]))
f.write('x0y5m ({} {} {}); \n'.format(x0y5m[0], x0y5m[1], x0y5m[2]))
f.write('x0y5p ({} {} {}); \n'.format(x0y5p[0], x0y5p[1], x0y5p[2]))
f.write('x0y5n ({} {} {}); \n'.format(x0y5n[0], x0y5n[1], x0y5n[2]))
f.write('x2y4m ({} {} {}); \n'.format(x2y4m[0], x2y4m[1], x2y4m[2]))
f.write('x2y4p ({} {} {}); \n'.format(x2y4p[0], x2y4p[1], x2y4p[2]))
f.write('x2y4n ({} {} {}); \n'.format(x2y4n[0], x2y4n[1], x2y4n[2]))
f.write('x2y5m ({} {} {}); \n'.format(x2y5m[0], x2y5m[1], x2y5m[2]))
f.write('x2y5p ({} {} {}); \n'.format(x2y5p[0], x2y5p[1], x2y5p[2]))
f.write('x2y5n ({} {} {}); \n'.format(x2y5n[0], x2y5n[1], x2y5n[2]))
f.write('xSy0m ({} {} {}); \n'.format(xSy0m[0], xSy0m[1], xSy0m[2]))
f.write('xSy1m ({} {} {}); \n'.format(xSy1m[0], xSy1m[1], xSy1m[2]))
f.write('xSy1p ({} {} {}); \n'.format(xSy1p[0], xSy1p[1], xSy1p[2]))
f.write('xSy1n ({} {} {}); \n'.format(xSy1n[0], xSy1n[1], xSy1n[2]))
f.write('xSy2m ({} {} {}); \n'.format(xSy2m[0], xSy2m[1], xSy2m[2]))
f.write('xSy2p ({} {} {}); \n'.format(xSy2p[0], xSy2p[1], xSy2p[2]))
f.write('xSy2n ({} {} {}); \n'.format(xSy2n[0], xSy2n[1], xSy2n[2]))
f.write('xSy4m ({} {} {}); \n'.format(xSy4m[0], xSy4m[1], xSy4m[2]))
f.write('xSy4p ({} {} {}); \n'.format(xSy4p[0], xSy4p[1], xSy4p[2]))
f.write('xSy4n ({} {} {}); \n'.format(xSy4n[0], xSy4n[1], xSy4n[2]))
f.write('xSy5m ({} {} {}); \n'.format(xSy5m[0], xSy5m[1], xSy5m[2]))
f.write('xSy5p ({} {} {}); \n'.format(xSy5p[0], xSy5p[1], xSy5p[2]))
f.write('xSy5n ({} {} {}); \n'.format(xSy5n[0], xSy5n[1], xSy5n[2]))
f.write('x3y2m ({} {} {}); \n'.format(x3y2m[0], x3y2m[1], x3y2m[2]))
f.write('x3y2p ({} {} {}); \n'.format(x3y2p[0], x3y2p[1], x3y2p[2]))
f.write('x3y2n ({} {} {}); \n'.format(x3y2n[0], x3y2n[1], x3y2n[2]))
f.write('x3y3m ({} {} {}); \n'.format(x3y3m[0], x3y3m[1], x3y3m[2]))
f.write('x3y3p ({} {} {}); \n'.format(x3y3p[0], x3y3p[1], x3y3p[2]))
f.write('x3y3n ({} {} {}); \n'.format(x3y3n[0], x3y3n[1], x3y3n[2]))
f.write('x3y4m ({} {} {}); \n'.format(x3y4m[0], x3y4m[1], x3y4m[2]))
f.write('x3y4p ({} {} {}); \n'.format(x3y4p[0], x3y4p[1], x3y4p[2]))
f.write('x3y4n ({} {} {}); \n'.format(x3y4n[0], x3y4n[1], x3y4n[2]))
f.write('x4y2m ({} {} {}); \n'.format(x4y2m[0], x4y2m[1], x4y2m[2]))
f.write('x4y2p ({} {} {}); \n'.format(x4y2p[0], x4y2p[1], x4y2p[2]))
f.write('x4y2n ({} {} {}); \n'.format(x4y2n[0], x4y2n[1], x4y2n[2]))
f.write('x4y3m ({} {} {}); \n'.format(x4y3m[0], x4y3m[1], x4y3m[2]))
f.write('x4y3p ({} {} {}); \n'.format(x4y3p[0], x4y3p[1], x4y3p[2]))
f.write('x4y3n ({} {} {}); \n'.format(x4y3n[0], x4y3n[1], x4y3n[2]))
f.write('x4y4m ({} {} {}); \n'.format(x4y4m[0], x4y4m[1], x4y4m[2]))
f.write('x4y4p ({} {} {}); \n'.format(x4y4p[0], x4y4p[1], x4y4p[2]))
f.write('x4y4n ({} {} {}); \n'.format(x4y4n[0], x4y4n[1], x4y4n[2]))
f.write('x3y5m ({} {} {}); \n'.format(x3y5m[0], x3y5m[1], x3y5m[2]))
f.write('x3y5p ({} {} {}); \n'.format(x3y5p[0], x3y5p[1], x3y5p[2]))
f.write('x3y5n ({} {} {}); \n'.format(x3y5n[0], x3y5n[1], x3y5n[2]))
f.write('x4y5m ({} {} {}); \n'.format(x4y5m[0], x4y5m[1], x4y5m[2]))
f.write('x4y5p ({} {} {}); \n'.format(x4y5p[0], x4y5p[1], x4y5p[2]))
f.write('x4y5n ({} {} {}); \n'.format(x4y5n[0], x4y5n[1], x4y5n[2]))
f.write('xr1y2m ({} {} {}); \n'.format(xr1y2m[0], xr1y2m[1], xr1y2m[2]))
f.write('xr1y2p ({} {} {}); \n'.format(xr1y2p[0], xr1y2p[1], xr1y2p[2]))
f.write('xr1y2n ({} {} {}); \n'.format(xr1y2n[0], xr1y2n[1], xr1y2n[2]))
f.write('xr1y3m ({} {} {}); \n'.format(xr1y3m[0], xr1y3m[1], xr1y3m[2]))
f.write('xr1y3p ({} {} {}); \n'.format(xr1y3p[0], xr1y3p[1], xr1y3p[2]))
f.write('xr1y3n ({} {} {}); \n'.format(xr1y3n[0], xr1y3n[1], xr1y3n[2]))
f.write('xr1y4m ({} {} {}); \n'.format(xr1y4m[0], xr1y4m[1], xr1y4m[2]))
f.write('xr1y4p ({} {} {}); \n'.format(xr1y4p[0], xr1y4p[1], xr1y4p[2]))
f.write('xr1y4n ({} {} {}); \n'.format(xr1y4n[0], xr1y4n[1], xr1y4n[2]))
f.write('xr1y5m ({} {} {}); \n'.format(xr1y5m[0], xr1y5m[1], xr1y5m[2]))
f.write('xr1y5p ({} {} {}); \n'.format(xr1y5p[0], xr1y5p[1], xr1y5p[2]))
f.write('xr1y5n ({} {} {}); \n'.format(xr1y5n[0], xr1y5n[1], xr1y5n[2]))
f.write('xr2y2m ({} {} {}); \n'.format(xr2y2m[0], xr2y2m[1], xr2y2m[2]))
f.write('xr2y2p ({} {} {}); \n'.format(xr2y2p[0], xr2y2p[1], xr2y2p[2]))
f.write('xr2y2n ({} {} {}); \n'.format(xr2y2n[0], xr2y2n[1], xr2y2n[2]))
f.write('xr2y3m ({} {} {}); \n'.format(xr2y3m[0], xr2y3m[1], xr2y3m[2]))
f.write('xr2y3p ({} {} {}); \n'.format(xr2y3p[0], xr2y3p[1], xr2y3p[2]))
f.write('xr2y3n ({} {} {}); \n'.format(xr2y3n[0], xr2y3n[1], xr2y3n[2]))
f.write('xr2y4m ({} {} {}); \n'.format(xr2y4m[0], xr2y4m[1], xr2y4m[2]))
f.write('xr2y4p ({} {} {}); \n'.format(xr2y4p[0], xr2y4p[1], xr2y4p[2]))
f.write('xr2y4n ({} {} {}); \n'.format(xr2y4n[0], xr2y4n[1], xr2y4n[2]))
f.write('xr2y5m ({} {} {}); \n'.format(xr2y5m[0], xr2y5m[1], xr2y5m[2]))
f.write('xr2y5p ({} {} {}); \n'.format(xr2y5p[0], xr2y5p[1], xr2y5p[2]))
f.write('xr2y5n ({} {} {}); \n'.format(xr2y5n[0], xr2y5n[1], xr2y5n[2]))
f.write('xs1y2m ({} {} {}); \n'.format(xs1y2m[0], xs1y2m[1], xs1y2m[2]))
f.write('xs1y2p ({} {} {}); \n'.format(xs1y2p[0], xs1y2p[1], xs1y2p[2]))
f.write('xs1y2n ({} {} {}); \n'.format(xs1y2n[0], xs1y2n[1], xs1y2n[2]))
f.write('xs1y3m ({} {} {}); \n'.format(xs1y3m[0], xs1y3m[1], xs1y3m[2]))
f.write('xs1y3p ({} {} {}); \n'.format(xs1y3p[0], xs1y3p[1], xs1y3p[2]))
f.write('xs1y3n ({} {} {}); \n'.format(xs1y3n[0], xs1y3n[1], xs1y3n[2]))
f.write('xs1y4m ({} {} {}); \n'.format(xs1y4m[0], xs1y4m[1], xs1y4m[2]))
f.write('xs1y4p ({} {} {}); \n'.format(xs1y4p[0], xs1y4p[1], xs1y4p[2]))
f.write('xs1y4n ({} {} {}); \n'.format(xs1y4n[0], xs1y4n[1], xs1y4n[2]))
f.write('xs1y5m ({} {} {}); \n'.format(xs1y5m[0], xs1y5m[1], xs1y5m[2]))
f.write('xs1y5p ({} {} {}); \n'.format(xs1y5p[0], xs1y5p[1], xs1y5p[2]))
f.write('xs1y5n ({} {} {}); \n'.format(xs1y5n[0], xs1y5n[1], xs1y5n[2]))
f.write('xs2y2m ({} {} {}); \n'.format(xs2y2m[0], xs2y2m[1], xs2y2m[2]))
f.write('xs2y2p ({} {} {}); \n'.format(xs2y2p[0], xs2y2p[1], xs2y2p[2]))
f.write('xs2y2n ({} {} {}); \n'.format(xs2y2n[0], xs2y2n[1], xs2y2n[2]))
f.write('xs2y3m ({} {} {}); \n'.format(xs2y3m[0], xs2y3m[1], xs2y3m[2]))
f.write('xs2y3p ({} {} {}); \n'.format(xs2y3p[0], xs2y3p[1], xs2y3p[2]))
f.write('xs2y3n ({} {} {}); \n'.format(xs2y3n[0], xs2y3n[1], xs2y3n[2]))
f.write('xs2y4m ({} {} {}); \n'.format(xs2y4m[0], xs2y4m[1], xs2y4m[2]))
f.write('xs2y4p ({} {} {}); \n'.format(xs2y4p[0], xs2y4p[1], xs2y4p[2]))
f.write('xs2y4n ({} {} {}); \n'.format(xs2y4n[0], xs2y4n[1], xs2y4n[2]))
f.write('xs2y5m ({} {} {}); \n'.format(xs2y5m[0], xs2y5m[1], xs2y5m[2]))
f.write('xs2y5p ({} {} {}); \n'.format(xs2y5p[0], xs2y5p[1], xs2y5p[2]))
f.write('xs2y5n ({} {} {}); \n'.format(xs2y5n[0], xs2y5n[1], xs2y5n[2]))
f.write('xCy2m ({} {} {}); \n'.format(xCy2m[0], xCy2m[1], xCy2m[2]))
f.write('xCy2p ({} {} {}); \n'.format(xCy2p[0], xCy2p[1], xCy2p[2]))
f.write('xCy2n ({} {} {}); \n'.format(xCy2n[0], xCy2n[1], xCy2n[2]))
f.write('xCy3m ({} {} {}); \n'.format(xCy3m[0], xCy3m[1], xCy3m[2]))
f.write('xCy3p ({} {} {}); \n'.format(xCy3p[0], xCy3p[1], xCy3p[2]))
f.write('xCy3n ({} {} {}); \n'.format(xCy3n[0], xCy3n[1], xCy3n[2]))
f.write('xCy4m ({} {} {}); \n'.format(xCy4m[0], xCy4m[1], xCy4m[2]))
f.write('xCy4p ({} {} {}); \n'.format(xCy4p[0], xCy4p[1], xCy4p[2]))
f.write('xCy4n ({} {} {}); \n'.format(xCy4n[0], xCy4n[1], xCy4n[2]))
f.write('xCy5m ({} {} {}); \n'.format(xCy5m[0], xCy5m[1], xCy5m[2]))
f.write('xCy5p ({} {} {}); \n'.format(xCy5p[0], xCy5p[1], xCy5p[2]))
f.write('xCy5n ({} {} {}); \n'.format(xCy5n[0], xCy5n[1], xCy5n[2]))
f.write('x0y2ArcP ({} {} {}); \n'.format(x0y2ArcP[0], x0y2ArcP[1], x0y2ArcP[2]))
f.write('x0y2ArcN ({} {} {}); \n'.format(x0y2ArcN[0], x0y2ArcN[1], x0y2ArcN[2]))
f.write('blRadial {}; \n'.format(blRadial))
f.write('coreRadial {}; \n'.format(coreRadial))
# f.write('xdMidG {}; \n'.format(lowMidGrad))
f.write('xdMidG {}; \n'.format(G))
f.write('xdBG {}; \n'.format(butterflyGrad))
f.write('xdP {}; \n'.format(xdP))
f.write('xdN {}; \n'.format(xdN))
f.write('xdupG {}; \n'.format(xdupG))
f.write('xddwG {}; \n'.format(xddwG))
f.write('ptX1 {}; \n'.format(ptX1))
f.write('ptXS {}; \n'.format(ptXS))
f.write('ptX2 {}; \n'.format(ptX2))
f.write('ptX3 {}; \n'.format(ptX3))
f.write('ptX4 {}; \n'.format(ptX4))
f.write('ptX5 {}; \n'.format(ptX5))
f.write('ptX6 {}; \n'.format(ptX6))
f.write('ptX7 {}; \n'.format(ptX7))
f.write('ptX8 {}; \n'.format(ptX8))
f.write('x2P {}; \n'.format(x2P))
f.write('x3P {}; \n'.format(x3P))
f.write('x4P {}; \n'.format(x4P))
f.write('x5P {}; \n'.format(x5P))
f.write('x6P {}; \n'.format(x6P))
f.write('x7P {}; \n'.format(x7P))
f.write('x8P {}; \n'.format(x8P))
f.write('x2N {}; \n'.format(x2N))
f.write('x3N {}; \n'.format(x3N))
f.write('x4N {}; \n'.format(x4N))
f.write('x5N {}; \n'.format(x5N))
f.write('x6N {}; \n'.format(x6N))
f.write('x7N {}; \n'.format(x7N))
f.write('x8N {}; \n'.format(x8N))
f.write('x2upG {}; \n'.format(x2upG))
f.write('x3upG {}; \n'.format(x3upG))
f.write('x4upG {}; \n'.format(x4upG))
f.write('x5upG {}; \n'.format(x5upG))
f.write('x6upG {}; \n'.format(x6upG))
f.write('x7upG {}; \n'.format(x7upG))
f.write('x8upG {}; \n'.format(x8upG))
f.write('x2dwG {}; \n'.format(x2dwG))
f.write('x3dwG {}; \n'.format(x3dwG))
f.write('x4dwG {}; \n'.format(x4dwG))
f.write('x5dwG {}; \n'.format(x5dwG))
f.write('x6dwG {}; \n'.format(x6dwG))
f.write('x7dwG {}; \n'.format(x7dwG))
f.write('x8dwG {}; \n'.format(x8dwG))
f.write('radial {}; \n'.format(radial))
f.write('x1dwG {}; \n'.format(x1dwG))
f.write('x1upG {}; \n'.format(x1upG))
f.write('xSdwG {}; \n'.format(xSdwG))
f.write('xSupG {}; \n'.format(xSupG))
f.write('tangential {}; \n'.format(int(Nq*0.5)))
f.write('midRadial {}; \n'.format(int(newNbcMid)))
f.close()















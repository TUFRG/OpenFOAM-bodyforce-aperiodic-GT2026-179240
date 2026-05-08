#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 23:31:46 2026

@author: adekola
"""
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

def rMatrix(points, R):
    r11,r12,r13,r21,r22,r23,r31,r32,r33 = R
    R = np.array([
        [ r11, r12, r13],
        [ r21, r22, r23],
        [r31, r32, r33]])
    return points @ R.T 

def continuous_f(thetaEval, c, M):
    """Vectorized Fourier series evaluation"""
    thetaEval = np.asarray(thetaEval) % (2 * np.pi)
    result = c[..., 0] * np.ones_like(thetaEval) 
    for k in range(1, M + 1):
        result += c[..., k] * np.cos(k * thetaEval) + c[..., M + k] * np.sin(k * thetaEval) #Evaluate the Field value at current thetaEval
    return result

def cMatrix(Ns, nSections, oldN, interpData, M=3):
    """Vectorized Fourier coefficient computation"""
    M = int(np.floor(M))
    cMat = np.zeros([Ns, nSections, oldN, M*2+1+3]) #The cMatrix also stores the 3 coordinates point, the constant term and An and Bn terms
    for d in range(Ns):
        for e in range(nSections):
            theta = interpData[d, e, :, 0] % (2 * np.pi)
            values = interpData[d, e, :, 3]
            # Build design matrix
            A = np.ones((oldN, 2 * M + 1))
            for k in range(1, M + 1):
                A[:, k] = np.cos(k * theta)
                A[:, M + k] = np.sin(k * theta)
            c, residuals, rank, s = np.linalg.lstsq(A, values, rcond=None) #Compute using least square method. 
            cMat[d, e, :, 3:M*2+1+3] = c
            cMat[d, e, :, 0] = theta
            cMat[d, e, :, 1] = interpData[d, e, :, 1]
            cMat[d, e, :, 2] = interpData[d, e, :, 2]  
    return cMat

def changeArrayStructure(Ns, nSections, N, data):
    """Vectorized array restructuring"""
    return np.transpose(data, (2, 1, 0, 3))

def origArrayStructure(Ns, nSections, N, data):
    """Vectorized array restructuring back"""
    return np.transpose(data, (2, 1, 0, 3))

#def read_coordinates(filePath):
#    """Read and parse coordinates from file"""
#    with open(filePath, 'r') as file:
#        cell_data = file.read()
#    pattern = r'\((-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\)'
#    coordinatesExtract = re.findall(pattern, cell_data)
#    coordinates = np.array(coordinatesExtract, dtype=float)
#    rho = np.sqrt(coordinates[:, 0]**2 + coordinates[:, 1]**2)
#    theta = np.arctan2(coordinates[:, 1], coordinates[:, 0])
#    return coordinates, rho, theta

def read_coordinates(file_path, R):
    """Read and parse coordinates from file"""
    with open(file_path, 'r') as file:
        cell_data = file.read()
    pattern = r'\((-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\)'
    coordinates_extract = re.findall(pattern, cell_data)
    coordinates = np.array(coordinates_extract, dtype=float)
    abcCord = rMatrix(coordinates, R)
    rho = np.sqrt(abcCord[:, 0]**2 + abcCord[:, 1]**2)
    theta = np.arctan2(abcCord[:, 1], abcCord[:, 0])
    return abcCord, rho, theta
    
def cellID(filename):
    """Extract cell IDs as numpy array"""
    cellIds = []
    with open(filename, 'r') as f:
        for line in f:
            match = re.search(r'Cell\s+(\d+)', line)
            if match:
                cellIds.append(int(match.group(1)))
    return np.array(cellIds)

def build_global_interpolators(origData, Nb, nSections, nPoints, M, scale=1.0):
    print(f"  Building global interpolators (combining {Nb} blades)...")
    # Stack coordinates from ALL blades
    zrList = []
    cList = [[] for _ in range(2*M + 1)]  # One list per coefficient
    for blade in range(Nb):
        # Extract (z, r) coordinates for this blade
        zrBlade = np.column_stack((
            origData[blade, :, :, 2].ravel(),  # z
            origData[blade, :, :, 1].ravel()   # r
        )) * scale
        # Extract all Fourier coefficients for this blade
        cBlade = origData[blade, :, :, 3:].reshape(-1, 2*M + 1)
        zrList.append(zrBlade)
        for j in range(2*M + 1):
            cList[j].append(cBlade[:, j])
    # Combine all blades into single arrays
    zrAll = np.vstack(zrList)  # Shape: (Nb * nSections * nPoints, 2)
    # Build one interpolator per coefficient
    globalInterpolators = []
    for j in range(2*M + 1):
        cAll = np.concatenate(cList[j])  
        interp = LinearNDInterpolator(zrAll, cAll, fill_value=np.nan)
        nearestInterp = NearestNDInterpolator(zrAll, cAll)
        globalInterpolators.append((interp, nearestInterp))
    print(f"  Global interpolators built with {len(zrAll)} points")
    return globalInterpolators

def hybrid_interpolate(interp_tuple, points):
    linear, nearest = interp_tuple
    result = linear(points)
    mask = np.isnan(result)
    if np.any(mask):
        result[mask] = nearest(points)[mask]
    return result

def signWithTolerance(d2, refMag=1.0, tol=1e-3):
    """
    Sign with tolerance relative to field magnitude.
    tol is a fraction of refMag — defaults to 0.1% of field magnitude.
    """
    if np.abs(d2) < tol * refMag:
        return 0
    return np.sign(d2)

def computeDiscreteD2(xVal, yVal):
    """
    Compute second derivative at each interior discrete point
    using central difference method
    """
    xMid = []
    d2 = []
    for i in range(1, int(len(xVal)-1)):
        h1 = xVal[i] - xVal[i-1]
        h2 = xVal[i+1] - xVal[i]
        xMid.append(xVal[i])
        y = (2/(h1+h2))*(((yVal[i+1]-yVal[i])/h2) - ((yVal[i] - yVal[i-1])/h1))
        d2.append(y)
    d2 = np.array(d2)
    xMid = np.array(xMid)
    return d2, xMid

def checkInflections(xVal, yVal, xFourier, yFourier, refMag=None):
    """
    For each interval between consecutive discrete points:
    - If discrete d2 does not change sign: fitted curve should have 0 sign changes
    - If discrete d2 changes sign once: fitted curve should have exactly 1 sign change
    refMag: reference magnitude for tolerance — passed in from optimizeHarmonics
            defaults to max of discrete d2 if not provided
    """
    d2Discrete, xMid = computeDiscreteD2(xVal, yVal)
    d2Fourier, xFourierD2 = computeDiscreteD2(xFourier, yFourier)
    
    # Fall back to local magnitude if not provided
    if refMag is None:
        refMag = np.max(np.abs(d2Discrete)) if np.max(np.abs(d2Discrete)) > 0 else 1.0
    violations = 0
    for i in range(len(d2Discrete) - 1):
        xLeft  = xMid[i]
        xRight = xMid[i+1]
        s1 = signWithTolerance(d2Discrete[i],   refMag=refMag)
        s2 = signWithTolerance(d2Discrete[i+1], refMag=refMag)
        if s1 == 0 or s2 == 0:
            continue
        expectedSignChanges = 0 if s1 == s2 else 1
        # Find fitted d2 points within this interval
        mask = (xFourierD2 >= xLeft) & (xFourierD2 <= xRight)
        d2Interval = d2Fourier[mask]
        if len(d2Interval) < 2:
            continue
        # Filter near-zero values before counting sign changes
        nonZero = d2Interval[np.abs(d2Interval) > 1e-3 * refMag]
        if len(nonZero) < 2:
            continue
        signChanges = np.sum(np.diff(np.sign(nonZero)) != 0)
        if signChanges > expectedSignChanges:
            violations += signChanges - expectedSignChanges
    return violations

def optimizeHarmonics(field, Nr, rSections, Nb, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale):
    mMax = int(np.floor(0.5*(Nb-1)))
    bestM = 2
    fieldCamber = changeArrayStructure(Nr, rSections, Nb, field)
    # Check if field is uniform across ALL sections and blades
    allDerivates = []
    for sec in range(fieldCamber.shape[1]): #axis 1 is the number of sections 
        for blade in range(fieldCamber.shape[2]): #axis 2 is the number of blades
            origVal = fieldCamber[:,sec, blade] #axis 0 is the number of points
            discrete_x = origVal[:, 0] % (2 * np.pi)
            discrete_y = origVal[:, 3]
            derivate, _ = computeDiscreteD2(discrete_x, discrete_y)
            allDerivates.append(np.mean(np.abs(derivate)))
    if np.max(allDerivates) < 1e-06:
        print(f"M={mMax} (uniform field)")
        return mMax
    # Find the most conservative M — lowest M with zero violations across ALL points
    for M in range(mMax, 1, -1):
        newFieldCamber = cMatrix(Nr, rSections, Nb, fieldCamber, M=M)
        angle = np.linspace(0, 2*np.pi, 1000) #Evaluate the signal at 1000 points for each z,r location
        num = np.linspace(1, M, M)
        Nl = int(2*M+1) + 3 
        totalViolations = 0
        for sec in range(fieldCamber.shape[1]): #axis 1 is the number of sections 
            for blade in range(fieldCamber.shape[2]): #axis 2 is the number of blades
                origVal = fieldCamber[:,sec, blade] #axis 0 is the number of points
                discrete_x = origVal[:, 0] % (2 * np.pi)
                discrete_y = origVal[:, 3]
                fieldData = newFieldCamber[:,sec, blade]           
                Fval = np.zeros(1000)
                for a in range(1000):
                    Fval[a] = (fieldData[0, 3] +
                              sum(fieldData[0, 4:4+M] * np.cos(angle[a] * num)) +
                              sum(fieldData[0, 4+M:Nl] * np.sin(angle[a] * num)))      
                refMag = np.max(np.abs(discrete_y)) if np.max(np.abs(discrete_y)) > 0 else 1.0
                violations = checkInflections(discrete_x, discrete_y, angle, Fval, refMag=refMag)
                totalViolations += violations 
        print(f"M={M:3d}  total violations across all sections/blades={totalViolations}")
        if totalViolations == 0:
            return M
    return bestM
    
def generateGridFields(field, Nr, rSections, Nb, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, filePath, fileName, plot=1):
    fieldCamber = changeArrayStructure(Nr, rSections, Nb, field)
    M = optimizeHarmonics(field, Nr, rSections, Nb, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale)
    newFieldCamber = cMatrix(Nr, rSections, Nb, fieldCamber, M=M)
    origFieldCamber = origArrayStructure(Nr, rSections, int(2*M+1), newFieldCamber)
    globalFieldInterp = build_global_interpolators(
        origFieldCamber, Nb, rSections, Nr, M, scale)
    zrCFDField = np.column_stack((rGridPoints[:, 2], rGridPoints[:, 1]))
    cCFDField = np.array([hybrid_interpolate(interp, zrCFDField) 
                            for interp in globalFieldInterp]).T
    fField = np.column_stack([rGrid[2], rGrid[1], rGrid[0][:, 2], 
                                continuous_f(rGrid[2], cCFDField, M)])
    if plot == 0:
        val = plotFourierFitting(fieldCamber, newFieldCamber, M, filePath, fileName)
    return fField

def writeFields(field, filePath, fieldName):
    fr = open(filePath + '/{}'.format(fieldName), 'w')
    fr.write('( \n')
    for j in range(NGrid):
        fr.write('{} \n'.format(field[j]))
    fr.write(') \n')
    fr.close() 
    return print("\n" + "="*70), print("Done! {} written successfully.".format(fieldName)), print("="*70)
def plotFourierFitting (origField, fittedField, M, filePath, fileName, N=1000):
    Nl = int(2*M+1)+3
    origVal = origField[-2,4]
    fieldData = fittedField[-2,4]
    Fval = np.zeros(N)
    num = np.linspace(1,M,M)
    angle = np.linspace(0,2*np.pi, N)
    for a in range(1000):
        Fval[a] = fieldData[0,3] + sum(fieldData[0,4:4+M]*np.cos(angle[a]*num)) + sum(fieldData[0,4+M:Nl]*np.sin(angle[a]*num)) 
    plt.figure()
    plt.plot(origVal[:,0]%(2 * np.pi), origVal[:,3], 'k.', label='geometric field value')
    plt.plot(angle, Fval, 'r', label='fitted Fourier Coefficient')
    plt.ylabel('{}'.format(fileName))
    plt.xlabel('Circumferential direction')
    plt.grid()
    plt.legend()
    # yVal = 0.5*(max(Fval) + min(Fval))
    # plt.text(3, yVal, 'M = {}'.format(M))    
    # Evaluate fitted curve at same points
    y_fitted = np.array([
        fieldData[0, 3] + 
        sum(fieldData[0, 4:4+M] * np.cos(t * num)) + 
        sum(fieldData[0, 4+M:Nl] * np.sin(t * num))
        for t in origVal[:,0]%(2 * np.pi)
    ])
    # RMS error
    rms = np.sqrt(np.mean((origVal[:,3] - y_fitted)**2))
    # Normalised RMS (as % of mean value)
    nrms = rms / np.mean(np.abs(origVal[:,3])) * 100
    # plt.text(3,newYval, f"Normalised RMS: {nrms:.2f}%")
    ax = plt.gca()
    ax.text(0.05, 0.95, 'M = {}'.format(M), transform=ax.transAxes, 
            verticalalignment='top')
    ax.text(0.05, 0.88, f"Normalised RMS: {nrms:.2f}%", transform=ax.transAxes,
            verticalalignment='top')
    ax.ticklabel_format(style='plain', useOffset=False)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.8f'))
    # print(f"Normalised RMS: {nrms:.2f}%")
    plt.savefig(filePath + '/' + fileName + '.png', dpi=500, bbox_inches='tight')
    plt.close()
    return Fval#,print('Fitted Data Plotted')

#%% Geometry parameters
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
    "scale", 
    type=float,
    help="scales"
)
parser.add_argument(
    "periodicOrAperiodic", 
    type=int,
    help="periodic (0) or nonPeriodic (1)"
)
parser.add_argument(
    "stage", 
    type=int,
    help="rotorAlone (0), statorAlone (1) or stage (2)"
)
parser.add_argument(
    "--RmatrixCoeff",
    type=float,
    nargs=9,
    default=[1,0,0,0,1,0,0,0,1],
    help="9 coefficients for rotation matrix (r11 r12 r13 r21 r22 r23 r31 r32 r33)"
)
parser.add_argument(
    "casePath", 
    type=str,
    help="The path to the caseFolder where the simulation will be run"
)
args = parser.parse_args()

# Use parsed arguments
Nbr = args.Nbr  # Total number of blades in the annulus for rotor
Nbs = args.Nbs  # Total number of blades in the annulus for stator
rSections = args.rSections  # Number of rotor blade profiles
sSections = args.sSections  # Number of stator blade profiles
scale = args.scale
casePath = args.casePath
periodicOrAperiodic = args.periodicOrAperiodic  # if periodic select 0 otherwise select 1
stage = args.stage
RmatrixCoeff = args.RmatrixCoeff  # Define based on your machine axis

#%%% Load all input Data
# filePath = '/home/adekola/Documents/New_PhD/RUN/bodyForceECL5/boundaryConditionTest/temp/'
Nr = 101
Ns = 101
NrB = Nr#201
NsB = Ns#241
Np = 10
filePath = '/{}/'.format(casePath)#'/fullWheel/'
if periodicOrAperiodic == 0:
    inputPath = '../inputData/periodic/'
    plotPath = '../plots/periodic/'
    if not os.path.exists(plotPath):
        os.makedirs(plotPath)
else:
    inputPath = '../inputData/nonPeriodic/'
    plotPath = '../plots/nonPeriodic/'
    if not os.path.exists(plotPath):
        os.makedirs(plotPath)    

#%%
if stage == 0:
    print("Loading grid data...")
    rGrid = read_coordinates(filePath + 'cellCoordinates_Rotor.txt', RmatrixCoeff)
    NrGrid = len(rGrid[1])
    rGridPoints = np.column_stack([rGrid[2], rGrid[1], rGrid[0][:, 2]])
    rCellID = cellID(filePath + 'cellCoordinates_Rotor.txt')
    grid = read_coordinates(filePath + 'cellCoordinates.txt', RmatrixCoeff)
    NGrid = len(grid[1])
    gridPoints = np.column_stack([grid[2], grid[1], grid[0][:, 2]])
    gridCellID = cellID(filePath + 'cellCoordinates.txt')
    
    print("Loading Blade Geometry data ...")
    rotorNr = np.loadtxt(inputPath + './rotorNr.txt', delimiter=',').reshape(Nbr, rSections, Nr+2, 4)
    rotorNc = np.loadtxt(inputPath + './rotorNc.txt', delimiter=',').reshape(Nbr, rSections, Nr+2, 4)
    rotorNth = np.loadtxt(inputPath + './rotorNth.txt', delimiter=',').reshape(Nbr, rSections, Nr+2, 4)
    xc_rotor = np.loadtxt(inputPath + './xc_rotor.txt', delimiter=',')
    xc_rotor[:, 3] *= scale
    xc_rotor = xc_rotor.reshape(Nbr, rSections, Nr, 4)
    Rt_hat = np.loadtxt(inputPath + './Rt_hat.txt', delimiter=',').reshape(Nbr, (rSections-1)*Np+1, NrB, 4)
    rotorBeta = np.loadtxt(inputPath + './rotorBeta.txt', delimiter=',').reshape(Nbr, rSections, Nr, 4)
elif stage == 1:
    print("Loading grid data...")
    sGrid = read_coordinates(filePath + 'cellCoordinates_Stator.txt', RmatrixCoeff)
    NsGrid = len(sGrid[1])
    sGridPoints = np.column_stack([sGrid[2], sGrid[1], sGrid[0][:, 2]])
    sCellID = cellID(filePath + 'cellCoordinates_Stator.txt')
    grid = read_coordinates(filePath + 'cellCoordinates.txt', RmatrixCoeff)
    NGrid = len(grid[1])
    gridPoints = np.column_stack([grid[2], grid[1], grid[0][:, 2]])
    gridCellID = cellID(filePath + 'cellCoordinates.txt')
    
    print("Loading blade geometry data...")
    statorNr = np.loadtxt(inputPath + './statorNr.txt', delimiter=',').reshape(Nbs, sSections, Ns+2, 4)
    statorNc = np.loadtxt(inputPath + './statorNc.txt', delimiter=',').reshape(Nbs, sSections, Ns+2, 4)
    statorNth = np.loadtxt(inputPath + './statorNth.txt', delimiter=',').reshape(Nbs, sSections, Ns+2, 4)
    xc_stator = np.loadtxt(inputPath + './xc_stator.txt', delimiter=',')
    xc_stator[:, 3] *= scale
    xc_stator = xc_stator.reshape(Nbs, sSections, Ns, 4)
    St_hat = np.loadtxt(inputPath + './St_hat.txt', delimiter=',').reshape(Nbs, (sSections-1)*Np+1, NsB, 4)
    statorBeta = np.loadtxt(inputPath + './statorBeta.txt', delimiter=',').reshape(Nbs, sSections, Ns, 4)
else:
    print("Loading grid data...")
    rGrid = read_coordinates(filePath + 'cellCoordinates_Rotor.txt', RmatrixCoeff)
    NrGrid = len(rGrid[1])
    rGridPoints = np.column_stack([rGrid[2], rGrid[1], rGrid[0][:, 2]])
    rCellID = cellID(filePath + 'cellCoordinates_Rotor.txt')
    sGrid = read_coordinates(filePath + 'cellCoordinates_Stator.txt', RmatrixCoeff)
    NsGrid = len(sGrid[1])
    sGridPoints = np.column_stack([sGrid[2], sGrid[1], sGrid[0][:, 2]])
    sCellID = cellID(filePath + 'cellCoordinates_Stator.txt')
    grid = read_coordinates(filePath + 'cellCoordinates.txt', RmatrixCoeff)
    NGrid = len(grid[1])
    gridPoints = np.column_stack([grid[2], grid[1], grid[0][:, 2]])
    gridCellID = cellID(filePath + 'cellCoordinates.txt')    
    
    print("Loading blade geometry data...")
    # Load and reshape blade data
    rotorNr = np.loadtxt(inputPath + './rotorNr.txt', delimiter=',').reshape(Nbr, rSections, Nr+2, 4)
    rotorNc = np.loadtxt(inputPath + './rotorNc.txt', delimiter=',').reshape(Nbr, rSections, Nr+2, 4)
    rotorNth = np.loadtxt(inputPath + './rotorNth.txt', delimiter=',').reshape(Nbr, rSections, Nr+2, 4)
    statorNr = np.loadtxt(inputPath + './statorNr.txt', delimiter=',').reshape(Nbs, sSections, Ns+2, 4)
    statorNc = np.loadtxt(inputPath + './statorNc.txt', delimiter=',').reshape(Nbs, sSections, Ns+2, 4)
    statorNth = np.loadtxt(inputPath + './statorNth.txt', delimiter=',').reshape(Nbs, sSections, Ns+2, 4)
    
    xc_rotor = np.loadtxt(inputPath + './xc_rotor.txt', delimiter=',')
    xc_rotor[:, 3] *= scale
    xc_rotor = xc_rotor.reshape(Nbr, rSections, Nr, 4)
    
    xc_stator = np.loadtxt(inputPath + './xc_stator.txt', delimiter=',')
    xc_stator[:, 3] *= scale
    xc_stator = xc_stator.reshape(Nbs, sSections, Ns, 4)
    
    St_hat = np.loadtxt(inputPath + './St_hat.txt', delimiter=',').reshape(Nbs, (sSections-1)*Np+1, NsB, 4)
    Rt_hat = np.loadtxt(inputPath + './Rt_hat.txt', delimiter=',').reshape(Nbr, (rSections-1)*Np+1, NrB, 4)
    
    rotorBeta = np.loadtxt(inputPath + './rotorBeta.txt', delimiter=',').reshape(Nbr, rSections, Nr, 4)
    statorBeta = np.loadtxt(inputPath + './statorBeta.txt', delimiter=',').reshape(Nbs, sSections, Ns, 4)

#%% Generate the fields
if stage == 0:
    rNr = generateGridFields(rotorNr, Nr+2, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_Nr', plot=0)
    rNc = generateGridFields(rotorNc, Nr+2, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale,  plotPath, 'rotor_Nc', plot=0)
    rNth = generateGridFields(rotorNth, Nr+2, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale,  plotPath, 'rotor_Nth', plot=0)
    rXc = generateGridFields(xc_rotor, Nr, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_arclength', plot=0)
    rSt = generateGridFields(Rt_hat, Nr, (rSections-1)*Np+1, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_thickness', plot=0)
    rBeta = generateGridFields(rotorBeta, Nr, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_pitch', plot=0)
    Nnr = np.zeros(NGrid)
    Nnc = np.zeros(NGrid)
    Nnth = np.zeros(NGrid)
    xC = np.zeros(NGrid)
    t_hat = np.zeros(NGrid)
    pitch = np.zeros(NGrid)
    rIndx = np.searchsorted(gridCellID, rCellID)
    rMask = (rIndx < len(gridCellID)) & (gridCellID[rIndx] == rCellID)
    Nnr[rIndx[rMask]] = rNr[rMask, 3]
    Nnc[rIndx[rMask]] = rNc[rMask, 3]
    Nnth[rIndx[rMask]] = rNth[rMask, 3]
    xC[rIndx[rMask]] = rXc[rMask, 3]
    t_hat[rIndx[rMask]] = rSt[rMask, 3]
    pitch[rIndx[rMask]] = rBeta[rMask, 3]
    writeNr = writeFields(Nnr, filePath, 'Nr_data')
    writeNc = writeFields(Nnc, filePath, 'Nz_data')
    writeNth = writeFields(Nnth, filePath, 'Nth_data')
    writeXc = writeFields(xC, filePath, 'xC_data')
    writet_hat = writeFields(t_hat, filePath, 't_hat')
    writePitch = writeFields(pitch, filePath, 'beta')
elif stage == 1:
    sNr = generateGridFields(statorNr, Ns+2, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_Nr', plot=0)
    sNc = generateGridFields(statorNc, Ns+2, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_Nc', plot=0)
    sNth = generateGridFields(statorNth, Ns+2, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_Nth', plot=0)
    sXc = generateGridFields(xc_stator, Ns, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_arclength', plot=0)
    sSt = generateGridFields(St_hat, Ns, (sSections-1)*Np+1, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale, plotPath, 'stator_thickness', plot=0)
    sBeta = generateGridFields(statorBeta, Ns, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_pitch', plot=0) 
    Nnr = np.zeros(NGrid)
    Nnc = np.zeros(NGrid)
    Nnth = np.zeros(NGrid)
    xC = np.zeros(NGrid)
    t_hat = np.zeros(NGrid)
    pitch = np.zeros(NGrid)
    sIndx = np.searchsorted(gridCellID, sCellID)
    sMask = (sIndx < len(gridCellID)) & (gridCellID[sIndx] == sCellID)
    Nnr[sIndx[sMask]] = sNr[sMask, 3]
    Nnc[sIndx[sMask]] = sNc[sMask, 3]
    Nnth[sIndx[sMask]] = sNth[sMask, 3]
    xC[sIndx[sMask]] = sXc[sMask, 3]
    t_hat[sIndx[sMask]] = sSt[sMask, 3]
    pitch[sIndx[sMask]] = sBeta[sMask, 3]
    writeNr = writeFields(Nnr, filePath, 'Nr_data')
    writeNc = writeFields(Nnc, filePath, 'Nz_data')
    writeNth = writeFields(Nnth, filePath, 'Nth_data')
    writeXc = writeFields(xC, filePath, 'xC_data')
    writet_hat = writeFields(t_hat, filePath, 't_hat')
    writePitch = writeFields(pitch, filePath, 'beta')   
else:    
    rNr = generateGridFields(rotorNr, Nr+2, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_Nr', plot=0)
    rNc = generateGridFields(rotorNc, Nr+2, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale,  plotPath, 'rotor_Nc', plot=0)
    rNth = generateGridFields(rotorNth, Nr+2, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale,  plotPath, 'rotor_Nth', plot=0)
    rXc = generateGridFields(xc_rotor, Nr, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_arclength', plot=0)
    rSt = generateGridFields(Rt_hat, Nr, (rSections-1)*Np+1, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_thickness', plot=0)
    rBeta = generateGridFields(rotorBeta, Nr, rSections, Nbr, rGridPoints, rGrid, NGrid, gridCellID, rCellID, scale, plotPath, 'rotor_pitch', plot=0)

    sNr = generateGridFields(statorNr, Ns+2, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_Nr', plot=0)
    sNc = generateGridFields(statorNc, Ns+2, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_Nc', plot=0)
    sNth = generateGridFields(statorNth, Ns+2, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_Nth', plot=0)
    sXc = generateGridFields(xc_stator, Ns, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_arclength', plot=0)
    sSt = generateGridFields(St_hat, Ns, (sSections-1)*Np+1, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale, plotPath, 'stator_thickness', plot=0)
    sBeta = generateGridFields(statorBeta, Ns, sSections, Nbs, sGridPoints, sGrid, NGrid, gridCellID, sCellID, scale,  plotPath, 'stator_pitch', plot=0)  
    
    Nnr = np.zeros(NGrid)
    Nnc = np.zeros(NGrid)
    Nnth = np.zeros(NGrid)
    xC = np.zeros(NGrid)
    t_hat = np.zeros(NGrid)
    pitch = np.zeros(NGrid)
    rIndx = np.searchsorted(gridCellID, rCellID)
    sIndx = np.searchsorted(gridCellID, sCellID)
    rMask = (rIndx < len(gridCellID)) & (gridCellID[rIndx] == rCellID)
    Nnr[rIndx[rMask]] = rNr[rMask, 3]
    Nnc[rIndx[rMask]] = rNc[rMask, 3]
    Nnth[rIndx[rMask]] = rNth[rMask, 3]
    xC[rIndx[rMask]] = rXc[rMask, 3]
    t_hat[rIndx[rMask]] = rSt[rMask, 3]
    pitch[rIndx[rMask]] = rBeta[rMask, 3]
    sMask = (sIndx < len(gridCellID)) & (gridCellID[sIndx] == sCellID)
    Nnr[sIndx[sMask]] = sNr[sMask, 3]
    Nnc[sIndx[sMask]] = sNc[sMask, 3]
    Nnth[sIndx[sMask]] = sNth[sMask, 3]
    xC[sIndx[sMask]] = sXc[sMask, 3]
    t_hat[sIndx[sMask]] = sSt[sMask, 3]
    pitch[sIndx[sMask]] = sBeta[sMask, 3]    
    writeNr = writeFields(Nnr, filePath, 'Nr_data')
    writeNc = writeFields(Nnc, filePath, 'Nz_data')
    writeNth = writeFields(Nnth, filePath, 'Nth_data')
    writeXc = writeFields(xC, filePath, 'xC_data')
    writet_hat = writeFields(t_hat, filePath, 't_hat')
    writePitch = writeFields(pitch, filePath, 'beta')









































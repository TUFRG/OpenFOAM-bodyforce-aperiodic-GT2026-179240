#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 22:23:11 2025

@author: adekola
"""

"""
This script is used for all the bodyforce cases. Instead of downloading the entire postProcessing folder, this file will extract the latest time information for further postProcessing
"""
import os
import glob
import pickle
import shutil
import numpy as np

#%% Names of folders in the postProcessing folder
inletMassFlow = 'flowRatePatch(name=inlet)'
outletMassFlow = 'flowRatePatch(name=outlet)'
outletData = 'Outlet'
massAvgInlet = 'massAvgInletpt'
massAvgRotIn = 'rotorInlet/'
massAvgRotOut = 'rotorOutlet/'
massAvgOutlet = 'massAvgOutletpt'
linesample = 'linesample'


filePath = './'
#%% Extract data for the inletMassFlow
allTimeFolderInletMdot = (glob.glob(os.path.join(filePath + 'postProcessing/'+inletMassFlow, '[0-9]*')))
timeFolderInletMdot = max(allTimeFolderInletMdot, key=lambda x: int(os.path.basename(x)))
numOfInletMdotFiles = (glob.glob(timeFolderInletMdot + '/*.dat'))
if len(numOfInletMdotFiles) > 1:
    inletMdot = os.path.join(timeFolderInletMdot + '/surfaceFieldValue_{}.dat'.format(int(os.path.basename(timeFolderInletMdot))))
else:
    inletMdot = os.path.join(timeFolderInletMdot + '/surfaceFieldValue.dat')
inletMassFlowRaw = np.loadtxt(inletMdot)
if inletMassFlowRaw.size < 3:
    surfaceInletMassFlow = {'iteration': inletMassFlowRaw[0], 'massFlow': inletMassFlowRaw[1]}
else:
    surfaceInletMassFlow = {'iteration': inletMassFlowRaw[:,0], 'massFlow': inletMassFlowRaw[:,1]}
#%% Extract data for the outletMassFlow
allTimeFolderOutletMdot = (glob.glob(os.path.join(filePath + 'postProcessing/'+outletMassFlow, '[0-9]*')))
timeFolderOutletMdot = max(allTimeFolderOutletMdot, key=lambda x: int(os.path.basename(x)))
numOfOutletMdotFiles = (glob.glob(timeFolderOutletMdot + '/*.dat'))
if len(numOfOutletMdotFiles) > 1:
    outletMdot = os.path.join(timeFolderOutletMdot + '/surfaceFieldValue_{}.dat'.format(int(os.path.basename(timeFolderOutletMdot))))
else:
    outletMdot = os.path.join(timeFolderOutletMdot + '/surfaceFieldValue.dat')
    
outletMassFlowRaw = np.loadtxt(outletMdot)
if outletMassFlowRaw.size < 3:
    surfaceOutletMassFlow = {'iteration': outletMassFlowRaw[0], 'massFlow': outletMassFlowRaw[1]}
else:
    surfaceOutletMassFlow = {'iteration': outletMassFlowRaw[:,0], 'massFlow': outletMassFlowRaw[:,1]}
#%% Extract data for the massAvgInlet Data
#Take note, you will have to copy the surface folder files and the timefolder files
dtypes = [('col1', 'float'),('col2', 'float'),('col3', 'float'), ('col4', 'U50'),('col5', 'float'),('col6', 'U50'),('col7', 'float'),('col8', 'float'),
          ('col9', 'float'),('col10', 'float'),('col11', 'U50'),('col12', 'float'),('col13', 'U50')]

allTimeFolderMassAvgInletpt = (glob.glob(os.path.join(filePath +  '/postProcessing/' + massAvgInlet, '[0-9]*')))
surfaceFolderMassAvgInletpt = (glob.glob(os.path.join(filePath +  '/postProcessing/' + massAvgInlet + '/surface/' +  '[0-9]*')))
timeFolderMassAvgeInlet = max(allTimeFolderMassAvgInletpt, key=lambda x: int(os.path.basename(x)))
timeSurfaceFolderMassAvgInlet = max(surfaceFolderMassAvgInletpt, key=lambda x:int(os.path.basename(x)))
numOfMassAvgInletFiles = (glob.glob(timeFolderMassAvgeInlet + '/*.dat'))
if len(numOfMassAvgInletFiles) > 1:
    massAvgInletpt = np.loadtxt(timeFolderMassAvgeInlet + '/surfaceFieldValue_{}.dat'.format(int(os.path.basename(timeFolderMassAvgeInlet))), dtype=dtypes)
else:
    massAvgInletpt = np.loadtxt(timeFolderMassAvgeInlet + '/surfaceFieldValue.dat', dtype=dtypes)
if massAvgInletpt.size == 1:
    inletMassAvg = np.zeros(len(dtypes))
    inletMassAvg[0], inletMassAvg[1], inletMassAvg[2] = massAvgInletpt.item()[0], massAvgInletpt.item()[1], massAvgInletpt.item()[2]
    inletMassAvg[3] = float(massAvgInletpt.item()[3].replace('(', ''))
    inletMassAvg[4] = massAvgInletpt.item()[4]
    inletMassAvg[5] = float(massAvgInletpt.item()[5].replace(')', ''))
    inletMassAvg[6], inletMassAvg[7], inletMassAvg[8], inletMassAvg[9] = massAvgInletpt.item()[6], massAvgInletpt.item()[7], massAvgInletpt.item()[8], massAvgInletpt.item()[9]
    inletMassAvg[10] = float(massAvgInletpt.item()[10].replace('(', ''))
    inletMassAvg[11] = massAvgInletpt.item()[11]
    inletMassAvg[12] = float(massAvgInletpt.item()[12].replace(')', ''))    
    massAvgInletSurface = {'iteration':inletMassAvg[0], 'p': inletMassAvg[1], 'totalP': inletMassAvg[2], 'Ux': inletMassAvg[3], 'Uy': inletMassAvg[4], 'Uz': inletMassAvg[5],
                           'alphaFlow': inletMassAvg[6], 'radialFlow': inletMassAvg[7], 'Urel': inletMassAvg[8], 'relFlow': inletMassAvg[9], 'Uradial': inletMassAvg[10], 
                           'Utangent': inletMassAvg[11], 'Uaxial': inletMassAvg[12]}

else:
    inletMassAvg = np.zeros([len(massAvgInletpt), len(dtypes)])
    index = 0
    for field in massAvgInletpt.dtype.names:
        # Handle string fields that contain parentheses
        if field in ['col4', 'col6', 'col11', 'col13']:  # These are U50 string types
            if field == 'col4':  # Has opening parenthesis
                col_data = np.array([float(str(b).replace('(', '')) for b in massAvgInletpt[field]])
            elif field == 'col6':  # Has closing parenthesis
                col_data = np.array([float(str(b).replace(')', '')) for b in massAvgInletpt[field]])
            elif field == 'col11':  # Has opening parenthesis
                col_data = np.array([float(str(b).replace('(', '')) for b in massAvgInletpt[field]])
            elif field == 'col13':  # Has closing parenthesis
                col_data = np.array([float(str(b).replace(')', '')) for b in massAvgInletpt[field]])
            inletMassAvg[:, index] = col_data
        else:
            # These are already floats, just assign directly
            inletMassAvg[:, index] = massAvgInletpt[field]
        
        index += 1
    
    massAvgInletSurface = {
        'iteration': inletMassAvg[:, 0], 
        'p': inletMassAvg[:, 1], 
        'totalP': inletMassAvg[:, 2], 
        'Ux': inletMassAvg[:, 3], 
        'Uy': inletMassAvg[:, 4], 
        'Uz': inletMassAvg[:, 5],
        'alphaFlow': inletMassAvg[:, 6], 
        'radialFlow': inletMassAvg[:, 7], 
        'Urel': inletMassAvg[:, 8], 
        'relFlow': inletMassAvg[:, 9], 
        'Uradial': inletMassAvg[:, 10], 
        'Utangent': inletMassAvg[:, 11], 
        'Uaxial': inletMassAvg[:, 12]
    }    

surfaceMassAvgInlet = {'x': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/alphaFlow_patch_inlet.raw')[:,0]),
                       'y': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/alphaFlow_patch_inlet.raw')[:,1]),
                       'z': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/alphaFlow_patch_inlet.raw')[:,2]),
                       'alphaFlow': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/alphaFlow_patch_inlet.raw')[:,3]),
                       'totalP':(np.loadtxt(timeSurfaceFolderMassAvgInlet + '/totalP_patch_inlet.raw')[:,3]),
                       'p': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/p_patch_inlet.raw')[:,3]),
                       'radialFlow': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/radialFlow_patch_inlet.raw')[:,3]),
                       'relFlow': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/relFlow_patch_inlet.raw')[:,3]),
                       'Urel': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/Urel_patch_inlet.raw')[:,3]),
                       'Uradial': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/U:Transformed_patch_inlet.raw')[:,3]),
                       'Utangent': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/U:Transformed_patch_inlet.raw')[:,4]),
                       'Uaxial': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/U:Transformed_patch_inlet.raw')[:,5]),
                       'Ux': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/U_patch_inlet.raw')[:,3]),
                       'Uy': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/U_patch_inlet.raw')[:,4]),
                       'Uz': (np.loadtxt(timeSurfaceFolderMassAvgInlet + '/U_patch_inlet.raw')[:,5])}

#%%Extract data for the massAvgRotIn Data
#Take note, you will have to copy the surface folder files and the timefolder files



surfaceFolderMassAvgRotInpt = (glob.glob(os.path.join(filePath +  '/postProcessing/' + massAvgRotIn  +  '[0-9]*')))
timeSurfaceFolderMassAvgRotIn = max(surfaceFolderMassAvgRotInpt, key=lambda x:int(os.path.basename(x)))

surfaceMassAvgRotIn = {'x': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/alphaFlow_zNormal.raw')[:,0]),
                       'y': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/alphaFlow_zNormal.raw')[:,1]),
                       'z': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/alphaFlow_zNormal.raw')[:,2]),
                       'alphaFlow': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/alphaFlow_zNormal.raw')[:,3]),
                       'totalP':(np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/totalP_zNormal.raw')[:,3]),
                       'p': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/p_zNormal.raw')[:,3]),
                       'radialFlow': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/radialFlow_zNormal.raw')[:,3]),
                       'relFlow': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/relFlow_zNormal.raw')[:,3]),
                       'Urel': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/Urel_zNormal.raw')[:,3]),
                       'Uradial': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/U:Transformed_zNormal.raw')[:,3]),
                       'Utangent': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/U:Transformed_zNormal.raw')[:,4]),
                       'Uaxial': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/U:Transformed_zNormal.raw')[:,5]),
                       'Ux': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/U_zNormal.raw')[:,3]),
                       'Uy': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/U_zNormal.raw')[:,4]),
                       'Uz': (np.loadtxt(timeSurfaceFolderMassAvgRotIn + '/U_zNormal.raw')[:,5])}

#%% Extract data for the massAvgRotOut Data
#Take note, you will have to copy the surface folder files and the timefolder files



surfaceFolderMassAvgRotOutpt = (glob.glob(os.path.join(filePath +  '/postProcessing/' + massAvgRotOut  +  '[0-9]*')))
timeSurfaceFolderMassAvgRotOut = max(surfaceFolderMassAvgRotOutpt, key=lambda x:int(os.path.basename(x)))


surfaceMassAvgRotOut = {'x': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/alphaFlow_zNormal.raw')[:,0]),
                       'y': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/alphaFlow_zNormal.raw')[:,1]),
                       'z': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/alphaFlow_zNormal.raw')[:,2]),
                       'alphaFlow': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/alphaFlow_zNormal.raw')[:,3]),
                       'totalP':(np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/totalP_zNormal.raw')[:,3]),
                       'p': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/p_zNormal.raw')[:,3]),
                       'radialFlow': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/radialFlow_zNormal.raw')[:,3]),
                       'relFlow': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/relFlow_zNormal.raw')[:,3]),
                       'Urel': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/Urel_zNormal.raw')[:,3]),
                       'Uradial': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/U:Transformed_zNormal.raw')[:,3]),
                       'Utangent': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/U:Transformed_zNormal.raw')[:,4]),
                       'Uaxial': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/U:Transformed_zNormal.raw')[:,5]),
                       'Ux': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/U_zNormal.raw')[:,3]),
                       'Uy': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/U_zNormal.raw')[:,4]),
                       'Uz': (np.loadtxt(timeSurfaceFolderMassAvgRotOut + '/U_zNormal.raw')[:,5])}


#%% Extract data for the massAvgOutlet Data
#Take note, you will have to copy the surface folder files and the timefolder files
dtypes = [('col1', 'float'),('col2', 'float'),('col3', 'float'), ('col4', 'U50'),('col5', 'float'),('col6', 'U50'),('col7', 'float'),('col8', 'float'),
          ('col9', 'float'),('col10', 'float'),('col11', 'U50'),('col12', 'float'),('col13', 'U50')]

allTimeFolderMassAvgOutletpt = (glob.glob(os.path.join(filePath +  '/postProcessing/' + massAvgOutlet, '[0-9]*')))
surfaceFolderMassAvgOutletpt = (glob.glob(os.path.join(filePath +  '/postProcessing/' + massAvgOutlet + '/surface/' +  '[0-9]*')))
timeFolderMassAvgeOutlet = max(allTimeFolderMassAvgOutletpt, key=lambda x: int(os.path.basename(x)))
timeSurfaceFolderMassAvgOutlet = max(surfaceFolderMassAvgOutletpt, key=lambda x:int(os.path.basename(x)))
numOfMassAvgOutletFiles = (glob.glob(timeFolderMassAvgeOutlet + '/*.dat'))
if len(numOfMassAvgOutletFiles) > 1:
    massAvgOutletpt = np.loadtxt(timeFolderMassAvgeOutlet + '/surfaceFieldValue_{}.dat'.format(int(os.path.basename(timeFolderMassAvgeOutlet))), dtype=dtypes)
else:
    massAvgOutletpt = np.loadtxt(timeFolderMassAvgeOutlet + '/surfaceFieldValue.dat', dtype=dtypes)
if massAvgOutletpt.size == 1:
    outletMassAvg = np.zeros(len(dtypes))
    outletMassAvg[0], outletMassAvg[1], outletMassAvg[2] = massAvgOutletpt.item()[0], massAvgOutletpt.item()[1], massAvgOutletpt.item()[2]
    outletMassAvg[3] = float(massAvgOutletpt.item()[3].replace('(', ''))
    outletMassAvg[4] = massAvgOutletpt.item()[4]
    outletMassAvg[5] = float(massAvgOutletpt.item()[5].replace(')', ''))
    outletMassAvg[6], outletMassAvg[7], outletMassAvg[8], outletMassAvg[9] = massAvgOutletpt.item()[6], massAvgOutletpt.item()[7], massAvgOutletpt.item()[8], massAvgOutletpt.item()[9]
    outletMassAvg[10] = float(massAvgOutletpt.item()[10].replace('(', ''))
    outletMassAvg[11] = massAvgOutletpt.item()[11]
    outletMassAvg[12] = float(massAvgOutletpt.item()[12].replace(')', ''))    
    massAvgOutletSurface = {'iteration':inletMassAvg[0], 'p': inletMassAvg[1], 'totalP': inletMassAvg[2], 'Ux': inletMassAvg[3], 'Uy': inletMassAvg[4], 'Uz': inletMassAvg[5],
                           'alphaFlow': inletMassAvg[6], 'radialFlow': inletMassAvg[7], 'Urel': inletMassAvg[8], 'relFlow': inletMassAvg[9], 'Uradial': inletMassAvg[10], 
                           'Utangent': inletMassAvg[11], 'Uaxial': inletMassAvg[12]}
 
else:
    outletMassAvg = np.zeros([len(massAvgOutletpt), len(dtypes)])
    index = 0
    for field in massAvgOutletpt.dtype.names:
        # Handle string fields that contain parentheses
        if field in ['col4', 'col6', 'col11', 'col13']:
            if field == 'col4':
                col_data = np.array([float(str(b).replace('(', '')) for b in massAvgOutletpt[field]])
            elif field == 'col6':
                col_data = np.array([float(str(b).replace(')', '')) for b in massAvgOutletpt[field]])
            elif field == 'col11':
                col_data = np.array([float(str(b).replace('(', '')) for b in massAvgOutletpt[field]])
            elif field == 'col13':
                col_data = np.array([float(str(b).replace(')', '')) for b in massAvgOutletpt[field]])
            outletMassAvg[:, index] = col_data
        else:
            outletMassAvg[:, index] = massAvgOutletpt[field]
        
        index += 1
    
    massAvgOutletSurface = {
        'iteration': outletMassAvg[:, 0], 
        'p': outletMassAvg[:, 1], 
        'totalP': outletMassAvg[:, 2], 
        'Ux': outletMassAvg[:, 3], 
        'Uy': outletMassAvg[:, 4], 
        'Uz': outletMassAvg[:, 5],
        'alphaFlow': outletMassAvg[:, 6], 
        'radialFlow': outletMassAvg[:, 7], 
        'Urel': outletMassAvg[:, 8], 
        'relFlow': outletMassAvg[:, 9], 
        'Uradial': outletMassAvg[:, 10], 
        'Utangent': outletMassAvg[:, 11], 
        'Uaxial': outletMassAvg[:, 12]
    }

surfaceMassAvgOutlet = {'x': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/alphaFlow_patch_outlet.raw')[:,0]),
                       'y': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/alphaFlow_patch_outlet.raw')[:,1]),
                       'z': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/alphaFlow_patch_outlet.raw')[:,2]),
                       'alphaFlow': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/alphaFlow_patch_outlet.raw')[:,3]),
                       'totalP':(np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/totalP_patch_outlet.raw')[:,3]),
                       'p': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/p_patch_outlet.raw')[:,3]),
                       'radialFlow': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/radialFlow_patch_outlet.raw')[:,3]),
                       'relFlow': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/relFlow_patch_outlet.raw')[:,3]),
                       'Urel': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/Urel_patch_outlet.raw')[:,3]),
                       'Uradial': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/U:Transformed_patch_outlet.raw')[:,3]),
                       'Utangent': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/U:Transformed_patch_outlet.raw')[:,4]),
                       'Uaxial': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/U:Transformed_patch_outlet.raw')[:,5]),
                       'Ux': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/U_patch_outlet.raw')[:,3]),
                       'Uy': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/U_patch_outlet.raw')[:,4]),
                       'Uz': (np.loadtxt(timeSurfaceFolderMassAvgOutlet + '/U_patch_outlet.raw')[:,5])}


    
#%% Combine all the dictionaries 

allDict = {
    'surfaceOutletMassFlow': surfaceOutletMassFlow,
    'surfaceInletMassFlow': surfaceInletMassFlow,
    'massAvgInletSurface': massAvgInletSurface,
    'surfaceMassAvgInlet': surfaceMassAvgInlet,
    'massAvgOutletSurface': massAvgOutletSurface,
    'surfaceMassAvgRotIn': surfaceMassAvgRotIn,
    'surfaceMassAvgRotOut': surfaceMassAvgRotOut,
    'massAvgOutletSurface': massAvgOutletSurface,
    'surfaceMassAvgOutlet': surfaceMassAvgOutlet,   
}

with open(filePath + 'bladedBaseline.pkl', 'wb') as file:
    pickle.dump(allDict, file)





























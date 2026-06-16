#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 21:39:06 2025

@author: adekola
"""

import argparse
import numpy as np

#%% Used function 
def cart2pol(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    z = z
    return theta, r, z

def cart2pol2(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    theta = np.unwrap(theta)
    z = z
    return theta, r, z

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
    "periodicOrAperiodic", 
    type=int,
    help="periodic (0) or nonPeriodic (1)"
)

# Parse arguments
args = parser.parse_args()

# Use parsed arguments
Nbr = args.Nbr  # Total number of blades in the annulus for rotor
Nbs = args.Nbs  # Total number of blades in the annulus for stator
rSections = args.rSections  # Number of rotor blade profiles
sSections = args.sSections  # Number of stator blade profiles
periodicOrAperiodic = args.periodicOrAperiodic  # if periodic select 0 otherwise select 1

Ns = 101 #number o
Nr = 101
if periodicOrAperiodic == 0:
    filePath = '../inputData/periodic/'
    dataPath = '../processedData/periodic/'
else:
    filePath = '../inputData/nonPeriodic/'
    dataPath = '../processedData/nonPeriodic/'

statorCamber = np.empty([Nbs, sSections, Ns, 3]) #x,y,z
rotorCamber = np.empty([Nbr, rSections, Nr, 3]) #x,y,z

for a in range(Nbs):
    for b in range(sSections):
        statorCamber[a,b,:] = np.loadtxt(dataPath + 'stator/blade{}/camber{}.txt'.format(a,b), delimiter=',')
       
                
for aa in range(Nbr):
    for bb in range(rSections):
        rotorCamber[aa,bb,:] = np.loadtxt(dataPath + 'rotor/blade{}/camber{}.txt'.format(aa,bb), delimiter=',')

        
statorCamberCyl = np.empty([Nbs, sSections, Ns, 3]) #theta, r, z
for c in range(Nbs):
    for d in range(sSections):
        for e in range(Ns):
            statorCamberCyl[c,d,e,:] = np.array(cart2pol(statorCamber[c,d,e,0], statorCamber[c,d,e,1], statorCamber[c,d,e,2])).T
            if max(statorCamberCyl[c,d,:,0]) - min(statorCamberCyl[c,d,:,0]) > np.pi:
                statorCamberCyl[c,d,:,:] = np.array(cart2pol2(statorCamber[c,d,:,0], statorCamber[c,d,:,1], statorCamber[c,d,:,2])).T
  
rotorCamberCyl = np.empty([Nbr, rSections, Nr, 3])
for cc in range(Nbr):
    for dd in range(rSections):
        for ee in range(Nr):
            rotorCamberCyl[cc,dd,ee,:] = np.array(cart2pol(rotorCamber[cc,dd,ee,0], rotorCamber[cc,dd,ee,1], rotorCamber[cc,dd,ee,2])).T    
            if max(rotorCamberCyl[cc,dd,:,0]) - min(rotorCamberCyl[cc,dd,:,0]) > np.pi:
                rotorCamberCyl[cc,dd,:,:] = np.array(cart2pol2(rotorCamber[cc,dd,:,0], rotorCamber[cc,dd,:,1], rotorCamber[cc,dd,:,2])).T            
#%%
chordWisePosStator = np.zeros((Nbs, sSections, Ns, 4))
chordWisePosRotor = np.zeros((Nbr, rSections, Nr, 4))
for f in range(Nbs):
    for g in range(sSections):
        chordWisePosStator[f,g,:,0:3] = statorCamberCyl[f,g,:]
        diffS = np.diff(statorCamberCyl[f,g,:][:,[2,1]], axis=0)        
        seg_lengthS = np.linalg.norm(diffS, axis=1)        
        chordWisePosStator[f,g,:,3] = np.concatenate([[0], np.cumsum(seg_lengthS)])

for h in range(Nbr):
    for i in range(rSections):
        chordWisePosRotor[h,i,:,0:3] = rotorCamberCyl[h,i,:]
        diffR = np.diff(rotorCamberCyl[h,i,:][:,[2,1]], axis=0)
        seg_lengthR = np.linalg.norm(diffR, axis=1)
        chordWisePosRotor[h,i,:,3] = np.concatenate([[0], np.cumsum(seg_lengthR)])
#%% create list to stack arc lengths
statorChordWisePos = chordWisePosStator.reshape(-1,4)
rotorChordWisePos = chordWisePosRotor.reshape(-1,4)

np.savetxt(filePath + "/xc_rotor.txt", rotorChordWisePos, delimiter=",")
np.savetxt(filePath + "/xc_stator.txt", statorChordWisePos, delimiter=",")

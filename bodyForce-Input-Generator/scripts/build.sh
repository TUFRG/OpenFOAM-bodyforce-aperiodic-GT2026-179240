#!/bin/bash

set -e

# Use Spyder's Python
#PYTHON=/home/adekola/.local/spyder-6/envs/spyder-runtime/bin/python3.11 
PYTHON=python
#if your python did not work, look for the correct location of the python version that works on your system, see example of mine above
#input Parameter
Nbr=16 #Total nuber of blades in the annulus for rotor
Nbs=31 #Total nuber of blades in the annulus for stator
Nr=201 #Number of points on each surface rotor blade profile
Ns=201 #Number of points on each surface stator blade profile
rSections=23 # Number of rotor blade profiles
sSections=23 # Number of stator blade profiles
RmatrixCoeff="1 0 0 0 1 0 0 0 1" #Rotation matrix for transformation
periodicOrAperiodic=0 #if periodic select 0 otherwise select 1
stage=2 #if rotorAlone select 0, if statorAlone select 1, if stage select 2
casePath='/caseSetup/fullWheel/'
#******************************************************************************************************#
#DO NOT TOUCH ANYTHING FROM HERE
echo "Running preprocessing..."
$PYTHON preProc.py $Nbr $Nbs $rSections $sSections $Nr $Ns $periodicOrAperiodic --RmatrixCoeff $RmatrixCoeff
echo "Preprocessing complete!"
echo "Running arcLenght Definition..."
$PYTHON thetaBasedChordWisePos.py $Nbr $Nbs $rSections $sSections $periodicOrAperiodic 
echo "arcLenght Definition complete!"
echo "Running pitch Definition..."
$PYTHON thetaBasedPitch.py $Nbr $Nbs $rSections $sSections $periodicOrAperiodic 
echo "pitch Definition complete!"
echo "Running camber surface normal Definition..."
$PYTHON thetaBasedCamber.py $Nbr $Nbs $rSections $sSections $periodicOrAperiodic 
echo "camber surface normal Definition complete!"
echo "Running blade thickness Definition..."
$PYTHON thetaBasedThickness.py  $Nbr $Nbs $rSections $sSections $Nr $Ns $periodicOrAperiodic 
echo "blade thickness Definition complete!"
echo "All input data generated and Stored"




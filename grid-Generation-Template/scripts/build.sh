#!/bin/bash

set -e


#PYTHON=/home/adekola/.local/spyder-6/envs/spyder-runtime/bin/python3.11 
PYTHON=python
#if your python did not work, look for the correct location of the python version that works on your system, see example of mine above
#input Parameter
Nbr=16 #Total nuber of blades in the annulus for rotor
Nbs=31 #Total nuber of blades in the annulus for stator
rSections=23 # Number of rotor blade profiles
sSections=21 # Number of stator blade profiles
RmatrixCoeff="1 0 0 0 1 0 0 0 1" #Rotation matrix for transformation
periodicOrAperiodic=0 ##if periodic select 0 otherwise select 1
##CFD grid Parameters
scales=0.001
res=360 #Resolution of stls, leave as is unless you want a different resolution
domainExtn=1.5 #This is how much you want to extend the domain in both inlet and outlet. 1.5 means 1.5*Diameter
hubZ=21.2 # Axial value of where rotor hub stops rotatiing
spinnerIdx=13 #This is the index on the hub curve where the spinner nose ends. 
velocity=166 # Farstream velocity for the definition of boundary layer cells
target_yPlus=50
#Please note, the value of NrCells,NzCells,NqCells does not mean the final cell count in the grid, 
#it is a way of controlling the cell count which is part of the solution
#Also, for a thin slice Nq is ignored as Nq is 1
NrCells=50 #Number of cells in the radial direction 
NzCells=25 # Total Number of cells in the axial direction 
NqCells=10 # Number of cells in circumferential direction
fullWheelOrThinSlice=0 # select 0 for fullwheel and 1 for thinSlice
# For the griz size setting select any of the following cases
#rotor and stator= 0
#rotor only = 1
#stator only = 2

gridSize=0 #Select the right number based on what size of grid to be built 
stage=2 #if rotorAlone select 0, if statorAlone select 1, if stage select 2
#******************************************************************************************************#
#DO NOT TOUCH ANYTHING FROM HERE
echo "Selecting the right case..."

# Display the size of grid to be built
if [[ $fullWheelOrThinSlice -eq 0 && $gridSize -eq 0  ]]; then
    echo "Building FullWheel: Both rotor and stator"
    $PYTHON fullWheelRotorStator.py  $Nbr $Nbs $rSections $sSections $res $domainExtn  $velocity $target_yPlus $NrCells $NzCells $NqCells $hubZ $spinnerIdx $scales 
elif [[ $fullWheelOrThinSlice -eq 0 && $gridSize -eq 1 ]]; then
    echo "Building FullWheel: Rotor Only"
    $PYTHON fullWheelRotor.py  $Nbr $Nbs $rSections $sSections $res $domainExtn  $velocity $target_yPlus $NrCells $NzCells $NqCells $hubZ $spinnerIdx $scales 
elif [[ $fullWheelOrThinSlice -eq 0 && $gridSize -eq 2 ]]; then
    echo "Building FullWheel: Stator Only"
    $PYTHON fullWheelStator.py $Nbr $Nbs $rSections $sSections $res $domainExtn  $velocity $target_yPlus $NrCells $NzCells $NqCells $hubZ $spinnerIdx $scales 
elif [[ $fullWheelOrThinSlice -eq 1 && $gridSize -eq 0 ]]; then
    echo "Building ThinSlice: Both rotor and stator"
    $PYTHON thinSliceRotorStator.py  $Nbr $Nbs $rSections $sSections $res $domainExtn  $velocity $target_yPlus $NrCells $NzCells $NqCells $hubZ $spinnerIdx $scales 
elif [[ $fullWheelOrThinSlice -eq 1 && $gridSize -eq 1 ]]; then
    echo "Building ThinSlice: Rotor Only"
    $PYTHON thinSliceRotor.py $Nbr $Nbs $rSections $sSections $res $domainExtn  $velocity $target_yPlus $NrCells $NzCells $NqCells $hubZ $spinnerIdx $scales
elif [[ $fullWheelOrThinSlice -eq 1 && $gridSize -eq 2 ]]; then
    echo "Building ThinSlice: Stator Only"
    $PYTHON thinSliceStator.py $Nbr $Nbs $rSections $sSections $res $domainExtn  $velocity $target_yPlus $NrCells $NzCells $NqCells $hubZ $spinnerIdx $scales        
else
    echo "Specify if you are working with stage or not or thinslice or full Wheel !"     
fi

cd ..
if [[ $fullWheelOrThinSlice -eq 0 && $gridSize -eq 0  ]]; then
    echo "FullWheel: Both rotor and stator"
    cd  fullWheelRotorStator
    bash geomUpdate1.sh
    blockMesh
    topoSet
    createPatch -overwrite 
    touch case.foam
    cd ..
    bash 90DegMeshRotorStator.sh
    cd fullWheel
    meshCoordinates
    cellZoneCoordinates
    cd ../../bodyForce-Input-Generator/scripts/
    casePath=$(realpath ../../grid-Generation-Template/fullWheelRotorStator)
    $PYTHON functionalizedInterpolation.py  $Nbr $Nbs $rSections $sSections $scales $periodicOrAperiodic $stage --RmatrixCoeff $RmatrixCoeff $casePath
elif [[ $fullWheelOrThinSlice -eq 0 && $gridSize -eq 1 ]]; then
    echo "FullWheel: Rotor Only"
    cd fullWheelRotorOnly
    bash geomUpdate1.sh
    blockMesh
    topoSet
    createPatch -overwrite
    touch case.foam
    cd ..
    bash 90DegMeshRotor.sh
    meshCoordinates
    cellZoneCoordinates    
    cd ../../bodyForce-Input-Generator/scripts/
    casePath=$(realpath ../../grid-Generation-Template/fullWheelRotorOnly)
    $PYTHON functionalizedInterpolation.py   $Nbr $Nbs $rSections $sSections $scales $periodicOrAperiodic $stage --RmatrixCoeff $RmatrixCoeff $casePath       
elif [[ $fullWheelOrThinSlice -eq 0 && $gridSize -eq 2 ]]; then
    echo "FullWheel: Stator Only"
    cd fullWheelStatorOnly 
    bash geomUpdate1.sh
	blockMesh
	topoSet
	createPatch -overwrite
    bash 90DegMeshStator.sh
    meshCoordinates
    cellZoneCoordinates
    cd ../../bodyForce-Input-Generator/scripts/
    casePath=$(realpath ../../grid-Generation-Template/fullWheelStatorOnly)
    $PYTHON functionalizedInterpolation.py   $Nbr $Nbs $rSections $sSections $scales $periodicOrAperiodic $stage --RmatrixCoeff $RmatrixCoeff $casePath
elif [[ $fullWheelOrThinSlice -eq 1 && $gridSize -eq 0 ]]; then
    echo "ThinSlice: Both rotor and stator"
    cd thinSliceRotorStator
    bash geomUpdate1.sh
    blockMesh
    topoSet
    createPatch -overwrite
    touch case.foam
    meshCoordinates
    cellZoneCoordinates    
    cd ../../bodyForce-Input-Generator/scripts/
    casePath=$(realpath ../../grid-Generation-Template/thinSliceRotorStator)
    $PYTHON functionalizedInterpolation.py   $Nbr $Nbs $rSections $sSections $scales $periodicOrAperiodic $stage --RmatrixCoeff $RmatrixCoeff $casePath       
elif [[ $fullWheelOrThinSlice -eq 1 && $gridSize -eq 1 ]]; then
    echo "ThinSlice: Rotor Only"
    cd thinSliceRotorOnly 
    bash geomUpdate1.sh
	blockMesh
	topoSet
	createPatch -overwrite
    meshCoordinates
    cellZoneCoordinates
    cd ../../bodyForce-Input-Generator/scripts/
    casePath=$(realpath ../../grid-Generation-Template/thinSliceRotorOnly)
    $PYTHON functionalizedInterpolation.py  $Nbr $Nbs $rSections $sSections $scales $periodicOrAperiodic $stage --RmatrixCoeff $RmatrixCoeff $casePath    
elif [[ $fullWheelOrThinSlice -eq 1 && $gridSize -eq 2 ]]; then
    echo "ThinSlice: Stator Only"
    cd thinSliceStatorOnly 
    bash geomUpdate1.sh
	blockMesh
	topoSet
	createPatch -overwrite
    meshCoordinates
    cellZoneCoordinates
    cd ../../bodyForce-Input-Generator/scripts/
    casePath=$(realpath ../../grid-Generation-Template/thinSliceStatorOnly)
    $PYTHON functionalizedInterpolation.py   $Nbr $Nbs $rSections $sSections $scales $periodicOrAperiodic $stage --RmatrixCoeff $RmatrixCoeff $casePath         
else
    echo "You need to specify values within range for fullWheelOrThinslice and gridSize"    
fi



#!/bin/bash 

set -e 


PYTHON=python
echo "building the grid generator"

cd ../caseSetup/bodyForceCase/
rm -rf case* fullWheel
mkdir case{1..60}
for i in {1..60};do
	cp -r mainCase/* case$i/
done
cd ../../scripts/
$PYTHON gridGenerator.py 
cd ../caseSetup/bodyForceCase/
bash buildCase.sh
mv case1 fullWheel 
rm -rf case*
cd fullWheel
meshCoordinates 
cellZoneCoordinates
cd ../
cd ../../scripts/
$PYTHON thetaBasedCamber.py
$PYTHON thetaBasedPitch.py
$PYTHON thetaBasedChordWisePos.py 
$PYTHON functionalizedInterpolation.py
cd ../caseSetup/bodyForceCase/fullWheel
bash serial.sh 


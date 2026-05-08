#!/bin/bash

set -e 

mkdir fullWheel
cd fullWheel
cp -r ../fullWheelRotorOnly/constant .
cp -r ../fullWheelRotorOnly/system .
blockMesh 
topoSet
createPatch -overwrite
#rm -rf constant/geometry
cd ..
#rm -rf rot*
cp -r fullWheel rot1
cd rot1
rm -rf constant/geometry 
transformPoints -rotate-angle '((1 0 0) 90)'
sed -i 's/nPer1/nPer2/g' constant/polyMesh/boundary
sed -i 's/pPer1/pPer2/g' constant/polyMesh/boundary
cd ..
cp -r rot1 rot2
cd rot2
transformPoints -rotate-angle '((1 0 0) 90)'
sed -i 's/nPer2/nPer3/g' constant/polyMesh/boundary
sed -i 's/pPer2/pPer3/g' constant/polyMesh/boundary
cd ..
cp -r rot2 rot3
cd rot3
transformPoints -rotate-angle '((1 0 0) 90)'
sed -i 's/nPer3/nPer4/g' constant/polyMesh/boundary
sed -i 's/pPer3/pPer4/g' constant/polyMesh/boundary
cd ..
mergeMeshes -overwrite fullWheel rot1
mergeMeshes -overwrite fullWheel rot2
mergeMeshes -overwrite fullWheel rot3
rm -rf rot*
cd fullWheel
#mv 0 0.old
stitchMesh -overwrite  nPer1 pPer4
stitchMesh -overwrite  nPer2 pPer1
stitchMesh -overwrite  nPer3 pPer2
stitchMesh -overwrite  nPer4 pPer3
createPatch -overwrite
rm -rf 0
#mv 0.old 0


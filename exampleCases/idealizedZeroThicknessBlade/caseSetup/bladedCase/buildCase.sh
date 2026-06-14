#!/bin/bash

set -e 

for i in {1..60};do 
	cd case$i;
	echo "working on case$i now"
	bash geomUpdate1.sh
	blockMesh 
	sed -i "s/nPer/nPer$i/g" constant/polyMesh/boundary;
	sed -i "s/pPer/pPer$i/g" constant/polyMesh/boundary;
	cd ..;
done	
for i in {2..60};do 
	mergeMeshes -overwrite case1 case$i
done
cd case1
rm -rf 0
for i in {1..59};do
	echo "working on case$i now"
	stitchMesh -overwrite  "nPer$((i+1))" "pPer$i"
done
stitchMesh -overwrite nPer1 pPer60
createPatch -overwrite
rm -rf 0
touch case.foam
cd ..

#!/bin/bash 

set -e 


PYTHON=python
echo "building the grid generator"
cd ../caseSetup/bladedCase/
rm -rf case*
mkdir case{1..60}
for i in {1..60};do
	cp -r mainCase/* case$i/
done
cd ../../scripts/
$PYTHON gridGenerator.py 
cd ../caseSetup/bladedCase/
bash buildCase.sh
mv case1 fullWheel 
rm -rf case*
cd fullWheel
bash run.sh 




rhoSimpleFoam
cd system
rm controlDict
rm fvOptions
rm fvSolution
ln -s controlDict.BFfirstOrder controlDict
ln -s fvOptions.BF fvOptions
ln -s fvSolution.BF fvSolution
cd ..

changeDictionary -latestTime >> log &
rhoSimpleFoamBF >> log &
#rhoSimpleFoam >> log &
#decomposePar -fields -latestTime > log &
#mpirun -np 5 changeDictionary -latestTime -parallel >> log &
#bash calibrate.sh -b 
#mpirun -np 5 rhoSimpleFoamBF -parallel >> log &
tail -f log 
#mv *.txt postProcessing/DevlocFiles/

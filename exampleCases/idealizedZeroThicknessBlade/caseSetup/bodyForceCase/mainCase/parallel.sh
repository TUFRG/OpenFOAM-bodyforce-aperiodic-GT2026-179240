#!/bin/bash 

set -e

decomposePar 
mpirun -np 40 simpleFoamBF -parallel >> log &
tail -f log 


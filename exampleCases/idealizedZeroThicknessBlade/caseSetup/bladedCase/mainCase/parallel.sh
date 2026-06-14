#!/bin/bash 

set -e

decomposePar 
mpirun -np 40 simpleFoam -parallel >> log &
tail -f log 


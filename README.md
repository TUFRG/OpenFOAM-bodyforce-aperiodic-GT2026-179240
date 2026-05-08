# OpenFOAM-bodyforce-aperiodic-GT2026-179240
Contains submodules sourced from repositories for body force input generation, custom OpenFOAM solvers, gird generation templates, and example cases

A sample Data has been included for a quick and easy test.

This readMe will make an attempt to explain what is in each folder and how to use them.
*****************************************************************************************
HOW TO USE 

- COMPILE ALL SOLVERS 

	In the OpenFOAM-Solvers folder, at each folder level, enter 

	wclean
	wmake

- GENERATE THE BODY FORCE INPUT DATA

	In the bodyForce-Input-Generator folder, 

	cd scripts
	bash build.sh 

	This will generate all the input data needed for the body force computation. 

- BUILD GRID / INTERPOLATE DATA UNTO GRID

	In the grid-Generation-Template, 

	cd scripts
	bash build.sh 

	This will build the grid and Interpolate the data unto the grid 

*****************************************************************************************
DETAILS OF FOLDER CONTENTS 


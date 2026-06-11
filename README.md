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
- OpenFOAM-Solvers

	This folder contains all the solvers needed to run the body force model. 
	- cellZoneCoordinates
		This folder contain the solver that extracts the cell center coordinates and cellID of the cellZones defined in the domain
		To run this solver on its own, at the case folder level run
		~$ cellZoneCoordinates
		This will generate the text (.txt) files containing the cell coordinates information

	- meshCoordinates 
		This folder contain the solver that extracts the cell center and cellID of the entire computational domain
		To run this solver on its own, at the case folder level run 
		~$ meshCoordinates
		This will generate the text (.txt) file containing the cell coordinate information
		
	- decomposeMeshCoordinates
		This is similar to the meshCoordinates solver, however this is to be used for decomposed cases. 
		To run this solver on its own, at the case folder level run 
		~$ mpirun -np 10 decomposeMeshCoordinates -parallel
		user can change 10 to the number of processors in the decomposed case
		This will generate the text (.txt) file containing the cell coordinate information for each processor folder

	- simpleFoamBF 
		This folder contains the solver for the body force model for incompressible flows. 
		To run this solver on its own, at the case folder level, run
		~$ simpleFoamBF
		Please note that to run this solver, either in serial or compressible mode, user needs to copy the body force inputs fields
		from the caliberation folder to the lastest time folder before running the solver

	- rhoSimpleFOamBF 
		This folder is similar with the simpleFoamBF folder, but for compressible flows,
		To run this solver on its own, at the case folder level, run 
		~$ rhoSimpleFOamBF

- bodyForce-Input-Generator

	This folder contains all the input 

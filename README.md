# OpenFOAM-bodyforce-aperiodic-GT2026-179240
Contains folders for body force input generation, custom OpenFOAM solvers, grid generation templates, and example cases.

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

- BUILD GRID / INTERPOLATE DATA ONTO GRID

	In the grid-Generation-Template, 

	cd scripts
	bash build.sh 

	This will build the grid and Interpolate the data unto the grid 

*****************************************************************************************
DETAILS OF FOLDER CONTENTS 
- OpenFOAM-Solvers

	This folder contains all the solvers needed to run the body force model. 
	- cellZoneCoordinates - 
		This folder contain the solver that extracts the cell center coordinates and cellID of the cellZones defined in the domain
		To run this solver on its own, at the case folder level run
		~$ cellZoneCoordinates
		This will generate the text (.txt) files containing the cell coordinates information

	- meshCoordinates - 
		This folder contain the solver that extracts the cell center and cellID of the entire computational domain
		To run this solver on its own, at the case folder level run 
		~$ meshCoordinates
		This will generate the text (.txt) file containing the cell coordinate information
		
	- decomposeMeshCoordinates - 
		This is similar to the meshCoordinates solver, however this is to be used for decomposed cases. 
		To run this solver on its own, at the case folder level run 
		~$ mpirun -np 10 decomposeMeshCoordinates -parallel
		user can change 10 to the number of processors in the decomposed case
		This will generate the text (.txt) file containing the cell coordinate information for each processor folder

	- simpleFoamBF -  
		This folder contains the solver for the body force model for incompressible flows. 
		To run this solver on its own, at the case folder level, run
		~$ simpleFoamBF
		Please note that to run this solver, either in serial or compressible mode, user needs to copy the body force inputs fields
		from the caliberation folder to the lastest time folder before running the solver

	- rhoSimpleFOamBF -  
		This folder is similar with the simpleFoamBF folder, but for compressible flows,
		To run this solver on its own, at the case folder level, run 
		~$ rhoSimpleFOamBF

- bodyForce-Input-Generator

	This folder contains all the input information for the body force model
	- calibration 
	This folder contains the fields files for the body force model and the calibration scripts to run the body force model off design 
	There are two subfolders here, the parallelEnrtyBF and the singleEntryBF, which are for parallel cases and serial cases respectively
	Both folders contain the same file content
	 - 0.BF -- This folder contains the field definition for the body force model, all the files in this folder should be copied to the latestTime
	 of the case folder being run. 
	 - constant/BFControlDict -- This file should be copied to the constant folder of the case folder being run. This is where the rotation matrix
	 is defined. It also contains the setting for calibration and uncalibration mode, when the calibrationMode is 1, the model runs uncalibrated and
	 when set to 0, the model runs in calibration mode. Everyother entries are controlled using the calibrate.sh bash script which is used to write
	 out the Devloc0 fields
	 - system/fvOptionsBF -- This file should be copied into the system folder of case folder being run. 
	 Additionally, during calibration to run the model off design, the calibrate.sh and Devloc.py files should be copied to the case folder being run
	 
	- inputData - 
	This folder contains all the processed data in the format that can be converted into the continous function and then evaluated on any
	CFD mesh. It contains a subfolder for periodic and nonPeriodic cases. 
	All the entries in this folder are populated using the thetaBased scripts in the scripts folder

	- processedData
	This folder contains the processed data from the raw data provided by the user. Here the camber curve from the user provided data is defined
	by the in house code we wrote. This folder contains subfolders for periodic and nonPeriodic cases. In each of the subfolders there are
	additional subfolders rotor and stator where the processed blade and camber data are stored. 
	The entries in this folder are populated with the preProc.py script in the scripts folder

	- rawData
	This is the folder that constains the user defined geometry data. It contains two subfolders periodic and nonPeriodic. In each folder there
	is a subfolder bladeSurface, this where the user geometry is copied to. In the bladeSurface folder, the user should put the rotor and stator 
	data as applicable. Please note, the name of these folders should not be changed, for example if the geometry is an IGV ,it should be kept 
	in the stator folder. For periodic blade row, the data should be stored in blade0 folder. Each profile should be named as blade0.txt, blade1.txt,
	blade3.txt ... bladeN.txt. The format of the input data is in Cartesian (x,y,z), however, the primary axis does not matter. The rotation
	matrix will take care of the primary axis. The data should be in .txt delimeted by a comma (,). For nonPeriodic cases, each blade is stored in
	blade0, blade1, blade2, blade3 ... bladeN and in each folder, the blade profiles are named as blade0.txt and so on. 
	An example case with the IGV data is placed in this repository as a guide and an example case. 

	- scripts
	This folder contains all the scripts used to compute the input data. A bash script (build.sh) is provided that automates the order which 
	the scripts run. The bash script also contains all the necessary input information required to run the scripts.
	However, if a user wants to run the scrip manually, here is the order which the scripts should be run 

	1. preProc.py -- This converts the data to Cylindrical coordinates and then generate the camber curves and place everything to the 
	processedData folder
	2. thetaBasedCamber.py -- This generates the camber surface normals, it takes the data from the processedData folder and put the normals
	in the inputData folder
	3. thetaBasedChordWisePos.py -- This generates the local arc length, it takes the data from the processedData folder and put the normals
		in the inputData folder
	4. thetaBasedPitch.py -- This generates the angular pitch, it takes the data from the processedData folder and put the normals
		in the inputData folder
	5. thetaBasedThichkness.py -- This generates the tangential blade thickness, it takes the data from the processedData folder and put the normals
		in the inputData folder
	6. functionalizedInterpolation.py -- This is the script that converts the discrete data to continuous field and evaluate the field on any
	CFD mesh

- grid-Generation-Template
	This folder contains the provided grid templates. There are 3 difference templates as detailled in the
	paper, however, there are 6 folders here for full wheel and thin slice of rotor and stator, rotor alone and 
	stator alone template. 
	The name of each folder is self explanatory for which template the case describes. 

	-	RawData
	This folder contains the geometry data that defines the grid. The grid generator requires the camber
	surface data and the gas path
	-	gasPath 
	The has path data should come in the cylindrical coordinate system (r,theta,z). The data is need for the
	hub and casing in text formart (.txt) and should be comma delimeted. The length of data for hub and 
	casing can be different but the unit must be consistent with the camber data. The names cannot be different 
	from the name shown in this folder
	-	camberData
	The camberData is the same camber data generated from our in house code. 

*****************************************************************************************
 

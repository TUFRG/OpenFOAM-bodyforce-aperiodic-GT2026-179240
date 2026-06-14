/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

/* Turbomachinery and Unsteady Flows Research Group | University of Windsor  *\
                    ░░▒▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒░░░                                 
                 ░▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒░                            
               ░▒▓▓▒▒░░░░▒▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒░                      
             ░▒▓▓░   ▒▓░    ░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓░░                  
            ▒▓▓░   ░▒▓▓▓░ ░▓▓▓▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒░               
           ▒▓▒░░░  ░▓▓▓▓▒░▓▓▓▓▒░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒░            
         ░▓▓▒ ░░░░ ░▓▓▓▓▓ ▒▓▓▓▒  ▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒          
         ▒▓▒  ░░░░░ ▒▓▓▓▓░░▓▓▓▒  ░▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒        
        ▓▓▒   ░░░░░░░▓▓▓▓░░▓▓▓▒                        ░░▒▒▒▓▓▓▓▓▓▓▓▓▓▓▓░      
       ▒▓▓░   ░░░░░░ ░▓▓▓▒ ▓▓▓▒       ░░░▒▒▒▒▒▒▒▒▒░░░░          ░░░▒▓▓▓▓▓░     
      ▒▓▓▒░    ░░░░░░ ▒▓▓▒ ▒▓▓▒ ░▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒░░     ▒░     
      ▓▓▓░▓▓▒    ░░░░ ░▓▓▓░▒▓▓▒  ░▓▓▓▓▒░  ▒▓▓▒                   ░▒▒▓▓▓▓▒      
     ▒▓▓░▒▓▓▓▓░   ░░░░ ▒▓▓░░▓▓▒           ▒▓▓▒    ░                      ░░▒▒░ 
    ░▓▓▒ ▒▓▓▓▓▓▒   ░░░░░▓▓░░▓▓▒           ▒▓▓▒    ▒▓  ▒▓░▒▓▓▓▒ ▓▓▓▓░  ░▓▓░     
    ▒▓▓░ ▒▓▓▓▓▓▓▓░  ░░░ ▒▓▒░▓▓░           ▒▓▓▒   ░▒▓  ▓▓░▒▓░   ▓▒░▒▓░▒▓▒░▓▒░   
   ░▒▓▒░ ▒▓▓▓▓▓▓▓▓░ ░░░░░▓▒░▓▓░           ▒▓▓▒   ░▒▓  ▓▓░▒▓▓▓░ ▓▒░▒▓ ▒▓░       
   ░▓▓▒    ░▒▓▓▓▓▓▓▒ ░░░ ▓░░▓▓            ▒▓▓▒   ░▒▓  ▓▒░▒▓░   ▓▓▓▓▒░▒▓░░▓▓░   
   ▒▓▓▒       ░▒▓▓▓▓▒ ░░ ▒░░▓░            ▒▓▓▒   ░▒▓ ░▓▒░▒▓░   ▓▓ ▒▓░▒▓▓▓▓▓░   
   ▒▓▓░▒▒▒░░░    ▒▓▓▓░                    ▒▓▓▒    ░▓▓▓▓░ ▒▓░   ▓▓ ░▓▓  ▓▓      
   ▓▓▓░▓▓▓▓▓▓▓▓▒  ░▒▓░   ░░               ▒▓▓▒                          ░▒▒▒▒░ 
  ░▓▓▓░▓▓▓▓▓▓▓▓▓▓▒░ ░  ░▓▓▓▓▓▓▒▒░░░       ░░░░               ░░░░▒▒▒▓▓▓▒▒░░    
  ░▓▓▒░▓▓▓▓▓▓▓▓▓▓▓▓░  ░▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒░            
  ░▒▒░                   ░░▒▒▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒░░░                   
\*---------------------------------------------------------------------------*/

FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      fvOptions;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Construct a temporary local dictionary
BFcontrolsLocal
{
	// Include the entries from the BFcontrolsDict
	#include			"$FOAM_CASE/constant/BFcontrolsDict"
}
BFZone				$BFcontrolsLocal.allBFZone; // pull out the value we want
rBFZone				$BFcontrolsLocal.rotorBFZone;
sBFZone 			$BFcontrolsLocal.statorBFZone;
	// This is the cellzone name for all BF regions
#remove				BFcontrolsLocal; // get rid of this temporary dictionary
	// Otherwise, fvOptions thinks we're trying to define a type of fvOption
	// and fails with an error.

limitT
{
    type       limitTemperature;
    min        150; //220;
    max        500; //310;
    selectionMode all;
}
velocityDamping
{
	type 				velocityDampingConstraint;
	selectionMode 		all;
	UMax				500;
}

rotorU // momentum source term for body force zones
{
	type			vectorCodedSource;
	name			RotorU;
	active 			yes;

	vectorCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$rBFZone;
		fields			(U);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volVectorField& tF = mesh_.lookupObject<volVectorField>("F");
			const IOdictionary& BFcontrols_ = mesh_.lookupObject<IOdictionary>("BFcontrolsDict");

			scalar rDevlocWrite;
			BFcontrols_.lookup("writeRotorDevloc0") >> rDevlocWrite;

			if (rDevlocWrite)
			{
				label zoneID = mesh_.cellZones().findZoneID("Rotor");
				const labelList& 	cellr	= mesh_.cellZones()[zoneID]; // list of cells in BF zone
				std::ofstream rDevlocFile;
    			std::string timeStr = std::to_string(mesh_.time().value());
    			const volScalarField& rDEVLOC = mesh_.lookupObject<volScalarField>("DEVLOC");
    			rDevlocFile.open((mesh_.time().path() / "rDevloc_" + timeStr + ".txt").c_str(), std::ios_base::app);

				forAll(cellr, i)
				{
    				rDevlocFile << "Cell: " << cellr[i] << " rDevloc: " << rDEVLOC[cellr[i]] << std::endl;
    			}
    			rDevlocFile.close();
			}
			
			eqn += tF;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}

statorU // momentum source term for body force zones
{
	type			vectorCodedSource;
	name			StatorU;
	active 			yes;

	vectorCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$sBFZone;
		fields			(U);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volVectorField& tF = mesh_.lookupObject<volVectorField>("F");
			const IOdictionary& BFcontrols_ = mesh_.lookupObject<IOdictionary>("BFcontrolsDict");

			scalar sDevlocWrite;
			BFcontrols_.lookup("writeStatorDevloc0") >> sDevlocWrite;

			if (sDevlocWrite)
			{
				label zoneID = mesh_.cellZones().findZoneID("Stator");
				const labelList& 	cells	= mesh_.cellZones()[zoneID]; // list of cells in BF zone
				std::ofstream sDevlocFile;
    			std::string timeStr = std::to_string(mesh_.time().value());
    			const volScalarField& sDEVLOC = mesh_.lookupObject<volScalarField>("DEVLOC");
    			sDevlocFile.open((mesh_.time().path() / "sDevloc_" + timeStr + ".txt").c_str(), std::ios_base::app);

				forAll(cells, i)
				{
    				sDevlocFile << "Cell: " << cells[i] << " sDevloc: " << sDEVLOC[cells[i]] << std::endl;
    			}
    			sDevlocFile.close();
			}
			
			eqn += tF;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}

rotorE // energy source term for body force zones
{
	type			scalarCodedSource;
	name			RotorE;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$rBFZone;
		fields			(e); // or h, have to select the active one here

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volVectorField& tF = mesh_.lookupObject<volVectorField>("F");
			const volVectorField& tabc = mesh_.lookupObject<volVectorField>("abc");
			const volVectorField& tOMEGA = mesh_.lookupObject<volVectorField>("OMEGA");
			
			eqn += (tOMEGA ^ tabc) & tF;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}

statorE // energy source term for body force zones
{
	type			scalarCodedSource;
	name			StatorE;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$sBFZone;
		fields			(e); // or h, have to select the active one here

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volVectorField& tF = mesh_.lookupObject<volVectorField>("F");
			const volVectorField& tabc = mesh_.lookupObject<volVectorField>("abc");
			const volVectorField& tOMEGA = mesh_.lookupObject<volVectorField>("OMEGA");
			
			eqn += (tOMEGA ^ tabc) & tF;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}

rotork // mass source related term for turbulence quantity "k"
{
	type			scalarCodedSource;
	name			RotorK;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$rBFZone;
		fields			(k);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volScalarField& tSb = mesh_.lookupObject<volScalarField>("Sb");
			const volScalarField& tk = mesh_.lookupObject<volScalarField>("k");

			eqn -= tSb * tk;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}

statork // mass source related term for turbulence quantity "k"
{
	type			scalarCodedSource;
	name			StatorK;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$sBFZone;
		fields			(k);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volScalarField& tSb = mesh_.lookupObject<volScalarField>("Sb");
			const volScalarField& tk = mesh_.lookupObject<volScalarField>("k");

			eqn -= tSb * tk;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}
rotoromega // mass source related term for turbulence quantity "omega"
{
	type			scalarCodedSource;
	name			RotorOmega;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$rBFZone;
		fields			(omega);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volScalarField& tSb = mesh_.lookupObject<volScalarField>("Sb");
			const volScalarField& tomega = mesh_.lookupObject<volScalarField>("omega");

			eqn -= tSb * tomega;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}

statoromega // mass source related term for turbulence quantity "omega"
{
	type			scalarCodedSource;
	name			StatorOmega;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$sBFZone;
		fields			(omega);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volScalarField& tSb = mesh_.lookupObject<volScalarField>("Sb");
			const volScalarField& tomega = mesh_.lookupObject<volScalarField>("omega");

			eqn -= tSb * tomega;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}
rotornuTilda // mass source related term for turbulence quantity "nuTilda"
{
	type			scalarCodedSource;
	name			RotorNuTilda;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$rBFZone;
		fields			(nuTilda);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volScalarField& tSb = mesh_.lookupObject<volScalarField>("Sb");
			const volScalarField& tnuTilda = mesh_.lookupObject<volScalarField>("nuTilda");

			eqn -= tSb * tnuTilda;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}
statornuTilda // mass source related term for turbulence quantity "nuTilda"
{
	type			scalarCodedSource;
	name			StatorNuTilda;
	active 			yes;

	scalarCodedSourceCoeffs
	{
		selectionMode	cellZone;
		cellZone 		$sBFZone;
		fields			(nuTilda);

	    codeInclude
	    #{   	
			#include "fvCFD.H"
			#include <IFstream.H>
			#include <OFstream.H>
			#include <ctime>
			#include "clock.H"
			#include "string.H"
			#include <sstream>
			#include <iomanip>
			#include <math.h>
			#include <volFieldsFwd.H>			
		#};
		
        codeCorrect
        #{
        #};

    	codeConstrain
    	#{
   	 	#};	
    
		codeAddSupRho
		#{
			const volScalarField& tSb = mesh_.lookupObject<volScalarField>("Sb");
			const volScalarField& tnuTilda = mesh_.lookupObject<volScalarField>("nuTilda");

			eqn -= tSb * tnuTilda;
		#};

		 codeSetValue
         #{
                //Pout<< "**codeSetValue**" << endl;
         #};

         // Dummy entry. Make dependent on above to trigger recompilation
         code
         #{
                $codeInclude
                $codeCorrect
                $codeAddSupRho
                $codeSetValue
         #};
            
	}
}

//************************************************************************ //

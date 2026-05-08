/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2006                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    "(p|rho|Phi)"
    {
        solver          		GAMG;
        smoother			symGaussSeidel;
        tolerance       		1e-08;
        relTol          		0.1;//0.001;
        minIter					1;
    }

    "(p|rho|Phi)Final"
    {
        $p;
        $Phi;
	tolerance	1e-8;
        relTol		0.1;
	minIter		1;
    }

    "(U|k|epsilon|omega|e|h|nuTilda)"
    {
        solver		PBiCGStab;
        preconditioner	DILU;
        tolerance       1e-08;
        relTol          0;
	minIter		1;
    }

    "(U|k|epsilon|omega|e|h|nuTilda)Final"
    {
        $e;
        tolerance       1e-08;
        relTol          0;
        minIter         1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 3;
    pMin 5000; //40000; 
    pMax 250000; //180000; 
    transonic	yes;
    consistent	yes;

    residualControl
    {
        p               	1e-6;
        U               	1e-6;
        "(e|h)"            	1e-6;
        "(k|epsilon|omega|nuTilda)"	1e-6;
    }

}



relaxationFactors
{
    fields
    {
        p               	1.0;
    }
    equations
    {
        U               	0.7; //0.9;
        "(k|epsilon|omega|nuTilda)"	0.7;
        "(e|h)"           	0.7; //0.9;
    }
}



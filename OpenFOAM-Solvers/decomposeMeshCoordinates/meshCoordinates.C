#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "meshTools.H"
#include <fstream>

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Reading decomposed case on processor " << Pstream::myProcNo() << endl;

    // Construct the filename for this processor
    word fileName = runTime.path()/"cellCoordinates.txt";

    // Open a file to write the coordinates for this processor
    std::ofstream outFile(fileName.c_str());
    if (!outFile.is_open())
    {
        FatalErrorInFunction
            << "Unable to open file for writing: " << fileName
            << exit(FatalError);
    }

    const volVectorField& C = mesh.C();

    // Write the coordinates to the file for this processor
    forAll(C, cellI)
    {
        outFile << "Cell " << cellI << ": (" 
                << C[cellI].x() << ", " 
                << C[cellI].y() << ", " 
                << C[cellI].z() << ")" 
                << std::endl;
    }

    // Close the file
    outFile.close();

    Info << "Cell coordinates written to " << fileName << " by processor " << Pstream::myProcNo() << endl;
    Info << "End\n" << endl;

    return 0;
}

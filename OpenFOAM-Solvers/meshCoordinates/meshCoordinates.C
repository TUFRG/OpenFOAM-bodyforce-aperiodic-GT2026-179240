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

    Info << "Reading mesh" << endl;

    const volVectorField& C = mesh.C();

    // Open a file to write the coordinates
    std::ofstream outFile("cellCoordinates.txt");
    if (!outFile.is_open())
    {
        FatalErrorInFunction
            << "Unable to open file for writing: cellCoordinates.txt"
            << exit(FatalError);
    }

    // Write the coordinates to the file
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

    Info << "Cell coordinates written to cellCoordinates.txt\n" << endl;

    Info << "End\n" << endl;

    return 0;
}

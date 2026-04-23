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

    // Access cell zones
    const cellZoneMesh& cellZones = mesh.cellZones();

    // Check if there are any cell zones
    if (cellZones.empty())
    {
        FatalErrorInFunction
            << "No cell zones found in the mesh."
            << exit(FatalError);
    }

    // Loop over each cell zone
    forAll(cellZones, zoneI)
    {
        const cellZone& cZone = cellZones[zoneI];
        const labelList& zoneCells = cZone;  // List of cell indices in this zone

        // Create a unique filename for each zone
        std::string zoneName = cZone.name();
        std::string outFileName = "cellCoordinates_" + zoneName + ".txt";
        std::ofstream outFile(outFileName);

        if (!outFile.is_open())
        {
            FatalErrorInFunction
                << "Unable to open file for writing: " << outFileName
                << exit(FatalError);
        }

        Info << "Extracting coordinates for cell zone: " << zoneName << endl;

        // Write the coordinates for cells in this zone
        forAll(zoneCells, cellI)
        {
            label cellIndex = zoneCells[cellI];
            outFile << "Cell " << cellIndex << " (Zone " << zoneI << "): (" 
                    << C[cellIndex].x() << ", " 
                    << C[cellIndex].y() << ", " 
                    << C[cellIndex].z() << ")" 
                    << std::endl;
        }

        // Close the file
        outFile.close();

        Info << "Cell coordinates written to " << outFileName << "\n" << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

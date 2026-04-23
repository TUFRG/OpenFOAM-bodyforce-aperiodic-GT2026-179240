#!/bin/bash
#created by Adekola Adeyemi
set -e 

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -r, --rotor-only     Calibrate rotor only"
    echo "  -s, --stator-only    Calibrate stator only"
    echo "  -b, --both          Calibrate both rotor and stator (default)"
    echo "  -h, --help          Show this help message"
    exit 1
}

# Default settings
calibrateRotor=1
calibrateStator=1

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--rotor-only)
            calibrateRotor=1
            calibrateStator=0
            echo "Mode: Calibrating ROTOR only"
            shift
            ;;
        -s|--stator-only)
            calibrateRotor=0
            calibrateStator=1
            echo "Mode: Calibrating STATOR only"
            shift
            ;;
        -b|--both)
            calibrateRotor=1
            calibrateStator=1
            echo "Mode: Calibrating BOTH rotor and stator"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Display selected calibration mode
if [[ $calibrateRotor -eq 1 && $calibrateStator -eq 1 ]]; then
    echo "Calibration mode: Both rotor and stator"
elif [[ $calibrateRotor -eq 1 ]]; then
    echo "Calibration mode: Rotor only"
elif [[ $calibrateStator -eq 1 ]]; then
    echo "Calibration mode: Stator only"
else
    echo "Error: At least one component must be calibrated"
    exit 1
fi

n=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
echo 'Setting BF solver to write out Devloc0 using current Devloc.'
echo 'Note! This assumes the flow field is converged.'

# Set calibration flags based on user selection
if [[ $calibrateRotor -eq 1 ]]; then
    echo 'Enabling rotor calibration...'
    foamDictionary -entry writeRotorDevloc0 -set 1 constant/BFcontrolsDict
    foamDictionary -entry rotorCalibrationMode -set 1 constant/BFcontrolsDict
else
    echo 'Disabling rotor calibration...'
    foamDictionary -entry writeRotorDevloc0 -set 0 constant/BFcontrolsDict
    foamDictionary -entry rotorCalibrationMode -set 0 constant/BFcontrolsDict
fi

if [[ $calibrateStator -eq 1 ]]; then
    echo 'Enabling stator calibration...'
    foamDictionary -entry writeStatorDevloc0 -set 1 constant/BFcontrolsDict
    foamDictionary -entry statorCalibrationMode -set 1 constant/BFcontrolsDict
else
    echo 'Disabling stator calibration...'
    foamDictionary -entry writeStatorDevloc0 -set 0 constant/BFcontrolsDict
    foamDictionary -entry statorCalibrationMode -set 0 constant/BFcontrolsDict
fi

foamDictionary -entry calibrationMode -set 1 constant/BFcontrolsDict 

echo 'Getting latest time...'
LT=$(foamListTimes -processor -latestTime)
echo 'Latest time = ' + $LT
priorEndTime=$(foamDictionary -entry endTime -value system/controlDict)
priorwriteInterval=$(foamDictionary -entry writeInterval -value system/controlDict)
newEndTime=$(echo "$LT + 1" | bc)
echo 'Updating controlDict to run 1 more iteration and write.'
foamDictionary -entry endTime -set $newEndTime system/controlDict
foamDictionary -entry writeInterval -set 1 system/controlDict

echo 'Writing out mesh coordinates...'
mpirun -np $n decomposeMeshCoordinates -parallel

echo 'Running solver 1 iteration to write out Devloc .txt file...'
mpirun -np $n rhoSimpleFoamBF -parallel 


echo 'Running Python script to process txt file into format used for reading in field...'
if [[ $calibrateRotor -eq 1 && $calibrateStator -eq 1 ]]; then
    python Devloc.py --numberOfSubDomains $n --both
elif [[ $calibrateRotor -eq 1 ]]; then
    python Devloc.py --numberOfSubDomains $n --rotor-only
elif [[ $calibrateStator -eq 1 ]]; then
    python Devloc.py --numberOfSubDomains $n --stator-only
fi

echo 'Devloc file written successfully.'
echo 'Updating BFcontrolsDict to get out of calibration / devloc writing mode...'
foamDictionary -entry rotorCalibrationMode -set 0 constant/BFcontrolsDict
foamDictionary -entry writeRotorDevloc0 -set 0 constant/BFcontrolsDict
foamDictionary -entry statorCalibrationMode -set 0 constant/BFcontrolsDict
foamDictionary -entry writeStatorDevloc0 -set 0 constant/BFcontrolsDict
foamDictionary -entry calibrationMode -set 0 constant/BFcontrolsDict

echo 'Restoring entries in controlDict...'
foamDictionary -entry endTime -set $priorEndTime system/controlDict
foamDictionary -entry writeInterval -set $priorwriteInterval system/controlDict

echo 'Copying codestream-based Devloc0 file to latest time folder...'
for dir in processor*/$newEndTime/; do
	cp 0/Devloc0 $dir	
done

# Summary message
echo 'Calibration setup complete. Running the solver again will proceed using'
if [[ $calibrateRotor -eq 1 && $calibrateStator -eq 1 ]]; then
    echo 'the calibrated body force model for BOTH rotor and stator.'
elif [[ $calibrateRotor -eq 1 ]]; then
    echo 'the calibrated body force model for ROTOR only.'
elif [[ $calibrateStator -eq 1 ]]; then
    echo 'the calibrated body force model for STATOR only.'
fi

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 10:22:19 2024
@author: Adekola Adeyemi
Modified to handle selective rotor/stator calibration
"""
import argparse
import numpy as np
import re
import subprocess
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description="Getting the number of subDomains.")
parser.add_argument(
    "-n", "--numberOfSubDomains", 
    type=int,
    required=True,
    help="numberOfSubDomains"
)
parser.add_argument(
    "--rotor-only", 
    action="store_true",
    help="Process rotor Devloc only"
)
parser.add_argument(
    "--stator-only", 
    action="store_true",
    help="Process stator Devloc only"
)
parser.add_argument(
    "--both", 
    action="store_true",
    help="Process both rotor and stator Devloc (default)"
)

args = parser.parse_args()
n = args.numberOfSubDomains

# Determine which components to process
if args.rotor_only and args.stator_only:
    print("Error: Cannot specify both --rotor-only and --stator-only")
    exit(1)
elif args.rotor_only:
    process_rotor = True
    process_stator = False
    print("Processing ROTOR Devloc only")
elif args.stator_only:
    process_rotor = False
    process_stator = True
    print("Processing STATOR Devloc only")
else:
    # Default: process both (maintains backward compatibility)
    process_rotor = True
    process_stator = True
    print("Processing BOTH rotor and stator Devloc")

time = subprocess.check_output(["foamListTimes -processor -latestTime"], shell=True).strip().decode()
print(time)

for k in range(n):
    # Initialize data containers
    rCell_devloc_pairs = []
    sCell_devloc_pairs = []
    
    # Process rotor file if needed
    if process_rotor:
        rotor_file_path = 'processor{}/rDevloc_'.format(k) + time + '.000000.txt'
        if os.path.exists(rotor_file_path):
            with open(rotor_file_path, 'r') as file:
                devlocFileR = file.read()
            rCell_devloc_pairs = re.findall(r"Cell:\s*(\d+)\s*rDevloc:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", devlocFileR)
        else:
            print("Warning: {} not found. Skipping rotor for processor{}".format(rotor_file_path, k))
    
    # Process stator file if needed
    if process_stator:
        stator_file_path = 'processor{}/sDevloc_'.format(k) + time + '.000000.txt'
        if os.path.exists(stator_file_path):
            with open(stator_file_path, 'r') as file:
                devlocFileS = file.read()
            sCell_devloc_pairs = re.findall(r"Cell:\s*(\d+)\s*sDevloc:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", devlocFileS)
        else:
            print("Warning: {} not found. Skipping stator for processor{}".format(stator_file_path, k))
                
    with open('processor{}/cellCoordinates.txt'.format(k), 'r') as file1:
        cellData = file1.read()
    
    cellNumExtract = re.findall(r"Cell\s*(\d+)", cellData)
    cellNum = [int(cell) for cell in cellNumExtract]
    
    # Combine rotor and stator data
    all_cell_devloc_pairs = rCell_devloc_pairs + sCell_devloc_pairs
    
    # Convert to list of tuples
    cell_devloc_pairs = [(int(cell), float(devloc)) for cell, devloc in all_cell_devloc_pairs]
    
    devlocData = np.empty([len(cell_devloc_pairs),2])

    index = 0
    for cell, devloc in cell_devloc_pairs:
        devlocData[index,0] = cell
        devlocData[index,1] = devloc
        index += 1
     
    newDevloc = np.zeros([len(cellNum)])

    devloc_dict = {int(cell): devloc for cell, devloc in cell_devloc_pairs}
    for i, cell in enumerate(cellNum):
        if cell in devloc_dict:
            newDevloc[i] = devloc_dict[cell]
            '''
    loc = np.where(np.isin(cellNum, devlocData[:,0]))[0]

    count = 0
    for i in loc:
        if devlocData[count,0] == loc[count]:
            newDevloc[i] = devlocData[count, 1]
        count += 1
    '''
    filename = 'processor{}/Devloc'.format(k)
    
    f = open(filename, 'w')
    f.write('( \n')
    for j in newDevloc:
        f.write('{} \n'.format(j))
    f.write(')')
    f.close()
    
    print("Devloc file written for processor{}".format(k))

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
import matplotlib.pyplot as plt
import subprocess


#%%
parser = argparse.ArgumentParser(description="Checking which Calibration Mode is Active")
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
    
time = subprocess.check_output(["foamListTimes -latestTime"], shell=True).strip().decode()
print(time)
rCell_devloc_pairs = []
sCell_devloc_pairs = []
if process_rotor:
    rotor_file_path = 'rDevloc_' + time + '.000000.txt'
    with open(rotor_file_path, 'r') as file:
        devlocFileR = file.read()
    rCell_devloc_pairs = re.findall(r"Cell:\s*(\d+)\s*rDevloc:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", devlocFileR)

if process_stator:
    stator_file_path = 'sDevloc_' + time + '.000000.txt'
    with open(stator_file_path, 'r') as file:
        devlocFileS = file.read()
    sCell_devloc_pairs = re.findall(r"Cell:\s*(\d+)\s*sDevloc:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", devlocFileS)
        
with open('cellCoordinates.txt', 'r') as file1:
    cellData = file1.read()



cellNumExtract = re.findall(r"Cell\s*(\d+)", cellData)  ##r"Cell\s*(\d+)"

# Combine rotor and stator data
all_cell_devloc_pairs = rCell_devloc_pairs + sCell_devloc_pairs

# Convert to a list of tuples
cell_devloc_pairs = [(int(cell), float(devloc)) for cell, devloc in all_cell_devloc_pairs]
cellNum = [int(cell) for cell in cellNumExtract]

devlocData = np.empty([len(cell_devloc_pairs),2])

index = 0
for cell, devloc in cell_devloc_pairs:
    devlocData[index,0] = cell
    devlocData[index,1] = devloc
    index += 1
 
newDevloc = np.zeros([len(cellNum)])

loc = np.where(np.isin(cellNum, devlocData[:,0]))[0]

count = 0
for i in loc:
    if devlocData[count,0] == loc[count]:
        newDevloc[i] = devlocData[count, 1]
    count += 1

filename = 'Devloc'

f = open(filename, 'w')
f.write('( \n')
for j in newDevloc:
    f.write('{} \n'.format(j))
f.write(')')
f.close()




















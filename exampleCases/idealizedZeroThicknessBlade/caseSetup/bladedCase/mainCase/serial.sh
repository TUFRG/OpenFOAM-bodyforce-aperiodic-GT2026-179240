#!/bin/bash 

set -e

simpleFoam >> log &
tail -f log 

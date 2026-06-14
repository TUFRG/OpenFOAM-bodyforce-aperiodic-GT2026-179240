#!/bin/bash 

set -e

simpleFoamBF >> log &
tail -f log 

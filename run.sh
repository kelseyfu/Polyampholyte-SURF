#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# Author     : Kelsey Fu 
# Date       : 2024-06-25
# Description: Script to modify the LAMMPS input file and run the simulation

# Define the current working directory
CWD_PATH=$(pwd)
LAMMPS_PATH="lmp"

# Define the input parameters (example values, you may replace them)
export TEMP="1"
export NCHAIN="1"

# export sequence="[-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1, 1]"
# export SEQUNAME="ABAB"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB"

 export sequence="[1,1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,-1,1,-1,1,-1,-1,1]"
 export SEQUNAME="RAND"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,1,1,1,1,0,0,0,0,0,0,0,0,-1,1,-1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="GRAD"

source "script/simulate.sh"

neg
pos
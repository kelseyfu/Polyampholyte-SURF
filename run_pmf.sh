#!/usr/bin/env bash
#SBATCH --job-name=a0.2_2_LBOUND_b
#SBATCH --time=7-00:00:00
# Slurm: Node configuration
#SBATCH --partition=cpu

#SBATCH --nodes=1 --ntasks-per-node=4 --mem=4G
#SBATCH --gres=gpu:0 --gpu-bind=closest

# -*- coding: utf-8 -*-
#
# Author     : Kelsey Fu 
# Date       : 2024-06-25
# Description: Script to modify the LAMMPS input file and run the simulation

# Define the current working directory
CWD_PATH=$(pwd)
LAMMPS_PATH="/home/pjwalker/Polyampholyte-SURF/software/lammps/build/lmp"


# source "script/simulate.sh" # 1 & 100 chain
source "script/simulate_pmf.sh"
# source "script/simulate_bulk.sh"
# source "script/simulate_pmf  Only 2 chain
# # source "script/analysis.sh" Only cluster analysis 4 100 chain







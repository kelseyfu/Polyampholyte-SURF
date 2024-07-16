#!/usr/bin/env bash
#SBATCH --job-name=GRADh12
#SBATCH --time=2-00:00:00
# Slurm: Node configuration
#SBATCH --partition=subset
#SBATCH --account=zgw
#SBATCH --qos=zgw
#SBATCH --exclude=node11

#SBATCH --nodes=1 --ntasks-per-node=2 --mem=4G
#SBATCH --gres=gpu:1 --gpu-bind=closest

#SBATCH --output=/home/kfu/slurm-reports/slurm-%j.out --error=/home/kfu/slurm-reports/slurm_error-%j.out

# -*- coding: utf-8 -*-
#
# Author     : Kelsey Fu 
# Date       : 2024-06-25
# Description: Script to modify the LAMMPS input file and run the simulation

# Define the current working directory
CWD_PATH=$(pwd)
LAMMPS_PATH="lmp"

source "inputs/GRAD_N_+12_nchain_100.sh"

source "script/simulate.sh" # 1 & 100 chain
# source "script/simulate_pmf.sh"
# source "script/simulate_pmf  Only 2 chain
# # source "script/analysis.sh" Only cluster analysis 4 100 chain







#!/usr/bin/env bash
#SBATCH --job-name=NEW_GRAD_2
#SBATCH --time=2-00:00:00
# Slurm: Node configuration
#SBATCH --partition=subset
#SBATCH --account=zgw
#SBATCH --qos=zgw

#SBATCH --nodes=1 --ntasks-per-node=2 --mem=4G
#SBATCH --gres=gpu:0 --gpu-bind=closest

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
export NCHAIN="2"
export LENGTH="50"

export NCPU="2"
export NGPU="0"

# export sequence="[-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1, 1]"
# export SEQUNAME="ABAB_N_2"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_N_2"

#  export sequence="[1,1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,-1,1,-1,1,-1,-1,1]"
#  export SEQUNAME="RAND_N_2"

export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,1]"
export SEQUNAME="NEWGRAD_N_2"




# # NET + 25% ------------------------------------------------
#  export sequence="[-1,1,-1,1,1,1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1]"
#  export SEQUNAME="ABAB_+25"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_+25"

#  export sequence="[1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1,-1,1,1]"
#  export SEQUNAME="RAND_+25"

# export sequence="[-1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="NEWGRAD_+25"


# # NET + 75% ------------------------------------------------
#  export sequence="[-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1]"
#  export SEQUNAME="ABAB_+75"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_+75"

#  export sequence="[1,1,-1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,-1,1,-1,1,1,1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,1]"
#  export SEQUNAME="RAND_+75"

# export sequence="[-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="NEWGRAD_+75"





# # NET - 25% ------------------------------------------------
#  export sequence="[-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1]"
#  export SEQUNAME="ABAB_-25"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_-25"

#  export sequence="[-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,-1,1]"
#  export SEQUNAME="RAND_-25"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,1,1]"
# export SEQUNAME="NEWGRAD_-25"




# # NET - 75% ------------------------------------------------
#  export sequence="[-1,-1,-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1]"
#  export SEQUNAME="ABAB_-75"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_-75"

#  export sequence="[-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1]"
#  export SEQUNAME="RAND_-75"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,1,1]"
# export SEQUNAME="NEWGRAD_-75"

#_________________________________________________________________

#reference neutral
# export sequence="[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]"
# export SEQUNAME="NEUTRAL_REF"

#Reference Charged
# export sequence="[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="CHARGED_REF"

# source "script/simulate.sh"
source "script/simulate_pmf.sh"
# source "script/analysis.sh"


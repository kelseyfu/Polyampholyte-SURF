#!/usr/bin/env bash
#SBATCH --job-name=GRAD_N
#SBATCH --time=2-00:00:00
# Slurm: Node configuration
#SBATCH --partition=subset
#SBATCH --account=zgw
#SBATCH --qos=zgw

#SBATCH --nodes=1 --ntasks-per-node=4 --mem=4G
#SBATCH --gres=gpu:1 --gpu-bind=closest

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
export NCHAIN="100"
export LENGTH="70"

export NCPU="2"
export NGPU="1"

# export sequence="[-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1, 1]"
# export SEQUNAME="ABAB"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_N"

#  export sequence="[1,1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,-1,1,-1,1,-1,-1,1]"
#  export SEQUNAME="RAND_N"

export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,1,1,1,1,0,0,0,0,0,0,0,0,-1,1,-1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,1]"
export SEQUNAME="GRAD_N"




# # NET + 25% ------------------------------------------------
#  export sequence="[-1,1,-1,1,1,1,-1,1,1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1]"
#  export SEQUNAME="ABAB_+25"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_+25"

#  export sequence="[1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,1,-1,1,-1,1,1]"
#  export SEQUNAME="RAND_+25"

# export sequence="[-1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,1,1,-1,1,0,0,0,0,0,0,0,0,-1,1,1,1,1,-1,1,-1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="GRAD_+25"


# # NET + 75% ------------------------------------------------
#  export sequence="[-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,1]"
#  export SEQUNAME="ABAB_+75"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_+75"

#  export sequence="[1,1,-1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,-1,1,-1,1,1,1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,1,-1,1,1,1]"
#  export SEQUNAME="RAND_+75"

# export sequence="[-1,-1,1,1,-1,-1,1,-1,1,1,1,-1,1,1,1,1,0,0,0,0,0,0,0,0,-1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="GRAD_+75"





# # NET - 25% ------------------------------------------------
#  export sequence="[-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1]"
#  export SEQUNAME="ABAB_-25"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_-25"

#  export sequence="[-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,-1,1]"
#  export SEQUNAME="RAND_-25"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,1,0,0,0,0,0,0,0,0,-1,1,-1,-1,1,-1,-1,-1,1,-1,1,1,1,1,1]"
# export SEQUNAME="GRAD_-25"




# # NET - 75% ------------------------------------------------
#  export sequence="[-1,-1,-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,1,1,-1,-1]"
#  export SEQUNAME="ABAB_-75"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="AABB_-75"

#  export sequence="[-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1]"
#  export SEQUNAME="RAND_-75"

# export sequence="[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,0,0,0,0,0,0,0,0,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,1,1]"
# export SEQUNAME="GRAD_-75"


#reference neutral
# export sequence="[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]"
# export SEQUNAME="NEUTRAL_REF"

#Reference Charged
# export sequence="[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]"
# export SEQUNAME="CHARGED_REF"

source "script/simulate.sh"


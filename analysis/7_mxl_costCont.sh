#!/bin/bash

#SBATCH --job-name=mxl_costCont
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/hdd/scratch/%u/logs/%x-%j.out
#SBATCH --error=/mnt/hdd/scratch/%u/logs/%x-%j.err

mkdir -p /mnt/hdd/scratch/$USER/logs

cd /mnt/hdd/scratch/$USER/jpn_dce

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE

Rscript analysis/7_mxl_costCont.r

exit $?
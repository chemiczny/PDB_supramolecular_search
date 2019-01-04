#!/bin/env bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=2GB
##### Nazwa kolejki
#SBATCH -p plgrid
#SBATCH -A supraml

module load plgrid/tools/pymol/2.2.0
module add plgrid/tools/python/3.6.5

python3 simpleRun.py

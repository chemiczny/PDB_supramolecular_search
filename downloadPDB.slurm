#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=70:00:00
##### Nazwa kolejki
#SBATCH -p plgrid
#SBATCH -A supraml

module add plgrid/tools/python/3.6.5

python3 downloadPDB.py
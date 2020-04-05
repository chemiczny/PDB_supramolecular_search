#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
##### Nazwa kolejki
#SBATCH -p plgrid-short

module add plgrid/tools/python/3.6.5

python3 updatePDB.py

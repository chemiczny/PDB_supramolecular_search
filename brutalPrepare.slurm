#!/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2GB
##### Nazwa kolejki
#SBATCH -p plgrid-short

module add plgrid/tools/python/3.6.5

python3 brutalPrepare.py

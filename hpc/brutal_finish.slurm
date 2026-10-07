#!/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2GB
##### Queue (partition) name
#SBATCH -p plgrid-testing
#SBATCH -A supraml

module add plgrid/tools/python/3.6.5

python3 scripts/brutal_finish.py

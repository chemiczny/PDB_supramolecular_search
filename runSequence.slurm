#!/bin/env bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=50GB
##### Nazwa kolejki
#SBATCH -p plgrid-testing

module add plgrid/tools/python/3.6.5

python3 sekwencja_new.py &> sequence.log

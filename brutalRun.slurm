#!/bin/env bash
#SBATCH -N 18
#SBATCH --ntasks-per-node=24
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=3GB
##### Nazwa kolejki
#SBATCH -p plgrid
#SBATCH -A supraml

module load plgrid/tools/pymol/2.2.0
module add plgrid/tools/python/3.6.5

cat cif2process.dat | xargs -t -d "\n" -P ${SLURM_NTASKS} -n 1 srun -n 1 -N 1 --mem=3gb python3 brutalRun.py

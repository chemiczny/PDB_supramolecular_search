#!/bin/env bash
#SBATCH -N 18
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2GB
##### Nazwa kolejki
#SBATCH -p plgrid-testing

module add plgrid/tools/python/3.6.5

python3 brutalPrepare.py

cat cif2process.dat | xargs -t -d "\n" -P ${SLURM_NTASKS} -n 1 srun -n 1 -N 1 --mem=2gb python3 brutalRun.py

python3 brutalFinish.py

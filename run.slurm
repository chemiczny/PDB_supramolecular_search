#!/bin/env bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=2GB
##### Nazwa kolejki
#SBATCH -p plgrid

module add plgrid/tools/python/3.6.5

python3 PDB_ligandPi_anion_search_full_screen.py

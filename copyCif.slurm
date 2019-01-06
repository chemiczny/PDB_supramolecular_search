#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=70:00:00
##### Nazwa kolejki
#SBATCH -p plgrid
#SBATCH -A supraml

cp /net/archive/groups/plggsuprm/cif/*cif cif/

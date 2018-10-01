#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:31:05 2018

@author: michal
"""

from cif_analyser import findSupramolecular
from os.path import basename, join
import sys


cif = sys.argv[1]
PDBcode = basename(cif).split(".")[0].upper()
cif = join( "/net/scratch/people/plgglanow/slurm_jobdir/cif/", cif)
arguments =  ( cif, PDBcode, True )
findSupramolecular(arguments)
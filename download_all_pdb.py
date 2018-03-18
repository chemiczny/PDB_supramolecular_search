#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 11:06:13 2018

@author: emilia
"""

from os.path import isdir
from os import remove, makedirs
import glob
from Bio.PDB import PDBList

if not isdir("cif"):
    makedirs("cif")
else:
    cifiles = glob.glob('/cif/*')
    for cifile in cifiles:
        remove(cifile)
        
pdb_dir = "cif/"
pl = PDBList(pdb=pdb_dir)
pl.flat_tree = True
pl.download_entire_pdb(file_format = "mmCif")

#uruchomiono: 18.03.2018 20:15
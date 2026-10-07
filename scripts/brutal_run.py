#!/usr/bin/env python3
"""
Created on Mon Oct  1 14:31:05 2018

@author: michal
"""

from supramolecular_search.cif_analyser import find_supramolecular
from os.path import basename
import sys
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(processName)s %(levelname)s %(name)s: %(message)s",
)

cif_files = []

cif_list_filename = sys.argv[1]
cif_list_file = open(cif_list_filename, "r")
line = cif_list_file.readline()

while line:
    cif_files.append(line.strip())
    line = cif_list_file.readline()

cif_list_file.close()

file_id = basename(cif_list_filename).split(".")[0]

for cif in cif_files:
    pdb_code = basename(cif).split(".")[0].upper()
    arguments = (cif, pdb_code, file_id)
    find_supramolecular(arguments)

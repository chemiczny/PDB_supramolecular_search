#!/usr/bin/env python3
"""
Created on Tue Jun  5 13:19:32 2018

@author: michal

Compares parsing speed of Biopython's MMCIF2Dict and primitiveCif2Dict.
Usage: python3 benchmarks/cif_parsing.py "cif/*.cif"
"""

import sys
import time
from glob import glob
from Bio.PDB import MMCIF2Dict
from supramolecular_search.primitive_cif2dict import PrimitiveCif2Dict


cifs = glob(sys.argv[1] if len(sys.argv) > 1 else "cif/*.cif")
time_start = time.time()

for cif in cifs:
    print(cif)
    test = MMCIF2Dict.MMCIF2Dict(cif)

time_stop = time.time()

print(time_stop - time_start)

time_start = time.time()

for cif in cifs:
    print(cif)
    test = PrimitiveCif2Dict(
        cif, ["_refine.ls_d_res_high", "_reflns_shell.d_res_high", "_exptl.method"]
    )
    print(test.result)

time_stop = time.time()

print(time_stop - time_start)

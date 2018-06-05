#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 13:19:32 2018

@author: michal
"""
import sys
if not "../" in sys.path:
    sys.path.append('../')
    
import time
from glob import glob
from Bio.PDB import MMCIF2Dict
from primitiveCif2Dict import primitiveCif2Dict


    
cifs = glob("../cif/*cif")
timeStart = time.time()

for cif in cifs:
    print(cif)
    test = MMCIF2Dict.MMCIF2Dict(cif)
    
timeStop = time.time()

print(timeStop-timeStart)

timeStart = time.time()

for cif in cifs:
    print(cif)
    test = primitiveCif2Dict(cif, ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" , "_exptl.method" ])
    print(test.result)
    
timeStop = time.time()

print(timeStop-timeStart)
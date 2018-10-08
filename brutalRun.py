#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:31:05 2018

@author: michal
"""

from cif_analyser import findSupramolecular
from os.path import basename
import sys

cifFiles = []

cifListFilename = sys.argv[1]
cifListFile = open(cifListFilename, 'r')
line = cifListFile.readline()

while line:
    cifFiles.append(line.strip())
    line = cifListFile.readline()

cifListFile.close()

fileId = basename(cifListFilename).split(".")[0]

for cif in cifFiles:
    PDBcode = basename(cif).split(".")[0].upper()
    arguments =  ( cif, PDBcode, fileId )
    findSupramolecular(arguments)
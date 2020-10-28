#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:26:23 2018

@author: michal
"""
from supramolecularLogging import writeAnionPiHeader, writeAnionCationHeader, writePiPiHeader, writeCationPiHeader, writeHbondsHeader, writeMetalLigandHeader
from supramolecularLogging import writeAnionPiLinearHeader, writeAnionPiPlanarHeader
from os.path import isdir, join
from os import makedirs, remove
import glob
import time
from configure import configure
config = configure()

cifFiles =  config["cif"]
        
if not isdir("logs"):
    makedirs("logs")
    
writeAnionPiHeader()
writeAnionCationHeader()
writePiPiHeader()
writeCationPiHeader()
writeHbondsHeader()
writeMetalLigandHeader()
writeAnionPiLinearHeader()
writeAnionPiPlanarHeader()

open("logs/additionalInfo.log", "w").close()

timeStart = time.time()
timeFile = open("logs/timeStart.log", 'w')
timeFile.write(str(timeStart))
timeFile.close()


cif_files = glob.glob(cifFiles)     
cifNoFile = open("logs/cif2process.log", 'w')
cifNoFile.write(str(len(cif_files)))
cifNoFile.close()

scratch = config["scratch"]
files2remove = glob.glob( join(scratch,"*"))

for f in files2remove:
    remove(f)

filesForStep = 70
actualId = 0
filesForCurrentStep = 0

inputFile = open("cif2process.dat", 'w')

actualFile = open(join(scratch, str(actualId)+".dat"), 'w')
inputFile.write(join(scratch, str(actualId)+".dat")+"\n")
for cif in cif_files:
    if filesForCurrentStep > filesForStep:
        actualId += 1
        actualFile.close()
        actualFile = open(join(scratch, str(actualId)+".dat"), 'w')
        inputFile.write(join(scratch, str(actualId)+".dat")+"\n")
        filesForCurrentStep = 0
        
    actualFile.write( cif + "\n" )
    filesForCurrentStep+=1
    
actualFile.close()
inputFile.close()
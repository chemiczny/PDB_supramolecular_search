#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:26:23 2018

@author: michal
"""
from supramolecularLogging import writeAnionPiHeader, writeAnionCationHeader, writePiPiHeader, writeCationPiHeader
from os.path import isdir, isfile
from os import makedirs, remove
import glob, json
import time

cifFiles =  "cif/*.cif" 

configurationFileName = "config.json"

if isfile(configurationFileName):
    configFile = open(configurationFileName)
    config = json.load(configFile)
    configFile.close()
    
    if "N" in config:
        numberOfProcesses = config["N"]
        
    if "cif" in config:
        cifFiles = config["cif"]
        
if not isdir("logs"):
        makedirs("logs")
        
log_files = glob.glob("logs/anionPi*.log")
log_files += glob.glob("logs/piPi*.log")
log_files += glob.glob("logs/cationPi*.log")
log_files += glob.glob("logs/anionCation*.log")
log_files += glob.glob("logs/partialProgress*")
log_files += glob.glob("logs/additionalInfo*.log")
log_files += glob.glob("logs/piPiPi.log")

for log_file in log_files:
    remove(log_file)
    
writeAnionPiHeader()
writeAnionCationHeader()
writePiPiHeader()
writeCationPiHeader()

timeStart = time.time()
timeFile = open("logs/timeStart.log", 'w')
timeFile.write(str(timeStart))
timeFile.close()


cif_files = glob.glob(cifFiles)     
cifNoFile = open("logs/cif2process.log", 'w')
cifNoFile.write(str(len(cif_files)))
cifNoFile.close()
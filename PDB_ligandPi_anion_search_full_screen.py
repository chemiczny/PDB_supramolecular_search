# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""
from cif_analyser import findSupramolecular
from supramolecularLogging import writeAnionPiHeader, writeAnionCationHeader, writePiPiHeader, writeCationPiHeader
from os.path import isdir, basename, isfile
from os import makedirs, remove
import glob
import time
import json
from multiprocessing import Pool

numberOfProcesses = 6
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


cif_files = glob.glob(cifFiles)                                                                                                                         

dataProcessed = 0
structure_saved = 0
dataLen = len(cif_files)

pool = Pool(numberOfProcesses)

def prepareArgumentsList(cifFiles ):
    arguments = []
    
    for cif in cifFiles:
        PDBcode = basename(cif).split(".")[0].upper()
        arguments.append( ( cif, PDBcode, True ) )
        
    return arguments

timeStart = time.time()
timeFile = open("logs/timeStart.log", 'w')
timeFile.write(str(timeStart))
timeFile.close()

cifNoFile = open("logs/cif2process.log", 'w')
cifNoFile.write(str(len(cif_files)))
cifNoFile.close()


argumentsList = prepareArgumentsList( cif_files)
pool.map(findSupramolecular, argumentsList)
##
#for arg in argumentsList:
#    findSupramolecular(arg)

def mergeLogs( logFinal, logs  ):
    final_log = open(logFinal, "a+")
    log_files = glob.glob(logs)
    for log_file in log_files:
        if log_file == logFinal:
            continue
        
        new_log = open(log_file, 'r')
        
        line = new_log.readline()
        while line:
            final_log.write(line)
            line = new_log.readline()
        
        new_log.close()
        remove(log_file)
    
    final_log.close()
#
mergeLogs("logs/anionPi.log", "logs/anionPi*.log"  )
mergeLogs("logs/cationPi.log", "logs/cationPi*.log"  )
mergeLogs("logs/piPi.log", "logs/piPi*.log"  )
mergeLogs("logs/anionCation.log", "logs/anionCation*.log"  )
mergeLogs("logs/additionalInfo.log", "logs/additionalInfo*.log"  )

timeStop = time.time()
timeFile = open("logs/timeStop.log", 'w')
timeFile.write(str(timeStop))
timeFile.close()
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""
from cif_analyser import findSupramolecular
from supramolecularLogging import writeAnionPiHeader, writeAnionCationHeader, writePiPiHeader, writeCationPiHeader
from os.path import isdir, basename
from os import makedirs, remove
import glob
import time
from multiprocessing import Pool

#write header
            
if not isdir("logs"):
        makedirs("logs")
        
log_files = glob.glob("logs/anionPi*.log")
log_files += glob.glob("logs/piPi*.log")
log_files += glob.glob("logs/cationPi*.log")
log_files += glob.glob("logs/anionCation*.log")
log_files += glob.glob("logs/partialProgress*")
log_files += glob.glob("logs/additionalInfo*.log")
for log_file in log_files:
    remove(log_file)
    
writeAnionPiHeader()
writeAnionCationHeader()
writePiPiHeader()
writeCationPiHeader()


cif_files = glob.glob( "cif/*.cif" )
    
dataProcessed = 0
structure_saved = 0
dataLen = len(cif_files)

numberOfProcesses = 6
pool = Pool(numberOfProcesses)

def prepareArgumentsList(cifFiles):
    arguments = []
    
    for cif in cifFiles:
        PDBcode = basename(cif).split(".")[0].upper()
        arguments.append( ( cif, PDBcode ) )
        
    return arguments

argumentsList = prepareArgumentsList( cif_files)

timeStart = time.time()
timeFile = open("logs/timeStart.log", 'w')
timeFile.write(str(timeStart))
timeFile.close()

cifNoFile = open("logs/cif2process.log", 'w')
cifNoFile.write(str(len(cif_files)))
cifNoFile.close()

pool.map(findSupramolecular, argumentsList)
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
    
    final_log.close()

mergeLogs("logs/anionPi.log", "logs/anionPi*.log"  )
mergeLogs("logs/cationPi.log", "logs/cationPi*.log"  )
mergeLogs("logs/piPi.log", "logs/piPi*.log"  )
mergeLogs("logs/anionCation.log", "logs/anionCation*.log"  )
mergeLogs("logs/additionalInfo.log", "logs/additionalInfo*.log"  )

timeStop = time.time()
timeFile = open("logs/timeStop.log", 'w')
timeFile.write(str(timeStop))
timeFile.close()
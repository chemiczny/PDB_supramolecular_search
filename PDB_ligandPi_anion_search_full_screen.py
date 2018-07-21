# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""
from cif_analyser import findSupramolecular
from supramolecularLogging import writeAnionPiHeader, writeAnionCationHeader, writePiPiHeader, writeCationPiHeader, writeAdditionalInfo
from os.path import isdir, basename
from os import makedirs, remove
import glob
import time
from multiprocessing import Pool

#write header
#restart = False

#if not restart:
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

#if restart:
#    log_files = glob.glob("logs/anionPi*.log")
#    print(log_files)
#    PDBprocessed = set()
#    for log in log_files:
#        if log != "logs/piPi.log":
#            continue
#        
#        logDataFrame = pd.read_csv(log, sep="\t")
#        newPDB = logDataFrame[0].unique()
#        print(newPDB)

dataProcessed = 0
structure_saved = 0
dataLen = len(cif_files)

numberOfProcesses = 6
pool = Pool(numberOfProcesses)


processedPDB = []

def prepareArgumentsList(cifFiles, PDB2doNotProcess ):
    arguments = []
    
    for cif in cifFiles:
        PDBcode = basename(cif).split(".")[0].upper()
        if not PDBcode in PDB2doNotProcess:
            arguments.append( ( cif, PDBcode ) )
            PDB2doNotProcess.append(PDBcode)
        
    return arguments, PDB2doNotProcess

timeStart = time.time()
timeFile = open("logs/timeStart.log", 'w')
timeFile.write(str(timeStart))
timeFile.close()

cifNoFile = open("logs/cif2process.log", 'w')
cifNoFile.write(str(len(cif_files)))
cifNoFile.close()


argumentsList, processedPDB = prepareArgumentsList( cif_files, processedPDB)
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
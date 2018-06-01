# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""
from cif_analyser import findSupramolecularAnionPiAllLigands
from supramolecularLogging import writeSupramolecularSearchHeader
from os.path import isdir, basename
from os import makedirs, remove
import glob
import time
from multiprocessing import Pool

#write header
            
if not isdir("logs"):
        makedirs("logs")
        
log_files = glob.glob("logs/MergeResultsFromLigprepOutput*.log")
log_files += glob.glob("logs/partialProgress*")
for log_file in log_files:
    remove(log_file)
    
writeSupramolecularSearchHeader()

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

pool.map(findSupramolecularAnionPiAllLigands, argumentsList)
#for arg in argumentsList:
#    findSupramolecularAnionPiAllLigandsMultiProcess(arg)

final_log = open("logs/MergeResultsFromLigprepOutput.log", "a+")
log_files = glob.glob("logs/MergeResultsFromLigprepOutput*.log")
for log_file in log_files:
    if log_file == "logs/MergeResultsFromLigprepOutput.log":
        continue
    
    new_log = open(log_file, 'r')
    
    line = new_log.readline()
    while line:
        final_log.write(line)
        line = new_log.readline()
    
    new_log.close()

final_log.close()

timeStop = time.time()
timeFile = open("logs/timeStop.log", 'w')
timeFile.write(str(timeStop))
timeFile.close()
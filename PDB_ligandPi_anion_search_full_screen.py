# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""

from cif_analyser import writeSupramolecularSearchHeader, findSupramolecularAnionPiAllLigandsMultiProcess
from PDB_requests import getLigandCodeFromSdf
from os.path import isdir, basename, isfile
from os import makedirs, remove
import glob
import time, datetime
from multiprocessing import Pool, Lock
    
def writeProgres(dataProcessed, allData, time_start):
    progressFile = open("logs/progress.log", "w")
    progressFile.write("Postep: "+ str(dataProcessed)+" / "+str(dataLen)+" "+ str(dataProcessed*100.0/dataLen)[0:5]+"%\n")
    
    timeActual = time.time()
    timeTaken = timeActual - timeStart
    timeEstimated = timeTaken/dataProcessed * (allData  - dataProcessed)
    prettyTimeTaken = str(datetime.timedelta(seconds = timeTaken))
    progressFile.write("Czas pracy: "+prettyTimeTaken+"\n")
    prettyTimeEstimated = str( datetime.timedelta(seconds = timeEstimated) )
    progressFile.write("Estymowany pozostaly czas: "+prettyTimeEstimated+"\n")
    progressFile.close()
    
#Prepare anions data
sdfFromLigprep = "sdf/ligprep_2-out_cutted.sdf"
ligprepData = {}

if isfile( sdfFromLigprep ):
    anionsCodes = getLigandCodeFromSdf(sdfFromLigprep)
    anionsCodes = list(set(anionsCodes))
    anionsCodes.append("CL")
    
    ligprepData = { "anionNames" : anionsCodes }

#write header
            
if not isdir("logs"):
        makedirs("logs")
        
log_files = glob.glob("logs/MergeResultsFromLigprepOutput*.log")
for log_file in log_files:
    remove(log_file)
    
writeSupramolecularSearchHeader()

cif_files = glob.glob( "cif/*.cif" )
    
dataProcessed = 0
structure_saved = 0
dataLen = len(cif_files)

numberOfProcesses = 6
pool = Pool(numberOfProcesses)
lock = Lock()

def prepareArgumentsList(cifFiles, lockObject):
    arguments = []
    
    for cif in cifFiles:
        PDBcode = basename(cif).split(".")[0].upper()
        arguments.append( ( cif, PDBcode ) )
        
    return arguments

argumentsList = prepareArgumentsList( cif_files, lock )

timeStart = time.time()
pool.map(findSupramolecularAnionPiAllLigandsMultiProcess, argumentsList)
        
writeProgres(dataLen, dataLen, timeStart)

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
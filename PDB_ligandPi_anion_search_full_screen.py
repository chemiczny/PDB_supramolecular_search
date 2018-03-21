# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""

from cif_analyser import writeSupramolecularSearchHeader, findSupramolecularAnionPiAllLigands
from PDB_requests import getLigandCodeFromSdf
from os.path import isdir, basename, isfile
from os import makedirs, remove
import glob
import time, datetime
    
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
log_file = "logs/MergeResultsFromLigprepOutput.log"
if isfile(log_file):
    remove(log_file)
    
writeSupramolecularSearchHeader()

cif_files = glob.glob( "cif/*.cif" )

necessary_dirs = [ "xyz", "xyz/ligands", "xyz/ligands_ENV", "xyz/ligands_anions", "xyz/aminoacids", "xyz/aminoacids_ENV", "xyz/aminoacids_anions" ]

for project_dir in necessary_dirs:
    if not isdir(project_dir):
        makedirs(project_dir)
    else:
        xyz2delete = glob.glob(project_dir+"/*.xyz")
        for xyz in xyz2delete:
            remove(xyz)
    

dataProcessed = 0
structure_saved = 0
dataLen = len(cif_files)

timeStart = time.time()
for cif_file in cif_files:
    PDBcode = basename(cif_file).split(".")[0].upper()

    dataProcessed += 1
    
    supramolecularFound = findSupramolecularAnionPiAllLigands( cif_file, PDBcode, ligprepData )
    
    if supramolecularFound:
        structure_saved += 1
        
    if structure_saved == 10:
        break      
    
    if dataProcessed % 10 == 0:
        writeProgres(dataProcessed, dataLen, timeStart)
        
        
writeProgres(dataLen, dataLen, timeStart)

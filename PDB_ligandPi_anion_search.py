# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""

from cif_analyser import writeSupramolecularSearchHeader, findSupramolecularAnionPiLigand
from PDB_requests import getLigandCodeFromSdf
from Bio.PDB import PDBList
from os.path import isfile, isdir
from os import remove, makedirs
import time

def getLigandCodePDBcodeFromLog( logFileName ):
    logFile = open(logFileName, 'r')
    results = []    
    line = logFile.readline()
    
    while line:
        data = line.split(":")
        ligandCode = data[0].strip()
        PDBcode = data[1].strip()
        results.append( { "ligandCode" : ligandCode, "PDBcode": PDBcode }  )
        line=logFile.readline()
    
    logFile.close()
    
    return results
    
def writeProgres(dataProcessed, allData):
    progressFile = open("logs/progress.log", "w")
    progressFile.write("Postep: "+ str(dataProcessed)+"/"+str(dataLen)+ str(dataProcessed*100.0/dataLen)[0:5]+"%\n")
    progressFile.close()
    

#logInput = "logs/anionPiSearchTest.log"
#logInput = "logs/aromaty_wiecej_niz_1_pierscien_podst_elektrofilowe_2Parsed.log"
print("Zaczynamy! Tralala")
logInput = "logs/MergeResultsFromLigprep.log"
sdfFromLigprep = "sdf/ligprep_2-out_cutted.sdf"

anionsCodes = getLigandCodeFromSdf(sdfFromLigprep)
anionsCodes = list(set(anionsCodes))
anionsCodes.append("CL")

ligprepData = { "anionNames" : anionsCodes }

writeSupramolecularSearchHeader()

data = getLigandCodePDBcodeFromLog(logInput)
pdbl = PDBList()
notFoundList = []

if not isdir("cif"):
    makedirs("cif")
    
if not isdir("xyz"):
    makedirs("xyz")
    
if not isdir("xyz/ligands"):
    makedirs("xyz/ligands")

if not isdir("xyz/ligands_ENV"):
    makedirs("xyz/ligands_ENV")

if not isdir("xyz/anions"):
    makedirs("xyz/anions")
    
dataLen = float(len(data))
dataProcessed = 0

timeStart = time.time()

previousPDBcode = ""
cifFile = ""

for record in data:
    ligandCode = record["ligandCode"]
    PDBcode = record["PDBcode"]
    
    if not PDBcode == previousPDBcode:
        if isfile(cifFile):
            remove(cifFile)
        cifFile = pdbl.retrieve_pdb_file( PDBcode, pdir="cif", file_format="mmCif" )
        
    dataProcessed += 1
    if not isfile(cifFile):
        notFoundList.append(PDBcode)
        continue
    
    supramolecularFound = findSupramolecularAnionPiLigand( ligandCode, cifFile, PDBcode, ligprepData )
    
    #if not supramolecularFound:
    
    previousPDBCode = PDBcode        
    
    if dataProcessed % 10 == 0:
        writeProgres(dataProcessed, dataLen)
        

timeStop = time.time()            
print("Kody PDB, ktorych nie udao sie pobrac: ", len(notFoundList))
print(notFoundList)
print("Czas analizy: ", timeStop-timeStart)
    

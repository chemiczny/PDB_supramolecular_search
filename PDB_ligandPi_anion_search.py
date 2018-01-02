# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""

from cif_analyser import writeSupramolecularSearchHeader, findSupramolecularAnionPiLigand
from Bio.PDB import PDBList
from os.path import isfile, isdir
from os import remove, makedirs
import time

def getLigandCodePDBcodeFromLog( logFileName ):
    logFile = open(logFileName, 'r')
    results = []    
    line = logFile.readline()
    
    while line:
        data = line.split()
        ligandCode = data[0].replace(":", "").strip()
        PDBcode = data[1].strip()
        results.append( { "ligandCode" : ligandCode, "PDBcode": PDBcode }  )
        line=logFile.readline()
    
    logFile.close()
    
    return results
    
def writeProgres(dataProcessed, allData):
    progressFile = open("logs/progress.log", "w")
    progressFile.write("Postep: ", str(dataProcessed)+"/"+str(dataLen), str(dataProcessed*100.0/dataLen)[0:5]+"%\n")
    progressFile.close()
    

logInput = "logs/anionPiSearchTest.log"
#logInput = "logs/aromaty_wiecej_niz_1_pierscien_podst_elektrofilowe_2Parsed.log"
writeSupramolecularSearchHeader()

data = getLigandCodePDBcodeFromLog(logInput)
pdbl = PDBList()
notFoundList = []

if not isdir("cif"):
    makedirs("cif")
    
dataLen = float(len(data))
dataProcessed = 0

timeStart = time.time()
for record in data:
    ligandCode = record["ligandCode"]
    PDBcode = record["PDBcode"]
    cifFile = pdbl.retrieve_pdb_file( PDBcode, pdir="cif", file_format="mmCif" )
    dataProcessed += 1
    if not isfile(cifFile):
        notFoundList.append(PDBcode)
        continue
    
    supramolecularFound = findSupramolecularAnionPiLigand( ligandCode, cifFile, PDBcode )
    
    #if not supramolecularFound:
    remove(cifFile)
        
    
    if dataProcessed % 10 == 0:
        writeProgres(dataProcessed, dataLen)

timeStop = time.time()            
print("Kody PDB, ktorych nie udao sie pobrac: ", len(notFoundList))
print(notFoundList)
print("Czas analizy: ", timeStop-timeStart)
    
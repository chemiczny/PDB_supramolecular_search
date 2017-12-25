# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 21:51:04 2017

@author: unknown
"""


def getPDBCodeFromLog(logFileName):
    logFile = open(logFileName, 'r')
    
    outputName = "allPDB.log"
    outputFile = open( outputName, "a+" )
    
    line = logFile.readline()
    pdbRubbish=[]
    pdbNo = 0
    while line:
        if "pdb code" in line:
            pdbRubbish = line.split("pdb code:")
            del pdbRubbish[0]
            for i in pdbRubbish:
                PDBCode = i.split("method")[0].strip()
                outputFile.write(PDBCode+", ")
                pdbNo += 1
    
        line = logFile.readline()
        
    outputFile.close()
    logFile.close()
    print("Emilko, skarbie moj, ty zawsze przy mnie stoj, rano, wieczor, we dnie w nocy jestem dla Ciebie do pomocy. Dzisiaj znalazlas: "+str(pdbNo)+"kodow PDB")

    

getPDBCodeFromLog("wiecej_niz_1_pierscien_obecny_aromat_i_metal.log")
getPDBCodeFromLog("aromaty_wiecej_niz_1_pierscien_podst_elektrofilowe_2.log")



# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 21:51:04 2017

@author: unknown
"""


def getPDBCodeFromLog(logFileName):
    logFile = open(logFileName, 'r')
    
    outputName = logFileName[:-4]+"Parsed.log"
    outputFile = open( outputName, "a+" )
    
    line = logFile.readline()
    pdbRubbish=[]
    pdbNo = 0
    LigandCodeNo = 0
    while line:
        if "ligand code" in line:
            pdbRubbish = line.split("ligand code: ")
            del pdbRubbish[0]
            for i in pdbRubbish:
                LigandCode = i.split(": ")[0].strip()
                #outputFile.write()
                LigandCodeNo += 1
                
                if "pdb code" in i:
                    pdbRubbish = i.split("pdb code:")
                    del pdbRubbish[0]
                    for i in pdbRubbish:
                        PDBCode = i.split("method")[0].strip()
                        outputFile.write(LigandCode + ": "+PDBCode+'\n')
                        pdbNo += 1       

        line = logFile.readline()
        
    outputFile.close()
    logFile.close()
    print("Znaleziono: "+str(pdbNo)+" kodow PDB" +" oraz " +str(LigandCodeNo)+" kodów liganda")
  
#getPDBCodeFromLog("logs/tescik.log")
#getPDBCodeFromLog("logs/wiecej_niz_1_pierscien_obecny_aromat_i_metal.log")
getPDBCodeFromLog("logs/aromaty_wiecej_niz_1_pierscien_podst_elektrofilowe_2.log")



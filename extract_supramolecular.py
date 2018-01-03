# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 14:54:53 2018

@author: cj4bdw
"""

def extract_data ( logFileName, badPDB = [] ):
    logFile = open(logFileName, 'r')
    extractedLog = open(logFileName[:-4]+"Extracted.log", 'w')
    line = logFile.readline()
    extractedLog.write(line)
    line = logFile.readline()
    
    under70 = 0
    over110 = 0
    underDict = {}
    overDict = {}
    
    pdbA = {}
    while line:
        lineSpl = line.split()        
#        distance = float(lineSpl[5])
        aType = lineSpl[2]
        angle  = float(lineSpl[6])
        pdb = lineSpl[0]
        
#        if pdb in [ "1I18", "1J5I", "2MGY", "2MGY", "2N4J", "2K31", "2KOT", "4KI7", "1MUX", "1I48", "3KNT", "2LRR", "1PUS", "1PUN", "2M53", "2ICZ",  ]:
        if pdb in badPDB:
            line= logFile.readline()
            continue
        
        extractedLog.write(line)
#        if aType == "ASP" or aType == "8XU":
        if not pdb in pdbA:
            pdbA[pdb] = 1
        else:
            pdbA[pdb] += 1
        
        if angle < 70 :
            under70 += 1
            
            if not aType in underDict.keys():
                underDict[aType] = 1
            else:
                underDict[aType] += 1
        elif angle > 110:
            over110 += 1
                
            if not aType in overDict.keys():
                overDict[aType] = 1
            else:
                overDict[aType] += 1
        
        line= logFile.readline()

    logFile.close()    
    extractedLog.close()
    print(under70)

    print(over110)

    #print(underDict)
    #print(overDict)
    print("common:")
    for underKey in sorted(underDict.keys()):
        if underKey in overDict:
            if underDict[underKey] > 100 or overDict[underKey] >100:
                print(underKey, underDict[underKey], overDict[underKey])
            
    print("unique under")
    for underKey in sorted(underDict.keys()):
        if not underKey in overDict and underDict[underKey] > 100:
            print(underKey, underDict[underKey])
            
    print("unique over")
    for underKey in sorted(overDict.keys()):
        if not underKey in underDict and overDict[underKey] > 100:
            print(underKey, overDict[underKey])
    
    print("mainPDB")
    averagePDBconn = 0
    for pdb in pdbA:
        averagePDBconn += pdbA[pdb]
        if pdbA[pdb] > 100:
            print(pdb, pdbA[pdb])
    
    print("All pdb: ", len(pdbA))
    print("Average connections in pdb: ", float(averagePDBconn)/len(pdbA))
    return pdbA
    
pdbDict = extract_data("logs/anionPiLigandBig.log")
badPDB = []
for pdb in pdbDict:
    if pdbDict[pdb] > 5:
        badPDB.append(pdb)

extract_data("logs/anionPiLigandBig.log", badPDB)
print("Wywalilem: ", len(badPDB))
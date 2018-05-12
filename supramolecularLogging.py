#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 13:44:59 2018

@author: michal
"""

from math import  sin, cos, radians, sqrt, acos, degrees
import numpy as np
from numpy_utilities import normalize

def writeSupramolecularSearchHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/MergeResultsFromLigprepOutput.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    resultsFile.write("Anion code\tAnion chain\tAnion id\tAnion type\t")
    resultsFile.write("Atom symbol\tDistance\tAngle\t")
    resultsFile.write("x\th\t")
    resultsFile.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    resultsFile.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    resultsFile.write("Model No\tDisordered\t")
    resultsFile.write("Ring size\t")
    resultsFile.write("Resolution\t")
    resultsFile.write("Metal cations\n")
    resultsFile.close()
    
            
def writeSupramolecularSearchResults( ligand, PDBcode, centroid, extractedAtoms, modelIndex, resolution, cationNear, fileId = None ):
    """
    Zapisz dane do pliku z wynikami
    """
    resultsFileName = "logs/MergeResultsFromLigprepOutput.log"
    if fileId != None:
        resultsFileName = "logs/MergeResultsFromLigprepOutput"+str(fileId)+".log"
    ligandCode = ligand.get_resname()
    ligandId = str(ligand.get_id()[1])
    ligandChain = ligand.get_parent().get_id()
    resultsFile = open(resultsFileName, "a+")
    newAtoms = []
    for atomData in extractedAtoms:
        distance = atomDistanceFromCentroid( atomData["Atom"], centroid )
        angle = atomAngleNomVecCentroid( atomData["Atom"], centroid )
        
        h = cos(radians( angle ))*distance
        x = sin(radians( angle ))*distance
        
#        angleOK = angle <= 45 or angle >= 135
#        xOK = x < 1.6
#        hOK = h >= 1.5 and h <= 4
#        if angleOK and xOK and hOK:
        newAtoms.append(atomData)
            
        atomCoords = atomData["Atom"].get_coord()
        centroidCoords = centroid["coords"]        
        
        anion = atomData["Atom"].get_parent()
        residueName = anion.get_resname()
        anionChain = anion.get_parent().get_id()
        anionId = str(anion.get_id()[1])
        resultsFile.write(PDBcode+"\t")
        resultsFile.write(ligandCode+"\t")
        resultsFile.write(ligandChain+"\t")
        resultsFile.write(ligandId+"\t")
        resultsFile.write(residueName+"\t")
        resultsFile.write(anionChain+"\t")
        resultsFile.write(anionId+"\t")
        resultsFile.write(atomData["AnionType"]+"\t")
        resultsFile.write(atomData["Atom"].element+"\t")
        
        resultsFile.write(str(distance)+"\t")
        resultsFile.write(str(angle)+"\t")
        
        resultsFile.write(str(x)+"\t")
        resultsFile.write(str(h)+"\t")
        
        resultsFile.write(str(centroidCoords[0])+"\t")
        resultsFile.write(str(centroidCoords[1])+"\t")
        resultsFile.write(str(centroidCoords[2])+"\t")
        
        resultsFile.write(str(atomCoords[0])+"\t")
        resultsFile.write(str(atomCoords[1])+"\t")
        resultsFile.write(str(atomCoords[2])+"\t")
        
        resultsFile.write(str(modelIndex)+"\t")
        resultsFile.write(str(atomData["Atom"].get_parent().is_disordered()) + "\t")
        resultsFile.write(str(centroid["ringSize"])+"\t")
        resultsFile.write(str(resolution)+"\t")
        if cationNear:
            resultsFile.write(str(cationNear)+"\n")
        else:
            resultsFile.write("None\n")
    
    resultsFile.close()
    
    return newAtoms

def atomDistanceFromCentroid( atom, centroid ):
    """
    Funkcja pomocnicza, oblicza odleglosc pomiedzy atomem a srodkiem 
    pierscienia
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    centroid - slownik, klucze: coords, normVec
    
    Wyjscie:
    odleglosc (float)
    """
    atomCoords = atom.get_coord()
    centoridCoords = centroid["coords"]
    
    dist = 0
    for atomCoord, centroidCoord in  zip( atomCoords, centoridCoords ):
        dist+= (atomCoord-centroidCoord)*(atomCoord-centroidCoord)
        
    return sqrt(dist)

def atomAngleNomVecCentroid( atom, centroid ):
    """
    Funkcja pomocnicza, oblicza kat pomiedzy kierunkiem od srodka 
    pierscienia do atomu a wektorem normalnym plaszczyzny pierscienia
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    centroid - slownik, klucze: coords, normVec
    
    Wyjscie:
    kat w stopniach
    """
    atomCoords = np.array(atom.get_coord())
    centroidCoords = np.array(centroid["coords"])
    normVec = centroid["normVec"]
    
    centrAtomVec= normalize( atomCoords - centroidCoords )
    inner_prod = np.inner( normVec, centrAtomVec )
    
    return degrees( acos(inner_prod) )
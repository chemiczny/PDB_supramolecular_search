#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 13:44:59 2018

@author: michal
"""

from math import  sin, cos, radians, sqrt, acos, degrees
import numpy as np
from numpy_utilities import normalize
from os.path import isfile, join
from ringDetection import getNormVec, getAverageCoords

def writeAnionPiHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/anionPi.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    resultsFile.write("Anion code\tAnion chain\tAnion id\tAnion type\tAnion group id\t")
    resultsFile.write("Atom symbol\tDistance\tAngle\t")
    resultsFile.write("x\th\t")
    resultsFile.write("CentroidId\t")
    resultsFile.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    resultsFile.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    resultsFile.write("Model No\tDisordered\t")
    resultsFile.write("Ring size\tRing elements\t")
    resultsFile.write("Resolution\t")
    resultsFile.write("Method\tStructure type\n")
    resultsFile.close()
    
def writeCationPiHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/cationPi.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    resultsFile.write("Cation code\tCation chain\tCation id\t")
    resultsFile.write("Atom symbol\tDistance\tAngle\t")
    resultsFile.write("x\th\t")
    resultsFile.write("RingChain\t")
    resultsFile.write("ChainFlat\t")
    resultsFile.write("Cation-Chain Distance\t")
    resultsFile.write("CentroidId\t")
    resultsFile.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    resultsFile.write("Cation x coord\tCation y coord\tCation z coord\t")
    resultsFile.write("Model No\n")
    resultsFile.close()
    
def writeMetalLigandHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/metalLigand.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\t")
    resultsFile.write("Cation code\tCation chain\tCation id\t")
    resultsFile.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    resultsFile.write("Cation element\t")
    resultsFile.write("Ligand element\t")
    resultsFile.write("isAnion\tanionType\tDistance\t")
    resultsFile.write("Cation x coord\tCation y coord\tCation z coord\t")
    resultsFile.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    resultsFile.write("Complex\tSummary\tCoordNo\t")
    resultsFile.write("Model No\n")
    resultsFile.close()
    
def writePiPiHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/piPi.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    resultsFile.write("Pi res code\tPi res chain\tPi res id\t")
    resultsFile.write("Distance\tAngle\t")
    resultsFile.write("x\th\t")
    resultsFile.write("theta\t")
    resultsFile.write("omega\t")
    resultsFile.write("CentroidId\t")
    resultsFile.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    resultsFile.write("Centroid 2 x coord\tCentroid 2 y coord\tCentroid 2 z coord\t")
    resultsFile.write("Model No\t")
    resultsFile.write("Ring size 2\n")
    resultsFile.close()
    
def writeAnionCationHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/anionCation.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tCation code\tCation chain\tCation id\t")
    resultsFile.write("Pi acid Code\tPi acid chain\tPiacid id\t")
    resultsFile.write("CentroidId\t")
    resultsFile.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    resultsFile.write("Anion symbol\tCation symbol\tDistance\t")
    resultsFile.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    resultsFile.write("Cation x coord\tCation y coord\tCation z coord\t")
    resultsFile.write("Same semisphere\t")
    resultsFile.write("Latitude diff\t")
    resultsFile.write("Model No\n")
    resultsFile.close()
    
def writeHbondsHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/hBonds.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tAnion code\tAnion chain\tAnion id\t")
    resultsFile.write("Donor code\tDonor chain\tDonor id\t")
    resultsFile.write("Acceptor group\tAcceptor atom\tAnion group id\t")
    resultsFile.write("Acceptor x coord\tAcceptor y coord\tAcceptor z coord\t")
    resultsFile.write("Donor group\tDonor atom\t")
    resultsFile.write("Donor x coord\tDonor y coord\tDonor z coord\t")
    resultsFile.write("Hydrogen x coord\tHydrogen y coord\tHydrogen z coord\t")
    resultsFile.write("H from Experm\tAngle\tDistance H Acc\t")
    resultsFile.write("Distance Don Acc\tModel No\n")
    resultsFile.close()
    
def writeAnionPiPlanarHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/planarAnionPi.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    resultsFile.write("CentroidId\t")
    resultsFile.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    resultsFile.write("Angle\t")
    resultsFile.write("DirectionalAngle\t")
    resultsFile.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    resultsFile.write("Anion group x coord\tAnion group y coord\tAnion group z coord\t")
    resultsFile.write("Model No\n")
    resultsFile.close()
    
def writeAnionPiLinearHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/linearAnionPi.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    resultsFile.write("CentroidId\t")
    resultsFile.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    resultsFile.write("Angle\t")
    resultsFile.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    resultsFile.write("Anion group x coord\tAnion group y coord\tAnion group z coord\t")
    resultsFile.write("Model No\n")
    resultsFile.close()

class SupramolecularLogger:
    def __init__(self, PDBcode, fileId = None, scratchDir = None):
        self.pdbCode = PDBcode
        self.fileId = fileId
        self.scratchDir = scratchDir
        self.finalLogDir = "logs"
        
        if fileId:
            self.additionalInfoLog = join( self.scratchDir, "additionalInfo{}.log".format(self.fileId) )
            self.anionPiLog = join( self.scratchDir, "anionPi{}.log".format(self.fileId) )
            self.planarAnionPiLog = join( self.scratchDir, "planarAnionPi{}.log".format(self.fileId) )
            self.linearAnionPiLog = join( self.scratchDir, "linearAnionPi{}.log".format(self.fileId) )
            self.cationPiLog = join( self.scratchDir, "cationPi{}.log".format(self.fileId) )
            self.metalLigandLog = join( self.scratchDir, "metalLigand{}.log".format(self.fileId) )
            self.piPiLog = join( self.scratchDir, "piPi{}.log".format(self.fileId) )
            self.anionCationLog = join( self.scratchDir, "anionCation{}.log".format(self.fileId) )
            self.hBondsLog = join( self.scratchDir, "hBonds{}.log".format(self.fileId) )
            
        else:
            self.additionalInfoLog = join( self.finalLogDir, "additionalInfo.log" )
            self.anionPiLog = join( self.finalLogDir, "anionPi.log" )
            self.planarAnionPiLog = join( self.finalLogDir, "planarAnionPi.log" )
            self.linearAnionPiLog = join( self.finalLogDir, "linearAnionPi.log" )
            self.cationPiLog = join( self.finalLogDir, "cationPi.log" )
            self.metalLigandLog = join( self.finalLogDir, "metalLigand.log" )
            self.piPiLog = join( self.finalLogDir, "piPi.log" )
            self.anionCationLog = join( self.finalLogDir, "anionCation.log" )
            self.hBondsLog = join( self.finalLogDir, "hBonds.log" )
            
        self.partialProgressLog = join( self.scratchDir, "partialProgress{}.log".format(self.fileId) )

    def writeAdditionalInfo(self, message):
        resultsFileName = self.additionalInfoLog
            
        results = open(resultsFileName, "a+")
        results.write(message+"\n")
        results.close()
    
    def writeAnionPiResults(self, ligand, centroid, extractedAtoms, modelIndex, resolution, method, structureType ):
        """
        Zapisz dane do pliku z wynikami
        """
        resultsFileName = self.anionPiLog
        ligandCode = ligand.get_resname()
        ligandId = str(ligand.get_id()[1])
        ligandChain = ligand.get_parent().get_id()
        resultsFile = open(resultsFileName, "a+")
        newAtoms = []
        for atomData in extractedAtoms:
            distance = atomDistanceFromCentroid( atomData["Atom"], centroid )
            angle = atomAngleNomVecCentroid( atomData["Atom"], centroid )
            
            h = abs(cos(radians( angle ))*distance)
            x = sin(radians( angle ))*distance
            if angle > 90.0 :
                angle = 180 - angle
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
            resultsFile.write(self.pdbCode+"\t")
            resultsFile.write(ligandCode+"\t")
            resultsFile.write(ligandChain+"\t")
            resultsFile.write(ligandId+"\t")
            resultsFile.write(residueName+"\t")
            resultsFile.write(anionChain+"\t")
            resultsFile.write(anionId+"\t")
            resultsFile.write(atomData["AnionType"]+"\t")
            resultsFile.write(str(atomData["AnionId"])+"\t")
            resultsFile.write(atomData["Atom"].element+"\t")
            
            resultsFile.write(str(distance)+"\t")
            resultsFile.write(str(angle)+"\t")
            
            resultsFile.write(str(x)+"\t")
            resultsFile.write(str(h)+"\t")
            
            resultsFile.write(str(centroid["cycleId"])+"\t")
            resultsFile.write(str(centroidCoords[0])+"\t")
            resultsFile.write(str(centroidCoords[1])+"\t")
            resultsFile.write(str(centroidCoords[2])+"\t")
            
            resultsFile.write(str(atomCoords[0])+"\t")
            resultsFile.write(str(atomCoords[1])+"\t")
            resultsFile.write(str(atomCoords[2])+"\t")
            
            resultsFile.write(str(modelIndex)+"\t")
            resultsFile.write(str(atomData["Atom"].get_parent().is_disordered()) + "\t")
            resultsFile.write(str(centroid["ringSize"])+"\t")
            resultsFile.write(str(centroid["ringElements"])+"\t")
            resultsFile.write(str(resolution)+"\t")
            resultsFile.write(str(method)+"\t")
            resultsFile.write(str(structureType)+"\n")
        
        resultsFile.close()
        
        return newAtoms

    def writeAnionPiPlanarResults( self, ligand, centroid, planeData, modelIndex, anionGroupId):
        """
        Zapisz naglowki do pliku z wynikami:
        """
        resultsFileName = self.planarAnionPiLog
            
        ligandCode = ligand.get_resname()
        ligandId = str(ligand.get_id()[1])
        ligandChain = ligand.get_parent().get_id()
        resultsFile = open(resultsFileName, "a+")
        
        normVec = getNormVec( planeData.atomsInvolved , list(range( len(planeData.atomsInvolved) )) )
        
        inner_prod = np.inner( normVec, centroid["normVec"] )
        if abs(inner_prod) > 1.0:
            if abs(inner_prod) < 1.1:
                angle = 0    
            else:
                angle = 666
        else:
            angle = degrees( acos(inner_prod) )
            
        if angle > 90.0 :
            angle = 180 - angle
        
        atom = planeData.atomsInvolved[0]
        anion = atom.get_parent()
        residueName = anion.get_resname()
        anionChain = anion.get_parent().get_id()
        anionId = str(anion.get_id()[1])
        centroidCoords = centroid["coords"]  
        anionGroupCoords = getAverageCoords( planeData.atomsInvolved , list(range( len(planeData.atomsInvolved) )) )
        
        directionalAngle = calcDirectionalVector(planeData, centroid)
        
        resultsFile.write(self.pdbCode+"\t")
        resultsFile.write(ligandCode+"\t")
        resultsFile.write(ligandChain+"\t")
        resultsFile.write(ligandId+"\t")
        resultsFile.write(str(centroid["cycleId"])+"\t")
        resultsFile.write(residueName+"\t")
        resultsFile.write(anionChain+"\t")
        resultsFile.write(anionId+"\t")
        resultsFile.write(str(anionGroupId)+"\t")
        resultsFile.write(str(angle)+"\t")
        resultsFile.write(str(directionalAngle)+"\t")
        
        resultsFile.write(str(centroidCoords[0])+"\t")
        resultsFile.write(str(centroidCoords[1])+"\t")
        resultsFile.write(str(centroidCoords[2])+"\t")
        
        resultsFile.write(str(anionGroupCoords[0])+"\t")
        resultsFile.write(str(anionGroupCoords[1])+"\t")
        resultsFile.write(str(anionGroupCoords[2])+"\t")
        
        resultsFile.write(str(modelIndex)+"\n")
        resultsFile.close()
    
    def writeAnionPiLinearResults( self, ligand, centroid, lineData, modelIndex, anionGroupId, symmetrizeAlpha = False ):
        resultsFileName = self.linearAnionPiLog
            
        ligandCode = ligand.get_resname()
        ligandId = str(ligand.get_id()[1])
        ligandChain = ligand.get_parent().get_id()
        resultsFile = open(resultsFileName, "a+")
        
        
        vector =  lineData.atomsInvolved[1].get_coord() - lineData.atomsInvolved[0].get_coord() 
        vector = normalize(vector)
        
        inner_prod = np.inner( vector, centroid["normVec"] )
        if abs(inner_prod) > 1.0:
            if abs(inner_prod) < 1.1:
                angle = 0    
            else:
                angle = 666
        else:
            angle = degrees( acos(inner_prod) )
            
        if symmetrizeAlpha and angle > 90.0:
            angle = 180 - angle
        
        atom = lineData.atomsInvolved[0]
        anion = atom.get_parent()
        residueName = anion.get_resname()
        anionChain = anion.get_parent().get_id()
        anionId = str(anion.get_id()[1])
        centroidCoords = centroid["coords"]  
        anionGroupCoords = getAverageCoords( lineData.atomsInvolved , list(range( len(lineData.atomsInvolved) )) )
        
        resultsFile.write(self.pdbCode+"\t")
        resultsFile.write(ligandCode+"\t")
        resultsFile.write(ligandChain+"\t")
        resultsFile.write(ligandId+"\t")
        resultsFile.write(str(centroid["cycleId"])+"\t")
        resultsFile.write(residueName+"\t")
        resultsFile.write(anionChain+"\t")
        resultsFile.write(anionId+"\t")
        resultsFile.write(str(anionGroupId)+"\t")
        resultsFile.write(str(angle)+"\t")
        
        resultsFile.write(str(centroidCoords[0])+"\t")
        resultsFile.write(str(centroidCoords[1])+"\t")
        resultsFile.write(str(centroidCoords[2])+"\t")
        
        resultsFile.write(str(anionGroupCoords[0])+"\t")
        resultsFile.write(str(anionGroupCoords[1])+"\t")
        resultsFile.write(str(anionGroupCoords[2])+"\t")
        
        resultsFile.write(str(modelIndex)+"\n")
        resultsFile.close()
            
    def writeCationPiResults(self, ligand, centroid, extractedAtoms, cationRingChainLens , modelIndex ):
        """
        Zapisz dane do pliku z wynikami
        """
        if not extractedAtoms:
            return
        
        resultsFileName = self.cationPiLog
        ligandCode = ligand.get_resname()
        ligandId = str(ligand.get_id()[1])
        ligandChain = ligand.get_parent().get_id()
        resultsFile = open(resultsFileName, "a+")
        
        for atom, chainLen in zip(extractedAtoms , cationRingChainLens):
            distance = atomDistanceFromCentroid( atom, centroid )
            angle = atomAngleNomVecCentroid( atom, centroid )
            
            h = abs(cos(radians( angle ))*distance)
            x = sin(radians( angle ))*distance
            if angle > 90.0 :
                angle = 180 - angle
                
            atomCoords = atom.get_coord()
            centroidCoords = centroid["coords"]        
            
            cation = atom.get_parent()
            residueName = cation.get_resname()
            cationChain = cation.get_parent().get_id()
            cationId = str(cation.get_id()[1])
            resultsFile.write(self.pdbCode+"\t")
            resultsFile.write(ligandCode+"\t")
            resultsFile.write(ligandChain+"\t")
            resultsFile.write(ligandId+"\t")
            resultsFile.write(residueName+"\t")
            resultsFile.write(cationChain+"\t")
            resultsFile.write(cationId+"\t")
            resultsFile.write(atom.element+"\t")
            
            resultsFile.write(str(distance)+"\t")
            resultsFile.write(str(angle)+"\t")
            
            resultsFile.write(str(x)+"\t")
            resultsFile.write(str(h)+"\t")
            
            resultsFile.write(str(chainLen[0])+"\t")
            resultsFile.write(str(chainLen[1])+"\t")
            resultsFile.write(str(chainLen[2])+"\t")
            
            resultsFile.write(str(centroid["cycleId"])+"\t")
            resultsFile.write(str(centroidCoords[0])+"\t")
            resultsFile.write(str(centroidCoords[1])+"\t")
            resultsFile.write(str(centroidCoords[2])+"\t")
            
            resultsFile.write(str(atomCoords[0])+"\t")
            resultsFile.write(str(atomCoords[1])+"\t")
            resultsFile.write(str(atomCoords[2])+"\t")
            
            resultsFile.write(str(modelIndex)+"\n")
        
        resultsFile.close()

    def writeMetalLigandResults( self,  extractedAtoms, complexData , modelIndex ):
        """
        Zapisz dane do pliku z wynikami
        """
        if not extractedAtoms:
            return
        
        resultsFileName = self.metalLigandLog
        resultsFile = open(resultsFileName, "a+")
        
        for atom, compData in zip(extractedAtoms , complexData):
            if not compData["ligands"]:
                continue
            
            cation = atom.get_parent()
            cationResidueName = cation.get_resname()
            cationChain = cation.get_parent().get_id()
            cationId = str(cation.get_id()[1])
            cationCoords = atom.get_coord()
            
            for ligandData in compData["ligands"]:
                distance = atom - ligandData["atom"]
                
                ligand = ligandData["atom"].get_parent()
                ligandResidueName = ligand.get_resname()
                ligandChain = ligand.get_parent().get_id()
                ligandId = str(ligand.get_id()[1])
                
                resultsFile.write(self.pdbCode+"\t")
                
                resultsFile.write(cationResidueName+"\t")
                resultsFile.write(cationChain+"\t")
                resultsFile.write(cationId+"\t")
                
                resultsFile.write(ligandResidueName+"\t")
                resultsFile.write(ligandChain+"\t")
                resultsFile.write(ligandId+"\t")
                resultsFile.write(str(ligandData["AnionId"])+"\t")
                
                resultsFile.write(atom.element+"\t")
                resultsFile.write(ligandData["atom"].element+"\t")
                
                resultsFile.write(str(ligandData["isAnion"])+"\t")
                resultsFile.write(ligandData["anionType"]+"\t")
                
                resultsFile.write(str(distance)+"\t")
                
                ligandCoords = ligandData["atom"].get_coord()
                
                resultsFile.write(str(cationCoords[0])+"\t")
                resultsFile.write(str(cationCoords[1])+"\t")
                resultsFile.write(str(cationCoords[2])+"\t")
                
                resultsFile.write(str(ligandCoords[0])+"\t")
                resultsFile.write(str(ligandCoords[1])+"\t")
                resultsFile.write(str(ligandCoords[2])+"\t")
                
                resultsFile.write(str(compData["complex"])+"\t")
                resultsFile.write(str(compData["summary"])+"\t")
                resultsFile.write(str(compData["coordNo"])+"\t")
                
                resultsFile.write(str(modelIndex)+"\n")
        
        resultsFile.close()
    
    def writePiPiResults(self, ligand, centroid, extractedRes, extractedCentroids, modelIndex ):
        """
        Zapisz dane do pliku z wynikami
        """
        resultsFileName = self.piPiLog
        ligandCode = ligand.get_resname()
        ligandId = str(ligand.get_id()[1])
        ligandChain = ligand.get_parent().get_id()
        resultsFile = open(resultsFileName, "a+")
        newAtoms = []
        for res, cent in zip(extractedRes, extractedCentroids):
            distance = cent["distance"]
            theta = angleBetweenNormVec(centroid, cent)
            angle = angleNormVecPoint(centroid, cent["coords"])
            omega = calcOmega(centroid, cent)
            
            h = abs(cos(radians( angle ))*distance)
            x = sin(radians( angle ))*distance
            
            if angle > 90.0 :
                angle = 180 - angle
                
            if theta > 90.0 :
                theta = 180 - theta
                
            if omega > 90.0:
                omega = 180 - omega
                
            centroid2Coords = cent["coords"]
            centroidCoords = centroid["coords"]        
            
            residueName = res.get_resname()
            resChain = res.get_parent().get_id()
            resId = str(res.get_id()[1])
            resultsFile.write(self.pdbCode+"\t")
            resultsFile.write(ligandCode+"\t")
            resultsFile.write(ligandChain+"\t")
            resultsFile.write(ligandId+"\t")
            resultsFile.write(residueName+"\t")
            resultsFile.write(resChain+"\t")
            resultsFile.write(resId+"\t")
            
            resultsFile.write(str(distance)+"\t")
            resultsFile.write(str(angle)+"\t")
            
            resultsFile.write(str(x)+"\t")
            resultsFile.write(str(h)+"\t")
            
            resultsFile.write(str(theta)+"\t")
            resultsFile.write(str(omega)+"\t")
            
            resultsFile.write(str(centroid["cycleId"])+"\t")
            resultsFile.write(str(centroidCoords[0])+"\t")
            resultsFile.write(str(centroidCoords[1])+"\t")
            resultsFile.write(str(centroidCoords[2])+"\t")
            
            resultsFile.write(str(centroid2Coords[0])+"\t")
            resultsFile.write(str(centroid2Coords[1])+"\t")
            resultsFile.write(str(centroid2Coords[2])+"\t")
            
            resultsFile.write(str(modelIndex)+"\t")
            resultsFile.write(str(cent["ringSize"])+"\n")
        
        resultsFile.close()
        
        return newAtoms

    def writeAnionCationResults(self, anionAtom, ligand, centroid, extractedCations, modelIndex):
        """
        Zapisz dane do pliku z wynikami
        """
        resultsFileName = self.anionCationLog
        anion = anionAtom.get_parent()
        anionCode = anion.get_resname()
        anionId = str(anion.get_id()[1])
        anionChain = anion.get_parent().get_id()
        anionCoord = anionAtom.get_coord()
        
        ligandCode = ligand.get_resname()
        ligandId = str(ligand.get_id()[1])
        ligandChain = ligand.get_parent().get_id()
        
        resultsFile = open(resultsFileName, "a+")
        newAtoms = []
        
        anionAngle = atomAngleNomVecCentroid( anionAtom, centroid )
        firstSemisphereAnion = anionAngle < 90
    
        for cat in extractedCations:
            distance = anionAtom - cat
            
            cationAngle = atomAngleNomVecCentroid(cat, centroid)
            firstSemisphereCation = cationAngle < 90
            sameSemisphere = firstSemisphereAnion == firstSemisphereCation
            
            catCoord = cat.get_coord()
            catRes = cat.get_parent()
            
            residueName = catRes.get_resname()
            resChain = catRes.get_parent().get_id()
            resId = str(catRes.get_id()[1])
            resultsFile.write(self.pdbCode+"\t")
            resultsFile.write(residueName+"\t")
            resultsFile.write(resChain+"\t")
            resultsFile.write(resId+"\t")
            
            resultsFile.write(ligandCode+"\t")
            resultsFile.write(ligandChain+"\t")
            resultsFile.write(ligandId+"\t")
            resultsFile.write(str(centroid["cycleId"])+"\t")
            
            resultsFile.write(anionCode+"\t")
            resultsFile.write(anionChain+"\t")
            resultsFile.write(anionId+"\t")
            resultsFile.write(str(anionAtom.anionData.anionId)+"\t")
            
            resultsFile.write(anionAtom.element+"\t")
            resultsFile.write(cat.element+"\t")
            
            resultsFile.write(str(distance)+"\t")
            
            resultsFile.write(str(anionCoord[0])+"\t")
            resultsFile.write(str(anionCoord[1])+"\t")
            resultsFile.write(str(anionCoord[2])+"\t")
            
            resultsFile.write(str(catCoord[0])+"\t")
            resultsFile.write(str(catCoord[1])+"\t")
            resultsFile.write(str(catCoord[2])+"\t")
            
            resultsFile.write(str(sameSemisphere)+"\t")
            resultsFile.write(str(abs(anionAngle - cationAngle))+"\t")
            resultsFile.write(str(modelIndex)+"\n")
        
        resultsFile.close()
        
        return newAtoms

    def writeHbondsResults(self ,hDonors, atom, modelIndex):
        resultsFileName = self.hBondsLog
        
        anionAtom = atom["Atom"]
        anion = anionAtom.get_parent()
        anionCode = anion.get_resname()
        anionId = str(anion.get_id()[1])
        anionChain = anion.get_parent().get_id()
        anionCoord = anionAtom.get_coord()
        
        resultsFile = open(resultsFileName, "a+")
    
        for hDonData in hDonors:
            hDon = hDonData["donor"]
            distance = anionAtom - hDon
            
            hDonCoords = hDon.get_coord()
            hDonRes = hDon.get_parent()
            
            residueName = hDonRes.get_resname()
            resChain = hDonRes.get_parent().get_id()
            resId = str(hDonRes.get_id()[1])
            
            resultsFile.write(self.pdbCode+"\t")
            
            resultsFile.write(anionCode+"\t")
            resultsFile.write(anionChain+"\t")
            resultsFile.write(anionId+"\t")
            
            resultsFile.write(residueName+"\t")
            resultsFile.write(resChain+"\t")
            resultsFile.write(resId+"\t")
            
            resultsFile.write(atom["AnionType"]+"\t")
            resultsFile.write(anionAtom.element+"\t")
            resultsFile.write(str(anionAtom.anionData.anionId)+"\t")
            
            resultsFile.write(str(anionCoord[0])+"\t")
            resultsFile.write(str(anionCoord[1])+"\t")
            resultsFile.write(str(anionCoord[2])+"\t")
            
            resultsFile.write(hDon.element+"\t"+hDon.element+"\t")
            
            resultsFile.write(str(hDonCoords[0])+"\t")
            resultsFile.write(str(hDonCoords[1])+"\t")
            resultsFile.write(str(hDonCoords[2])+"\t")
            
            hydrogenAtom = hDonData["hydrogen"]
            hydrogenCoords = hydrogenAtom.get_coord()
            
            resultsFile.write(str(hydrogenCoords[0])+"\t")
            resultsFile.write(str(hydrogenCoords[1])+"\t")
            resultsFile.write(str(hydrogenCoords[2])+"\t")
            
            resultsFile.write(str(hDonData["HFromExp"])+"\t")
            
            vec1 = normalize(anionCoord - hydrogenCoords)
            vec2 = normalize(hDonCoords - hydrogenCoords)
            angle = degrees( acos( np.inner(vec1, vec2 ) ))
            hDistance = anionAtom - hydrogenAtom
            
            resultsFile.write(str(angle)+"\t")
            resultsFile.write(str(hDistance)+"\t")
            
            resultsFile.write(str(distance)+"\t")
            
            resultsFile.write(str(modelIndex)+"\n")
        
        resultsFile.close()
    
    def incrementPartialProgress(self ):
        fileName = self.partialProgressLog
        if not isfile(fileName):
            partialProgressFile = open(fileName, 'w')
            partialProgressFile.write("1")
            partialProgressFile.close()
            
        else:
            partialProgressFile = open(fileName, 'r')
            actualNo = int(partialProgressFile.readline())
            partialProgressFile.close()
            
            actualNo += 1
            partialProgressFile = open(fileName, 'w')
            partialProgressFile.write(str(actualNo))
            partialProgressFile.close()

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

def calcOmega( piAcidCentroid, piResCentroid ):
    cent1cent2Vec = np.array(piResCentroid["coords"]) - np.array(piAcidCentroid["coords"])
    cent1cent2Vec = normalize(cent1cent2Vec)
    
    planeNormVec = np.cross(cent1cent2Vec, piAcidCentroid["normVec"])
    planeNormVec = normalize(planeNormVec)
    
    return degrees( acos( np.inner(planeNormVec, piResCentroid["normVec"] ) ) )
    

def angleBetweenNormVec( centroid1, centroid2):
    normVec1 = centroid1["normVec"]
    normVec2 = centroid2["normVec"]
    
    inner_prod = np.inner( normVec1, normVec2 )
    if abs(inner_prod) > 1.0:
        if abs(inner_prod) < 1.1:
            return 0.0
        else:
            return 666.0
    
    return degrees( acos(inner_prod) )

def angleNormVecPoint( centroid, point):
    coords = np.array(point)
    centroidCoords = np.array(centroid["coords"])
    normVec = centroid["normVec"]
    
    centrAtomVec= normalize( coords - centroidCoords )
    inner_prod = np.inner( normVec, centrAtomVec )
    
    return degrees( acos(inner_prod) )
        
def calcDirectionalVector(planeData, centroid):
    vecCoords = []
    for pointData in planeData.directionalVector:
        kind = list(pointData.keys())[0]
        if kind == "atom":
            vecCoords.append( pointData[kind].get_coord() )
        elif kind == "center":
            atomsList = pointData[kind]
            vecCoords.append( getAverageCoords( atomsList, list( range(len(atomsList)))) )
        elif kind == "closest":
            atomsList = pointData[kind]
            
            minDist = 1000
            closestAtom = None
            
            for atom in atomsList:
                dist = atomDistanceFromCentroid(atom, centroid)
                
                if dist < minDist:
                    minDist = dist
                    closestAtom = atom
                    
            vecCoords.append(closestAtom.get_coord())

            
    vec = normalize( vecCoords[1] - vecCoords[0] )
    vec2 = normalize( np.array(centroid["coords"]) - vecCoords[1] )
    
    inner_prod = np.inner( vec , vec2)
    if abs(inner_prod) > 1.0:
        if abs(inner_prod) < 1.1:
            return 0.0
        else:
            return 666.0
        
    return degrees( acos(inner_prod) )
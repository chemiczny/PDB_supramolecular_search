# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 13:56:41 2017

Skrypt powstał we wspolpracy z najpiekniejsza kobieta na swiecie, Zrodzona
w Olkuszu, Pania mojego serca Emilia Kuzniak.

Sluzy do analizy pojedynczego pliku cif celem znalezienia potencjalnych struktur
zawierajacych oddzialywania supramolekularne anion-pi oraz obliczenie
wielkosci geometrycznych tychze czasteczek na potrzeby dalszej analizy
"""
#fast ang ugly code especially for prometheus
import sys
from os.path import isdir, getsize

if isdir("/net/archive/groups/plggsuprm/pythonPackages") and not "/net/archive/groups/plggsuprm/pythonPackages" in sys.path :
    sys.path.insert(0, "/net/archive/groups/plggsuprm/pythonPackages" )

from Bio.PDB import FastMMCIFParser, NeighborSearch, Selection
from primitiveCif2Dict import primitiveCif2Dict
import numpy as np                
from supramolecularLogging import writeAnionPiResults, incrementPartialProgress, writeAdditionalInfo
from supramolecularLogging import writeCationPiResults, writePiPiResults, writeAnionCationResults, writeHbondsResults, writeMetalLigandResults
from ringDetection import getRingsCentroids, findInGraph, isFlatPrimitive, normalize, molecule2graph
from protonate import Protonate
from anionRecogniser import AnionRecogniser, createResId
from multiprocessing import current_process
import networkx as nx
from time import time

def findSupramolecular( cifData):
    """
    Przeanalizuj pojedynczy plik cif pod katem oddzialywan suporamolekularnych
    zwiazanych z konkretnym ligandem
    
    Wejscie:
    ligandCode - kod liganda
    cifFile    - plik cif pobrany z bazy PDB
    
    Wyjscie:
    Hehehe, czas pokaze...
    """
    cifFile = cifData[0]
    PDBcode = cifData[1]
    logId = cifData[2]
    
    cifAnalyser = CifAnalyser(cifFile, PDBcode, logId)
    cifAnalyser.analyseCif()
    

class CifAnalyser:
    def __init__(self, cifFile, PDBcode, logId = "default"):
        self.cifFile  = cifFile
        self.PDBcode = PDBcode
        
        if logId == "default":
            self.fileId = current_process()
            self.fileId = str(self.fileId)
            self.fileId = self.fileId.replace(" ", "")
            self.fileId = self.fileId.replace("(", "")
            self.fileId = self.fileId.replace(")", "")
            self.fileId = self.fileId.replace(",", "")
            self.fileId = self.fileId.replace("<", "")
            self.fileId = self.fileId.replace(">", "")
        else:
            self.fileId = logId
    
        self.ns = None
        self.resolution = None
        self.method = None
        self.hAtomsPresent = False
        
        self.anionRecogniser = AnionRecogniser()
    
    def analyseCif(self):
        parser = FastMMCIFParser(QUIET=True)
        
        timeStart = time()
        try:
            time1 = time()
            structure = parser.get_structure('temp', self.cifFile)
            timeBioParsing = time()-time1
        except:
            errorMessage = "FastMMCIFParser cannot handle with: " + self.cifFile
    #        fileId = current_process()
            incrementPartialProgress(self.fileId)
            writeAdditionalInfo(errorMessage, self.fileId)
            #Zeby zobaczyc co sie dzieje
            return True        
            
        writeAdditionalInfo("Zaczynam analize: "+self.PDBcode, self.fileId)
        writeAdditionalInfo("Rozmiar pliku: " + str(getsize(self.cifFile)), self.fileId)
        
    #    supramolecularFound = False 
        notPiacids = [ "HOH", "DOD", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "ILE", "LEU", "LYS", "MET", "PRO", "SER", "THR", "VAL" ] 
            
        time2 = time()
        self.resolution, self.method = self.readResolutionAndMethod()
        timeResolutionReading = time() - time2
        
        self.hAtomsPresent = hydrogensPresent(structure)
        
        times = { "getCentroidTime" : 0,
            "extractAnionsAroundRingTime" : 0,
            "extractCationsTime" : 0,
            "cationAnalysisTime" : 0,
            "extractPiPiTime" : 0,
            "extractHBondsTime" :  0,
            "timeBioParsing" : timeBioParsing,
            "timeResolutionReading" : timeResolutionReading}
        
    #    try:
        for modelIndex, model in enumerate(structure):
    #        print("model "+str(modelIndex))
            atoms = Selection.unfold_entities(model, 'A')  
            not_disordered_atoms = []
            for atom in atoms:
                if not atom.is_disordered() or atom.get_altloc() == 'A':
    
                    not_disordered_atoms.append(atom)
                
            if len(not_disordered_atoms) == 0:
                continue
    #        print("only atoms A selected")
            self.ns = NeighborSearch(not_disordered_atoms)
    #        print("neighbor search utworzony")
            for residue in model.get_residues():
                residueName = residue.get_resname().upper()
                if not residueName in notPiacids  :
                    supramolecularFound , newTimes = self.analysePiacid(residue, modelIndex, structure)
                    for key in newTimes:
                        times[key] += newTimes[key]
                
    #    fileId = current_process()
        incrementPartialProgress(self.fileId)
    #    except:
    #        fileId = current_process()
    #        writeAdditionalInfo("UNKNOWN ERROR!!!! PDB: "+PDBcode, fileId)
        timeStop = time()
        timeTaken = timeStop - timeStart
        
        for key in times:
            writeAdditionalInfo("time "+key+" "+str(times[key]), self.fileId)
        
        writeAdditionalInfo("Analiza skonczona "+self.PDBcode+ " czas: "+str(timeTaken), self.fileId)
#    return supramolecularFound

    def readResolutionAndMethod( self ):
        try:
            mmcif_dict_parser = primitiveCif2Dict(self.cifFile, ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" , "_exptl.method"] )
            mmcif_dict = mmcif_dict_parser.result
        except:
            errorMessage = "primitiveCif2Dict cannot handle with: "+self.cifFile
            writeAdditionalInfo(errorMessage, self.fileId)
            return -666, -666
        
        method = "Unknown"
        if "_exptl.method" in mmcif_dict:
            method = mmcif_dict["_exptl.method"]
            if len(method) == 1:
                method = method[0]
            else:
                method = sorted(method)
    #    res_keys = ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" ]
        res_key = "_refine.ls_d_res_high"
        resolution = []
        if res_key in mmcif_dict:
            resolution = mmcif_dict[res_key]
    
                    
        if not resolution:
            return -1, method
                    
        if len(resolution) == 1  :
            if isfloat( resolution[0] ):
                return resolution[0], method
            else:
                return -1, method
        else:
            resFloats = []
            for res in resolution:
                if isfloat(res):
                    resFloats.append(float(res))
                    
            if not resFloats:
                return -1, method
            
            resFloats = sorted(resFloats)
            return resFloats[0], method

    def analysePiacid(self, ligand, modelIndex, structure):
    #    print("analysePiacid - start ")
        firstAtom = list(ligand.get_atoms())[0]
        if firstAtom.is_disordered() and firstAtom.get_altloc() != 'A':
            return False, {}
        
        getCentroidTime = 0
        extractAnionsAroundRingTime = 0
        extractCationsTime = 0
        cationAnalysisTime = 0
        extractPiPiTime = 0
        extractHBondsTime = 0
        
        time0 = time()
        centroids, ligandGraph = getRingsCentroids( ligand, True )
        getCentroidTime += time() - time0
        
        ligandWithAnions = False
        
        for centroid in centroids:
            bigDistance =12
            distance = 4.5
            
            atoms = self.ns.search(np.array(centroid["coords"]), bigDistance, 'A')
            
            nsSmall = NeighborSearch(atoms)
            neighbors = nsSmall.search(np.array(centroid["coords"]), distance, 'A')
            
            time1 = time()
            extractedAnionAtoms = self.anionRecogniser.extractAnionAtoms( neighbors, ligand, nsSmall )
            extractAnionsAroundRingTime += time()-time1
    
            if len(extractedAnionAtoms) > 0:
                time2 = time()
                extractedMetalCations = extractMetalCations( centroid["coords"], nsSmall, 10 )
                extractedAAcations = extractAACations( centroid["coords"], nsSmall, 4.5 )
                extractCationsTime += time() - time2
                
                cationRingLenChains = []
                cationComplexData = []
                
                time3 = time()
                for cat in extractedMetalCations:
                    cationRingLenChains.append( findChainLenCationRing(cat, ligand, centroid , self.ns, ligandGraph, self.fileId) )
                    cationComplexData.append( self.findCationComplex( cat, self.ns, ligand ) )
                cationAnalysisTime += time() - time3
                        
                extractedCations = extractedMetalCations + extractedAAcations
                cationRingLenChainsFull = cationRingLenChains + len(extractedAAcations)* [ (0, False) ]
                
                writeCationPiResults(ligand, self.PDBcode, centroid, extractedCations, cationRingLenChainsFull, modelIndex, self.fileId )
                writeMetalLigandResults( self.PDBcode , extractedMetalCations, cationComplexData, modelIndex, self.fileId)
                
                time4 = time()
                extractedCentroids, ringMolecules = extractRingCentroids(centroid["coords"], ligand, nsSmall)
                extractPiPiTime += time() - time4
                
                writePiPiResults(ligand, self.PDBcode, centroid, ringMolecules, extractedCentroids, modelIndex, self.fileId)
                
                for atom in extractedAnionAtoms:
                    writeAnionCationResults(atom["Atom"], self.PDBcode, ligand, centroid, extractedCations, modelIndex, self.fileId)
                    
                    time5 = time()
                    hDonors = extractHbonds( atom , nsSmall, 3.5, self.hAtomsPresent, self.fileId, structure)
                    extractHBondsTime += time() - time5
                    
                    writeHbondsResults( self.PDBcode,hDonors, atom, modelIndex, self.fileId)
                    
            extractedAtoms =  writeAnionPiResults(ligand, self.PDBcode, centroid, extractedAnionAtoms, modelIndex, self.resolution, self.method, self.fileId)
    
            if len(extractedAtoms) > 0:
                ligandWithAnions = True
                
        times = { "getCentroidTime" : getCentroidTime,
            "extractAnionsAroundRingTime" : extractAnionsAroundRingTime,
            "extractCationsTime" : extractCationsTime,
            "cationAnalysisTime" : cationAnalysisTime,
            "extractPiPiTime" : extractPiPiTime,
            "extractHBondsTime" :  extractHBondsTime }
            
        return ligandWithAnions, times

    def findCationComplex(self, cation, ns, ligand):
        potentialLigands = ns.search(cation.get_coord(), 2.6 , 'A')
        anionSpaceWithCation = ns.search(cation.get_coord(), 4.5 , 'A') 
        metals =  [
            "LI", "BE", 
            "NA", "MG", "AL", 
            "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", 
            "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", 
            "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", 
            "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG"]
        
        anionSpace = []
        for a in anionSpaceWithCation:
            if not a.element.upper() in metals:
                anionSpace.append(a)
                
        if len(anionSpace) == 0:
            return  { "complex" : False, "coordNo" : 0, "ligands" : [] }
        
        nsSmall = NeighborSearch(anionSpace)
    #    if not potentialLigands:
    #        return  { "complex" : False, "coordNo" : 0 }
        
        resultantVector = np.array([0.0 , 0.0, 0.0])
        coordNo = 0
        ligands = []
        for atom in potentialLigands:
            if atom.element == "H":
                continue
            if atom != cation:
                newVector = normalize( atom.get_coord() - cation.get_coord() )
                resultantVector += newVector
                coordNo += 1
                isAnion, anionType = self.anionRecogniser.searchInAnionTemplates(atom, anionSpace, nsSmall)
                anionId = -1
                if isAnion:
                    anionId = atom.anionData.anionId
                    
                ligands.append({ "isAnion" : isAnion, "anionType" : anionType, "atom" :atom, "AnionId" : anionId })
                
        if coordNo == 0 :
            return  { "complex" : False, "coordNo" : 0, "ligands" : [] }
                
        vectorLen = np.linalg.norm(resultantVector)
         
        if vectorLen < 0.2:
             return { "complex" : True, "coordNo" : coordNo , "ligands" : ligands}
        else:
             return  { "complex" : False, "coordNo" : coordNo , "ligands" : ligands }

def extractMetalCations ( point,  ns, distance  ):
    neighbors = ns.search(np.array(point), distance, 'A')
    
    metalCationsFound = []
    
    metalCations =  [
        "Li", "Be", 
        "Na", "Mg", "Al", 
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",  
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", 
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", 
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]
    
    for atom in neighbors:
        if atom.element.upper() in [ element.upper() for element in metalCations ]:
            metalCationsFound.append( atom )
    
    return metalCationsFound

def extractAACations ( point,  ns, distance  ):
    neighbors = ns.search(np.array(point), distance, 'A')
    
    aaCationsFound = []
    
    for atom in neighbors:
            
        if atom.get_parent().get_resname().upper() in [ "ARG" , "LYS" ] and atom.get_name() != "N" and atom.element.upper() == "N":
            aaCationsFound.append(atom)
    
    return  aaCationsFound

def extractHbonds( atom , nsSmall, distance, hAtomsPresent, fileId, structure):
    
    neighbors = nsSmall.search( np.array(atom["Atom"].get_coord()) , distance + 2, 'A' )
    acceptorResidue = atom["Atom"].get_parent()
    
    if not hAtomsPresent:       
        protonationWorker = Protonate()
        protonationWorker.protonate(neighbors, atom["Atom"])
        graph = protonationWorker.moleculeGraph
        neighbors = protonationWorker.atomList
        anionAtomInd = protonationWorker.anionId
    else:
        graph, anionAtomInd = molecule2graph( neighbors, atom["Atom"], False, False )
    
    donors = []
    
    for potentialDonorInd in graph.nodes():
        element = neighbors[potentialDonorInd].element
        
        if element in [ "O" , "N" ] :
            connected = list(graph.neighbors(potentialDonorInd))
            
            if atom["Atom"] - neighbors[potentialDonorInd] > distance:
                continue
                
            if neighbors[potentialDonorInd].get_parent() == acceptorResidue:
                continue
            
            connectedElements = [ neighbors[c].element for c in connected ]
            connectedElements = list(set(connectedElements))
            
            for c in connected:
                element = neighbors[c].element
                if element != "H":
                    continue
                else:
                    donors.append( { "donor" : neighbors[potentialDonorInd], "hydrogen" : neighbors[c] , "HFromExp" : hAtomsPresent} )
            
        
    return donors
    
    

def findChainLenCationRing( cation, piAcid, centroidData, ns, ligandGraph, fileId ):
    piAcidAtoms = list(piAcid.get_atoms())
    
    cationNeighbours = ns.search( cation.get_coord(), 2.6, 'A' )
    firstAtomInRing = centroidData["cycleAtoms"][0]
    shortestPath = []
    somethingFound = False
#    print("Analizuje :", piAcid.get_resname())
    for catN in cationNeighbours:
        if catN.element == "H":
            continue
        if catN == cation:
            continue
        if catN.get_parent() == piAcid:
            tempG, catNIndex = findInGraph( ligandGraph, catN, piAcidAtoms )
            if catNIndex == None:
                print(piAcid.get_resname(), piAcid.get_id(), piAcid.get_parent().get_id())
                continue
            
            if not  nx.has_path(ligandGraph, firstAtomInRing, catNIndex):
                continue 
            
            newPath = nx.shortest_path(ligandGraph, firstAtomInRing, catNIndex)
#            print(catN.get_name())
#            print("przed obieciem: ", len(newPath))
            newPath = [ node for node in newPath if not node in centroidData["cycleAtoms"] ]
#            print("znalazlem sciezke dlugosci", len(newPath))
            if len(newPath) < len(shortestPath) or not somethingFound:
#                print("znalazlem krotsza sciezke")
#                print("stara dlugosc: ", len(shortestPath))
                shortestPath = newPath
                
            somethingFound = True
                
    if not somethingFound:
        return 0, False
    
    flat = True
    for node in shortestPath:
        neighbors = list(nx.neighbors(ligandGraph, node))
        if len( neighbors ) > 2:
            verdict = isFlatPrimitive(piAcidAtoms, neighbors + [ node ], 0.25)
            if not verdict["isFlat"]:
                writeAdditionalInfo( "Nieplaski lancuch! "+cation.element, fileId)
                flat = False
            
    return len(shortestPath)+1, flat
        

def extractRingCentroids(point, residue, ns):
    distAtom = 4.7
    distCent = 4.5
    neighbors = ns.search(np.array(point), distAtom, 'A')
    residues = Selection.unfold_entities(neighbors, 'R') 
    point = np.array(point)
    
    notAromatic = [ "HOH", "DOD", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "ILE", "LEU", "LYS", "MET", "PRO", "SER", "THR", "VAL" ] 
    
    res1id = createResId(residue)
    centroidsFound = []
    ringMolecules = []
    
    for res in residues:
        if res.get_resname() in notAromatic:
            continue
        
        if createResId(res) == res1id:
            continue
        
        centroids = getRingsCentroids(res)
        for centroid in centroids:
            dist = np.linalg.norm( point - np.array(centroid["coords"]) )
            if dist < distCent:
                centroid["distance"] = dist
                centroidsFound.append(centroid)
                ringMolecules.append(res)
                
    return centroidsFound, ringMolecules
    
def getResiduesListFromAtomData( atomDataList ):
    atomsList = []
    for atomData in atomDataList:
        atomsList.append(atomData["Atom"])
    return Selection.unfold_entities(atomsList, 'R')  
    
            
def anionScreening( atoms, ligprepData ):
    selectedAtoms = []
    anionsNames = ligprepData["anionNames"]
    
    for atomData in atoms:
        parentName = atomData["Atom"].get_parent().get_resname()        
        if parentName in anionsNames:
            selectedAtoms.append(atomData)
        
    return selectedAtoms
    
def hydrogensPresent(structure):
    for a in structure.get_atoms():
        if a.element == "H":
            return True
        
    return False
    
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

if __name__ == "__main__":
    pass
#    writeSupramolecularSearchHeader( )
#    timeStart = time.time()
##    findSupramolecularAnionPiLigand( "MCY", "cif/106d.cif", "106D" )
#    #findSupramolecularAnionPiLigand( "NCZ", "cif/1j5i.cif", "1J5I" )
#    findSupramolecularAnionPiAllLigands( "cif/1bp0.cif", "1BP0")
#    findSupramolecularAnionPiAllLigands( "cif/3bdj.cif", "3BDJ")
##    findSupramolecularAnionPiLigand( "7NC", "cif/5wqk.cif", "5WQK" )
##    findSupramolecularAnionPiLigand( "HPA", "cif/3nrz.cif", "3NRZ" )
##    findSupramolecularAnionPiLigand( "LUM", "cif/1he5.cif", "1HE5" )
##    findSupramolecularAnionPiLigand( "NAP", "cif/3bcj.cif", "3BCJ" )
##    findSupramolecularAnionPiLigand( "NAP", "cif/4lbs.cif", "4LBS" )
#    timeStop = time.time()
#    print("Calosc: ", timeStop-timeStart)
#    print(readResolutionAndMethod("cif/1bp0.cif"))
#    print(readResolutionAndMethod("cif/1j5i.cif"))
#    print(readResolutionAndMethod("cif/2eet.cif"))
#    print(readResolutionAndMethod("cif/4c3y.cif"))
#findSupramolecular( "LUM", 666 )

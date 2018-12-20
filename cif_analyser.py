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
from os.path import isdir
from os import system
#if isdir("/net/people/plgglanow/pythonPackages") and not "/net/people/plgglanow/pythonPackages" in sys.path :
#    sys.path.insert(0, "/net/people/plgglanow/pythonPackages" )
if isdir("/net/archive/groups/plggsuprm/pythonPackages") and not "/net/people/plgglanow/pythonPackages" in sys.path :
    sys.path.insert(0, "/net/people/plgglanow/pythonPackages" )

from Bio.PDB import FastMMCIFParser, NeighborSearch, Selection, PDBIO, PDBParser
from primitiveCif2Dict import primitiveCif2Dict
import numpy as np                
from supramolecularLogging import writeAnionPiResults, incrementPartialProgress, writeAdditionalInfo
from supramolecularLogging import writeCationPiResults, writePiPiResults, writeAnionCationResults, writeHbondsResults
from ringDetection import getRingsCentroids, findInGraph, isFlatPrimitive, normalize, molecule2graph
from anionRecogniser import extractAnionAtoms, createResId
from multiprocessing import current_process
import networkx as nx
from buildStructure import primitiveBuildStructure

#from time import time

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
    parser = FastMMCIFParser()
    
    if logId == "default":
        fileId = current_process()
        fileId = str(fileId)
        fileId = fileId.replace(" ", "")
        fileId = fileId.replace("(", "")
        fileId = fileId.replace(")", "")
        fileId = fileId.replace(",", "")
    else:
        fileId = logId
    
    try:
        structure = parser.get_structure('temp', cifFile)
    except:
        errorMessage = "FastMMCIFParser cannot handle with: " + cifFile
#        fileId = current_process()
        incrementPartialProgress(fileId)
        writeAdditionalInfo(errorMessage, fileId)
        #Zeby zobaczyc co sie dzieje
        return True        
        
    supramolecularFound = False 
    notPiacids = [ "HOH", "DOD", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "ILE", "LEU", "LYS", "MET", "PRO", "SER", "THR", "VAL" ] 
        
    resolution, method = readResolutionAndMethod(cifFile, fileId)
    hAtomsPresent = hydrogensPresent(structure)
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
        ns = NeighborSearch(not_disordered_atoms)
#        print("neighbor search utworzony")
        for residue in model.get_residues():
            residueName = residue.get_resname().upper()
            if not residueName in notPiacids  :
                if analysePiacid(residue, PDBcode, modelIndex, ns, resolution, method, fileId, hAtomsPresent):
                    supramolecularFound = True
            
#    fileId = current_process()
    incrementPartialProgress(fileId)
#    except:
#        fileId = current_process()
#        writeAdditionalInfo("UNKNOWN ERROR!!!! PDB: "+PDBcode, fileId)
    return supramolecularFound

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

def readResolutionAndMethod( cifFile, fileId ):
    try:
        mmcif_dict_parser = primitiveCif2Dict(cifFile, ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" , "_exptl.method"] )
        mmcif_dict = mmcif_dict_parser.result
    except:
        errorMessage = "primitiveCif2Dict cannot handle with: "+cifFile
        writeAdditionalInfo(errorMessage, fileId)
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

def analysePiacid(ligand, PDBcode, modelIndex, ns, resolution, method, fileId, hAtomsPresent):
#    print("analysePiacid - start ")
    firstAtom = list(ligand.get_atoms())[0]
    if firstAtom.is_disordered() and firstAtom.get_altloc() != 'A':
        return False
    
    
    centroids, ligandGraph = getRingsCentroids( ligand, True )
    
    ligandWithAnions = False
    
    for centroid in centroids:
        bigDistance =12
        distance = 4.5
        
        atoms = ns.search(np.array(centroid["coords"]), bigDistance, 'A')
        
        nsSmall = NeighborSearch(atoms)
        neighbors = nsSmall.search(np.array(centroid["coords"]), distance, 'A')
        
        extractedAnionAtoms = extractAnionAtoms( neighbors, ligand, nsSmall )

        if len(extractedAnionAtoms) > 0:
            extractedCationAtoms = extractCationAtoms( centroid["coords"], nsSmall, 10 )
            
            cationRingLenChains = []
            cationComplexData = []
            for cat in extractedCationAtoms:
                cationRingLenChains.append( findChainLenCationRing(cat, ligand, centroid , ns, ligandGraph, fileId) )
                cationComplexData.append( findCationComplex( cat, ns ) )
                    
            writeCationPiResults(ligand, PDBcode, centroid, extractedCationAtoms, cationRingLenChains, cationComplexData, modelIndex, fileId )
            
            extractedCentroids, ringMolecules = extractRingCentroids(centroid["coords"], ligand, nsSmall)
            writePiPiResults(ligand, PDBcode, centroid, ringMolecules, extractedCentroids, modelIndex, fileId)
            
            for atom in extractedAnionAtoms:
                cationNearAnion = extractCationAtoms( atom["Atom"].get_coord(), nsSmall, 4.5  )
                writeAnionCationResults(atom["Atom"], PDBcode, ligand, centroid, cationNearAnion, modelIndex, fileId)
                
                hDonors = extractHbonds( atom , nsSmall, 3.5, hAtomsPresent, fileId)
                writeHbondsResults( PDBcode,hDonors, atom, modelIndex, fileId)
                
            

        extractedAtoms =  writeAnionPiResults(ligand, PDBcode, centroid, extractedAnionAtoms, modelIndex, resolution, method, fileId)

        if len(extractedAtoms) > 0:
            ligandWithAnions = True
        
    return ligandWithAnions

def findCationComplex(cation, ns):
    potentialLigands = ns.search(cation.get_coord(), 2.6 , 'A')
    
#    if not potentialLigands:
#        return  { "complex" : False, "coordNo" : 0 }
    
    resultantVector = np.array([0.0 , 0.0, 0.0])
    coordNo = 0
    for atom in potentialLigands:
        if atom.element == "H":
            continue
        if atom != cation:
            newVector = normalize( atom.get_coord() - cation.get_coord() )
            resultantVector += newVector
            coordNo += 1
            
    if coordNo == 0 :
        return  { "complex" : False, "coordNo" : 0 }
            
    vectorLen = np.linalg.norm(resultantVector)
     
    if vectorLen < 0.2:
         return { "complex" : True, "coordNo" : coordNo }
    else:
         return  { "complex" : False, "coordNo" : coordNo }

def extractCationAtoms ( point,  ns, distance  ):
    neighbors = ns.search(np.array(point), distance, 'A')
    
    metalCationsFound = []
    metalCations =  [
        "Li", "Be", 
        "Na", "Mg", "Al", 
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", 
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", 
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", 
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]
    
    for atom in neighbors:
        if atom.element.upper() in [ element.upper() for element in metalCations ]:
            metalCationsFound.append( atom )
    
    return metalCationsFound

def extractHbonds( atom , nsSmall, distance, hAtomsPresent, fileId):
    
    neighbors = nsSmall.search( np.array(atom["Atom"].get_coord()) , distance + 2, 'A' )
    acceptorResidue = atom["Atom"].get_parent()
    
    if not hAtomsPresent:
#        return []
        structure = primitiveBuildStructure(neighbors, acceptorResidue)
        io = PDBIO()
        io.set_structure(structure)
        pdbFileName = "hBondsScr/"+str(fileId)+".pdb"
        io.save(pdbFileName)
        system("python primitiveAddHydrogens.py "+pdbFileName)
        parser = PDBParser()
        structure = parser.get_structure('hBondTemp', pdbFileName)
        structureAtoms = structure.get_atoms()
        hydrogensAtoms = [ ]
        for sa in structureAtoms:
            if sa.element == "H":
                hydrogensAtoms.append(sa)
        
        neighbors += hydrogensAtoms
    
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
    for catN in cationNeighbours:
        if catN.element == "H":
            continue
        if catN.get_parent() == piAcid:
            tempG, catNIndex = findInGraph( ligandGraph, catN, piAcidAtoms )
            if catNIndex == None:
                print(piAcid.get_resname(), piAcid.get_id(), piAcid.get_parent().get_id())
                continue
            
            if not  nx.has_path(ligandGraph, firstAtomInRing, catNIndex):
                continue 
            
            somethingFound = True
            newPath = nx.shortest_path(ligandGraph, firstAtomInRing, catNIndex)
            newPath = [ node for node in newPath if not node in centroidData["cycleAtoms"] ]
            
            if len(newPath) < len(shortestPath) or not shortestPath:
                shortestPath = newPath
                
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
    

if __name__ == "__main__":
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
    print(readResolutionAndMethod("cif/1bp0.cif"))
    print(readResolutionAndMethod("cif/1j5i.cif"))
    print(readResolutionAndMethod("cif/2eet.cif"))
    print(readResolutionAndMethod("cif/4c3y.cif"))
#findSupramolecular( "LUM", 666 )

# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 13:56:41 2017

Skrypt powstał we wspolpracy z najpiekniejsza kobieta na swiecie, Zrodzona
w Olkuszu, Pania mojego serca Emilia Kuzniak.

Sluzy do analizy pojedynczego pliku cif celem znalezienia potencjalnych struktur
zawierajacych oddzialywania supramolekularne anion-pi oraz obliczenie
wielkosci geometrycznych tychze czasteczek na potrzeby dalszej analizy
"""


from Bio.PDB import FastMMCIFParser, NeighborSearch, Selection
from primitiveCif2Dict import primitiveCif2Dict
import numpy as np                
from supramolecularLogging import writeAnionPiResults, incrementPartialProgress, writeAdditionalInfo
from supramolecularLogging import writeCationPiResults, writePiPiResults, writeAnionCationResults
from ringDetection import getRingsCentroids
from anionRecogniser import extractAnionAtoms, createResId
from multiprocessing import current_process

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
    parser = FastMMCIFParser()
    
    try:
        structure = parser.get_structure('temp', cifFile)
    except:
        errorMessage = "FastMMCIFParser cannot handle with: " + cifFile
        fileId = current_process()
        incrementPartialProgress(fileId)
        writeAdditionalInfo(errorMessage, fileId)
        #Zeby zobaczyc co sie dzieje
        return True        
        
    supramolecularFound = False 
    notPiacids = [ "HOH", "DOD", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "ILE", "LEU", "LYS", "MET", "PRO", "SER", "THR", "VAL" ] 
        
    resolution, method = readResolutionAndMethod(cifFile)
#    print("jade z pliku: ", cifFile)
    for modelIndex, model in enumerate(structure):
        atoms = Selection.unfold_entities(model, 'A')  
        ns = NeighborSearch(atoms)
           
        for residue in model.get_residues():
            residueName = residue.get_resname().upper()
            if not residueName in notPiacids  :
#                print("Analizuje: ", residueName)
                if analysePiacid(residue, PDBcode, modelIndex, ns, resolution, method):
                    supramolecularFound = True
            
    fileId = current_process()
    incrementPartialProgress(fileId)
    return supramolecularFound
    
def readResolutionAndMethod( cifFile ):
    try:
        mmcif_dict_parser = primitiveCif2Dict(cifFile, ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" , "_exptl.method"] )
        mmcif_dict = mmcif_dict_parser.result
    except:
        fileId = current_process()
        errorMessage = "primitiveCif2Dict cannot handle with: "+cifFile
        writeAdditionalInfo(errorMessage, fileId)
        return -666, -666
    
    method = "Unknown"
    if "_exptl.method" in mmcif_dict:
        method = mmcif_dict["_exptl.method"]
    res_keys = ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" ]
    res_list = []
    for key in res_keys:
        if key in mmcif_dict:
            resolution = mmcif_dict[key]
            
            try:
                resolution = float(resolution)
            except:
                return -2, method

            res_list.append(resolution)

                
    if not res_list:
        return -1, method
                
    variance = np.var(res_list)
    if variance < 0.001:
        return res_list[0], method
    else:
        return -3, method

def analysePiacid(ligand, PDBcode, modelIndex, ns, resolution, method):
    centroids = getRingsCentroids( ligand )
#    print("Znalazlem pierscienie w ilosci: ", len(centroids))
    
    ligandWithAnions = False
    fileId = current_process()
    for centroid in centroids:
        distance = 4.5
        neighbors = ns.search(np.array(centroid["coords"]), distance, 'A')
        extractedAnionAtoms = extractAnionAtoms( neighbors, ligand, ns )
            
        if len(extractedAnionAtoms) > 0:
            extractedCationAtoms = extractCationAtoms( centroid["coords"], ns, 10 )
            writeCationPiResults(ligand, PDBcode, centroid, extractedCationAtoms, modelIndex, fileId )
            
            extractedCentroids, ringMolecules = extractRingCentroids(centroid["coords"], ligand, ns)
            writePiPiResults(ligand, PDBcode, centroid, ringMolecules, extractedCentroids, modelIndex, fileId)
            
            for atom in extractedAnionAtoms:
                cationNearAnion = extractCationAtoms( atom["Atom"].get_coord(), ns, 4.5  )
                writeAnionCationResults(atom["Atom"], PDBcode, cationNearAnion, modelIndex, fileId)
            
        
        extractedAtoms =  writeAnionPiResults(ligand, PDBcode, centroid, extractedAnionAtoms, modelIndex, resolution, [], method, fileId)

        
        if len(extractedAtoms) > 0:
            ligandWithAnions = True
        
    return ligandWithAnions

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
